#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/Info.h"
#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/util/KokkosSession.h"
#include "polyhedralGravity/model/Polyhedron.h"


namespace py = pybind11;

namespace {

    using namespace polyhedralGravity;

    /**
     * Keeps a Python object alive for as long as an unmanaged Kokkos view aliases its buffer.
     *
     * A {@link Polyhedron} built from a NumPy array or a PyTorch tensor does not copy it, so the array has
     * to outlive the polyhedron and every GravityEvaluable sharing its mesh. Holding a reference to the
     * object achieves that without the user having to think about it.
     *
     * @param object the array whose buffer is handed over
     * @return an owner which releases the reference once the last view of the buffer is gone
     */
    std::shared_ptr<void> keepAlive(py::object object) {
        auto *reference = new py::object{std::move(object)};
        return std::shared_ptr<void>{reference, [](void *pointer) {
                                         // The last mesh may be released from C++ code without the GIL held
                                         const py::gil_scoped_acquire gil{};
                                         delete static_cast<py::object *>(pointer);
                                     }};
    }

    /**
     * A mesh buffer of an array library which lives on an accelerator, as described by its
     * @code __cuda_array_interface__ @endcode.
     */
    struct DeviceBuffer {
        /** The address of the first element in the device's memory */
        const void *data;
        /** The number of rows, i.e. the number of vertices or faces */
        size_t count;
        /** The NumPy type string of the elements, e.g. @code "<f4" @endcode */
        std::string typestr;
    };

    /**
     * Reads the @code __cuda_array_interface__ @endcode of an array which lives on an accelerator.
     *
     * This is the protocol PyTorch, JAX, CuPy, and Numba expose their device pointers through. Only
     * C-contiguous @f$(N, 3)@f$ arrays are accepted, since anything else would have to be copied and can be
     * made contiguous by the caller much more cheaply.
     *
     * @param array the array exposing the interface
     * @param name the name of the argument, for the error messages
     * @return the device pointer together with its shape and element type
     *
     * @throws std::invalid_argument if the array is not a C-contiguous (N, 3) array
     *
     * @note The caller has to make sure that the work producing the array has finished, e.g. by calling
     * torch.cuda.synchronize(), since the buffer is read without synchronizing on the producing stream.
     */
    DeviceBuffer readDeviceBuffer(const py::object &array, const std::string &name) {
        const auto interface = array.attr("__cuda_array_interface__").cast<py::dict>();
        const auto shape = interface["shape"].cast<std::vector<size_t>>();
        if (shape.size() != 2 || shape[1] != 3) {
            throw std::invalid_argument{"The " + name + " must be an (N, 3) array, but the given array has " +
                                        std::to_string(shape.size()) + " dimensions!"};
        }
        if (interface.contains("strides") && !interface["strides"].is_none()) {
            throw std::invalid_argument{"The " + name +
                                        " on the device must be C-contiguous, i.e. not a strided view. "
                                        "Call .contiguous() on it before handing it over!"};
        }
        const auto data = interface["data"].cast<py::tuple>();
        return {reinterpret_cast<const void *>(data[0].cast<uintptr_t>()), shape[0],
                interface["typestr"].cast<std::string>()};
    }

    /**
     * Builds a mesh from two device buffers, resolving the element type of the faces at runtime.
     *
     * @tparam VertexType the element type the vertex buffer has already been resolved to
     * @param vertices the vertex buffer
     * @param faces the face buffer
     * @param owner an owner of both buffers
     * @return the mesh aliasing both buffers
     *
     * @throws std::invalid_argument if the face buffer's element type is not a 32 or 64 bit integer
     */
    template<typename VertexType>
    PolyhedralMesh meshFromDeviceBuffers(const DeviceBuffer &vertices, const DeviceBuffer &faces,
                                         std::shared_ptr<void> owner) {
        const auto build = [&](const auto *typedFaces) {
            return PolyhedralMesh::fromBuffers(static_cast<const VertexType *>(vertices.data), vertices.count,
                                               typedFaces, faces.count, MemoryLocation::DEVICE, std::move(owner));
        };
        // The type string's first character is the byte order, which is always the platform's own for a
        // device buffer, so only the kind and the width matter
        const std::string kind = faces.typestr.substr(1);
        if (kind == "i4") {
            return build(static_cast<const int32_t *>(faces.data));
        }
        if (kind == "i8") {
            return build(static_cast<const int64_t *>(faces.data));
        }
        if (kind == "u4") {
            return build(static_cast<const uint32_t *>(faces.data));
        }
        if (kind == "u8") {
            return build(static_cast<const uint64_t *>(faces.data));
        }
        throw std::invalid_argument{"The faces on the device have the unsupported element type '" +
                                    faces.typestr + "'. Use a 32 or 64 bit integer type!"};
    }

    /**
     * Builds a mesh from the vertices and faces given by the user, copying as little as possible.
     *
     * An array which lives on an accelerator is recognized by its @code __cuda_array_interface__ @endcode
     * and its buffer is used where it is. Everything else is read through the buffer protocol, which does
     * not copy a C-contiguous NumPy array of the matching element type either.
     *
     * @param vertices the vertices as an @f$(N, 3)@f$ array-like
     * @param faces the faces as an @f$(M, 3)@f$ array-like
     * @return the mesh, aliasing the given arrays wherever possible
     */
    PolyhedralMesh makeMesh(const py::object &vertices, const py::object &faces) {
        const bool verticesOnDevice = py::hasattr(vertices, "__cuda_array_interface__");
        const bool facesOnDevice = py::hasattr(faces, "__cuda_array_interface__");
        if (verticesOnDevice != facesOnDevice) {
            throw std::invalid_argument{
                    "The vertices and the faces of a polyhedron must live in the same memory, but only one of "
                    "them was given as an array on the device!"};
        }
        if (verticesOnDevice) {
            const DeviceBuffer vertexBuffer = readDeviceBuffer(vertices, "vertices");
            const DeviceBuffer faceBuffer = readDeviceBuffer(faces, "faces");
            auto owner = keepAlive(py::make_tuple(vertices, faces));
            const std::string vertexKind = vertexBuffer.typestr.substr(1);
            if (vertexKind == "f4") {
                return meshFromDeviceBuffers<float>(vertexBuffer, faceBuffer, std::move(owner));
            }
            if (vertexKind == "f8") {
                return meshFromDeviceBuffers<double>(vertexBuffer, faceBuffer, std::move(owner));
            }
            throw std::invalid_argument{"The vertices on the device have the unsupported element type '" +
                                        vertexBuffer.typestr + "'. Use float32 or float64!"};
        }
        // pybind11 only copies here if the array is not already a C-contiguous array of the target type,
        // i.e. a float64 NumPy array and a uint64 index array are handed over as they are
        const py::array_t<double, py::array::c_style | py::array::forcecast> vertexArray{vertices};
        const py::array_t<size_t, py::array::c_style | py::array::forcecast> faceArray{faces};
        if (vertexArray.ndim() != 2 || vertexArray.shape(1) != 3) {
            throw std::invalid_argument{"The vertices must be an (N, 3) array-like!"};
        }
        if (faceArray.ndim() != 2 || faceArray.shape(1) != 3) {
            throw std::invalid_argument{"The faces must be an (M, 3) array-like!"};
        }
        return PolyhedralMesh::fromBuffers(vertexArray.data(), static_cast<size_t>(vertexArray.shape(0)),
                                           faceArray.data(), static_cast<size_t>(faceArray.shape(0)),
                                           MemoryLocation::HOST,
                                           keepAlive(py::make_tuple(vertexArray, faceArray)));
    }

    /**
     * Turns whatever the user passed as polyhedral_source into a mesh.
     *
     * @param polyhedralSource either a list of mesh file names or a pair of vertices and faces
     * @return the mesh of the polyhedron
     *
     * @throws std::invalid_argument if the source is neither of the two
     */
    PolyhedralMesh readPolyhedralSource(const py::object &polyhedralSource) {
        if (py::isinstance<py::str>(polyhedralSource)) {
            throw std::invalid_argument{
                    "The polyhedral_source must be a list of mesh file names or a pair of vertices and faces, "
                    "not a single string!"};
        }
        const auto source = py::cast<py::sequence>(polyhedralSource);
        if (py::len(source) > 0 && py::isinstance<py::str>(source[0])) {
            const auto [vertices, faces] = MeshReader::getPolyhedralSource(polyhedralSource.cast<PolyhedralFiles>());
            return PolyhedralMesh{vertices, faces};
        }
        if (py::len(source) != 2) {
            throw std::invalid_argument{
                    "The polyhedral_source must be a pair of vertices and faces, but a sequence of " +
                    std::to_string(py::len(source)) + " elements was given!"};
        }
        return makeMesh(source[0], source[1]);
    }

    /**
     * Exposes a host resident @f$(N, 3)@f$ view as a read-only NumPy array without copying it.
     *
     * @tparam Scalar the element type of the view
     * @param view the view to expose
     * @param owner the Python object keeping the view alive, i.e. the polyhedron
     * @return the NumPy array aliasing the view's memory
     */
    template<typename Scalar>
    py::array_t<Scalar> asReadOnlyArray(const kokkos::Vector3View<Scalar, kokkos::HostMemory> &view,
                                        const py::object &owner) {
        py::array_t<Scalar> array{{view.extent(0), view.extent(1)},
                                  {sizeof(Scalar) * 3, sizeof(Scalar)},
                                  view.data(),
                                  owner};
        array.attr("setflags")(py::arg("write") = false);
        return array;
    }

}// namespace

PYBIND11_MODULE(_core, m, py::mod_gil_not_used()) {
    using namespace polyhedralGravity;
    m.doc() = R"mydelimiter(
    The evaluation of the polyhedral gravity model requires the following parameters:

    +------------------------------------------------------------------------------+
    | Name                                                                         |
    +==============================================================================+
    | Polyhedral Mesh (either as vertices & faces or as polyhedral source files)   |
    +------------------------------------------------------------------------------+
    | Constant Density :math:`\rho`                                                |
    +------------------------------------------------------------------------------+

    In the Python Interface, you define these parameters as :py:class:`polyhedral_gravity.Polyhedron`

    .. code-block:: python

        from polyhedral_gravity import Polyhedron, GravityEvaluable, evaluate, PolyhedronIntegrity, NormalOrientation, ComputeBackend

        polyhedron = Polyhedron(
            polyhedral_source=(vertices, faces),           # (N,3) and (M,3) array-like
            density=density,                               # Density of the Polyhedron, Unit must match to the mesh's scale
            normal_orientation=NormalOrientation.OUTWARDS, # Possible values OUTWARDS (default) or INWARDS
            integrity_check=PolyhedronIntegrity.VERIFY,    # Possible values AUTOMATIC (default), VERIFY, HEAL, DISABLE
        )

    .. note::

        *Tsoulis et al.*'s formulation requires that the normals point :code:`OUTWARDS`.
        The implementation **can handle both cases and also can automatically determine the property** if initially set wrong.
        Using :code:`AUTOMATIC` (default for first-time-user) or :code:`VERIFY` raises a :code:`ValueError` if the :py:class:`polyhedral_gravity.NormalOrientation` is wrong.
        Using :code:`HEAL` will re-order the vertex sorting to fix errors.
        Using :code:`DISABLE` will turn this check off and avoid :math:`O(n^2)` runtime complexity of this check! Highly recommended, when you "know your mesh"!

    The polyhedron's mesh's units must match with the constant density!
    For example, if the mesh is in :math:`[m]`, then the constant density should be in :math:`[\frac{kg}{m^3}]`.

    Afterwards one can use :py:func:`polyhedral_gravity.evaluate` the gravity at a single point *P* via:

    .. code-block:: python

        potential, acceleration, tensor = evaluate(
            polyhedron=polyhedron,
            computation_points=P,
            backend=ComputeBackend.CPU_PARALLEL,
        )

    or via use the cached approach :py:class:`polyhedral_gravity.GravityEvaluable` (desirable for subsequent evaluations using the same :py:class:`polyhedral_gravity.Polyhedron`)

    .. code-block:: python

        evaluable = GravityEvaluable(
            polyhedron=polyhedron,
            backend=ComputeBackend.CPU_PARALLEL,
        )
        potential, acceleration, tensor = evaluable(computation_points=P)

    .. note::

        If :code:`P` would be an array of points, the return value would be a :code:`List[Tuple[potential, acceleration, tensor]]`!

    The calculation outputs the following parameters for every Computation Point *P*.
    The units of the respective output depend on the units of the input parameters (mesh and density)!

    Hence, if e.g. your mesh is in :math:`km`, the density must match. Further, the output units will match the input units.

    +------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+
    |         Name                                                                                   | If mesh :math:`[m]` and density :math:`[\frac{kg}{m^3}]`                   |                             Comment                             |
    +================================================================================================+============================================================================+=================================================================+
    |         :math:`V`                                                                              |  :math:`\frac{m^2}{s^2}` or :math:`\frac{J}{kg}`                           |           The potential or also called specific energy          |
    +------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+
    |     :math:`V_x`, :math:`V_y`, :math:`V_z`                                                      |   :math:`\frac{m}{s^2}`                                                    |The gravitational acceleration in the three cartesian directions |
    +------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+
    | :math:`V_{xx}`, :math:`V_{yy}`, :math:`V_{zz}`, :math:`V_{xy}`, :math:`V_{xz}`, :math:`V_{yz}` |   :math:`\frac{1}{s^2}`                                                    |The spatial rate of change of the gravitational acceleration     |
    +------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+

    This model's output obeys to the geodesy and geophysics sign conventions.
    Hence, the potential :math:`V` for a polyhedron with a mass :math:`m > 0` is defined as **positive**.

    The accelerations :math:`V_x`, :math:`V_y`, :math:`V_z` are defined as

    .. math::

        \textbf{g} = + \nabla V = \left( \frac{\partial V}{\partial x}, \frac{\partial V}{\partial y}, \frac{\partial V}{\partial z} \right)

    Accordingly, the second derivative tensor is defined as the derivative of :math:`\textbf{g}`.
    )mydelimiter";

    // We embedded the version and compilation information into the Python Interface
    m.attr("__version__") = POLYHEDRAL_GRAVITY_VERSION;
    m.attr("__parallelization__") = kokkos::getEnabledExecutionSpaces();
    m.attr("__commit__") = POLYHEDRAL_GRAVITY_COMMIT_HASH;
    m.attr("__logging__") = POLYHEDRAL_GRAVITY_LOGGING_LEVEL;
    m.attr("__host_compiler__") = POLYHEDRAL_GRAVITY_HOST_COMPILER;
    m.attr("__device_compiler__") = POLYHEDRAL_GRAVITY_DEVICE_COMPILER;

    py::enum_<NormalOrientation>(m, "NormalOrientation", R"mydelimiter(
        The orientation of the plane unit normals of the polyhedron.
        *Tsoulis et al.* equations require the normals to point outwards of the polyhedron.
        If the opposite hold, the result is negated.
        The implementation can handle both cases.
        )mydelimiter")
        .value("OUTWARDS", NormalOrientation::OUTWARDS, "Outwards pointing plane unit normals")
        .value("INWARDS", NormalOrientation::INWARDS, "Inwards pointing plane unit normals");

    py::enum_<PolyhedronIntegrity>(m, "PolyhedronIntegrity", R"mydelimiter(
        The pointing direction of the normals of a Polyhedron.
        They can either point outwards or inwards the polyhedron.
        )mydelimiter")
        .value("DISABLE", PolyhedronIntegrity::DISABLE,
               "All activities regarding MeshChecking are disabled. No runtime overhead!")
        .value("VERIFY", PolyhedronIntegrity::VERIFY,
               "Only verification of the NormalOrientation. "
               "A misalignment (e.g. specified OUTWARDS, but is not) leads to a runtime_error. Runtime Cost :math:`O(n^2)`")
        .value("AUTOMATIC", PolyhedronIntegrity::AUTOMATIC,
               "Like :code:`VERIFY`, but also informs the user about the option in any case on the runtime costs. "
               "This is the implicit default option. Runtime Cost: :math:`O(n^2)` and output to stdout in every case!")
        .value("HEAL", PolyhedronIntegrity::HEAL,
               "Verification and Automatic Healing of the NormalOrientation. "
               "A misalignment does not lead to a runtime_error, but to an internal correction of vertices ordering. Runtime Cost: :math:`O(n^2)`");

    py::enum_<ComputeBackend>(m, "ComputeBackend", R"mydelimiter(
        The compute backend on which the polyhedral gravity model is evaluated.
        Every backend is a `Kokkos <https://kokkos.org/>`__ execution space; which of them are compiled in
        is listed in :code:`polyhedral_gravity.__parallelization__`.
        )mydelimiter")
        .value("CPU_SERIAL", ComputeBackend::CPU_SERIAL,
               "Evaluation on a single CPU thread, i.e. the Kokkos Serial execution space")
        .value("CPU_PARALLEL", ComputeBackend::CPU_PARALLEL,
               "Evaluation on all CPU threads, i.e. the Kokkos OpenMP execution space. This is the default.")
        .value("GPU_PARALLEL", ComputeBackend::GPU_PARALLEL,
               "Evaluation on the GPU using the vendor's native paradigm, i.e. the Kokkos CUDA (NVIDIA), "
               "HIP (AMD), or SYCL (Intel) execution space. "
               "Raises a :code:`RuntimeError` if this build has no GPU backend.");

    py::enum_<ComputePrecision>(m, "ComputePrecision", R"mydelimiter(
        The floating point precision in which the polyhedral gravity model is evaluated.
        The polyhedron's mesh and the results are always given in double precision;
        this only selects the precision of the evaluation itself.
        )mydelimiter")
        .value("FLOAT32", ComputePrecision::FLOAT32,
               "Single precision :math:`[32\\ bit]`. Roughly halves the memory traffic of the evaluation, but the "
               "polyhedral gravity model cancels large terms against each other, so expect a relative accuracy "
               "of only about :math:`10^{-4}`.")
        .value("FLOAT64", ComputePrecision::FLOAT64, "Double precision :math:`[64\\ bit]`. This is the default.");

    py::enum_<MetricUnit>(m, "MetricUnit", R"mydelimiter(
        The metric unit of for example a polyhedral mesh source.
        )mydelimiter")
    .value("METER", MetricUnit::METER, "Representing meter :math:`[m]`")
    .value("KILOMETER", MetricUnit::KILOMETER, "Representing kilometer :math:`[km]`")
    .value("UNITLESS", MetricUnit::UNITLESS, "Representing no unit :math:`[1]`");

    py::class_<Polyhedron>(m, "Polyhedron", R"mydelimiter(
            A constant density Polyhedron stores the mesh data consisting of vertices and triangular faces.

            The density and the coordinate system in which vertices and faces are defined need to have the same scale/ units.
            The vertices are indexed starting with zero, not one. If the polyhedral source starts indexing with one, the counting is shifted by -1.

            Tsoulis et al.'s polyhedral gravity model requires that the plane unit normals of every face are pointing outwards
            of the polyhedron. Otherwise the results are negated.
            The class by default enforces this constraints and offers utility to (automatically) make the input data obey to this constraint.
            )mydelimiter")
            .def(py::init([](const py::object &polyhedralSource, const double density,
                             const NormalOrientation &normalOrientation, const PolyhedronIntegrity &integrityCheck,
                             const MetricUnit &metricUnit) {
                     return Polyhedron{readPolyhedralSource(polyhedralSource), density, normalOrientation,
                                       integrityCheck, metricUnit};
                 }),
                 R"mydelimiter(
            Creates a new Polyhedron from vertices and faces and a constant density.
            If the integrity_check is not set to DISABLE, the mesh integrity is checked
            (so that it fits the specification of the polyhedral model by *Tsoulis et al.*)

            Args:
                polyhedral_source:  The vertices (:math:`(N, 3)`-array-like) and faces (:math:`(M, 3)`-array-like) of the polyhedron as pair or
                                    The filenames of the files containing the vertices & faces as list of strings
                density:            The constant density of the polyhedron, it must match the mesh's units, e.g. mesh in :math:`[m]` then density in :math:`[kg/m^3]`
                                    or mesh in :math:`[km]` then density in :math:`[kg/km^3]`
                normal_orientation: The pointing direction of the mesh's plane unit normals, i.e., either :code:`OUTWARDS` or :code:`INWARDS` of the polyhedron.
                                    One of :py:class:`polyhedral_gravity.NormalOrientation`.
                                    (default: :code:`OUTWARDS`)
                integrity_check:    Conducts an Integrity Check (degenerated faces/ vertex ordering) depending on the values. One of :py:class:`polyhedral_gravity.PolyhedronIntegrity`:

                                        * :code:`AUTOMATIC` (Default): Prints to stdout and throws ValueError if normal_orientation is wrong/ inconsistent
                                        * :code:`VERIFY`: Like :code:`AUTOMATIC`, but does not print to stdout
                                        * :code:`DISABLE`: Recommend, when you are familiar with the mesh to avoid :math:`O(n^2)` runtime cost. Disables ALL checks
                                        * :code:`HEAL`: Automatically fixes the normal_orientation and vertex ordering to the correct values
                metric_unit:        The metric unit of the mesh. Can be either :code:`METER`, :code:`KILOMETER`, or :code:`UNITLESS`.
                                    (default: :code:`METER`)

            Raises:
                ValueError: If :code:`integrity_check` is set to :code:`AUTOMATIC` or :code:`VERIFY` and the mesh is inconsistent
                RuntimeError: If files given as :code:`polyhedral_source` do not exist

            Note:
                The :code:`integrity_check` is automatically enabled to avoid wrong results due to the wrong vertex ordering.
                The check requires :math:`O(n^2)` operations. You want to turn this off, when you know you mesh!
                The faces array's indexing is shifted by -1 if the indexing started previously from vertex one (i.e., the first index is referred to as one).
                In other words, the first vertex is always referred to as vertex zero not one!

            Note:
                The vertices and faces are handed over through the buffer protocol, i.e. they are **not copied**
                if they already are C-contiguous arrays of :code:`float64` (vertices) and :code:`uint64` (faces).
                Any other :code:`dtype` or a non-contiguous array costs one conversion.
                The polyhedron keeps a reference to the arrays alive for as long as it needs them.

                An array which lives on an accelerator, i.e. a CUDA :code:`torch.Tensor` or a GPU
                :code:`jax.Array`, is recognized by its :code:`__cuda_array_interface__` and is used **where it
                is**, so its mesh never travels through the host. Make sure the work producing such an array has
                finished (e.g. :code:`torch.cuda.synchronize()`) before handing it over, and note that a build
                without a GPU backend raises a :code:`RuntimeError` for it.
            )mydelimiter",
                 py::arg("polyhedral_source"),
                 py::arg("density"),
                 py::arg("normal_orientation") = NormalOrientation::OUTWARDS,
                 py::arg("integrity_check") = PolyhedronIntegrity::AUTOMATIC,
                 py::arg("metric_unit") = MetricUnit::METER
                    )
            .def("check_normal_orientation", &Polyhedron::checkPlaneUnitNormalOrientation, R"mydelimiter(
            Returns a tuple consisting of majority plane unit normal orientation,
            i.e. the direction in which at least more than half of the plane unit normals point,
            and the indices of the faces violating this orientation, i.e. the faces whose plane unit normals point in the other direction.
            The set of indices violating the property is empty if the mesh has a clear ordering.
            The set contains values if the mesh is inconsistent.

            Returns:
                Tuple consisting consisting of majority plane unit normal orientation and the indices of the faces violating this orientation.

            Note:
                This utility is mainly for diagnostics and debugging purposes. If the polyhedron is constructed with `integrity_check`
                set to :code:`AUTOMATIC` or :code:`VERIFY`, the construction fails anyways.
                If set to :code:`HEAL`, this method should return an empty set (but maybe a different ordering than initially specified)
                Only if set to :code:`DISABLE`, then this method might actually return a set with faulty indices.
                Hence, if you want to know your mesh error. Construct the polyhedron with :code:`integrity_check=DISABLE` and call this method.
            )mydelimiter")
            .def("__getitem__", static_cast<Array3Triplet (Polyhedron::*)(size_t) const>(&Polyhedron::getResolvedFace), R"mydelimiter(
            Returns the the three coordinates of the vertices making the face at the requested index.
            This does not return the face as list of vertex indices, but resolved with the actual coordinates.

            Args:
                index:  The index of the face

            Returns:
                :math:`(3, 3)`-array-like: The resolved face

            Raises:
                IndexError if face index is out-of-bounds
            )mydelimiter", py::arg("index"))
            .def("__repr__", &Polyhedron::toString, R"mydelimiter(
            :py:class:`str`: A string representation of this polyhedron
            )mydelimiter")
            .def_property_readonly("vertices", [](const py::object &self) {
                        return asReadOnlyArray(self.cast<const Polyhedron &>().getMesh().getHostMesh().vertices, self);
                    }, R"mydelimiter(
            (N, 3) :py:class:`numpy.ndarray` of :py:class:`numpy.float64`: The vertices of the polyhedron.
            Coordinates in the unit of the mesh (Read-Only).

            This is a view of the polyhedron's own memory and therefore not writeable. If the polyhedron's mesh
            lives on an accelerator, reading this property downloads it to the host once.
            )mydelimiter")
            .def_property_readonly("faces", [](const py::object &self) {
                        return asReadOnlyArray(self.cast<const Polyhedron &>().getMesh().getHostMesh().faces, self);
                    }, R"mydelimiter(
            (M, 3) :py:class:`numpy.ndarray` of :py:class:`numpy.uint64`: The faces of the polyhedron (Read-Only).

            This is a view of the polyhedron's own memory and therefore not writeable. If the polyhedron's mesh
            lives on an accelerator, reading this property downloads it to the host once.
            )mydelimiter")
            .def_property("density", &Polyhedron::getDensity, &Polyhedron::setDensity, R"mydelimiter(
            :py:class:`float`: The density of the polyhedron in :math:`[kg/X^3]` with X being the unit of the mesh (Read/ Write).
            )mydelimiter")
            .def_property_readonly("normal_orientation", &Polyhedron::getOrientation, R"mydelimiter(
            :py:class:`polyhedral_gravity.NormalOrientation`: The orientation of the plane unit normals (Read-Only).
            )mydelimiter")
            .def_property_readonly("mesh_unit", &Polyhedron::getMeshUnitAsString, R"mydelimiter(
            :py:class:`str`: The metric unit of the polyhedral mesh (Read-Only).
            )mydelimiter")
            .def_property_readonly("density_unit", &Polyhedron::getDensityUnit, R"mydelimiter(
            :py:class:`str`: The metric unit of the density (Read-Only).
            )mydelimiter")
            .def(py::pickle(
                    [](const Polyhedron &polyhedron) {
                        const auto &[vertices, faces, density, orientation, metricUnit] = polyhedron.getState();
                        return py::make_tuple(vertices, faces, density, orientation, metricUnit);
                    },
                    [](const py::tuple &tuple) {
                        constexpr static size_t POLYHEDRON_STATE_SIZE = 5;
                        if (tuple.size() != POLYHEDRON_STATE_SIZE) {
                            throw std::runtime_error("Invalid state!");
                        }
                        Polyhedron polyhedron{
                                tuple[0].cast<std::vector<Array3>>(), tuple[1].cast<std::vector<IndexArray3>>(),
                                tuple[2].cast<double>(), tuple[3].cast<NormalOrientation>(), PolyhedronIntegrity::DISABLE,
                                tuple[4].cast<MetricUnit>()
                        };
                        return polyhedron;
                    }
                    ));

    py::class_<GravityEvaluable>(m, "GravityEvaluable", R"mydelimiter(
             A class to evaluate the polyhedral gravity model for a given constant density polyhedron at a given computation point.
             It provides a :py:meth:`polyhedral_gravity.GravityEvaluable.__call__` method to evaluate the polyhedral gravity model for computation points while
             also caching the polyhedron & intermediate results over the lifetime of the object.
             )mydelimiter")
            .def(py::init<const Polyhedron &, const ComputeBackend &, const ComputePrecision &>(),R"mydelimiter(
             Creates a new GravityEvaluable for a given constant density polyhedron.
             It provides a :py:meth:`polyhedral_gravity.GravityEvaluable.__call__` method to evaluate the polyhedral gravity model for computation points while
             also caching the polyhedron & intermediate results over the lifetime of the object.

             The compute backend is fixed here rather than per call, so that the cached properties stay in the
             memory of that backend for the whole lifetime of this GravityEvaluable.

             Args:
                 polyhedron: The polyhedron for which to evaluate the gravity model
                 backend:    The compute backend on which every evaluation runs (default: :code:`ComputeBackend.CPU_PARALLEL`)
                 precision:  The floating point precision of the evaluation (default: :code:`ComputePrecision.FLOAT64`)

             Raises:
                 RuntimeError: If :code:`ComputeBackend.GPU_PARALLEL` is requested, but this build has no GPU backend
             )mydelimiter", py::arg("polyhedron"), py::arg("backend") = ComputeBackend::CPU_PARALLEL,
                 py::arg("precision") = ComputePrecision::FLOAT64)
            .def_property_readonly("backend", &GravityEvaluable::getComputeBackend,R"mydelimiter(
            :py:class:`polyhedral_gravity.ComputeBackend`: The compute backend every evaluation of this GravityEvaluable runs on (Read-Only).
            )mydelimiter")
            .def_property_readonly("precision", &GravityEvaluable::getComputePrecision,R"mydelimiter(
            :py:class:`polyhedral_gravity.ComputePrecision`: The floating point precision of this GravityEvaluable (Read-Only).
            )mydelimiter")
            .def_property_readonly("output_units", &GravityEvaluable::getOutputMetricUnit,R"mydelimiter(
            (3)-array-like of :py:class:`str`: A human-readable string representation of the output units. This depends on the polyhedron's definition (Read-Only).
            )mydelimiter")
            .def("__repr__", &GravityEvaluable::toString,R"mydelimiter(
            :py:class:`str`: A string representation of this GravityEvaluable.
            )mydelimiter")
            .def("__call__", &GravityEvaluable::operator(),
             R"mydelimiter(
             Evaluates the polyhedral gravity model for a given constant density polyhedron at a given computation point.

             The results' units depend on the polyhedron's input units.
             For example, if the polyhedral mesh is in :math:`[m]` and the density in :math:`[kg/m^3]`, then the potential is in :math:`[m^2/s^2]`.
             In case the polyhedron is unitless, the results are **not** multiplied with the Gravitational Constant :math:`G`, but returned raw.

             The evaluation runs on the compute backend this GravityEvaluable was created for.

             Args:
                 computation_points: The computation points as tuple or list of points

             Returns:
                 Either a triplet of potential :math:`V`, acceleration :math:`[V_x, V_y, V_z]`
                 and second derivatives :math:`[V_{xx}, V_{yy}, V_{zz}, V_{xy},V_{xz}, V_{yz}]` at the computation points or
                 if multiple computation points are given a list of these triplets
             )mydelimiter", py::arg("computation_points"))
            .def(py::pickle(
                    [](const GravityEvaluable &evaluable) {
                        const auto &[polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals] = evaluable.getState();
                        return py::make_tuple(polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals,
                                              evaluable.getComputeBackend(), evaluable.getComputePrecision());
                    },
                    [](const py::tuple &tuple) {
                        constexpr size_t GRAVITY_EVALUABLE_STATE_SIZE = 6;
                        if (tuple.size() != GRAVITY_EVALUABLE_STATE_SIZE) {
                            throw std::runtime_error("Invalid state!");
                        }
                        GravityEvaluable evaluable{
                                tuple[0].cast<Polyhedron>(), tuple[1].cast<std::vector<Array3Triplet>>(),
                                tuple[2].cast<std::vector<Array3>>(), tuple[3].cast<std::vector<Array3Triplet>>(),
                                tuple[4].cast<ComputeBackend>(), tuple[5].cast<ComputePrecision>()
                        };
                        return evaluable;
                    }
                    ));

    m.def("evaluate", [](const Polyhedron &polyhedron,
                         const std::variant<Array3, std::vector<Array3>> &computationPoints,
                         ComputeBackend backend,
                         ComputePrecision precision) -> std::variant<GravityModelResult, std::vector<GravityModelResult>> {
                    return std::visit(util::overloaded{
                            [&](const Array3 &point) {
                                return std::variant<GravityModelResult, std::vector<GravityModelResult>>(GravityModel::evaluate(polyhedron, point, backend, precision));
                            },
                            [&](const std::vector<Array3> &points) {
                                return std::variant<GravityModelResult, std::vector<GravityModelResult>>(GravityModel::evaluate(polyhedron, points, backend, precision));
                            }
                        }, computationPoints);
          }, R"mydelimiter(
             Evaluates the polyhedral gravity model for a given constant density polyhedron at a given computation point.

             The results' units depend on the polyhedron's input units.
             For example, if the polyhedral mesh is in :math:`[m]` and the density in :math:`[kg/m^3]`, then the potential is in :math:`[m^2/s^2]`.
             In case the polyhedron is unitless, the results are **not** multiplied with the Gravitational Constant :math:`G`, but returned raw.

             Args:
                 polyhedron:            The polyhedron for which to evaluate the gravity model
                 computation_points:    The computation points as tuple or list of points
                 backend:               The compute backend on which to evaluate
                                        (default: :code:`ComputeBackend.CPU_PARALLEL`)
                 precision:             The floating point precision of the evaluation
                                        (default: :code:`ComputePrecision.FLOAT64`)

             Returns:
                 Either a triplet of potential :math:`V`, acceleration :math:`[V_x, V_y, V_z]`
                 and second derivatives :math:`[V_{xx}, V_{yy}, V_{zz}, V_{xy},V_{xz}, V_{yz}]` at the computation points or
                 if multiple computation points are given a list of these triplets

             Raises:
                 RuntimeError: If :code:`ComputeBackend.GPU_PARALLEL` is requested, but this build has no GPU backend
             )mydelimiter", py::arg("polyhedron"), py::arg("computation_points"),
             py::arg("backend") = ComputeBackend::CPU_PARALLEL,
             py::arg("precision") = ComputePrecision::FLOAT64);

}