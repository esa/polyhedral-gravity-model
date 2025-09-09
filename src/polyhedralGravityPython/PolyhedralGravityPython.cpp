#include <tuple>
#include <variant>
#include <string>
#include <array>
#include <vector>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/Info.h"
#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/Polyhedron.h"


namespace py = pybind11;

PYBIND11_MODULE(polyhedral_gravity, m) {
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

        from polyhedral_gravity import Polyhedron, GravityEvaluable, evaluate, PolyhedronIntegrity, NormalOrientation

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
            parallel=True,
        )

    or via use the cached approach :py:class:`polyhedral_gravity.GravityEvaluable` (desirable for subsequent evaluations using the same :py:class:`polyhedral_gravity.Polyhedron`)

    .. code-block:: python

        evaluable = GravityEvaluable(polyhedron=polyhedron)
        potential, acceleration, tensor = evaluable(
            computation_points=P,
            parallel=True,
        )

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
    m.attr("__parallelization__") = POLYHEDRAL_GRAVITY_PARALLELIZATION;
    m.attr("__commit__") = POLYHEDRAL_GRAVITY_COMMIT_HASH;
    m.attr("__logging__") = POLYHEDRAL_GRAVITY_LOGGING_LEVEL;

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
            .def(py::init<const std::variant<PolyhedralSource, PolyhedralFiles> &, double, const NormalOrientation &, const PolyhedronIntegrity &, const MetricUnit &>(), R"mydelimiter(
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
            .def("__getitem__", &Polyhedron::getResolvedFace, R"mydelimiter(
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
            .def_property_readonly("vertices", &Polyhedron::getVertices, R"mydelimiter(
            (N, 3)-array-like of :py:class:`float`: The vertices of the polyhedron. Coordinates in the unit of the mesh (Read-Only).
            )mydelimiter")
            .def_property_readonly("faces", &Polyhedron::getFaces, R"mydelimiter(
            (M, 3)-array-like of :py:class:`int`: The faces of the polyhedron (Read-Only).
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
            .def(py::init<const Polyhedron &>(),R"mydelimiter(
             Creates a new GravityEvaluable for a given constant density polyhedron.
             It provides a :py:meth:`polyhedral_gravity.GravityEvaluable.__call__` method to evaluate the polyhedral gravity model for computation points while
             also caching the polyhedron & intermediate results over the lifetime of the object.

             Args:
                 polyhedron: The polyhedron for which to evaluate the gravity model
             )mydelimiter", py::arg("polyhedron"))
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

             Args:
                 computation_points: The computation points as tuple or list of points
                 parallel:           If :code:`True`, the computation is done in parallel on the CPU using the technology specified by
                                     :code:`polyhedral_gravity.__parallelization__` (default: :code:`True`)

             Returns:
                 Either a triplet of potential :math:`V`, acceleration :math:`[V_x, V_y, V_z]`
                 and second derivatives :math:`[V_{xx}, V_{yy}, V_{zz}, V_{xy},V_{xz}, V_{yz}]` at the computation points or
                 if multiple computation points are given a list of these triplets
             )mydelimiter", py::arg("computation_points"), py::arg("parallel") = true)
            .def(py::pickle(
                    [](const GravityEvaluable &evaluable) {
                        const auto &[polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals] = evaluable.getState();
                        return py::make_tuple(polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals);
                    },
                    [](const py::tuple &tuple) {
                        constexpr size_t GRAVITY_EVALUABLE_STATE_SIZE = 4;
                        if (tuple.size() != GRAVITY_EVALUABLE_STATE_SIZE) {
                            throw std::runtime_error("Invalid state!");
                        }
                        GravityEvaluable evaluable{
                                tuple[0].cast<Polyhedron>(), tuple[1].cast<std::vector<Array3Triplet>>(),
                                tuple[2].cast<std::vector<Array3>>(), tuple[3].cast<std::vector<Array3Triplet>>()
                        };
                        return evaluable;
                    }
                    ));

    m.def("evaluate", [](const Polyhedron &polyhedron,
                         const std::variant<Array3, std::vector<Array3>> &computationPoints,
                         bool parallel) -> std::variant<GravityModelResult, std::vector<GravityModelResult>> {
                    return std::visit(util::overloaded{
                            [&](const Array3 &point) {
                                return std::variant<GravityModelResult, std::vector<GravityModelResult>>(GravityModel::evaluate(polyhedron, point, parallel));
                            },
                            [&](const std::vector<Array3> &points) {
                                return std::variant<GravityModelResult, std::vector<GravityModelResult>>(GravityModel::evaluate(polyhedron, points, parallel));
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
                 parallel:              If :code:`True`, the computation is done in parallel on the CPU using the technology specified by
                                        :code:`polyhedral_gravity.__parallelization__` (default: :code:`True`)

             Returns:
                 Either a triplet of potential :math:`V`, acceleration :math:`[V_x, V_y, V_z]`
                 and second derivatives :math:`[V_{xx}, V_{yy}, V_{zz}, V_{xy},V_{xz}, V_{yz}]` at the computation points or
                 if multiple computation points are given a list of these triplets
             )mydelimiter", py::arg("polyhedron"), py::arg("computation_points"), py::arg("parallel") = true);

}