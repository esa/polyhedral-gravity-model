#include <tuple>
#include <variant>
#include <string>
#include <array>
#include <vector>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/input/TetgenAdapter.h"


namespace py = pybind11;

PYBIND11_MODULE(polyhedral_gravity, m) {
    using namespace polyhedralGravity;
    m.doc() = "Computes the full gravity tensor for a given constant density polyhedron which consists of some "
            "vertices and triangular faces at a given computation points accoring to the geodetic convention";

    py::enum_<NormalOrientation>(m, "NormalOrientation", R"mydelimiter(
        The orientation of the plane unit normals of the polyhedron.
        Tsoulis et al. equations require the normals to point outwards of the polyhedron.
        If the opposite hold, the result is negated.
        The implementation can handle both cases.
        )mydelimiter")
            .value("OUTWARDS", NormalOrientation::OUTWARDS, "Outwards pointing plane unit normals")
            .value("INWARDS", NormalOrientation::INWWARDS, "Inwards pointing plane unit normals");

    py::enum_<PolyhedronIntegrity>(m, "PolyhedronIntegrity", R"mydelimiter(
        The pointing direction of the normals of a Polyhedron.
        They can either point outwards or inwards the polyhedron.
        )mydelimiter")
            .value("DISABLE", PolyhedronIntegrity::DISABLE,
                   "All activities regarding MeshChecking are disabled. No runtime overhead!")
            .value("VERIFY", PolyhedronIntegrity::VERIFY,
                   "Only verification of the NormalOrientation. "
                   "A misalignment (e.g. specified OUTWARDS, but is not) leads to a runtime_error. Runtime Cost O(n^2)")
            .value("AUTOMATIC", PolyhedronIntegrity::AUTOMATIC,
                   "Like VERIFY, but also informs the user about the option in any case on the runtime costs. "
                   "This is the implicit default option. Runtime Cost: O(n^2) and output to stdout in every case!")
            .value("HEAL", PolyhedronIntegrity::HEAL,
                   "Verification and Autmatioc Healing of the NormalOrientation. "
                   "A misalignemt does not lead to a runtime_error, but to an internal correction of vertices ordering. Runtime Cost: @f$O(n^2)$@f");

    py::class_<Polyhedron>(m, "Polyhedron", R"mydelimiter(
            A constant density Polyhedron stores the mesh data consisting of vertices and triangular faces.
            The density and the coordinate system in which vertices and faces are defined need to have the same scale/ units.
            Tsoulis et al.'s polyhedral gravity model requires that the plane unit normals of every face are pointing outwards
            of the polyhedron. Otherwise the results are negated.
            The class by default enforces this constraints and offers utility to (automatically) make the input data obey to this constraint.
            )mydelimiter")
            .def(py::init<const PolyhedralSource &, double, const NormalOrientation &, const PolyhedronIntegrity &>(), R"mydelimiter(
            Creates a new Polyhedron from vertices and faces and a constant density.
            If the integrity_check is not set to DISABLE, the mesh integrity is checked
            (so that it fits the specification of the polyhedral model by Tsoulis et al.)

            Args:
                polyhedral_source:  The vertices & faces of the polyhedron as pair of (N, 3)-arrays
                density:            The constant density of the polyhedron, it must match the mesh's units, e.g. mesh in [m] then density in [kg/m^3]
                normal_orientation: The pointing direction of the mesh's plane unit normals, i.e., either OUTWARDS or INWARDS of the polyhedron
                                    (default: OUTWARDS)
                integrity_check:    Conducts an Integrity Check (degenerated faces/ vertex ordering) depending on the values. Options
                                    - AUTOMATIC (Default): Prints to stdout and throws ValueError if normal_orientation is wrong/ inconsisten
                                    - VERIFY: Like AUTOMATIC, but does not print to stdout
                                    - DISABLE: Recommened, when you know the mesh to avoid to pay O(n^2) runtime. Disables ALL checks
                                    - HEAL: Automatically fixes the normal_orientation and vertex ordering to the correct values
            Raises:
                ValueError if `integrity_check` is set to AUTOMATIC or VERIFY and the mesh is inconsistent

            Note:
                The integrity_check is automatically enabled to avoid wrong results due to the wrong vertex ordering.
                The check requires O(n^2) operations. You want to turn this off, when you know you mesh!
            )mydelimiter",
                 py::arg("polyhedral_source"),
                 py::arg("density"),
                 py::arg("normal_orientation") = NormalOrientation::OUTWARDS,
                 py::arg("integrity_check") = PolyhedronIntegrity::AUTOMATIC
                    )
            .def(py::init<const PolyhedralFiles &, double, const NormalOrientation &, const PolyhedronIntegrity &>(), R"mydelimiter(
            Creates a new Polyhedron from vertices and faces and a constant density.
            If the integrity_check is not set to DISABLE, the mesh integrity is checked
            (so that it fits the specification of the polyhedral model by Tsoulis et al.)

            Args:
                polyhedral_source:  The filenames of the files containing the vertices & faces as list of strings
                density:            The constant density of the polyhedron, it must match the mesh's units, e.g. mesh in [m] then density in [kg/m^3]
                normal_orientation: The pointing direction of the mesh's plane unit normals, i.e., either OUTWARDS or INWARDS of the polyhedron
                                    (default: OUTWARDS)
                integrity_check:    Conducts an Integrity Check (degenerated faces/ vertex ordering) depending on the values. Options
                                    - AUTOMATIC (Default): Prints to stdout and throws ValueError if normal_orientation is wrong/ inconsisten
                                    - VERIFY: Like AUTOMATIC, but does not print to stdout
                                    - DISABLE: Recommened, when you know the mesh to avoid to pay O(n^2) runtime. Disables ALL checks
                                    - HEAL: Automatically fixes the normal_orientation and vertex ordering to the correct values
            Raises:
                ValueError if `integrity_check` is set to AUTOMATIC or VERIFY and the mesh is inconsistent

            Note:
                The integrity_check is automatically enabled to avoid wrong results due to the wrong vertex ordering.
                The check requires O(n^2) operations. You want to turn this off, when you know you mesh!
            )mydelimiter",
                 py::arg("polyhedral_source"),
                 py::arg("density"),
                 py::arg("normal_orientation") = NormalOrientation::OUTWARDS,
                 py::arg("integrity_check") = PolyhedronIntegrity::AUTOMATIC
                    )
            .def("check_normal_orientation", &Polyhedron::checkPlaneUnitNormalOrientation, R"mydelimiter(
            Returns a tuple consisting of majority plane unit normal orientation,
            i.e. the direction in which at least more than half of the plane unit normals point,
            and the indices of the faces violating this orientation, i.e. the faces whose plane unit normals point in the other direction.
            The set of incides vioalting the property is empty if the mesh has a clear ordering.
            The set contains values if the mesh is incosistent.

            Returns:
                Tuple consisting consisting of majority plane unit normal orientation and the indices of the faces violating this orientation.

            Note:
                This utility is mainly for diagnostics and debugging purposes. If the polyhedron is constrcuted with `integrity_check`
                set to AUTOMATIC or VERIFY, the construction fails anyways.
                If set to HEAL, this method should return an empty set (but maybe a different ordering than initially specified)
                Only if set to DISABLE, then this method might actually return a set with faulty indices.
                Hence, if you want to know your mesh error. Construct the polyhedron with `integrity_check=DISABLE` and call this method.
            )mydelimiter")
            .def("__getitem__", &Polyhedron::getResolvedFace, R"mydelimiter(
            Returns the the three cooridnates of the vertices making the face at the requested index.

            Args:
                index:  The index of the face

            Returns:
                The resolved face (so the vertices as coordinates) as (3, 3)-array

            Raises:
                IndexError if face index is out-of-bounds
            )mydelimiter", py::arg("index"))
            .def("__repr__", &Polyhedron::toString, "Returns a string representation of this polyhedron")
            .def_property_readonly("vertices", &Polyhedron::getVertices, "The vertices (type: floats) of the polyhedron as (N, 3)-array")
            .def_property_readonly("faces", &Polyhedron::getFace, "The faces (type: int) of the polyhedron as (M, 3) array")
            .def_property("density", &Polyhedron::getDensity, &Polyhedron::setDensity, "The density (mutable) of the polyhedron")
            .def_property_readonly("normal_orientation", &Polyhedron::getOrientation, "The plane unit normal orientation of the polyhedron")
            .def(py::pickle(
                    [](const Polyhedron &polyhedron) {
                        const auto &[vertices, faces, density, orientation] = polyhedron.getState();
                        return py::make_tuple(vertices, faces, density, orientation);
                    },
                    [](const py::tuple &tuple) {
                        constexpr size_t tupleSize = 4;
                        if (tuple.size() != tupleSize) {
                            throw std::runtime_error("Invalid state!");
                        }
                        Polyhedron polyhedron{
                                tuple[0].cast<std::vector<Array3>>(), tuple[1].cast<std::vector<IndexArray3>>(),
                                tuple[2].cast<double>(), tuple[3].cast<NormalOrientation>()
                        };
                        return polyhedron;
                    }
                    ));

    py::class_<GravityEvaluable>(m, "GravityEvaluable", R"mydelimiter(
             A class to evaluate the polyhedral gravity model for a given constant density polyhedron at a given computation point.
             It provides a __call__ method to evaluate the polyhedral gravity model for computation points while
             also caching the polyhedron & intermediate results over the lifetime of the object.
             )mydelimiter")
            .def(py::init<const Polyhedron &>(),R"mydelimiter(
             Creates a new GravityEvaluable for a given constant density polyhedron.
             It provides a __call__ method to evaluate the polyhedral gravity model for computation points while
             also caching the polyhedron & intermediate results over the lifetime of the object.

             Args:
                 polyhedron: The polyhedron for which to evaluate the gravity model
             )mydelimiter", py::arg("polyhedron"))
            .def("__repr__", &GravityEvaluable::toString, "Returns a string representation of this GravityEvaluable")
            .def("__call__", &GravityEvaluable::operator(),
             R"mydelimiter(
             Evaluates the polyhedral gravity model for a given constant density polyhedron at a given computation point.

             Args:
                 computation_points: The computation points as tuple or list of points
                 parallel:           If true, the computation is done in parallel (default: true)

             Returns:
                 Either a tuple of potential, acceleration and second derivatives at the computation points or
                 if multiple computation points are given, the result is a list of tuples
             )mydelimiter", py::arg("computation_points"), py::arg("parallel") = true)
            .def(py::pickle(
                    [](const GravityEvaluable &evaluable) {
                        const auto &[polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals] = evaluable.getState();
                        return py::make_tuple(polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals);
                    },
                    [](const py::tuple &tuple) {
                        constexpr size_t tupleSize = 4;
                        if (tuple.size() != tupleSize) {
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

             Args:
                 polyhedron:            The polyhedron for which to evaluate the gravity model
                 computation_points:    The computation points as tuple or list of points
                 parallel:              If true, the computation is done in parallel (default: true)

             Returns:
                 Either a tuple of potential, acceleration and second derivatives at the computation points or
                 if multiple computation points are given, the result is a list of tuples
             )mydelimiter", py::arg("polyhedron"), py::arg("computation_points"), py::arg("parallel") = true);

}