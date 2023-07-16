#include <tuple>
#include <variant>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/calculation/GravityEvaluable.h"
#include "polyhedralGravity/calculation/MeshChecking.h"
#include "polyhedralGravity/input/TetgenAdapter.h"


namespace py = pybind11;

PYBIND11_MODULE(polyhedral_gravity, m) {
    using namespace polyhedralGravity;
    m.doc() = "Computes the full gravity tensor for a given constant density polyhedron which consists of some "
              "vertices and triangular faces at a given computation points";

    m.def("evaluate", [](const PolyhedralSource &polyhedralSource, double density,
                         const std::variant<Array3, std::vector<Array3>> &computationPoints,
                         bool parallel) -> std::variant<GravityModelResult, std::vector<GravityModelResult>> {
              if (std::holds_alternative<Array3>(computationPoints)) {
                  return std::visit([&density, &computationPoints, parallel](const auto &source) {
                      return GravityModel::evaluate(source, density, std::get<Array3>(computationPoints), parallel);
                  }, polyhedralSource);
              } else {
                  return std::visit([&density, &computationPoints, parallel](const auto &source) {
                      return GravityModel::evaluate(source, density, std::get<std::vector<Array3>>(computationPoints),
                                                    parallel);
                  }, polyhedralSource);
              }
          }, R"mydelimiter(
            Evaluates the polyhedral gravity model for a given constant density polyhedron at a given computation point.

            Args:
                polyhedral_source:  The vertices & faces of the polyhedron as tuple or
                                    the filenames of the files containing the vertices & faces as list of strings
                density:            The constant density of the polyhedron in [kg/m^3]
                computation_points: The computation points as tuple or list of points
                parallel:           If true, the computation is done in parallel (default: true)

            Returns:
                Either a tuple of potential, acceleration and second derivatives at the computation points or
                if multiple computation points are given, the result is a list of tuples
        )mydelimiter", py::arg("polyhedral_source"), py::arg("density"), py::arg("computation_points"),
          py::arg("parallel") = true);

    pybind11::class_<GravityEvaluable>(m, "GravityEvaluable", R"mydelimiter(
                A class to evaluate the polyhedral gravity model for a given constant density polyhedron at a given computation point.
                It provides a __call__ method to evaluate the polyhedral gravity model for computation points while
                also caching the polyhedron data over the lifetime of the object.
            )mydelimiter")
            .def(pybind11::init<const PolyhedralSource &, double>(),
                 R"mydelimiter(
                    Creates a new GravityEvaluable for a given constant density polyhedron.
                    It provides a __call__ method to evaluate the polyhedral gravity model for computation points while
                    also caching the polyhedron data over the lifetime of the object.

                    Args:
                        polyhedral_source:  The vertices & faces of the polyhedron as tuple or
                                            the filenames of the files containing the vertices & faces as list of strings
                        density:            The constant density of the polyhedron in [kg/m^3]
                    )mydelimiter",
                 pybind11::arg("polyhedral_source"),
                 pybind11::arg("density"))
            .def("__repr__", &GravityEvaluable::toString)
            .def("__call__", &GravityEvaluable::operator(),
                 R"mydelimiter(
                    Evaluates the polyhedral gravity model for a given constant density polyhedron at a given computation point.

                    Args:
                        computation_points: The computation points as tuple or list of points
                        parallel:           If true, the computation is done in parallel (default: true)

                    Returns:
                        Either a tuple of potential, acceleration and second derivatives at the computation points or
                        if multiple computation points are given, the result is a list of tuples
                )mydelimiter", py::arg("computation_points"), py::arg("parallel") = true);

    py::module_ utility = m.def_submodule("utility",
                                          "This submodule contains useful utility functions like parsing meshes "
                                          "or checking if the polyhedron's mesh plane unit normals point outwards "
                                          "like it is required by the polyhedral-gravity model.");

    utility.def("read", [](const std::vector<std::string> &filenames) {
        TetgenAdapter tetgen{filenames};
        auto polyhedron = tetgen.getPolyhedron();
        return std::make_tuple(polyhedron.getVertices(), polyhedron.getFaces());
    }, R"mydelimiter(
            Reads a polyhedron from a mesh file. The vertices and faces are read from input
            files (either .node/.face, mesh, .ply, .off, .stl). File-Order matters in case of the first option!

            Args:
                input_files (List[str]): list of input files.
            Returns:
                tuple of vertices (N, 3) (floats) and faces (N, 3) (ints).
          )mydelimiter", py::arg("input_files"));

    utility.def("check_mesh", [](const std::vector<Array3> &vertices, const std::vector<IndexArray3> &faces) {
        Polyhedron poly{vertices, faces};
        return MeshChecking::checkTrianglesNotDegenerated(poly) && MeshChecking::checkNormalsOutwardPointing(poly);
    }, R"mydelimiter(
                Checks if no triangles of the polyhedral mesh are degenerated by checking that their surface area
                is greater zero.
                Further, Checks if all the polyhedron's plane unit normals are pointing outwards.

                Args:
                    vertices (2-D array-like): (N, 3) array of vertex coordinates (floats).
                    faces (2-D array-like): (N, 3) array of faces, vertex-indices (ints).
                Returns:
                    True if no triangle is degenerate and the polyhedron's plane unit normals are all pointing outwards.
          )mydelimiter", py::arg("vertices"), py::arg("faces"));


    utility.def("check_mesh", [](const std::vector<std::string> &filenames) {
        TetgenAdapter tetgen{filenames};
        Polyhedron poly = tetgen.getPolyhedron();
        return MeshChecking::checkTrianglesNotDegenerated(poly) && MeshChecking::checkNormalsOutwardPointing(poly);
    }, R"mydelimiter(
                Checks if no triangles of the polyhedral mesh are degenerated by checking that their surface area
                is greater zero.
                Further, Checks if all the polyhedron's plane unit normals are pointing outwards.
                Reads a polyhedron from a mesh file. The vertices and faces are read from input
                files (either .node/.face, mesh, .ply, .off, .stl). File-Order matters in case of the first option!

                Args:
                    input_files (List[str]): list of input files.
                Returns:
                    True if no triangle is degenerate and the polyhedron's plane unit normals are all pointing outwards.
          )mydelimiter", py::arg("input_files"));
}
