#include <tuple>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/input/TetgenAdapter.h"


namespace py = pybind11;

/**
 * This method converts a named tuples GravityModelResult to a tuple
 * @param result a named tuple of type GravityModelResult
 * @return tuple of potential, acceleration and gradiometric tensor
 */
std::tuple<double, std::array<double, 3>, std::array<double, 6>> convertToTuple(
        const polyhedralGravity::GravityModelResult &result) {
    return std::make_tuple(result.gravitationalPotential,
                           result.acceleration,
                           result.gradiometricTensor);
}

/**
 * This method converts a vector of named tuples GravityModelResult to a vector of tuples
 * @param resultVector a vector of named tuples of type GravityModelResult
 * @return vector of tuples of potential, acceleration and gradiometric tensor
 */
std::vector<std::tuple<double, std::array<double, 3>, std::array<double, 6>>> convertToTupleVector(
        const std::vector<polyhedralGravity::GravityModelResult> &resultVector) {
    std::vector<std::tuple<double, std::array<double, 3>, std::array<double, 6>>> resVectorTuple{resultVector.size()};
    std::transform(resultVector.cbegin(), resultVector.cend(), resVectorTuple.begin(), &convertToTuple);
    return resVectorTuple;
}

PYBIND11_MODULE(polyhedral_gravity, m) {
    using namespace polyhedralGravity;
    m.doc() = "Computes the full gravity tensor for a given constant density polyhedron which consists of some "
              "vertices and triangular faces at a given computation point P";

    /*
     * Methods for vertices and faces from the script context
     */

    m.def("evaluate",
          [](const std::vector<std::array<double, 3>> &vertices, const std::vector<std::array<size_t, 3>> &faces,
             double density, const std::array<double, 3> &computationPoint) {
              return convertToTuple(GravityModel::evaluate({vertices, faces}, density, computationPoint));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices "
          "and triangular faces at a given computation point P",
          py::arg("vertices"), py::arg("faces"), py::arg("density"), py::arg("computation_point"));

    m.def("evaluate",
          [](const std::vector<std::array<double, 3>> &vertices, const std::vector<std::array<size_t, 3>> &faces,
             double density, const std::vector<std::array<double, 3>> &computationPoints) {
              return convertToTupleVector(GravityModel::evaluate({vertices, faces}, density, computationPoints));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices "
          "and triangular faces at multiple given computation points",
          py::arg("vertices"), py::arg("faces"), py::arg("density"), py::arg("computation_points"));

    /*
     * Methods for vertices and faces from files via TetgenAdapter
     */

    m.def("evaluate",
          [](const std::vector<std::string> &filenames,
             double density, const std::array<double, 3> &computationPoint) {
              TetgenAdapter tetgen{filenames};
              return convertToTuple(GravityModel::evaluate(tetgen.getPolyhedron(), density, computationPoint));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices"
          "and triangular faces at a given computation point P. The vertices and faces are read from input "
          "files (either .node/.face, mesh, .ply, .off, .stl). File-Order matters in case of the first option!",
          py::arg("input_files"), py::arg("density"), py::arg("computation_point"));

    m.def("evaluate",
          [](const std::vector<std::string> &filenames,
             double density, const std::vector<std::array<double, 3>> &computationPoints) {
              TetgenAdapter tetgen{filenames};
              return convertToTupleVector(GravityModel::evaluate(tetgen.getPolyhedron(), density, computationPoints));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices "
          "and triangular faces at multiple given computation points. The vertices and faces are read from "
          "input files (either .node/.face, mesh, .ply, .off, .stl). File-Order matters in case of the first option!",
          py::arg("input_files"), py::arg("density"), py::arg("computation_points"));
}
