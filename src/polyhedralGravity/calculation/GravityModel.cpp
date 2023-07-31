#include "GravityModel.h"

namespace polyhedralGravity::GravityModel {

    GravityModelResult evaluate(
            const PolyhedronOrSource &polyhedron,
            double density, const Array3 &computationPoint, bool parallel) {
        return std::get<GravityModelResult>(
                std::visit([&parallel, &density, &computationPoint](const auto &polyhedron) {
                    GravityEvaluable evaluable{polyhedron, density};
                    return evaluable(computationPoint, parallel);
                }, polyhedron));
    }

    std::vector<GravityModelResult>
    evaluate(const PolyhedronOrSource &polyhedron, double density, const std::vector<Array3> &computationPoints,
             bool parallel) {
        return std::get<std::vector<GravityModelResult>>(std::visit([&parallel, &density, &computationPoints](const auto &polyhedron) {
            GravityEvaluable evaluable{polyhedron, density};
            return evaluable(computationPoints, parallel);
        }, polyhedron));
    }

}
