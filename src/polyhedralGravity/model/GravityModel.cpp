#include "GravityModel.h"

namespace polyhedralGravity::GravityModel {

    GravityModelResult evaluate(
            const PolyhedralSource &polyhedron,
            double density, const Array3 &computationPoint, bool parallel, const NormalOrientation &orientation,
            const std::optional<bool> check) {
        return std::get<GravityModelResult>(
                std::visit([&](const auto &polyhedron) {
                    GravityEvaluable evaluable{polyhedron, density, orientation, check};
                    return evaluable(computationPoint, parallel);
                }, polyhedron));
    }

    std::vector<GravityModelResult>
    evaluate(const PolyhedralSource &polyhedron, double density, const std::vector<Array3> &computationPoints,
             bool parallel, const NormalOrientation &orientation,
             const std::optional<bool> check) {
        return std::get<std::vector<GravityModelResult>>(std::visit([&](const auto &polyhedron) {
            GravityEvaluable evaluable{polyhedron, density, orientation, check};
            return evaluable(computationPoints, parallel);
        }, polyhedron));
    }

}
