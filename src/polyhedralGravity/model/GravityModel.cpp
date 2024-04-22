#include "GravityModel.h"

namespace polyhedralGravity::GravityModel {

    GravityModelResult evaluate(const Polyhedron &polyhedron, const Array3 &computationPoint, bool parallel) {
        GravityEvaluable evaluable{polyhedron};
        return std::get<GravityModelResult>(evaluable(computationPoint, parallel));
    }

    std::vector<GravityModelResult> evaluate(const Polyhedron &polyhedron, const std::vector<Array3> &computationPoints, bool parallel) {
        GravityEvaluable evaluable{polyhedron};
        return std::get<std::vector<GravityModelResult>>(evaluable(computationPoints, parallel));
    }

}
