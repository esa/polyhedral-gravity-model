#include "GravityModel.h"

namespace polyhedralGravity::GravityModel {

    GravityModelResult evaluate(const Polyhedron &polyhedron, const Array3 &computationPoint,
                                const ComputeBackend backend, const ComputePrecision precision) {
        const GravityEvaluable evaluable{polyhedron, precision};
        return std::get<GravityModelResult>(evaluable(computationPoint, backend));
    }

    std::vector<GravityModelResult> evaluate(const Polyhedron &polyhedron, const std::vector<Array3> &computationPoints,
                                             const ComputeBackend backend, const ComputePrecision precision) {
        const GravityEvaluable evaluable{polyhedron, precision};
        return std::get<std::vector<GravityModelResult>>(evaluable(computationPoints, backend));
    }

}
