#include "GravityModel.h"

namespace polyhedralGravity::GravityModel {

    GravityModelResult evaluate(const Polyhedron &polyhedron, const Array3 &computationPoint, bool parallel,
                                const ComputeBackend backend, const ComputePrecision precision) {
        const GravityEvaluable evaluable{polyhedron, backend, precision};
        return std::get<GravityModelResult>(evaluable(computationPoint, parallel));
    }

    std::vector<GravityModelResult> evaluate(const Polyhedron &polyhedron, const std::vector<Array3> &computationPoints,
                                             bool parallel, const ComputeBackend backend,
                                             const ComputePrecision precision) {
        const GravityEvaluable evaluable{polyhedron, backend, precision};
        return std::get<std::vector<GravityModelResult>>(evaluable(computationPoints, parallel));
    }

}
