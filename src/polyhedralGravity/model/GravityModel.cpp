#include "GravityModel.h"

namespace polyhedralGravity::GravityModel {

    GravityModelResult evaluate(const Polyhedron &polyhedron, const Array3 &computationPoint,
                                const ComputeBackend backend, const ComputePrecision precision) {
        const GravityEvaluable evaluable{polyhedron, backend, precision};
        return std::get<GravityModelResult>(evaluable(computationPoint));
    }

    std::vector<GravityModelResult> evaluate(const Polyhedron &polyhedron, const std::vector<Array3> &computationPoints,
                                             const ComputeBackend backend, const ComputePrecision precision) {
        const GravityEvaluable evaluable{polyhedron, backend, precision};
        return std::get<std::vector<GravityModelResult>>(evaluable(computationPoints));
    }

}
