#include "GravityModel.h"

namespace polyhedralGravity {

    GravityModelResult
    GravityModel::evaluate(const Polyhedron &polyhedron, double density, const Array3 &computationPoint) {
        GravityEvaluable evaluable{polyhedron, density};
        return evaluable.evaluate(computationPoint);
    }

    std::vector<GravityModelResult>
    GravityModel::evaluate(const Polyhedron &polyhedron, double density, const std::vector<Array3> &computationPoints) {
        std::vector<GravityModelResult> result{computationPoints.size()};
        GravityEvaluable evaluable{polyhedron, density};
        thrust::transform(computationPoints.begin(), computationPoints.end(), result.begin(),
                          [&evaluable](const Array3 &computationPoint) {
                              return evaluable.evaluate(computationPoint);
                          });
        return result;
    }


}
