#pragma once

#include <array>
#include <vector>

#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/Polyhedron.h"

/**
 * Namespace containing the methods used to evaluate the polyhedral Gravity Model
 * @note Naming scheme corresponds to the following:
 * evaluate()           --> main Method for evaluating the gravity model
 * *()                  --> Methods calculating one property for the evaluation
 */
namespace polyhedralGravity::GravityModel {

    /**
     * Evaluates the polyhedral gravity model for a given constant density polyhedron at computation
     * point P.
     *
     * The results' units depend on the polyhedron's input units.
     * For example, if the polyhedral mesh is in @f$[m]@f$ and the density in @f$[kg/m^3]@f$, then the potential is in @f$[m^2/s^2]@f$.
     * In case the polyhedron is unitless, the results are not multiplied with the Gravitational Constant @f$G@f$, but returned raw.
     *
     * @param polyhedron the polyhedron consisting of vertices and triangular faces
     * @param computationPoint the computation Point P
     * @param backend the compute backend the evaluation runs on (default: CPU_PARALLEL)
     * @param precision the floating point precision the evaluation computes in (default: FLOAT64)
     * @return the GravityModelResult containing the potential, the acceleration, and the change of acceleration
     * at computation Point P
     *
     * @throws std::runtime_error if GPU_PARALLEL is requested, but this build has no GPU backend
     *
     * @note Every call sets up its own {@link GravityEvaluable}, which means uploading the polyhedron to the
     * compute device again. Prefer a {@link GravityEvaluable} for repeated evaluations.
     */
    GravityModelResult evaluate(const Polyhedron &polyhedron, const Array3 &computationPoint,
                                ComputeBackend backend = ComputeBackend::CPU_PARALLEL,
                                ComputePrecision precision = ComputePrecision::FLOAT64);

    /**
     * Evaluates the polyhedral gravity model for a given constant density polyhedron at multiple computation
     * points.
     *
     * The results' units depend on the polyhedron's input units.
     * For example, if the polyhedral mesh is in @f$[m]@f$ and the density in @f$[kg/m^3]@f$, then the potential is in @f$[m^2/s^2]@f$.
     * In case the polyhedron is unitless, the results are not multiplied with the Gravitational Constant @f$G@f$, but returned raw.
     *
     * @param polyhedron the polyhedron consisting of vertices and triangular faces
     * @param computationPoints vector of computation points
     * @param backend the compute backend the evaluation runs on (default: CPU_PARALLEL)
     * @param precision the floating point precision the evaluation computes in (default: FLOAT64)
     * @return the GravityModelResult containing the potential, the acceleration, and the change of acceleration
     * foreach computation Point P
     *
     * @throws std::runtime_error if GPU_PARALLEL is requested, but this build has no GPU backend
     *
     * @note Every call sets up its own {@link GravityEvaluable}, which means uploading the polyhedron to the
     * compute device again. Prefer a {@link GravityEvaluable} for repeated evaluations.
     */
    std::vector<GravityModelResult>
    evaluate(const Polyhedron &polyhedron, const std::vector<Array3> &computationPoints,
             ComputeBackend backend = ComputeBackend::CPU_PARALLEL,
             ComputePrecision precision = ComputePrecision::FLOAT64);

}