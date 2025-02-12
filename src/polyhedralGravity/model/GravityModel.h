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
     * @param parallel whether to evaluate in parallel or serial
     * @return the GravityModelResult containing the potential, the acceleration, and the change of acceleration
     * at computation Point P
     */
    GravityModelResult evaluate(const Polyhedron &polyhedron, const Array3 &computationPoint, bool parallel = true);

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
     * @param parallel whether to evaluate in parallel or serial
     * @return the GravityModelResult containing the potential, the acceleration, and the change of acceleration
     * foreach computation Point P
     */
    std::vector<GravityModelResult>
    evaluate(const Polyhedron &polyhedron, const std::vector<Array3> &computationPoints, bool parallel = true);

}