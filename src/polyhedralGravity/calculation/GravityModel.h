#pragma once

#include <utility>
#include <array>
#include <vector>
#include <algorithm>
#include <variant>

#include "polyhedralGravity/calculation/GravityEvaluable.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/Polyhedron.h"

/**
 * Namespace containing the methods used to evaluate the polyhedrale Gravity Model
 * @note Naming scheme corresponds to the following:
 * evaluate()           --> main Method for evaluating the gravity model
 * *()                  --> Methods calculating one property for the evaluation
 */
namespace polyhedralGravity::GravityModel {

    /**
     * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
     * point P.
     * @param polyhedron - the polyhedron consisting of vertices and triangular faces
     * @param density - the constant density in [kg/m^3]
     * @param computationPoint - the computation Point P
     * @param parallel - whether to evaluate in parallel or serial
     * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration
     * at computation Point P
     */
    GravityModelResult evaluate(const PolyhedronOrSource &polyhedron, double density, const Array3 &computationPoint, bool parallel = true);

    /**
     * Evaluates the polyhedral gravity model for a given constant density polyhedron at multiple computation
     * points.
     * @param polyhedron - the polyhedron consisting of vertices and triangular faces
     * @param density - the constant density in [kg/m^3]
     * @param computationPoints - vector of computation points
     * @param parallel - whether to evaluate in parallel or serial
     * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration
     * foreach computation Point P
     */
    std::vector<GravityModelResult>
    evaluate(const PolyhedronOrSource &polyhedron, double density, const std::vector<Array3> &computationPoints, bool parallel = true);

}
