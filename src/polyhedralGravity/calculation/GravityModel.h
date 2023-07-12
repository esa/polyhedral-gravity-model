#pragma once

#include <utility>
#include <array>
#include <vector>
#include <algorithm>
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/calculation/GravityEvaluable.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "spdlog/spdlog.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"
#include "thrust/transform_reduce.h"
#include "thrust/execution_policy.h"
#include "polyhedralGravity/util/UtilityThrust.h"
#include "polyhedralGravity/output/Logging.h"
#include "xsimd/xsimd.hpp"

namespace polyhedralGravity {

    /**
     * Namespace containing the methods used to evaluate the polyhedrale Gravity Model
     * @note Naming scheme corresponds to the following:
     * evaluate()           --> main Method for evaluating the gravity model
     * *()                  --> Methods calculating one property for the evaluation
     */
    namespace GravityModel {

        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
         * point P.
         * @param polyhedron - the polyhedron consisting of vertices and triangular faces
         * @param density - the constant density in [kg/m^3]
         * @param computationPoint - the computation Point P
         * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration
         * at computation Point P
         */
        GravityModelResult evaluate(
                const Polyhedron &polyhedron,
                double density,
                const Array3 &computationPoint);

        /**
         * Evaluates the polyhedral gravity model for a given constant density polyhedron at multiple computation
         * points.
         * @param polyhedron - the polyhedron consisting of vertices and triangular faces
         * @param density - the constant density in [kg/m^3]
         * @param computationPoints - vector of computation points
         * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration
         * foreach computation Point P
         */
        std::vector<GravityModelResult> evaluate(
                const Polyhedron &polyhedron,
                double density,
                const std::vector<Array3> &computationPoints);

        /**
         * An iterator transforming the polyhedron's coordinates on demand by a given offset.
         * This function returns a pair of transform iterators (first --> begin(), second --> end()).
         * @param polyhedron - reference to the polyhedron
         * @param offset - the offset to apply
         * @return pair of transform iterators
         */
        inline auto transformPolyhedron(const Polyhedron &polyhedron, const Array3 &offset = {0.0, 0.0, 0.0}) {
            //The offset is captured by value to ensure its lifetime!
            const auto lambdaOffsetApplication = [&polyhedron, offset](
                    const std::array<size_t, 3> &face) -> Array3Triplet {
                using namespace util;
                return {polyhedron.getVertex(face[0]) - offset,
                        polyhedron.getVertex(face[1]) - offset,
                        polyhedron.getVertex(face[2]) - offset};
            };
            auto first = thrust::make_transform_iterator(polyhedron.getFaces().begin(), lambdaOffsetApplication);
            auto last = thrust::make_transform_iterator(polyhedron.getFaces().end(), lambdaOffsetApplication);
            return std::make_pair(first, last);
        }

    }

}