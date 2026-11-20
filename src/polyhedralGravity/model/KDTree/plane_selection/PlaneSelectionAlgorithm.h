#pragma once

#include <cstddef>
#include <iterator>
#include <limits>
#include <tuple>
#include <utility>
#include <variant>

#include "polyhedralGravity/model/KDTree/KdDefinitions.h"

namespace polyhedralGravity {

    //forward declaration
    struct SplitParam;

    class PlaneSelectionAlgorithm {
    public:
        virtual ~PlaneSelectionAlgorithm() = default;
        /**
        * Finds the optimal split plane to split a provided rectangle section optimally.
        * @param splitParam specifies the polyhedron section to be split @link SplitParam.
        * @return Tuple of the optimal plane to split the specified bounding box, its cost as double and a list of triangle sets with respective positions to the found plane. Refer to {@link TriangleIndexVectors<2>} for more information.
        */
        virtual std::tuple<Plane, double, std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2>>> findPlane(const SplitParam &splitParam) = 0;

        enum class Algorithm {
            NOTREE,
            QUADRATIC,
            LOGSQUARED,
            LOG
        };


        /**
        * Constant that describes the cost of traversing the KDTree by one step.
        */
        constexpr static double traverseStepCost{1.0};

        /**
        * Constant that describes the cost of intersecting a ray and a single object.
        */
        constexpr static double triangleIntersectionCost{1.0};

    protected:
        /**
       * Evaluates the cost function should the specified bounding box and it's faces be divided by the specified plane. Used to evaluate possible split planes.
       * @param boundingBox the bounding box encompassing the scene to be split.
       * @param plane the candidate split plane to be evaluated.
       * @param trianglesMin the number of triangles overlapping with the min side of the bounding box.
       * @param trianglesMax the number of triangles overlapping with the max side of the bounding box.
       * @param trianglesPlanar the number of triangles lying in the plane.
       * @return A pair of: 1. the cost for performing intersection operations on the finalized tree later, should the KDTree be built using the specified split plane and the triangle sets resulting through division by the plane.
       * 2. true if the planar triangles should be added to the min side of the bounding box.
       */
        static std::pair<const double, bool> costForPlane(Box boundingBox, Plane plane, size_t trianglesMin, size_t trianglesMax, size_t trianglesPlanar);
    };

}// namespace polyhedralGravity