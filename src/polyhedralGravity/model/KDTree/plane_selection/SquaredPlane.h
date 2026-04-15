#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/KDTree/KdDefinitions.h"
#include "polyhedralGravity/model/KDTree/SplitParam.h"
#include "polyhedralGravity/model/KDTree/plane_selection/PlaneSelectionAlgorithm.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "thrust/detail/execution_policy.h"
#include "thrust/execution_policy.h"
#include "thrust/iterator/iterator_facade.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/system/detail/sequential/for_each.h"

namespace polyhedralGravity {
struct SplitParam;

    /**
* O(N^2) implementation to finding optimal split planes.
*/
    class SquaredPlane final : public PlaneSelectionAlgorithm {
        /**
        * Finds the optimal split plane to split a provided rectangle section optimally.
        * @param splitParam specifies the polyhedron section to be split @link SplitParam.
        * @return Tuple of the optimal plane to split the specified bounding box, its cost as double and a list of triangle sets with respective positions to the found plane. Refer to {@link TriangleIndexVectors<2>} for more information.
        */
        std::tuple<Plane, double, std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2> > > findPlane(
            const SplitParam &splitParam) override;

        /**
        * Splits a section of a polyhedron into two bounding boxes and calculates the triangle face sets contained in the new bounding boxes.
        * @param splitParam specifies the polyhedron section to be split.
        * @param split the plane by which to split the polyhedron section.
        * @return Three triangle lists contained in an array. Those being the set of triangles with non-zero area in the bounding box closer to the origin with respect to the split plane,
        * the set of triangles with non-zero area in the bounding box further away from the origin with respect to the split plane.
        * The set of triangles that lies on the plane.
        */
        static TriangleIndexVectors<3> containedTriangles(const SplitParam &splitParam, const Plane &split);
    };
} // namespace polyhedralGravity
