#pragma once

#include "polyhedralGravity/model/KDTree/SplitParam.h"
#include "polyhedralGravity/model/KDTree/plane_selection/PlaneSelectionAlgorithm.h"


namespace polyhedralGravity {

    class NoTreePlane final : public PlaneSelectionAlgorithm {
        /**
        * Builds a plane with infinite cost. That way splitting is considered unpractical and no tree is built -> intersection is performed directly on the trinagle faces.
        * @param splitParam specifies the polyhedron section to be split @link SplitParam.
        * @return Tuple of the default plane to split the specified bounding box, infinite cost as double and a list of empty triangle sets. Refer to {@link TriangleIndexVectors<2>} for more information.
        */
        std::tuple<Plane, double, std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2>>> findPlane(const SplitParam &splitParam) override {
            return {Plane{}, std::numeric_limits<double>::infinity(), TriangleIndexVectors<2>{std::make_unique<TriangleIndexVector>(), std::make_unique<TriangleIndexVector>()}};
        }
    };
}// namespace polyhedralGravity