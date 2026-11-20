#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <sstream>
#include <utility>
#include <variant>
#include <vector>

#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/KDTree/KdDefinitions.h"
#include "polyhedralGravity/model/KDTree/SplitParam.h"
#include "polyhedralGravity/model/KDTree/TreeNode.h"
#include "polyhedralGravity/model/KDTree/TreeNodeFactory.h"
#include "polyhedralGravity/util/UtilityContainer.h"

namespace polyhedralGravity {
struct SplitParam;

    /**
     * A TreeNode contained in a KDTree that splits the spatial hierarchy into two new sub boxes. Intersection tests are delegated to the child nodes.
     */
    class SplitNode final : public TreeNode {
        /**
         * friend declaration for testing purposes.
         */
        friend class KDTreeTest_AlgorithmRegressionTest_Test;

        /**
        * The SplitNode that contains the bounding box closer to the origin with respect to the split plane
        */
        std::shared_ptr<TreeNode> _lesser;
        /**
        * The SplitNode that contains the bounding box further away from the origin with respect to the split plane
        */
        std::shared_ptr<TreeNode> _greater;
        /**
         * Flags set when child node is created. Index 0 for lesser and 1 for greater.
         */
        std::array<std::once_flag, 2> childNodeCreated;
        /**
        * The plane splitting the two TreeNodes contained in this SplitNode
        */
        Plane _plane;
        /**
         * the bounding box for all triangles contained in this node.
         */
        Box _boundingBox;
        /**
         * Contains the triangle lists for the lesser and greater bounding boxes. {@link TriangleIndexVectors}
        */
        std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2>> _triangleLists;

    public:
        /**
         * Takes parameters from the parent node and stores them for lazy child node creation.
         * @param splitParam Parameters produced during the split that resulted in the creation of this node.
         * @param plane The plane that splits this node's bounding box into two sub boxes. The child nodes are created based on these boxes.
         * @param triangleIndexLists Index sets of the triangles contained in the lesser and greater child nodes. {@link TriangleIndexVector}
         * @param nodeId Unique Id given by the TreeNodeFactory.
         */
        SplitNode(const SplitParam &splitParam, const Plane &plane, std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2>> &triangleIndexLists, size_t nodeId);
        /**
         * Computes the child node decided by the given index (0 for lesser, 1 for greater) if not present already and returns it to the caller.
         * @param index Specifies which node to build. 0 or LESSER for _lesser, 1 or GREATER for _greater.
         * @return the built TreeNode.
        */
        std::shared_ptr<TreeNode> getChildNode(size_t index);
        /**
         * Gets the children of this node whose bounding boxes are hit by the ray.
         * @param origin The point where the ray originates from.
         * @param ray Specifies the ray direction.
         * @param inverseRay The inverse of the ray (1/ray), used to speed up calculations where it is divided by ray. Instead, we multiply with inverseRay.
         * @return the child nodes that intersect with the ray
         */
        [[nodiscard]] std::vector<std::shared_ptr<TreeNode>> getChildrenForIntersection(const Array3 &origin, const Array3 &ray, const Array3 &inverseRay);

        [[nodiscard]] std::string toString() const override;

        friend std::ostream &operator<<(std::ostream &os, const SplitNode &node);
    };

}// namespace polyhedralGravity