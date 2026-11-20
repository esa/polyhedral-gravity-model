#pragma once

#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/KDTree/KdDefinitions.h"
#include "polyhedralGravity/model/KDTree/SplitParam.h"
#include "polyhedralGravity/model/KDTree/TreeNode.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/util/UtilityFloatArithmetic.h"
#include "thrust/detail/execution_policy.h"
#include "thrust/execution_policy.h"
#include "thrust/system/detail/sequential/for_each.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <set>
#include <string>
#include <sstream>
#include <variant>
#include <vector>

namespace polyhedralGravity {
struct SplitParam;

    /**
     * A TreeNode contained in a KDTree that doesn't split the spatial hierarchy any further. Intersection tests are directly performed on the contained triangles here.
     */
    class LeafNode final : public TreeNode {
    public:
        /**
         * Takes parameters from the parent node and stores them for later intersection tests.
         * @param splitParam Parameters produced during the split that resulted in the creation of this node.
         * @param nodeId Unique Id given by the TreeNodeFactory.
         */
        explicit LeafNode(const SplitParam &splitParam, size_t nodeId);

        /**
        * Used to calculated intersections of a ray and the polyhedron's faces contained in this node.
        * @param origin The point where the ray originates from.
        * @param ray Specifies the ray direction.
        * @param intersections The set intersection points are added to.
        */
        void getFaceIntersections(const Array3 &origin, const Array3 &ray, std::set<Array3> &intersections);

        [[nodiscard]] std::string toString() const override;

        friend std::ostream &operator<<(std::ostream &os, const LeafNode &node);

    private:
        /**
         * Möller-Trumbore Algorithm for Ray-Face intersection.
         * @param rayOrigin The point where the ray originates from.
         * @param rayVector Specifies the ray direction.
         * @param triangleVertexIndex the face to test against, described by the vertices that comprise it (passed by index reference).
         * @return
         */
        [[nodiscard]] std::optional<Array3> rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector,
                                                                  const IndexArray3 &triangleVertexIndex) const;

        /**
         * Möller-Trumbore Algorithm for Ray-Face intersection.
         * @param rayOrigin The point where the ray originates from.
         * @param rayVector Specifies the ray direction.
         * @param triangleVertices the face to test against, described by the vertices that comprise it (passed by value).
         * @return
         */
        static std::optional<Array3> rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector,
                                                           const Array3Triplet &triangleVertices);

        /**
         * Flags set when _splitParam boundFaces are converted from PlaneEvents to faces
         */
        std::once_flag convertedToFace;
    };
} // namespace polyhedralGravity
