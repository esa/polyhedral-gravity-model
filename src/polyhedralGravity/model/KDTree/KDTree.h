#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <deque>
#include <iterator>
#include <memory>
#include <mutex>
#include <ostream>
#include <set>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <utility>
#include <vector>

#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/KDTree/KdDefinitions.h"
#include "polyhedralGravity/model/KDTree/LeafNode.h"
#include "polyhedralGravity/model/KDTree/SplitNode.h"
#include "polyhedralGravity/model/KDTree/SplitParam.h"
#include "polyhedralGravity/model/KDTree/TreeNode.h"
#include "polyhedralGravity/model/KDTree/TreeNodeFactory.h"
#include "polyhedralGravity/model/KDTree/plane_selection/PlaneSelectionAlgorithm.h"
#include "polyhedralGravity/model/KDTree/plane_selection/PlaneSelectionAlgorithmFactory.h"
#include "polyhedralGravity/util/UtilityContainer.h"

namespace polyhedralGravity {
    /**
     * A KDTree for a given polyhedron to speed up ray intersections with the polyhedron. It is thread safe.
     */
    class KDTree {
        /**
         * friend declaration for testing purposes.
         */
        friend class KDTreeTest_AlgorithmRegressionTest_Test;

        /**
        * The entry node of the KDTree. Only access using getter.
        */
        std::shared_ptr<TreeNode> _rootNode;

        /**
         * The polyhedron's vertices.
         */
        const std::vector<Array3> _vertices;
        /**
         * The polyhedron's faces: A face is a triplet of vertex indices.
         */
        const std::vector<IndexArray3> _faces;

        /**
         * Set when the root node has been created.
         */
        std::once_flag _rootNodeCreated;

        /**
        * Parameters for lazily building the root node {@link SplitParam}
        */
        std::unique_ptr<SplitParam> _splitParam;

    public:
        /**
        * Call to build a KDTree to speed up intersections of rays with a polyhedron's faces.
        * @param vertices The vertex coordinates of the polyhedron
        * @param faces The faces of the polyhedron with a face being a triplet of vertex indices
        * @param algorithm Specifies which algorithm to use for finding optimal split planes.
        * @return the lazily built KDTree.
        */
        KDTree(const std::vector<Array3> &vertices, const std::vector<IndexArray3> &faces,
               PlaneSelectionAlgorithm::Algorithm algorithm = PlaneSelectionAlgorithm::Algorithm::LOG);

        /**
        * Creates the root tree node if not initialized and returns it.
        * @return the root tree Node.
        */
        std::shared_ptr<TreeNode> getRootNode();

        /**
        * Used to calculate intersections of a ray and the polyhedron's faces contained in this node.
        * @param origin The point where the ray originates from.
        * @param ray Specifies the ray direction.
        * @param intersections The set found intersection points are added to.
        */
        void getFaceIntersections(const Array3 &origin, const Array3 &ray, std::set<Array3> &intersections);

        /**
         * Calculates the number of intersections of a ray with the polyhedron.
         * @param origin The origin point of the ray.
         * @param ray The ray direction vector.
         * @return the number of intersections.
         */
        size_t countIntersections(const Array3 &origin, const Array3 &ray);

        /**
         * Prebuilds the whole KDTree bypassing lazy loading entirely.
         */
        KDTree &prebuildTree();

        friend std::ostream &operator<<(std::ostream &os, const KDTree &kdTree);
    };
} // namespace polyhedralGravity
