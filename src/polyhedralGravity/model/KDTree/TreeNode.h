#pragma once

#include <cstddef>
#include <memory>

#include "polyhedralGravity/model/KDTree/KdDefinitions.h"
#include "polyhedralGravity/model/KDTree/SplitParam.h"

namespace polyhedralGravity {

    /**
    * Abstract super class for nodes in the {@link KDTree}.
    */
    class TreeNode {
    public:
        virtual ~TreeNode() = default;

        /**
        * The current node's id. Follows the convention that the left child gets the id 2 * <current_id> + 1 and
        * the right child 2 * <currrent_id> + 2. Starts with 0 at the root node.
         */
        const size_t nodeId;

        [[nodiscard]] virtual std::string toString() const = 0;

        friend std::ostream& operator<<(std::ostream& os, const TreeNode& node);

    protected:
        /**
        * Protected constructor intended only for child classes. Please use {@link TreeNodeFactory} instead.
        */
        explicit TreeNode(const SplitParam &splitParam, size_t nodeId);
        /**
        * Stores parameters required for building child nodes lazily. Gets freed if the Node is an inner node and after both children are built.
        */
        std::unique_ptr<SplitParam> _splitParam;
    };

}// namespace polyhedralGravity
