#include "polyhedralGravity/model/KDTree/SplitNode.h"

namespace polyhedralGravity {
    SplitNode::SplitNode(const SplitParam &splitParam, const Plane &plane,
                         std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2> > &triangleIndexLists,
                         const size_t nodeId)
        : TreeNode(splitParam, nodeId), _plane{plane}, _boundingBox{splitParam.boundingBox},
          _triangleLists{std::move(triangleIndexLists)} {
    }

    std::shared_ptr<TreeNode> SplitNode::getChildNode(const size_t index) {
        //create a reference to store the built node in
        std::shared_ptr<TreeNode> &node = index == 0 ? _lesser : _greater;
        //node is not yet built
        std::call_once(childNodeCreated[index], [this, &node, &index] {
            //copy parent param and modify to fit new node
            SplitParam childParam{*_splitParam};
            //get the bounding box after splitting;
            auto [lesserBox, greaterBox] = this->_boundingBox.splitBox(this->_plane);
            childParam.boundingBox = index == 0 ? lesserBox : greaterBox;
            //get the triangles of the box
            std::visit([&childParam, index](auto &typeLists) -> void {
                childParam.boundFaces = *std::move(typeLists[index]);
            }, _triangleLists);
            childParam.splitDirection = static_cast<Direction>(
                (static_cast<int>(_splitParam->splitDirection) + 1) % DIMENSIONS);
            //increase the recursion depth of the direct child by 1
            node = TreeNodeFactory::createTreeNode(childParam, 2 * nodeId + 1 + index);
            if (_lesser != nullptr && _greater != nullptr) {
                _splitParam.reset();
            }
        });
        return node;
    }

    std::vector<std::shared_ptr<TreeNode> > SplitNode::getChildrenForIntersection(
        const Array3 &origin, const Array3 &ray, const Array3 &inverseRay) {
        using namespace polyhedralGravity::util;
        std::vector<std::shared_ptr<TreeNode> > delegates{};
        //a SplitNode has max two children, so no more space needed.
        delegates.reserve(2);
        //calculate entry and exit points of the ray hitting the bounding box
        auto [t_enter, t_exit] = _boundingBox.rayBoxIntersection(origin, inverseRay);
        // bounding box was not hit because the ray passed the box or is moving into the opposite direction of it,
        if (t_exit < t_enter || t_exit < 0) {
            //empty
            return delegates;
        }
        //calculate point where plane was hit
        const double t_split{_plane.rayPlaneIntersection(origin, inverseRay)};
        //the split plane is hit inside of the bounding box -> both child boxes need to be checked
        const bool isParallel = std::isinf(t_split);
        bool planeIsHitInsideBox = 0 <= t_split && t_enter <= t_split && t_split <= t_exit;
        if (!isParallel && planeIsHitInsideBox) {
            delegates.push_back(getChildNode(0));
            delegates.push_back(getChildNode(1));
            return delegates;
        }
        // the split plane is behind the ray origin
        if (t_split < 0) {
            //check in which point the origin lies in order to continue intersection in that box
            delegates.push_back(origin[static_cast<int>(_plane.orientation)] < _plane.axisCoordinate
                                    ? getChildNode(0)
                                    : getChildNode(1));
            return delegates;
        }
        //intersection point of the ray and the bounding box
        const double intersectionCoord{
            ray[static_cast<int>(_plane.orientation)] * t_enter + origin[static_cast<int>(_plane.orientation)]
        };
        // the entry point of the ray to the bounding box is nearer to the origin than the split plane -> ray hits lesser box
        if (intersectionCoord < _plane.axisCoordinate) {
            delegates.push_back(getChildNode(0));
        }
        // only the greater box is hit by the ray
        else {
            delegates.push_back(getChildNode(1));
        }
        return delegates;
    }

    std::string SplitNode::toString() const {
        std::stringstream sstream{};
        sstream << "SplitNode ID:  " << this->nodeId << ", Depth: " << recursionDepth(this->nodeId) << ", Plane: " <<
                this->_plane << std::endl;
        sstream << "Children; Lesser: " << (this->_lesser != nullptr ? std::to_string(this->_lesser->nodeId) : "None")
                << "; Greater: " << (this->_greater != nullptr ? std::to_string(this->_greater->nodeId) : "None") <<
                std::endl;
        if (this->_lesser != nullptr) {
            sstream << *(this->_lesser);
        }
        if (this->_greater != nullptr) {
            sstream << *(this->_greater);
        }
        return sstream.str();
    }

    std::ostream &operator<<(std::ostream &os, const SplitNode &node) {
        std::cout << node.toString();

        return os;
    }
} // namespace polyhedralGravity
