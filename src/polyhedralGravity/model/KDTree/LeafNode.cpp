#include "polyhedralGravity/model/KDTree/LeafNode.h"

namespace polyhedralGravity {
    LeafNode::LeafNode(const SplitParam &splitParam, const size_t nodeId)
        : TreeNode(splitParam, nodeId) {
    }

    void LeafNode::getFaceIntersections(const Array3 &origin, const Array3 &ray,
                                        std::set<Array3> &intersections) {
        if (std::holds_alternative<PlaneEventVector>(_splitParam->boundFaces)) {
            std::call_once(convertedToFace, [this]() {
                _splitParam->boundFaces = convertEventsToFaces(std::get<PlaneEventVector>(_splitParam->boundFaces));
            });
        }
        std::mutex writeLock{};
        const TriangleIndexVector &boundTriangles{std::get<TriangleIndexVector>(_splitParam->boundFaces)};
        std::vector<Array3> results(boundTriangles.size());
        //traverses all contained faces and performs intersection tests with them -> store results in the buffer passed in the arguments
        thrust::for_each(thrust::device, boundTriangles.cbegin(), boundTriangles.cend(),
                         [this, &ray, &origin, &intersections, &writeLock](const size_t faceIndex) {
                             const std::optional<Array3> intersection = rayIntersectsTriangle(
                                 origin, ray, _splitParam->faces[faceIndex]);
                             if (intersection.has_value()) {
                                 std::unique_lock lock(writeLock);
                                 intersections.insert(intersection.value());
                             }
                         });
    }

    std::optional<Array3> LeafNode::rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector,
                                                          const IndexArray3 &triangleVertexIndex) const {
        Array3Triplet edgeVertices{};
        //transforms a face to its vertices
        std::transform(triangleVertexIndex.cbegin(), triangleVertexIndex.cend(), edgeVertices.begin(),
                       [this](const size_t vertexIndex) {
                           return _splitParam->vertices[vertexIndex];
                       });
        return rayIntersectsTriangle(rayOrigin, rayVector, edgeVertices);
    }

    std::optional<Array3> LeafNode::rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector,
                                                          const Array3Triplet &triangleVertices) {
        // Adapted Möller–Trumbore intersection algorithm
        // see https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
        using namespace polyhedralGravity;
        using namespace util;
        const Array3 edge1 = triangleVertices[1] - triangleVertices[0];
        const Array3 edge2 = triangleVertices[2] - triangleVertices[0];
        const Array3 h = cross(rayVector, edge2);
        const double a = dot(edge1, h);
        if (a > -EPSILON_ZERO_OFFSET && a < EPSILON_ZERO_OFFSET) {
            return std::nullopt;
        }

        const double f = 1.0 / a;
        const Array3 s = rayOrigin - triangleVertices[0];
        const double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0) {
            return std::nullopt;
        }

        const Array3 q = cross(s, edge1);
        const double v = f * dot(rayVector, q);
        if (v < 0.0 || u + v > 1.0) {
            return std::nullopt;
        }

        const double t = f * dot(edge2, q);
        if (t > EPSILON_ZERO_OFFSET) {
            return rayOrigin + rayVector * t;
        }
        return std::nullopt;
    }

    std::string LeafNode::toString() const {
        std::stringstream sstream{};
        sstream << "LeafNode ID: " << this->nodeId << ", Depth: " << recursionDepth(this->nodeId) << std::endl;
        return sstream.str();
    }

    std::ostream &operator<<(std::ostream &os, const LeafNode &node) {
        os << node.toString();
        return os;
    }
}
