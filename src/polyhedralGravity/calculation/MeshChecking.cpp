#include "MeshChecking.h"

namespace polyhedralGravity::MeshChecking {

    bool checkPlaneUnitNormalOrientation(const Polyhedron &polyhedron, const NormalOrientation &normalOrientation) {
        auto it = transformPolyhedron(polyhedron);
        // All normals have to point outwards (intersect the polyhedron even times)
        return thrust::transform_reduce(
                thrust::device,
                it.first, it.second, [&polyhedron, &normalOrientation](const Array3Triplet &face) {
                    // If the ray intersects the polyhedron an even number of times then the normal points outwards
                    // If the ray intersects the polyhedron odd number if times the normal points inwards
                    const size_t intersects = detail::countRayPolyhedronIntersections(face, polyhedron);
                    return (intersects % 2 == 0) == (normalOrientation == NormalOrientation::OUTWARDS);
                }, true, thrust::logical_and<bool>());
    }

    std::pair<NormalOrientation, std::set<size_t>> getPlaneUnitNormalOrientation(const Polyhedron &polyhedron) {
        // 1. Step: Find all indices of normals which vioate the constraint outwards pointing
        auto indexIterator = thrust::make_counting_iterator<size_t>(0);
        size_t numberOfFaces = polyhedron.countFaces();
        std::set<size_t> violatingOutwards{};
        thrust::copy_if(
                thrust::device,
                indexIterator,
                indexIterator + numberOfFaces,
                std::inserter(violatingOutwards, violatingOutwards.end()),
                [&polyhedron](const size_t index) {
                    // If the ray intersects the polyhedron odd number if times the normal points inwards
                    // Hence, violating the OUTWARDS constraint
                    const size_t intersects = detail::countRayPolyhedronIntersections(polyhedron.getFace(index), polyhedron);
                    return intersects % 2 != 0;
                });
        // 2. Step: If we have more than n/2 violations, the majority is inwards pointing --> Build complemntary set!
        if (violatingOutwards.size() > numberOfFaces / 2) {
            std::set<size_t> fullSet{};
            std::copy_n(indexIterator, numberOfFaces, std::inserter(fullSet, fullSet.end()));
            std::set<size_t> violatingInwards{};
            std::set_difference(
                    fullSet.begin(), fullSet.end(),
                    violatingOutwards.begin(), violatingOutwards.end(),
                    std::inserter(violatingInwards, violatingInwards.end())
                    );
            // 3a. Step: Return the inwards pointing as major orientation and
            // the violating faces, i.e. which have outwards pointing normals
            return std::make_pair(NormalOrientation::INWWARDS, violatingInwards);
        }
        // 3b. Step: Return the outwards pointing as major orientation
        // and the violating faces, i.e. which have inwards pointing normals
        return std::make_pair(NormalOrientation::OUTWARDS, violatingOutwards);
    }

    bool checkTrianglesNotDegenerated(const Polyhedron &polyhedron) {
        auto it = transformPolyhedron(polyhedron);
        // All triangles surface area needs to be greater than zero
        return thrust::transform_reduce(
                thrust::device,
                it.first, it.second, [](const Array3Triplet &face) {
                    return util::surfaceArea(face) > 0.0;
                }, true, thrust::logical_and<bool>());
    }

    size_t
    detail::countRayPolyhedronIntersections(const Array3Triplet &face, const Polyhedron &polyhedron) {
        using namespace util;
        // The centroid of the triangular face
        const Array3 centroid = (face[0] + face[1] + face[2]) / 3.0;

        // The normal of the plane calculated with two segments of the triangle
        // The normal is the rayVector starting at the rayOrigin
        const Array3 segmentVector1 = face[1] - face[0];
        const Array3 segmentVector2 = face[2] - face[1];
        const Array3 rayVector = util::normal(segmentVector1, segmentVector2);

        // The origin of the array has a slight offset in direction of the normal
        const Array3 rayOrigin = centroid + (rayVector * EPSILON);

        // Count every triangular face which is intersected by the ray
        auto it = transformPolyhedron(polyhedron);
        std::set<Array3> intersections{};
        std::for_each(it.first, it.second, [&rayOrigin, &rayVector, &intersections](const Array3Triplet &otherFace) {
            const std::unique_ptr<Array3> intersection = rayIntersectsTriangle(rayOrigin, rayVector, otherFace);
            if (intersection != nullptr) {
                intersections.insert(*intersection);
            }
        });
        return intersections.size();
    }

    std::unique_ptr<Array3>
    detail::rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector, const Array3Triplet &triangle) {
        // Adapted Möller–Trumbore intersection algorithm
        // see https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
        using namespace util;
        const Array3 edge1 = triangle[1] - triangle[0];
        const Array3 edge2 = triangle[2] - triangle[0];
        const Array3 h = cross(rayVector, edge2);
        const double a = dot(edge1, h);
        if (a > -EPSILON && a < EPSILON) {
            return nullptr;
        }

        const double f = 1.0 / a;
        const Array3 s = rayOrigin - triangle[0];
        const double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0) {
            return nullptr;
        }

        const Array3 q = cross(s, edge1);
        const double v = f * dot(rayVector, q);
        if (v < 0.0 || u + v > 1.0) {
            return nullptr;
        }

        const double t = f * dot(edge2, q);
        if (t > EPSILON) {
            return std::make_unique<Array3>(rayOrigin + rayVector * t);
        } else {
            return nullptr;
        }
    }
}