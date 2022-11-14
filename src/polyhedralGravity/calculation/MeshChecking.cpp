#include "MeshChecking.h"

namespace polyhedralGravity::MeshChecking {

    bool checkNormalsOutwardPointing(const Polyhedron &polyhedron) {
        using namespace util;
        auto it = GravityModel::transformPolyhedron(polyhedron);
        // All normals have to point outwards (intersect the polyhedron even times)
        return thrust::transform_reduce(
                thrust::device,
                it.first, it.second, [&polyhedron](const Array3Triplet &face) {
                    // The centroid of the triangular face
                    const Array3 centroid = (face[0] + face[1] + face[2]) / 3.0;

                    // The normal of the plane calculated with two segments of the triangle
                    const Array3 segmentVector1 = face[1] - face[0];
                    const Array3 segmentVector2 = face[2] - face[1];
                    const Array3 normal = util::normal(segmentVector1, segmentVector2);

                    // The origin of the array has a slight offset in direction of the normal
                    const Array3 rayOrigin = centroid + (normal * EPSILON);

                    // If the ray intersects the polyhedron even-times than the normal points outwards
                    const size_t intersects = detail::countRayPolyhedronIntersections(rayOrigin, normal, polyhedron);
                    return intersects % 2 == 0;
                }, true, thrust::logical_and<bool>());
    }

    size_t
    detail::countRayPolyhedronIntersections(const Array3 &rayOrigin, const Array3 &rayVector, const Polyhedron &polyhedron) {
        auto it = GravityModel::transformPolyhedron(polyhedron);
        std::set<Array3> intersections{};
        // Count every triangular face which is intersected by the ray
        std::for_each(it.first, it.second, [&rayOrigin, &rayVector, &intersections](const Array3Triplet &triangle) {
            const auto intersection = rayIntersectsTriangle(rayOrigin, rayVector, triangle);
            if (intersection.has_value()) {
                intersections.insert(intersection.value());
            }
        });
        return intersections.size();
    }

    std::optional<Array3>
    detail::rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector, const Array3Triplet &triangle) {
        // Adapted Möller–Trumbore intersection algorithm
        // see https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
        using namespace util;
        const Array3 edge1 = triangle[1] - triangle[0];
        const Array3 edge2 = triangle[2] - triangle[0];
        const Array3 h = cross(rayVector, edge2);
        const double a = dot(edge1, h);
        if (a > -EPSILON && a < EPSILON) {
            return std::nullopt;
        }

        const double f = 1.0 / a;
        const Array3 s = rayOrigin - triangle[0];
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
        if (t > EPSILON) {
            return rayOrigin + rayVector * t;
        } else {
            return std::nullopt;
        }
    }
}