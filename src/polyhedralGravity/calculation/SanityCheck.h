#pragma once

#include <set>
#include <optional>
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/calculation/GravityModel.h"

namespace polyhedralGravity::SanityCheck {

    /**
     * Checks if the vertices are in such a order that the unit normals of each plane point outwards the polyhedron
     * @param polyhedron - the polyhedron consisting of vertices and triangular faces
     * @return true if all the unit normals are pointing outwards
     *
     * @note This has quadratic complexity O(n^2)! For bigger problem sizes, K-D trees/ Octrees could improve
     * the complexity to determine the intersections
     * (see https://stackoverflow.com/questions/45603469/how-to-calculate-the-normals-of-a-box)
     */
    bool checkNormalsOutwardPointing(const Polyhedron &polyhedron);


    namespace detail {

        /**
         * Calculates how often a vector starting at a specific origin intersects a polyhedron's mesh's triangles.
         * @param rayOrigin - the origin of the ray
         * @param rayVector - the vector describing the ray
         * @param polyhedron - the polyhedron consisting of vertices and triangular faces
         * @return true if the ray intersects the triangle
         */
        size_t rayIntersectsPolyhedron(const Array3 &rayOrigin, const Array3 &rayVector, const Polyhedron &polyhedron);

        /**
         * Calculates how often a vector starting at a specific origin intersects a triangular face.
         * Uses the Möller–Trumbore intersection algorithm.
         * @param rayOrigin - the origin of the ray
         * @param rayVector - the vector describing the ray
         * @param triangle - a triangular face
         * @return true if the ray intersects the triangle
         *
         * @related Adapted from https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
         */
        std::optional<Array3>
        rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector, const Array3Triplet &triangle);

    }

}