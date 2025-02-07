#pragma once

#include <array>
#include <vector>
#include <string>
#include <tuple>

namespace polyhedralGravity {

    /**
     * Alias for an array of size 3 (double) for x, y, z coordinates.
     */
    using Array3 = std::array<double, 3>;

    /**
     * Alias for an array of size 3 (size_t) for the vertex indices in a triangular face.
     */
    using IndexArray3 = std::array<size_t, 3>;

    /**
     * Alias for an array of size 6 for xx, yy, zz, xy, xz, yz second derivatives.
     */
    using Array6 = std::array<double, 6>;

    /**
     * Alias for a triplet of arrays of size 3 for the segment of a triangular face
     */
    using Array3Triplet = std::array<Array3, 3>;

    /**
     * Contains in the order of the tuple:
     * The gravitational potential in [m^2/s^2] <--> [J/kg] at point P.
     * Related are Equation (1) and (11) of Tsoulis Paper, here referred as V.
     *
     * The first order derivatives of the gravitational potential in [m/s^2].
     * The array contains the derivatives depending on the coordinates x-y-z in this order.
     * Related are Equation (2) and (12) of Tsoulis Paper, here referred as Vx, Vy, Vz.
     *
     * The second order derivatives or also called gradiometric Tensor in [1/s^2].
     * The array contains the second order derivatives in the following order xx, yy, zz, xy, xz, yz.
     * Related are Equation (3) and (13) of Tsoulis Paper, here referred as Vxx, Vyy, Vzz, Vxy, Vxz, Vyz.
     */
    using GravityModelResult = std::tuple<double, Array3, Array6>;

    /** A polyhedron defined by a set of filenames, as available in {@link TetgenAdapter} */
    using PolyhedralFiles = std::vector<std::string>;

    /** A polyhedron defined by a set of vertices and face indices */
    using PolyhedralSource = std::tuple<std::vector<Array3>, std::vector<IndexArray3>>;

}