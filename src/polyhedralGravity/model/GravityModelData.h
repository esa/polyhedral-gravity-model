#pragma once

#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/util/UtilityFloatArithmetic.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <ostream>
#include <tuple>

namespace polyhedralGravity {

    /**
     * Contains the 3D distances l1_pq and l2_pq between P and the endpoints of segment pq and
     * the 1D distances s1_pq and s2_pq between P'' and the segment endpoints.
     * @note This struct is basically a named tuple
     */
    struct Distance {
        /**
         * the 3D distance between computation point P and the first endpoint of line segment pq
         */
        double l1;
        /**
         * the 3D distance between computation point P and the second endpoint of line segment pq
         */
        double l2;
        /**
         * the 1D distance between projection of the computation point on line segment pq and the first endpoint of
         * line segment pq
         */
        double s1;
        /**
         * the 1D distance between projection of the computation point on line segment pq and the second endpoint of
         * line segment pq
         */
        double s2;

        /**
         * Checks two Distance structs for equality with another one by ensuring that the members are
         * almost equal.
         * @param rhs the other Distance struct
         * @return true if equal
         *
         * @note Just used for testing purpose
         */
        bool operator==(const Distance &rhs) const {
            return util::almostEqualRelative(l1, rhs.l1) &&
                   util::almostEqualRelative(l2, rhs.l2) &&
                   util::almostEqualRelative(s1, rhs.s1) &&
                   util::almostEqualRelative(s2, rhs.s2);
        }

        /**
         * Checks two Distance structs for inequality with another one by ensuring that the members are
         * not almost equal.
         * @param rhs the other Distance struct
         * @return false if unequal
         *
         * @note Just used for testing purpose
         */
        bool operator!=(const Distance &rhs) const {
            return !(rhs == *this);
        }

        /**
         * Pretty prints this struct on the given ostream.
         * @param os ostream
         * @param distance a Distance struct
         * @return os
         */
        friend std::ostream &operator<<(std::ostream &os, const Distance &distance) {
            os << "l1: " << distance.l1 << " l2: " << distance.l2 << " s1: " << distance.s1 << " s2: " << distance.s2;
            return os;
        }
    };

    /**
     * Contains the Transcendental Expressions LN_pq and AN_pq for a given line segment pq of the polyhedron.
     * @note This struct is basically a named tuple
     */
    struct TranscendentalExpression {
        /**
         * The LN values for plane p and segment q of this plane is calculated in the following way:
         * LN_pq = ln ((s_2_pq + l_2_pq) / (s_1_pq + l_1_pq))
         * @note see Tsoulis Paper Equation (14)
         */
        double ln;
        /**
         * The AN values for plane p and segment q of this plane is calculated in the following way:
         * AN_pq = arctan ((h_p * s_2_pq) / (h_pq * l_2_pq)) - arctan ((h_pq * s_1_pq) / (h_pq * l_1_pq))
         * @note see Tsoulis Paper Equation (15)
         */
        double an;

        /**
         * Checks two TranscendentalExpressions for equality with another one by ensuring that the members are
         * almost equal.
         * @param rhs the other TranscendentalExpressions
         * @return true if equal
         *
         * @note Just used for testing purpose
         */
        bool operator==(const TranscendentalExpression &rhs) const {
            return util::almostEqualRelative(ln, rhs.ln) && util::almostEqualRelative(an, rhs.an);
        }

        /**
         * Checks two TranscendentalExpressions for inequality with another one by ensuring that the members are
         * not almost equal.
         * @param rhs the other TranscendentalExpressions
         * @return false if unequal
         *
         * @note Just used for testing purpose
         */
        bool operator!=(const TranscendentalExpression &rhs) const {
            return !(rhs == *this);
        }

        /**
         * Pretty output of this struct on the given ostream.
         * @param os the ostream
         * @param expression a TranscendentalExpression
         * @return os
         */
        friend std::ostream &operator<<(std::ostream &os, const TranscendentalExpression &expression) {
            os << "ln: " << expression.ln << " an: " << expression.an;
            return os;
        }
    };


    /**
     * A struct describing a plane in Hessian Normal Form:
     * ax + by + cz + d = 0
     * where a,b,c are the plane's normal
     * and d as the signed distance to the plane from the origin along the normal.
     */
    struct HessianPlane {
        /**
         * part of the planes normal [a, b, c]
         */
        double a;
        /**
         * part of the panes normal [a, b, c]
         */
        double b;
        /**
         * part of the planes normal [a, b, c]
         */
        double c;
        /**
         * the signed distance to the plane from the origin along the normal
         */
        double d;

        /**
         * Checking the equality of two this Hessian Plane with another one by ensuring that the members are
         * almost equal.
         * @param rhs other HessianPlane
         * @return true if equal
         *
         * @note Just used for testing purpose
         */
        bool operator==(const HessianPlane &rhs) const {
            return util::almostEqualRelative(a, rhs.a) &&
                   util::almostEqualRelative(b, rhs.b) &&
                   util::almostEqualRelative(c, rhs.c) &&
                   util::almostEqualRelative(d, rhs.d);
        }

        /**
         * Checking the inequality of two this Hessian Plane with another one by ensuring that the members are
         * not almost equal.
         * @param rhs other HessianPlane
         * @return true if unequal
         *
         * @note Just used for testing purpose
         */
        bool operator!=(const HessianPlane &rhs) const {
            return !(rhs == *this);
        }

         /**
         * Pretty output of this struct on the given ostream.
         * @param os the ostream
         * @param hessianPlane a HessianPlane
         * @return os
         */
        friend std::ostream &operator<<(std::ostream &os, const HessianPlane &hessianPlane) {
            os << "a: " << hessianPlane.a << " b: " << hessianPlane.b << " c: " << hessianPlane.c << " d: " << hessianPlane.d;
            return os;
        }
    };

}
