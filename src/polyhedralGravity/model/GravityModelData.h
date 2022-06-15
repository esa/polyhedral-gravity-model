#pragma once

#include <array>
#include <ostream>
#include <cmath>
#include <algorithm>
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/util/UtilityConstants.h"

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
         * Checks two Distance structs for equality.
         * @warning This method compares doubles! So only exact copies will evaluate to true.
         * @param rhs - the other Distance struct
         * @return true if equal
         *
         * @note Just used for testing purpose
         */
        bool operator==(const Distance &rhs) const {
            return l1 == rhs.l1 &&
                   l2 == rhs.l2 &&
                   s1 == rhs.s1 &&
                   s2 == rhs.s2;
        }

        /**
         * Checks two Distance structs for inequality.
         * @warning This method compares doubles! So only exact copies will evaluate to false.
         * @param rhs - the other Distance struct
         * @return false if unequal
         *
         * @note Just used for testing purpose
         */
        bool operator!=(const Distance &rhs) const {
            return !(rhs == *this);
        }

        /**
         * Pretty prints this struct on the given ostream.
         * @param os - ostream
         * @param distance - a Distance struct
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
         * Checks two TranscendentalExpressions for equality.
         * @warning This method compares doubles! So only exact copies will evaluate to true.
         * @param rhs - the other TranscendentalExpressions
         * @return true if equal
         *
         * @note Just used for testing purpose
         */
        bool operator==(const TranscendentalExpression &rhs) const {
            return ln == rhs.ln &&
                   an == rhs.an;
        }

        /**
         * Checks two TranscendentalExpressions for inequality.
         * @warning This method compares doubles! So only exact copies will evaluate to false.
         * @param rhs - the other TranscendentalExpressions
         * @return false if unequal
         *
         * @note Just used for testing purpose
         */
        bool operator!=(const TranscendentalExpression &rhs) const {
            return !(rhs == *this);
        }

        /**
         * Pretty output of this struct on the given ostream.
         * @param os - the ostream
         * @param expression - a TranscendentalExpression
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
         * Checking the equality of two this Hessian Plane with another one by comparing their members.
         * @warning This method compares doubles! So only exact copies will evaluate to true.
         * @param rhs - other HessianPlane
         * @return true if equal
         *
         * @note Just used for testing purpose
         */
        bool operator==(const HessianPlane &rhs) const {
            return a == rhs.a &&
                   b == rhs.b &&
                   c == rhs.c &&
                   d == rhs.d;
        }

        /**
         * Checking the inequality of two this Hessian Plane with another one by comparing their members.
         * @warning This method compares doubles! So only exact copies will evaluate to false.
         * @param rhs - other HessianPlane
         * @return true if unequal
         *
         * @note Just used for testing purpose
         */
        bool operator!=(const HessianPlane &rhs) const {
            return !(rhs == *this);
        }
    };

    /**
     * A data structure containing the result of the polyhedral gravity model's evaluation.
    */
    struct GravityModelResult {

        /**
         * The gravitational potential in [m^2/s^2] <--> [J/kg] at point P.
         * @related Equation (1) and (11) of Tsoulis Paper, here referred as V
         */
        double gravitationalPotential{};

        /**
         * The first order derivatives of the gravitational potential in [m/s^2].
         * The array contains the derivatives depending on the coordinates x-y-z in this order.
         * @related Equation (2) and (12) of Tsoulis Paper, here referred as Vx, Vy, Vz
         */
        std::array<double, 3> acceleration{};

        /**
         * The second order derivatives or also called gradiometric Tensor in [1/s^2].
         * The array contains the second order derivatives in the following order xx, yy, zz, xy, xz, yz.
         * @related Equation (3) and (13) of Tsoulis Paper, here referred as Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
         */
        std::array<double, 6> gradiometricTensor{};

        /**
         * Creates new empty Result
         */
        GravityModelResult() = default;

        /**
         * Creates a new Result with the given values
         * @param gravitationalPotential - the potential
         * @param gravitationalPotentialDerivative - the acceleration
         * @param gradiometricTensor - the gradiometric tensor i.e. the second derivatives
         */
        GravityModelResult(double gravitationalPotential, const std::array<double, 3> &gravitationalPotentialDerivative,
                           const std::array<double, 6> &gradiometricTensor)
                : gravitationalPotential(gravitationalPotential),
                  acceleration(gravitationalPotentialDerivative),
                  gradiometricTensor(gradiometricTensor) {}

        /**
         * Pretty output of this struct on the given ostream.
         * @param os - the ostream
         * @param result - the GravityModelResult to output
         * @return os
         */
        friend std::ostream &operator<<(std::ostream &os, const GravityModelResult &result) {
            using util::operator<<;
            os << "gravitationalPotential: " << result.gravitationalPotential << " acceleration: "
               << result.acceleration << " gradiometricTensor: " << result.gradiometricTensor;
            return os;
        }

        /**
         * Compares the results to used epsilon and if one is smaller resets it to zero.
         * This step is performed due to rounding errors when calculating with doubles.
         */
        inline void eliminateRoundingErrors(){
            if (std::abs(gravitationalPotential) < util::EPSILON) {
                gravitationalPotential = 0.0;
            }
            for (double &acc : acceleration) {
                if (std::abs(acc) < util::EPSILON) {
                    acc = 0.0;
                }
            }
            for (double &secondDev : gradiometricTensor) {
                if (std::abs(secondDev) < util::EPSILON) {
                    secondDev = 0.0;
                }
            }
        }

    };

}
