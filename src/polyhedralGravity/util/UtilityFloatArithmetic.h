#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace polyhedralGravity::util {

    /**
     * The maximal allowed ULP distance as computed by {@ref polyhedralGravity::util::floatNear}
     * utilized for FloatingPoint comparisons in this implementation
     *
     * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
     */
    constexpr int MAX_ULP_DISTANCE = 4;


    /**
     * Function for comparing closeness of two floating point numbers using ULP (Units in the Last Place) method.
     *
     * @param lhs The left hand side floating point number to compare.
     * @param rhs The right hand side floating point number to compare.
     * @param ulpDistance The maximum acceptable ULP distance between the two floating points
     *  for which they would be considered near each other. This is optional and by default, it will be MAX_ULP_DISTANCE.
     *
     * @return true if the ULP distance between lhs and rhs is less than or equal to the provided ulpDistance value, otherwise, false.
     *  Returns true if both numbers are exactly the same. Returns false if the signs do not match.
     *
     * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
     */
    bool floatNear(double lhs, double rhs, int ulpDistance = MAX_ULP_DISTANCE);

    /**
     * Function for comparing closeness of two floating point numbers using ULP (Units in the Last Place) method.
     *
     * @param lhs The left hand side floating point number to compare.
     * @param rhs The right hand side floating point number to compare.
     * @param ulpDistance The maximum acceptable ULP distance between the two floating points
     *  for which they would be considered near each other. This is optional and by default, it will be MAX_ULP_DISTANCE.
     *
     * @return true if the ULP distance between lhs and rhs is less than or equal to the provided ulpDistance value, otherwise, false.
     *  Returns true if both numbers are exactly the same. Returns false if the signs do not match.
     *
     * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
     */
    bool floatNear(float lhs, float rhs, int ulpDistance = MAX_ULP_DISTANCE);

}
