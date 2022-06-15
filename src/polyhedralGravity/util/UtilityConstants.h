#pragma once

namespace polyhedralGravity::util {

    /**
     * The EPSILON used in the polyhedral gravity model.
     * @related Used to determine if a floating point number is equal to zero as threshold for rounding errors
     * @related Used for the sgn() function to determine the sign of a double value. Different compilers
     * produce different results if no EPSILON is applied for the comparison!
     */
    constexpr double EPSILON = 1e-14;


    /**
     * PI with enough precision
     */
    constexpr double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;


    /**
    * PI multiplied by 2 with enough precision
    */
    constexpr double PI2 = 6.2831853071795864769252867665590057683943387987502116419498891846;


    /**
    * PI divided by 2 with enough precision
    */
    constexpr double PI_2 = 1.5707963267948966192313216916397514420985846996875529104874722961;

    /**
     * The gravitational constant G in [m^3/(kg*s^2)].
     * @related in his paper above Equation (4)
     */
    constexpr double GRAVITATIONAL_CONSTANT = 6.67430e-11;

    /**
     * The assumed constant density rho for a polyhedron after Tsoulis paper in [kg/m^3].
     * @related in his paper above Equation (4)
     */
    constexpr double DEFAULT_CONSTANT_DENSITY = 2670.0;

}
