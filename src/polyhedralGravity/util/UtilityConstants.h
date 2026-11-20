#pragma once

namespace polyhedralGravity::util {

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
     * The gravitational constant G in @f$[m^3/(kg*s^2)]@f$.
     * In Tsoulis et al. (2012) above Equation (4)
     */
    constexpr double GRAVITATIONAL_CONSTANT = 6.67430e-11;

    /**
     * The assumed constant density rho for a polyhedron after Tsoulis paper in @f$[kg/m^3]@f$.
     * In Tsoulis et al. (2012) above Equation (4)
     */
    constexpr double DEFAULT_CONSTANT_DENSITY = 2670.0;

}
