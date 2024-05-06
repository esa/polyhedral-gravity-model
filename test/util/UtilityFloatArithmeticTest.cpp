#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <algorithm>
#include <cmath>

#include "polyhedralGravity/util/UtilityFloatArithmetic.h"


TEST(UtilityFloatArithmeticTest, Test01) {
    ASSERT_FALSE(polyhedralGravity::util::floatNear(3.0, 4.0));
    ASSERT_FALSE(polyhedralGravity::util::floatNear(-3.0, -4.0));
    ASSERT_FALSE(polyhedralGravity::util::floatNear(-3.0, 4.0));
    ASSERT_FALSE(polyhedralGravity::util::floatNear(3.0, -4.0));
    ASSERT_TRUE(polyhedralGravity::util::floatNear(1.0, 1.0));
}

TEST(UtilityFloatArithmeticTest, Test02) {
    ASSERT_TRUE(polyhedralGravity::util::floatNear(9.40569e-05, 9.40569e-05));
    ASSERT_TRUE(polyhedralGravity::util::floatNear(-0.000150712, -0.000150712));
    ASSERT_TRUE(polyhedralGravity::util::floatNear(0.000135291, 0.000135291));
    ASSERT_TRUE(polyhedralGravity::util::floatNear(-8.63978e-05, -8.63978e-05));
}

TEST(UtilityFloatArithmeticTest, Test03) {
    double x = 300.3;
    double y = std::nextafter(x, 400.0);
    double z = std::nextafter(y, 400.0);
    ASSERT_EQ(1, reinterpret_cast<int64_t&>(y) - reinterpret_cast<int64_t&>(x));
    ASSERT_EQ(2, reinterpret_cast<int64_t&>(z) - reinterpret_cast<int64_t&>(x));
}