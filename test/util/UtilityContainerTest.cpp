#include "gtest/gtest.h"

#include <utility>
#include <array>
#include "polyhedralGravity/util/UtilityContainer.h"


TEST(UtilityContainer, VectorContainerPlus) {
    using namespace ::polyhedralGravity::util;
    std::array<int, 3> a{3, 4, 5};
    std::array<int, 3> b{6, 8, 10};

    std::array<int, 3> expectedResult{9, 12, 15};
    auto actualResult = a + b;
    ASSERT_EQ(actualResult, expectedResult);
}

TEST(UtilityContainer, VectorScalarPlus) {
    using namespace ::polyhedralGravity::util;
    std::array<int, 3> a{3, 4, 5};
    int b = 100;

    std::array<int, 3> expectedResult{103, 104, 105};
    auto actualResult = a + b;
    ASSERT_EQ(actualResult, expectedResult);
}

TEST(UtilityContainer, VectorContainerMul) {
    using namespace ::polyhedralGravity::util;
    std::array<int, 3> a{3, 4, 5};
    std::array<int, 3> b{6, 8, 10};

    std::array<int, 3> expectedResult{18, 32, 50};
    auto actualResult = a * b;
    ASSERT_EQ(actualResult, expectedResult);
}

TEST(UtilityContainer, VectorScalarMul) {
    using namespace ::polyhedralGravity::util;
    std::array<int, 3> a{3, 4, 5};
    int b = 100;

    std::array<int, 3> expectedResult{300, 400, 500};
    auto actualResult = a * b;
    ASSERT_EQ(actualResult, expectedResult);
}

TEST(UtilityContainer, Determinant_1) {
    using namespace ::polyhedralGravity::util;
    Matrix<double, 3, 3> matrix{{{3, 0, 1}, {1, 2, 5}, {-1, 4, 2}}};

    double expectedDeterminant = -42.0;
    double actualDeterminant = det(matrix);
    ASSERT_DOUBLE_EQ(actualDeterminant, expectedDeterminant);
}

TEST(UtilityContainer, Determinant_2) {
    using namespace ::polyhedralGravity::util;
    Matrix<int, 3, 3> matrix{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};

    double expectedDeterminant = 0.0;
    double actualDeterminant = det(matrix);
    ASSERT_DOUBLE_EQ(actualDeterminant, expectedDeterminant);
}