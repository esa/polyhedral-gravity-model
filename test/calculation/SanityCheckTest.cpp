#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <tuple>
#include <fstream>
#include <sstream>
#include <string>
#include "polyhedralGravity/calculation/SanityCheck.h"
#include "polyhedralGravity/model/Polyhedron.h"


/**
 * Contains Tests if the sanity checks works
 */
class SanityCheckTest : public ::testing::Test {

protected:

    const polyhedralGravity::Polyhedron _correctCube{{
                                                      {-1.0, -1.0, -1.0},
                                                      {1.0, -1.0, -1.0},
                                                      {1.0, 1.0, -1.0},
                                                      {-1.0, 1.0, -1.0},
                                                      {-1.0, -1.0, 1.0},
                                                      {1.0, -1.0, 1.0},
                                                      {1.0, 1.0, 1.0},
                                                      {-1.0, 1.0, 1.0}},
                                              {
                                                      {1,    3,    2},
                                                      {0,   3,    1},
                                                      {0,   1,   5},
                                                      {0,    5,   4},
                                                      {0,    7,    3},
                                                      {0,   4,    7},
                                                      {1,   2,   6},
                                                      {1,    6,   5},
                                                      {2, 3, 6},
                                                      {3, 7, 6},
                                                      {4, 5, 6},
                                                      {4, 6, 7}}
    };

    const polyhedralGravity::Polyhedron _wrongCube{{
                                                             {-1.0, -1.0, -1.0},
                                                             {1.0, -1.0, -1.0},
                                                             {1.0, 1.0, -1.0},
                                                             {-1.0, 1.0, -1.0},
                                                             {-1.0, -1.0, 1.0},
                                                             {1.0, -1.0, 1.0},
                                                             {1.0, 1.0, 1.0},
                                                             {-1.0, 1.0, 1.0}},
                                                     {
                                                             {3,    1,    2},
                                                             {3,   0,    1},
                                                             {1,   0,   5},
                                                             {5,    0,   4},
                                                             {7,    0,    3},
                                                             {4,   0,    7},
                                                             {2,   1,   6},
                                                             {6,    1,   5},
                                                             {3, 2, 6},
                                                             {7, 3, 6},
                                                             {5, 4, 6},
                                                             {6, 4, 7}}
    };
};

TEST_F(SanityCheckTest, CorrectCube) {
    using namespace testing;
    ASSERT_TRUE(polyhedralGravity::SanityCheck::checkNormalsOutwardPointing(_correctCube));
}

TEST_F(SanityCheckTest, WrongCube) {
    using namespace testing;
    ASSERT_FALSE(polyhedralGravity::SanityCheck::checkNormalsOutwardPointing(_wrongCube));
}
