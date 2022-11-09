#include "gtest/gtest.h"

#include <array>
#include <utility>
#include <tuple>
#include "polyhedralGravity/calculation/SanityCheck.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/input/TetgenAdapter.h"


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
                                                             {1, 3, 2},
                                                             {0, 3, 1},
                                                             {0, 1, 5},
                                                             {0, 5, 4},
                                                             {0, 7, 3},
                                                             {0, 4, 7},
                                                             {1, 2, 6},
                                                             {1, 6, 5},
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
                                                           {3, 1, 2},
                                                           {3, 0, 1},
                                                           {1, 0, 5},
                                                           {5, 0, 4},
                                                           {7, 0, 3},
                                                           {4, 0, 7},
                                                           {2, 1, 6},
                                                           {6, 1, 5},
                                                           {3, 2, 6},
                                                           {7, 3, 6},
                                                           {5, 4, 6},
                                                           {6, 4, 7}}
    };


    const polyhedralGravity::Polyhedron _correctPrism{{
                                                              {-20.0, 0.0, 25.0},
                                                              {0.0, 0.0, 25.0},
                                                              {0.0, 10.0, 25.0},
                                                              {-20.0, 10.0, 25.0},
                                                              {-20.0, 0.0, 15.0},
                                                              {0.0, 0.0, 15.0},
                                                              {0.0, 10.0, 15.0},
                                                              {-20.0, 10.0, 15.0}},
                                                      {
                                                              {0, 4, 5},
                                                              {0, 5, 1},
                                                              {0, 1, 3},
                                                              {1, 2, 3},
                                                              {1, 5, 6},
                                                              {1, 6, 2},
                                                              {0, 7, 4},
                                                              {0, 3, 7},
                                                              {4, 7, 5},
                                                              {5, 7, 6},
                                                              {2, 7, 3},
                                                              {2, 6, 7}
                                                      }};

    const polyhedralGravity::Polyhedron _wrongPrism{{
                                                            {-20.0, 0.0, 25.0},
                                                            {0.0, 0.0, 25.0},
                                                            {0.0, 10.0, 25.0},
                                                            {-20.0, 10.0, 25.0},
                                                            {-20.0, 0.0, 15.0},
                                                            {0.0, 0.0, 15.0},
                                                            {0.0, 10.0, 15.0},
                                                            {-20.0, 10.0, 15.0}},
                                                    {
                                                            {4, 0, 5},
                                                            {5, 0, 1},
                                                            {1, 0, 3},
                                                            {2, 1, 3},
                                                            {5, 1, 6},
                                                            {6, 1, 2},
                                                            {7, 0, 4},
                                                            {3, 0, 7},
                                                            {7, 4, 5},
                                                            {7, 5, 6},
                                                            {7, 2, 3},
                                                            {6, 2, 7}
                                                    }};
};

TEST_F(SanityCheckTest, CorrectCube) {
    using namespace testing;
    ASSERT_TRUE(polyhedralGravity::SanityCheck::checkNormalsOutwardPointing(_correctCube))
    << "The vertices of the cube are correctly sorted, however the Sanity Check came to the wrong conclusion!";
}

TEST_F(SanityCheckTest, WrongCube) {
    using namespace testing;
    ASSERT_FALSE(polyhedralGravity::SanityCheck::checkNormalsOutwardPointing(_wrongCube))
    << "The vertices of the cube are incorrectly sorted, however the Sanity Check failed to notice that!";
}

TEST_F(SanityCheckTest, CorrectPrism) {
    using namespace testing;
    ASSERT_TRUE(polyhedralGravity::SanityCheck::checkNormalsOutwardPointing(_correctPrism))
    << "The vertices of the prism are correctly sorted, however the Sanity Check came to the wrong conclusion!";
}

TEST_F(SanityCheckTest, WrongPrism) {
    using namespace testing;
    ASSERT_FALSE(polyhedralGravity::SanityCheck::checkNormalsOutwardPointing(_wrongPrism))
    << "The vertices of the prism are incorrectly sorted, however the Sanity Check failed to notice that!";
}

TEST_F(SanityCheckTest, CorrectBigPolyhedron) {
    using namespace testing;
    using namespace polyhedralGravity;
    Polyhedron polyhedron{
            TetgenAdapter{
                    {"resources/GravityModelBigTest.node", "resources/GravityModelBigTest.face"}}.getPolyhedron()};
    ASSERT_TRUE(polyhedralGravity::SanityCheck::checkNormalsOutwardPointing(polyhedron))
    << "The vertices of the polyhedron are correctly sorted, however the Sanity Check came to the wrong conclusion!";
}

