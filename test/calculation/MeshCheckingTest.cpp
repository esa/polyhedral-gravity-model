#include "gtest/gtest.h"

#include <array>
#include <utility>
#include <tuple>
#include "polyhedralGravity/calculation/MeshChecking.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/input/TetgenAdapter.h"


/**
 * Contains Tests if the sanity checks works
 */
class MeshCheckingTest : public ::testing::Test {

protected:
    const std::vector<polyhedralGravity::Array3> _cubeVertices{
            {-1.0, -1.0, -1.0},
            {1.0, -1.0, -1.0},
            {1.0, 1.0, -1.0},
            {-1.0, 1.0, -1.0},
            {-1.0, -1.0, 1.0},
            {1.0, -1.0, 1.0},
            {1.0, 1.0, 1.0},
            {-1.0, 1.0, 1.0}
    };

    const polyhedralGravity::Polyhedron _cubeOutwardNormals{_cubeVertices,
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

    const polyhedralGravity::Polyhedron _cubeInwardNormals{_cubeVertices,
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

    const polyhedralGravity::Polyhedron _cubeMajorOutwardNormals{_cubeVertices,
                                                                 {
                                                                         {3, 1, 2},
                                                                         {0, 3, 1},
                                                                         {0, 1, 5},
                                                                         {0, 5, 4},
                                                                         {7, 0, 3},
                                                                         {0, 4, 7},
                                                                         {1, 2, 6},
                                                                         {1, 6, 5},
                                                                         {2, 3, 6},
                                                                         {3, 7, 6},
                                                                         {4, 5, 6},
                                                                         {4, 6, 7}}
    };

    const polyhedralGravity::Polyhedron _cubeMajorInwardNormals{_cubeVertices,
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
                                                                        {3, 7, 6},
                                                                        {4, 5, 6},
                                                                        {4, 6, 7}}
    };

    const polyhedralGravity::Polyhedron _degeneratedCube{_cubeVertices,
                                                         {
                                                                 {1, 3, 2},
                                                                 {0, 3, 1},
                                                                 {0, 1, 5},
                                                                 {0, 5, 4},
                                                                 {7, 7, 3},
                                                                 {0, 4, 7},
                                                                 {1, 2, 6},
                                                                 {1, 6, 5},
                                                                 {2, 3, 6},
                                                                 {3, 7, 6},
                                                                 {4, 5, 6},
                                                                 {4, 6, 7}}
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

TEST_F(MeshCheckingTest, CubeOutwardNormals) {
    using namespace polyhedralGravity;
    using namespace testing;
    // All normals are pointing outwards
    ASSERT_TRUE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeOutwardNormals, NormalOrientation::OUTWARDS));
    ASSERT_FALSE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeOutwardNormals, NormalOrientation::INWWARDS));

    const auto &[actualOrientation, actualViolatingFaces] = MeshChecking::getPlaneUnitNormalOrientation(_cubeOutwardNormals);
    ASSERT_EQ(actualOrientation, NormalOrientation::OUTWARDS);
    ASSERT_EQ(actualViolatingFaces, std::set<size_t>{});

    ASSERT_TRUE(MeshChecking::checkTrianglesNotDegenerated(_cubeOutwardNormals));
}

TEST_F(MeshCheckingTest, CubeInwardNormals) {
    using namespace polyhedralGravity;
    using namespace testing;
    // All normals are pointing inwards (reversed order)
    ASSERT_FALSE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeInwardNormals, NormalOrientation::OUTWARDS));
    ASSERT_TRUE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeInwardNormals, NormalOrientation::INWWARDS));

    const auto &[actualOrientation, actualViolatingFaces] = MeshChecking::getPlaneUnitNormalOrientation(_cubeInwardNormals);
    ASSERT_EQ(actualOrientation, NormalOrientation::INWWARDS);
    ASSERT_EQ(actualViolatingFaces, std::set<size_t>{});

    ASSERT_TRUE(MeshChecking::checkTrianglesNotDegenerated(_cubeInwardNormals));
}

TEST_F(MeshCheckingTest, CubeMajorOutwardNormals) {
    using namespace polyhedralGravity;
    using namespace testing;
    // All normals are pointing outwards expect two which points inwards (at indices 0 and 4)
    ASSERT_FALSE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeMajorOutwardNormals, NormalOrientation::INWWARDS));
    ASSERT_FALSE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeMajorOutwardNormals, NormalOrientation::OUTWARDS));

    const auto &[actualOrientation, actualViolatingFaces] = MeshChecking::getPlaneUnitNormalOrientation(_cubeMajorOutwardNormals);
    ASSERT_EQ(actualOrientation, NormalOrientation::OUTWARDS);
    ASSERT_EQ(actualViolatingFaces, std::set<size_t>({0, 4}));

    ASSERT_TRUE(MeshChecking::checkTrianglesNotDegenerated(_cubeMajorOutwardNormals));
}

TEST_F(MeshCheckingTest, CubeMajorInwardNormals) {
    using namespace polyhedralGravity;
    using namespace testing;
    // All normals are pointing outwards expect three which points inwards (at indices 9, 10, 11)
    ASSERT_FALSE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeMajorInwardNormals, NormalOrientation::INWWARDS));
    ASSERT_FALSE(MeshChecking::checkPlaneUnitNormalOrientation(_cubeMajorInwardNormals, NormalOrientation::OUTWARDS));

    const auto &[actualOrientation, actualViolatingFaces] = MeshChecking::getPlaneUnitNormalOrientation(_cubeMajorInwardNormals);
    ASSERT_EQ(actualOrientation, NormalOrientation::INWWARDS);
    ASSERT_EQ(actualViolatingFaces, std::set<size_t>({9, 10, 11}));

    ASSERT_TRUE(MeshChecking::checkTrianglesNotDegenerated(_cubeMajorInwardNormals));
}

TEST_F(MeshCheckingTest, DegeneratedCube) {
    using namespace testing;
    // One surface has the same point twice as vertex
    ASSERT_FALSE(polyhedralGravity::MeshChecking::checkTrianglesNotDegenerated(_degeneratedCube));
}

TEST_F(MeshCheckingTest, CorrectPrism) {
    using namespace polyhedralGravity;
    using namespace testing;
    // All normals are pointing outwards
    ASSERT_TRUE(polyhedralGravity::MeshChecking::checkPlaneUnitNormalOrientation(_correctPrism))
    << "The vertices of the prism are correctly sorted, however the Sanity Check came to the wrong conclusion!";

    const auto &[actualOrientation, actualViolatingFaces] = MeshChecking::getPlaneUnitNormalOrientation(_correctPrism);
    ASSERT_EQ(actualOrientation, NormalOrientation::OUTWARDS);
    ASSERT_EQ(actualViolatingFaces, std::set<size_t>{});

    ASSERT_TRUE(MeshChecking::checkTrianglesNotDegenerated(_correctPrism));
}

TEST_F(MeshCheckingTest, WrongPrism) {
    using namespace polyhedralGravity;
    using namespace testing;
    // All normals are pointing inwards (reversed order)
    ASSERT_FALSE(MeshChecking::checkPlaneUnitNormalOrientation(_wrongPrism))
    << "The vertices of the prism are incorrectly sorted, however the Sanity Check failed to notice that!";

    const auto &[actualOrientation, actualViolatingFaces] = MeshChecking::getPlaneUnitNormalOrientation(_wrongPrism);
    ASSERT_EQ(actualOrientation, NormalOrientation::INWWARDS);
    ASSERT_EQ(actualViolatingFaces, std::set<size_t>{});

    ASSERT_TRUE(MeshChecking::checkTrianglesNotDegenerated(_wrongPrism));
}

TEST_F(MeshCheckingTest, CorrectBigPolyhedron) {
    using namespace testing;
    using namespace polyhedralGravity;
    // All normals are pointing outwards, extensive Eros example
    Polyhedron polyhedron{
            TetgenAdapter{
                    {"resources/GravityModelBigTest.node", "resources/GravityModelBigTest.face"}}.getPolyhedron()};
    ASSERT_TRUE(MeshChecking::checkPlaneUnitNormalOrientation(polyhedron))
    << "The vertices of the polyhedron are correctly sorted, however the Sanity Check came to the wrong conclusion!";
    ASSERT_TRUE(MeshChecking::checkTrianglesNotDegenerated(polyhedron));

    const auto &[actualOrientation, actualViolatingFaces] = MeshChecking::getPlaneUnitNormalOrientation(polyhedron);
    ASSERT_EQ(actualOrientation, NormalOrientation::OUTWARDS);
    ASSERT_EQ(actualViolatingFaces, std::set<size_t>{});
}

