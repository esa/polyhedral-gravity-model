#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <array>
#include "polyhedralGravity/model/Polyhedron.h"


/**
 * Contains Tests if the mesh sanity checks of the Polyhedron work
 * Parameterization not possible since the patterns what and where the error exactly is varies leading
 * to different assertions being required in the single test case.
 */
class PolyhedronTest : public ::testing::Test {

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

    // The indexing starts at 1 rather than at zero
    const std::vector<polyhedralGravity::IndexArray3> _facesCorrection{
                {2, 4, 3},
                {1, 4, 2},
                {1, 2, 6},
                {1, 6, 5},
                {1, 8, 4},
                {1, 5, 8},
                {2, 3, 7},
                {2, 7, 6},
                {3, 4, 7},
                {4, 8, 7},
                {5, 6, 7},
                {5, 7, 8}
    };

    const std::vector<polyhedralGravity::IndexArray3> _facesOutwards{
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
            {4, 6, 7}
    };


    const std::vector<polyhedralGravity::IndexArray3> _facesInwards{
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
            {6, 4, 7}
    };

    const std::vector<polyhedralGravity::IndexArray3> _facesOutwardsMajority{
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
            {4, 6, 7}
    };

    const std::vector<polyhedralGravity::IndexArray3> _facesInwardsMajority{
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
            {4, 6, 7}
    };

    const std::vector<polyhedralGravity::IndexArray3> _degeneratedFaces{
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
            {4, 6, 7}
    };


    const std::vector<polyhedralGravity::Array3> _prismVertices{
            {-20.0, 0.0, 25.0},
            {0.0, 0.0, 25.0},
            {0.0, 10.0, 25.0},
            {-20.0, 10.0, 25.0},
            {-20.0, 0.0, 15.0},
            {0.0, 0.0, 15.0},
            {0.0, 10.0, 15.0},
            {-20.0, 10.0, 15.0}
    };

    const std::vector<polyhedralGravity::IndexArray3> _prismOutwards{
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
    };

    const std::vector<polyhedralGravity::IndexArray3> _prismInwards{

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
    };
};

TEST_F(PolyhedronTest, FaceCorrection) {
    using namespace polyhedralGravity;
    using namespace testing;

    const auto polyhedron = Polyhedron(_cubeVertices, _facesCorrection, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE);
    // The faces indices need to be shifted by -1, i.e., the indexing starts with vertex 0 not 1
    ASSERT_THAT(polyhedron.getFaces(), ContainerEq(_facesOutwards));
}

TEST_F(PolyhedronTest, CubeOutwardNormals) {
    using namespace polyhedralGravity;
    using namespace testing;
    // Correct Set-Up, No-throw for all options
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::AUTOMATIC));
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY));
    EXPECT_NO_THROW(Polyhedron (_cubeVertices, _facesOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL));

    // Wrong Set-Up, Throws in case of AUTOMATIC and VERIFY, DISABLE and HEAL do not throw but respectivley ignore or repair
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL));

    // Checking that the healing leads to the correct result (More than n/2 normals (actuall all) are outwards. Hence, only change of the Orientation)
    Polyhedron healedPolyhedron(_cubeVertices, _facesOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL);
    EXPECT_EQ(healedPolyhedron.getOrientation(), NormalOrientation::OUTWARDS);
    EXPECT_THAT(healedPolyhedron.getFaces(), ContainerEq(_facesOutwards));
}

TEST_F(PolyhedronTest, CubeInwardsNormals) {
    using namespace polyhedralGravity;
    using namespace testing;
    // Correct Set-Up, No-throw for all options
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::AUTOMATIC));
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::VERIFY));
    EXPECT_NO_THROW(Polyhedron (_cubeVertices, _facesInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL));

    // Wrong Set-Up, Throws in case of AUTOMATIC and VERIFY, DISABLE and HEAL do not throw but respectivley ignore or repair
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL));

    // Checking that checkPlaneUnitNormalOrientation() returns the correct set of indices & majority orientation
    Polyhedron intenionalDisabledCheckingPolyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE);
    const auto &[majorityOrientation, violatingIndices] = intenionalDisabledCheckingPolyhedron.checkPlaneUnitNormalOrientation();
    EXPECT_EQ(majorityOrientation, NormalOrientation::OUTWARDS);
    EXPECT_THAT(violatingIndices, ContainerEq(std::set<size_t>({0, 4})));

    // Checking that the healing leads to the correct result (More than n/2 normals (actuall all) are outwards. Hence, only change of the Orientation)
    Polyhedron healedPolyhedron(_cubeVertices, _facesInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL);
    EXPECT_EQ(healedPolyhedron.getOrientation(), NormalOrientation::INWARDS);
    EXPECT_THAT(healedPolyhedron.getFaces(), ContainerEq(_facesInwards));
}

TEST_F(PolyhedronTest, CubeOutwardNormalsMajor) {
    using namespace polyhedralGravity;
    using namespace testing;
    // The majority is correct. However, 2 faces have inward pointing normals (Index 0 and 4) --> Will throw but not in the healing scenario
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron (_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL));

    // Wrong Set-Up, Throws in case of AUTOMATIC and VERIFY, DISABLE and HEAL do not throw but respectivley ignore or repair
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL));

    // Checking that checkPlaneUnitNormalOrientation() returns the correct set of indices & majority orientation
    Polyhedron intenionalDisabledCheckingPolyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE);
    const auto &[majorityOrientation, violatingIndices] = intenionalDisabledCheckingPolyhedron.checkPlaneUnitNormalOrientation();
    EXPECT_EQ(majorityOrientation, NormalOrientation::OUTWARDS);
    EXPECT_THAT(violatingIndices, ContainerEq(std::set<size_t>({0, 4})));

    // Checking that the healing leads to the correct result (More than n/2 normals are outwards. Here orientation & faces' vertex ordering changes)
    Polyhedron healedPolyhedron(_cubeVertices, _facesOutwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL);
    EXPECT_EQ(healedPolyhedron.getOrientation(), NormalOrientation::OUTWARDS);
    EXPECT_THAT(healedPolyhedron.getFaces(), ContainerEq(_facesOutwards));
}

TEST_F(PolyhedronTest, CubeOutwardInwardsMajor) {
    using namespace polyhedralGravity;
    using namespace testing;
    // The majority is correct. However, 2 faces have inward pointing normals (Index 9, 10, and 11) --> Will throw but not in the healing scenario
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron (_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL));

    // Wrong Set-Up, Throws in case of AUTOMATIC and VERIFY, DISABLE and HEAL do not throw but respectivley ignore or repair
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL));

    // Checking that checkPlaneUnitNormalOrientation() returns the correct set of indices & majority orientation
    Polyhedron intenionalDisabledCheckingPolyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE);
    const auto &[majorityOrientation, violatingIndices] = intenionalDisabledCheckingPolyhedron.checkPlaneUnitNormalOrientation();
    EXPECT_EQ(majorityOrientation, NormalOrientation::INWARDS);
    EXPECT_THAT(violatingIndices, ContainerEq(std::set<size_t>({9, 10, 11})));

    // Checking that the healing leads to the correct result (More than n/2 normals are outwards. Here orientation & faces' vertex ordering changes)
    Polyhedron healedPolyhedron(_cubeVertices, _facesInwardsMajority, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL);
    EXPECT_EQ(healedPolyhedron.getOrientation(), NormalOrientation::INWARDS);
    EXPECT_THAT(healedPolyhedron.getFaces(), ContainerEq(_facesInwards));
}

TEST_F(PolyhedronTest, CubeDegenerated) {
    using namespace polyhedralGravity;
    using namespace testing;
    // A degenarted mesh is not repaired, throw always expect DISABLE is checked
    for (auto orientation: std::set({NormalOrientation::INWARDS, NormalOrientation::OUTWARDS})) {
        EXPECT_NO_THROW(Polyhedron(_cubeVertices, _degeneratedFaces, 1.0, orientation, PolyhedronIntegrity::DISABLE));
        EXPECT_THROW(Polyhedron(_cubeVertices, _degeneratedFaces, 1.0, orientation, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
        EXPECT_THROW(Polyhedron(_cubeVertices, _degeneratedFaces, 1.0, orientation, PolyhedronIntegrity::VERIFY), std::invalid_argument);
        EXPECT_THROW(Polyhedron(_cubeVertices, _degeneratedFaces, 1.0, orientation, PolyhedronIntegrity::HEAL), std::invalid_argument);
    }
}

TEST_F(PolyhedronTest, PrsimOutwards) {
    using namespace polyhedralGravity;
    using namespace testing;
    // Correct Set-Up, No-throw for all options
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::AUTOMATIC));
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY));
    EXPECT_NO_THROW(Polyhedron (_prismVertices, _prismOutwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL));

    // Wrong Set-Up, Throws in case of AUTOMATIC and VERIFY, DISABLE and HEAL do not throw but respectivley ignore or repair
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL));

    // Checking that the healing leads to the correct result (More than n/2 normals (actuall all) are outwards. Hence, only change of the Orientation)
    Polyhedron healedPolyhedron(_prismVertices, _prismOutwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL);
    EXPECT_EQ(healedPolyhedron.getOrientation(), NormalOrientation::OUTWARDS);
    EXPECT_THAT(healedPolyhedron.getFaces(), ContainerEq(_prismOutwards));
}

TEST_F(PolyhedronTest, PrsimInwards) {
    using namespace polyhedralGravity;
    using namespace testing;
    // Correct Set-Up, No-throw for all options
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::AUTOMATIC));
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::VERIFY));
    EXPECT_NO_THROW(Polyhedron (_prismVertices, _prismInwards, 1.0, NormalOrientation::INWARDS, PolyhedronIntegrity::HEAL));

    // Wrong Set-Up, Throws in case of AUTOMATIC and VERIFY, DISABLE and HEAL do not throw but respectivley ignore or repair
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE));
    EXPECT_THROW(Polyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::AUTOMATIC), std::invalid_argument);
    EXPECT_THROW(Polyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY), std::invalid_argument);
    EXPECT_NO_THROW(Polyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL));

    // Checking that the healing leads to the correct result (More than n/2 normals (actuall all) are outwards. Hence, only change of the Orientation)
    Polyhedron healedPolyhedron(_prismVertices, _prismInwards, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL);
    EXPECT_EQ(healedPolyhedron.getOrientation(), NormalOrientation::INWARDS);
    EXPECT_THAT(healedPolyhedron.getFaces(), ContainerEq(_prismInwards));
}

TEST_F(PolyhedronTest, CorrectBigPolyhedron) {
    using namespace testing;
    using namespace polyhedralGravity;
    // All normals are pointing outwards, extensive Eros example
    ASSERT_NO_THROW(Polyhedron(
                std::vector<std::string>({"resources/GravityModelBigTest.node", "resources/GravityModelBigTest.face"}),
                1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY)
            );
}

