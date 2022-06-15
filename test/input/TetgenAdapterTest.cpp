#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include "polyhedralGravity/input/TetgenAdapter.h"

class TetgenAdapterTest : public ::testing::Test {

protected:

    std::vector<std::array<double, 3>> _expectedNodes = {
            {-20, 0,  25},
            {0,   0,  25},
            {0,   10, 25},
            {-20, 10, 25},
            {-20, 0,  15},
            {0,   0,  15},
            {0,   10, 15},
            {-20, 10, 15}
    };

    std::vector<std::array<size_t, 3>> _expectedFaces = {
            {0, 1, 3},
            {1, 2, 3},
            {0, 4, 5},
            {0, 5, 1},
            {0, 7, 4},
            {0, 3, 7},
            {1, 5, 6},
            {1, 6, 2},
            {3, 6, 7},
            {2, 6, 3},
            {4, 6, 5},
            {4, 7, 6}
    };

};

TEST_F(TetgenAdapterTest, readSimpleNode) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    std::vector<std::string> simpleFiles{
            "resources/TetgenAdapterTestReadSimple.node",
            "resources/TetgenAdapterTestReadSimple.face",
    };

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    ASSERT_THAT(actualPolyhedron.getVertices(), ContainerEq(_expectedNodes));
}


TEST_F(TetgenAdapterTest, readSimpleFace) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    std::vector<std::string> simpleFiles{
            "resources/TetgenAdapterTestReadSimple.node",
            "resources/TetgenAdapterTestReadSimple.face"
    };

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    ASSERT_THAT(actualPolyhedron.getFaces(), ContainerEq(_expectedFaces));
}

TEST_F(TetgenAdapterTest, readSimpleMesh) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    std::vector<std::string> simpleFiles{"resources/TetgenAdapterTestReadSimple.mesh"};

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    for (const auto &actualVertice: actualPolyhedron.getVertices()) {
        ASSERT_THAT(_expectedNodes, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualPolyhedron.countFaces());
}

TEST_F(TetgenAdapterTest, readSimpleOff) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    std::vector<std::string> simpleFiles{"resources/TetgenAdapterTestReadSimple.off"};

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    for (const auto &actualVertice: actualPolyhedron.getVertices()) {
        ASSERT_THAT(_expectedNodes, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualPolyhedron.countFaces());
}

TEST_F(TetgenAdapterTest, readSimplePly) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    std::vector<std::string> simpleFiles{"resources/TetgenAdapterTestReadSimple.ply"};

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    for (const auto &actualVertice: actualPolyhedron.getVertices()) {
        ASSERT_THAT(_expectedNodes, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualPolyhedron.countFaces());
}

TEST_F(TetgenAdapterTest, readSimpleStl) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    std::vector<std::string> simpleFiles{"resources/TetgenAdapterTestReadSimple.stl"};

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    for (const auto &actualVertice: actualPolyhedron.getVertices()) {
        ASSERT_THAT(_expectedNodes, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualPolyhedron.countFaces());
}