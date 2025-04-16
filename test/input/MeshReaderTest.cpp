#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include "polyhedralGravity/input/MeshReader.h"

class MeshReaderTest : public ::testing::Test {

protected:

    std::vector<std::array<double, 3>> _expectedVertices = {
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

TEST_F(MeshReaderTest, readSimpleNode) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{
            "resources/MeshReaderTestReadSimple.node",
            "resources/MeshReaderTestReadSimple.face",
    };
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    ASSERT_THAT(actualVertices, ContainerEq(_expectedVertices));
}


TEST_F(MeshReaderTest, readSimpleFace) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{
            "resources/MeshReaderTestReadSimple.node",
            "resources/MeshReaderTestReadSimple.face"
    };
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    ASSERT_THAT(actualFaces, ContainerEq(_expectedFaces));
}

TEST_F(MeshReaderTest, readSimpleMesh) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{"resources/MeshReaderTestReadSimple.mesh"};
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    for (const auto &actualVertice: actualVertices) {
        ASSERT_THAT(_expectedVertices, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualFaces.size());
}

TEST_F(MeshReaderTest, readSimpleOff) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{"resources/MeshReaderTestReadSimple.off"};
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    for (const auto &actualVertice: actualVertices) {
        ASSERT_THAT(_expectedVertices, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualFaces.size());
}

TEST_F(MeshReaderTest, readSimplePly) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{"resources/MeshReaderTestReadSimple.ply"};
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    for (const auto &actualVertice: actualVertices) {
        ASSERT_THAT(_expectedVertices, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualFaces.size());
}

TEST_F(MeshReaderTest, readSimpleStl) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{"resources/MeshReaderTestReadSimple.stl"};
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    for (const auto &actualVertice: actualVertices) {
        ASSERT_THAT(_expectedVertices, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualFaces.size());
}

TEST_F(MeshReaderTest, readSimpleObj) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{"resources/MeshReaderTestReadSimple.obj"};
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    for (const auto &actualVertice: actualVertices) {
        ASSERT_THAT(_expectedVertices, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualFaces.size());
}

TEST_F(MeshReaderTest, readSimpleTab) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    const std::vector<std::string> simpleFiles{"resources/MeshReaderTestReadSimple.tab"};
    const auto&[actualVertices, actualFaces] = MeshReader::getPolyhedralSource(simpleFiles);

    for (const auto &actualVertice: actualVertices) {
        ASSERT_THAT(_expectedVertices, Contains(actualVertice));
    }
    ASSERT_EQ(_expectedFaces.size(), actualFaces.size());
}