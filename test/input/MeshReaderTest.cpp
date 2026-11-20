#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
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
/**
 * Writes the given content into a temporary .obj file which is removed again on destruction.
 */
class TemporaryObjFile {
public:
    explicit TemporaryObjFile(const std::string &content)
        : _path{std::filesystem::temp_directory_path() /
                ("MeshReaderTest_" + std::to_string(std::hash<std::string>{}(content)) + ".obj")} {
        std::ofstream file{_path};
        file << content;
    }

    ~TemporaryObjFile() {
        std::error_code errorCode{};
        std::filesystem::remove(_path, errorCode);
    }

    [[nodiscard]] std::string string() const {
        return _path.string();
    }

private:
    const std::filesystem::path _path;
};

TEST_F(MeshReaderTest, readObjIgnoresUnrelatedKeywords) {
    using namespace testing;
    using namespace ::polyhedralGravity;
    // Normals, texture coordinates and grouping must not be mistaken for vertices
    const TemporaryObjFile objFile{
            "o cube\n"
            "v 0.0 0.0 0.0\n"
            "vn 0.0 0.0 1.0\n"
            "vt 0.5 0.5\n"
            "vp 0.1 0.2\n"
            "v 1.0 0.0 0.0\n"
            "g group\n"
            "s off\n"
            "v 0.0 1.0 0.0\n"
            "f 0 1 2\n"};

    const auto &[actualVertices, actualFaces] = MeshReader::readObj(objFile.string());

    ASSERT_EQ(actualVertices.size(), 3);
    ASSERT_EQ(actualFaces.size(), 1);
}

TEST_F(MeshReaderTest, readObjSupportsVertexTextureNormalNotation) {
    using namespace testing;
    using namespace ::polyhedralGravity;
    // Only the vertex part of the v/vt/vn notation is relevant for the gravity model
    const TemporaryObjFile objFile{
            "v 0.0 0.0 0.0\n"
            "v 1.0 0.0 0.0\n"
            "v 0.0 1.0 0.0\n"
            "f 0/1/1 1/2/2 2//3\n"};

    const auto &[actualVertices, actualFaces] = MeshReader::readObj(objFile.string());

    ASSERT_THAT(actualFaces, ContainerEq(std::vector<std::array<size_t, 3>>{{0, 1, 2}}));
}

TEST_F(MeshReaderTest, readObjThrowsOnMalformedInput) {
    using namespace testing;
    using namespace ::polyhedralGravity;
    const std::string validVertices{
            "v 0.0 0.0 0.0\n"
            "v 1.0 0.0 0.0\n"
            "v 0.0 1.0 0.0\n"};
    // A truncated face must not silently reuse a previously parsed index
    const TemporaryObjFile truncatedFace{validVertices + "f 0 1\n"};
    EXPECT_THROW(MeshReader::readObj(truncatedFace.string()), std::runtime_error);
    // A quadrilateral must not silently be truncated to a triangle
    const TemporaryObjFile quadFace{validVertices + "f 0 1 2 3\n"};
    EXPECT_THROW(MeshReader::readObj(quadFace.string()), std::runtime_error);
    // A negative index must not wrap around when converted to an unsigned type
    const TemporaryObjFile negativeIndex{validVertices + "f 0 1 -1\n"};
    EXPECT_THROW(MeshReader::readObj(negativeIndex.string()), std::runtime_error);
    // A non-numeric index must be rejected
    const TemporaryObjFile nonNumericIndex{validVertices + "f 0 1 two\n"};
    EXPECT_THROW(MeshReader::readObj(nonNumericIndex.string()), std::runtime_error);
    // A truncated vertex must not be padded with stale coordinates
    const TemporaryObjFile truncatedVertex{"v 0.0 0.0 0.0\nv 1.0 0.0\n"};
    EXPECT_THROW(MeshReader::readObj(truncatedVertex.string()), std::runtime_error);
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