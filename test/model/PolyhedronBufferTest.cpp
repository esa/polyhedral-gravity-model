#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cstdint>
#include <vector>

#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/util/KokkosSession.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"

/**
 * Checks that a Polyhedron can be built on top of a foreign buffer, i.e. on the memory of an array
 * library, and that doing so neither copies nor modifies that buffer.
 */
class PolyhedronBufferTest : public ::testing::Test {

protected:
    /** A unitless cube with an edge length of two, centered around the origin, as a flat buffer */
    std::vector<double> _vertices{
            -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0,
            -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0};

    /** The cube's twelve triangular faces, with OUTWARDS pointing plane unit normals */
    std::vector<size_t> _faces{1, 3, 2, 0, 3, 1, 0, 1, 5, 0, 5, 4, 0, 7, 3, 0, 4, 7,
                               1, 2, 6, 1, 6, 5, 2, 3, 6, 3, 7, 6, 4, 5, 6, 4, 6, 7};

    /** A computation point which is neither symmetric nor on the cube's surface */
    polyhedralGravity::Array3 _computationPoint{0.5, 0.25, 0.75};

    /**
     * Builds the cube from the flat buffers above.
     * @tparam VertexType the element type the vertices are handed over in
     * @tparam IndexType the element type the face indices are handed over in
     * @param vertices the vertex buffer
     * @param faces the face buffer
     * @param integrity the integrity check to run
     * @return the cube
     */
    template<typename VertexType, typename IndexType>
    static polyhedralGravity::Polyhedron makeCube(
            const std::vector<VertexType> &vertices, const std::vector<IndexType> &faces,
            const polyhedralGravity::PolyhedronIntegrity integrity = polyhedralGravity::PolyhedronIntegrity::DISABLE) {
        using namespace polyhedralGravity;
        return Polyhedron{vertices.data(),  vertices.size() / 3,
                          faces.data(),     faces.size() / 3,
                          1.0,              MemoryLocation::HOST,
                          NormalOrientation::OUTWARDS, integrity, MetricUnit::UNITLESS};
    }
};

/**
 * A polyhedron built from a buffer must describe the very same mesh as one built from vectors.
 */
TEST_F(PolyhedronBufferTest, MatchesTheVectorConstructor) {
    using namespace polyhedralGravity;
    const Polyhedron fromBuffer = makeCube(_vertices, _faces);
    const Polyhedron fromVectors{fromBuffer.getVertices(), fromBuffer.getFaces(), 1.0,
                                 NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE,
                                 MetricUnit::UNITLESS};

    ASSERT_EQ(fromBuffer.countVertices(), 8);
    ASSERT_EQ(fromBuffer.countFaces(), 12);
    ASSERT_THAT(fromBuffer.getVertices(), ::testing::ContainerEq(fromVectors.getVertices()));
    ASSERT_THAT(fromBuffer.getFaces(), ::testing::ContainerEq(fromVectors.getFaces()));

    const auto [expectedPotential, expectedAcceleration, expectedTensor] =
            GravityModel::evaluate(fromVectors, _computationPoint);
    const auto [potential, acceleration, tensor] = GravityModel::evaluate(fromBuffer, _computationPoint);
    EXPECT_NEAR(potential, expectedPotential, 1e-13);
    for (size_t index = 0; index < 3; ++index) {
        EXPECT_NEAR(acceleration[index], expectedAcceleration[index], 1e-13);
    }
}

/**
 * A buffer whose element type is the library's own is used where it is, i.e. it is not copied.
 */
TEST_F(PolyhedronBufferTest, MatchingElementTypesAreNotCopied) {
    const polyhedralGravity::Polyhedron cube = makeCube(_vertices, _faces);
    EXPECT_EQ(cube.getMesh().getHostMesh().vertices.data(), _vertices.data());
    EXPECT_EQ(cube.getMesh().getHostMesh().faces.data(), _faces.data());
}

/**
 * A buffer of any other element type is converted, which must not change the evaluation's result.
 */
TEST_F(PolyhedronBufferTest, ForeignElementTypesAreConverted) {
    using namespace polyhedralGravity;
    const std::vector<float> narrowVertices{_vertices.begin(), _vertices.end()};
    const std::vector<int32_t> narrowFaces{_faces.begin(), _faces.end()};

    const Polyhedron cube = makeCube(narrowVertices, narrowFaces);
    EXPECT_NE(static_cast<const void *>(cube.getMesh().getHostMesh().vertices.data()),
              static_cast<const void *>(narrowVertices.data()));
    ASSERT_THAT(cube.getFaces(), ::testing::ContainerEq(makeCube(_vertices, _faces).getFaces()));

    const auto [expectedPotential, expectedAcceleration, expectedTensor] =
            GravityModel::evaluate(makeCube(_vertices, _faces), _computationPoint);
    const auto [potential, acceleration, tensor] = GravityModel::evaluate(cube, _computationPoint);
    // The vertices went through single precision, so only their significant digits survived
    EXPECT_NEAR(potential, expectedPotential, 1e-6);
}

/**
 * Healing the vertex ordering allocates and must never write into the caller's buffer.
 */
TEST_F(PolyhedronBufferTest, HealingDoesNotModifyTheBuffer) {
    using namespace polyhedralGravity;
    std::vector<size_t> violatingFaces = _faces;
    std::swap(violatingFaces[0], violatingFaces[1]);
    const std::vector<size_t> unhealedFaces = violatingFaces;

    const Polyhedron healed = makeCube(_vertices, violatingFaces, PolyhedronIntegrity::HEAL);

    EXPECT_THAT(violatingFaces, ::testing::ContainerEq(unhealedFaces));
    EXPECT_THAT(healed.getFaces(), ::testing::ContainerEq(makeCube(_vertices, _faces).getFaces()));
}

/**
 * Asking for a device pointer on a build without a GPU backend has to fail with a clear message.
 */
TEST_F(PolyhedronBufferTest, DevicePointerNeedsAGpuBackend) {
    using namespace polyhedralGravity;
    if (kokkos::GPU_AVAILABLE) {
        GTEST_SKIP() << "This build has a GPU backend, so a device pointer is accepted";
    }
    EXPECT_THROW(
            (Polyhedron{_vertices.data(), _vertices.size() / 3, _faces.data(), _faces.size() / 3, 1.0,
                        MemoryLocation::DEVICE, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE,
                        MetricUnit::UNITLESS}),
            std::runtime_error);
}

/**
 * A null pointer is never a mesh.
 */
TEST_F(PolyhedronBufferTest, NullPointerIsRejected) {
    using namespace polyhedralGravity;
    EXPECT_THROW((Polyhedron{static_cast<const double *>(nullptr), 8, _faces.data(), 12, 1.0}),
                 std::invalid_argument);
}
