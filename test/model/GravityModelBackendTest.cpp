#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <tuple>
#include <vector>

#include "polyhedralGravity/model/KokkosSession.h"
#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"

/**
 * Skips the surrounding test if this build has no GPU backend.
 * This is a macro and not a helper function since GTEST_SKIP() returns from the test body.
 */
#define SKIP_WITHOUT_GPU()                                                                    \
    if (!polyhedralGravity::kokkos::GPU_AVAILABLE) {                                          \
        GTEST_SKIP() << "This build has no GPU backend, the enabled execution spaces are "     \
                     << polyhedralGravity::kokkos::getEnabledExecutionSpaces();                \
    }

/**
 * Checks that all three compute backends and both floating point precisions agree with each other
 * on a unitless cube, whose computation points cover every singularity case of Tsoulis' algorithm.
 */
class GravityModelBackendTest : public ::testing::Test {

protected:
    /** A unitless cube with an edge length of two, centered around the origin */
    polyhedralGravity::Polyhedron _cube{
            std::vector<polyhedralGravity::Array3>{
                    {-1.0, -1.0, -1.0},
                    {1.0, -1.0, -1.0},
                    {1.0, 1.0, -1.0},
                    {-1.0, 1.0, -1.0},
                    {-1.0, -1.0, 1.0},
                    {1.0, -1.0, 1.0},
                    {1.0, 1.0, 1.0},
                    {-1.0, 1.0, 1.0}},
            std::vector<polyhedralGravity::IndexArray3>{
                    {1, 3, 2}, {0, 3, 1}, {0, 1, 5}, {0, 5, 4}, {0, 7, 3}, {0, 4, 7},
                    {1, 2, 6}, {1, 6, 5}, {2, 3, 6}, {3, 7, 6}, {4, 5, 6}, {4, 6, 7}},
            1.0,
            polyhedralGravity::NormalOrientation::OUTWARDS,
            polyhedralGravity::PolyhedronIntegrity::DISABLE,
            polyhedralGravity::MetricUnit::UNITLESS};

    /**
     * Computation points covering the interesting positions relative to the cube: far outside, close to a
     * face, on a face, on an edge, on a vertex, and inside the polyhedron.
     */
    std::vector<polyhedralGravity::Array3> _computationPoints{
            {0.0, 0.0, 0.0},   // the center, i.e. P' is inside every face
            {0.5, 0.25, 0.75}, // inside, but not symmetric
            {1.0, 0.0, 0.0},   // on a face
            {1.0, 1.0, 0.0},   // on an edge
            {1.0, 1.0, 1.0},   // on a vertex
            {2.0, 0.0, 0.0},   // just outside a face
            {-3.0, 2.0, 1.5},  // outside
            {10.0, 10.0, 10.0},// far outside
            {0.0, 0.0, 1.0},   // on a face, at its center
            {-1.0, -1.0, 0.5}, // on an edge
    };

    /**
     * Asserts that two results are equal up to a given relative epsilon.
     * @param actual the result to check
     * @param expected the reference result
     * @param epsilon the allowed relative deviation
     * @param message added to every failure message
     */
    static void expectResultsNear(const polyhedralGravity::GravityModelResult &actual,
                                  const polyhedralGravity::GravityModelResult &expected, const double epsilon,
                                  const std::string &message) {
        const auto &[actualPotential, actualAcceleration, actualTensor] = actual;
        const auto &[expectedPotential, expectedAcceleration, expectedTensor] = expected;
        // The values of this cube are of order 1, so an absolute tolerance is a relative one here
        EXPECT_NEAR(actualPotential, expectedPotential, epsilon) << "The potential differed " << message;
        for (size_t index = 0; index < 3; ++index) {
            EXPECT_NEAR(actualAcceleration[index], expectedAcceleration[index], epsilon)
                    << "The acceleration's component " << index << " differed " << message;
        }
        for (size_t index = 0; index < 6; ++index) {
            EXPECT_NEAR(actualTensor[index], expectedTensor[index], epsilon)
                    << "The tensor's component " << index << " differed " << message;
        }
    }
};

/**
 * The Serial and the OpenMP backend run the identical kernel and only differ in how the faces are
 * distributed over the threads, so they may only disagree by the reassociation of the reduction.
 */
TEST_F(GravityModelBackendTest, CpuParallelMatchesCpuSerial) {
    using namespace polyhedralGravity;
    const GravityEvaluable evaluable{_cube};

    for (const auto &computationPoint: _computationPoints) {
        const auto serial = std::get<GravityModelResult>(evaluable(computationPoint, ComputeBackend::CPU_SERIAL));
        const auto parallel = std::get<GravityModelResult>(evaluable(computationPoint, ComputeBackend::CPU_PARALLEL));
        expectResultsNear(parallel, serial, 1e-13,
                          "between CPU_SERIAL and CPU_PARALLEL at the computation point [" +
                                  std::to_string(computationPoint[0]) + ", " + std::to_string(computationPoint[1]) +
                                  ", " + std::to_string(computationPoint[2]) + "]");
    }
}

/**
 * The GPU runs the identical kernel as the CPU, so it must agree with it up to the reassociation the
 * different reduction order causes.
 */
TEST_F(GravityModelBackendTest, GpuParallelMatchesCpu) {
    using namespace polyhedralGravity;
    SKIP_WITHOUT_GPU()
    const GravityEvaluable evaluable{_cube};

    for (const auto &computationPoint: _computationPoints) {
        const auto host = std::get<GravityModelResult>(evaluable(computationPoint, ComputeBackend::CPU_PARALLEL));
        const auto device = std::get<GravityModelResult>(evaluable(computationPoint, ComputeBackend::GPU_PARALLEL));
        expectResultsNear(device, host, 1e-10, "between CPU_PARALLEL and GPU_PARALLEL");
    }
}

/**
 * Requesting the GPU on a build without a GPU backend must fail loudly instead of silently computing
 * somewhere else.
 */
TEST_F(GravityModelBackendTest, GpuParallelThrowsWithoutGpuBackend) {
    using namespace polyhedralGravity;
    if (kokkos::GPU_AVAILABLE) {
        GTEST_SKIP() << "This build has a GPU backend, so requesting it must not throw";
    }
    const GravityEvaluable evaluable{_cube};
    EXPECT_THROW(std::ignore = evaluable(Array3{0.0, 0.0, 0.0}, ComputeBackend::GPU_PARALLEL), std::runtime_error);
    EXPECT_THROW(std::ignore = evaluable(_computationPoints, ComputeBackend::GPU_PARALLEL), std::runtime_error);
    EXPECT_THROW(std::ignore = GravityModel::evaluate(_cube, Array3{0.0, 0.0, 0.0}, ComputeBackend::GPU_PARALLEL),
                 std::runtime_error);
}

/**
 * Single precision cancels large terms against each other, so it only reproduces the double precision
 * result to a few significant digits. This test pins down that it is still in the right ballpark.
 */
TEST_F(GravityModelBackendTest, Float32ApproximatesFloat64) {
    using namespace polyhedralGravity;
    const GravityEvaluable float64{_cube, ComputePrecision::FLOAT64};
    const GravityEvaluable float32{_cube, ComputePrecision::FLOAT32};
    ASSERT_EQ(float32.getComputePrecision(), ComputePrecision::FLOAT32);
    ASSERT_EQ(float64.getComputePrecision(), ComputePrecision::FLOAT64);

    for (const auto &computationPoint: _computationPoints) {
        const auto expected = std::get<GravityModelResult>(float64(computationPoint));
        const auto actual = std::get<GravityModelResult>(float32(computationPoint));
        expectResultsNear(actual, expected, 1e-4, "between FLOAT32 and FLOAT64");
    }
}

/**
 * The multi point kernel uses a team policy instead of a range policy, so it is a genuinely different
 * code path which has to agree with the single point one.
 */
TEST_F(GravityModelBackendTest, MultiPointMatchesSinglePoint) {
    using namespace polyhedralGravity;
    const GravityEvaluable evaluable{_cube};

    for (const auto backend: {ComputeBackend::CPU_SERIAL, ComputeBackend::CPU_PARALLEL}) {
        const auto actual = std::get<std::vector<GravityModelResult>>(evaluable(_computationPoints, backend));
        ASSERT_EQ(actual.size(), _computationPoints.size());
        for (size_t index = 0; index < _computationPoints.size(); ++index) {
            const auto expected = std::get<GravityModelResult>(evaluable(_computationPoints[index], backend));
            expectResultsNear(actual[index], expected, 1e-13, "between the multi and the single point evaluation");
        }
    }
}

/**
 * Restoring a GravityEvaluable from its state, which is what unpickling does, must yield an evaluable
 * computing exactly the same values.
 */
TEST_F(GravityModelBackendTest, RestoredFromStateMatchesOriginal) {
    using namespace polyhedralGravity;
    const GravityEvaluable original{_cube};
    const auto &[polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals] = original.getState();
    const GravityEvaluable restored{polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals};

    for (const auto &computationPoint: _computationPoints) {
        const auto expected = std::get<GravityModelResult>(original(computationPoint));
        const auto actual = std::get<GravityModelResult>(restored(computationPoint));
        expectResultsNear(actual, expected, 1e-15, "between the original and the restored GravityEvaluable");
    }
}

/**
 * The free evaluate function is a thin wrapper around a GravityEvaluable and must not change the result.
 */
TEST_F(GravityModelBackendTest, FreeFunctionMatchesEvaluable) {
    using namespace polyhedralGravity;
    const GravityEvaluable evaluable{_cube};

    for (const auto &computationPoint: _computationPoints) {
        const auto expected = std::get<GravityModelResult>(evaluable(computationPoint, ComputeBackend::CPU_SERIAL));
        const auto actual = GravityModel::evaluate(_cube, computationPoint, ComputeBackend::CPU_SERIAL);
        expectResultsNear(actual, expected, 1e-15, "between GravityModel::evaluate and GravityEvaluable");
    }
}

/**
 * A polyhedron without faces cannot be evaluated and must be rejected while constructing the evaluable
 * rather than producing a zero result.
 */
TEST_F(GravityModelBackendTest, EmptyPolyhedronIsRejected) {
    using namespace polyhedralGravity;
    const Polyhedron empty{std::vector<Array3>{}, std::vector<IndexArray3>{}, 1.0, NormalOrientation::OUTWARDS,
                           PolyhedronIntegrity::DISABLE, MetricUnit::UNITLESS};
    EXPECT_THROW(GravityEvaluable{empty}, std::invalid_argument);
}
