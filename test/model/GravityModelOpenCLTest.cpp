#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include <array>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

/**
 * Verifies that the OpenCL compute backend agrees with the CPU backend, which the remaining test
 * suite pins against Tsoulis' reference values.
 *
 * The tests skip themselves where the requested backend is unavailable -- a library compiled without
 * OpenCL, a machine without an OpenCL device, or (for FLOAT64) a device without cl_khr_fp64, which
 * notably includes every Apple Silicon GPU.
 */

/**
 * Skips the running test if the given evaluable did not end up on the OpenCL backend.
 * This cannot be a helper function since GTEST_SKIP() returns from the function it appears in.
 */
#define SKIP_UNLESS_OPENCL(evaluable)                                                            \
    if ((evaluable).getComputeBackend() != polyhedralGravity::ComputeBackend::OPENCL) {          \
        GTEST_SKIP() << "No OpenCL device supporting the requested precision is available.";     \
    }                                                                                            \
    static_assert(true, "swallow the trailing semicolon")

class GravityModelOpenCLTest : public ::testing::Test {

protected:
    /**
     * A unit cube. Being unitless, the results are not scaled by the gravitational constant, which
     * keeps the compared magnitudes in a range where a relative comparison is meaningful.
     */
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
     * Computation points covering the ordinary case as well as every singularity case of
     * Tsoulis' algorithm, since those are the ones the device kernel reimplements.
     */
    static std::vector<polyhedralGravity::Array3> computationPoints() {
        return {
                {0.0, 0.0, 0.0},   // inside the cube
                {2.0, 0.0, 0.0},   // outside, facing a face
                {5.0, 5.0, 5.0},   // outside, facing a vertex
                {-3.0, 1.5, 0.25}, // outside, arbitrary
                {1.0, 0.0, 0.0},   // 1. case: P' lies inside a face
                {1.0, 0.5, 0.0},   // 1. case: P' lies inside a face, off-centre
                {1.0, 1.0, 0.0},   // 2. case: P' lies on an edge, but not on a vertex
                {0.0, 1.0, 1.0},   // 2. case: P' lies on an edge, but not on a vertex
                {1.0, 1.0, 1.0},   // 3. case: P' lies on a vertex
                {-1.0, -1.0, -1.0} // 3. case: P' lies on a vertex
        };
    }

    /**
     * Asserts that two results agree relative to their magnitude.
     * @param actual the result of the OpenCL backend
     * @param expected the result of the CPU backend
     * @param epsilon the tolerated relative deviation
     * @param context a description of the computation point, printed on failure
     */
    static void expectResultsNear(const polyhedralGravity::GravityModelResult &actual,
                                  const polyhedralGravity::GravityModelResult &expected,
                                  const double epsilon, const std::string &context) {
        const auto &[actualPotential, actualAcceleration, actualTensor] = actual;
        const auto &[expectedPotential, expectedAcceleration, expectedTensor] = expected;

        // A relative comparison is meaningless around zero, so the tolerance also has an absolute floor
        const auto near = [epsilon](const double lhs, const double rhs) {
            return std::abs(lhs - rhs) <= epsilon * std::max({std::abs(lhs), std::abs(rhs), 1.0});
        };

        EXPECT_TRUE(near(actualPotential, expectedPotential))
                << context << ": potential " << actualPotential << " != " << expectedPotential;
        for (size_t index = 0; index < actualAcceleration.size(); ++index) {
            EXPECT_TRUE(near(actualAcceleration[index], expectedAcceleration[index]))
                    << context << ": acceleration[" << index << "] " << actualAcceleration[index]
                    << " != " << expectedAcceleration[index];
        }
        for (size_t index = 0; index < actualTensor.size(); ++index) {
            EXPECT_TRUE(near(actualTensor[index], expectedTensor[index]))
                    << context << ": tensor[" << index << "] " << actualTensor[index]
                    << " != " << expectedTensor[index];
        }
    }
};

/**
 * The double precision OpenCL backend has to reproduce the CPU backend closely, since both evaluate
 * the very same algorithm in the very same precision and only differ in the summation order.
 */
TEST_F(GravityModelOpenCLTest, Float64MatchesCpuBackend) {
    using namespace polyhedralGravity;
    constexpr double EPSILON = 1e-10;

    const GravityEvaluable openCLEvaluable{_cube, ComputeBackend::OPENCL, ComputePrecision::FLOAT64};
    SKIP_UNLESS_OPENCL(openCLEvaluable);
    const GravityEvaluable cpuEvaluable{_cube, ComputeBackend::CPU};

    for (const Array3 &point: computationPoints()) {
        const auto actual = std::get<GravityModelResult>(openCLEvaluable(point));
        const auto expected = std::get<GravityModelResult>(cpuEvaluable(point));
        expectResultsNear(actual, expected, EPSILON,
                          "P = [" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ", " +
                                  std::to_string(point[2]) + "]");
    }
}

/**
 * The single precision backend follows the same algorithm but accumulates considerably more rounding
 * error, so it is only held to a correspondingly looser tolerance.
 */
TEST_F(GravityModelOpenCLTest, Float32MatchesCpuBackendApproximately) {
    using namespace polyhedralGravity;
    constexpr double EPSILON = 1e-4;

    const GravityEvaluable openCLEvaluable{_cube, ComputeBackend::OPENCL, ComputePrecision::FLOAT32};
    SKIP_UNLESS_OPENCL(openCLEvaluable);
    const GravityEvaluable cpuEvaluable{_cube, ComputeBackend::CPU};

    for (const Array3 &point: computationPoints()) {
        const auto actual = std::get<GravityModelResult>(openCLEvaluable(point));
        const auto expected = std::get<GravityModelResult>(cpuEvaluable(point));
        expectResultsNear(actual, expected, EPSILON,
                          "P = [" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ", " +
                                  std::to_string(point[2]) + "]");
    }
}

/**
 * Evaluating multiple points at once has to give the same results as evaluating them one by one.
 */
TEST_F(GravityModelOpenCLTest, MultiPointEvaluationMatchesSinglePoint) {
    using namespace polyhedralGravity;

    const GravityEvaluable evaluable{_cube, ComputeBackend::OPENCL, ComputePrecision::FLOAT32};
    SKIP_UNLESS_OPENCL(evaluable);
    const std::vector<Array3> points = computationPoints();

    const auto batched = std::get<std::vector<GravityModelResult>>(evaluable(points));
    ASSERT_EQ(batched.size(), points.size());
    for (size_t index = 0; index < points.size(); ++index) {
        const auto single = std::get<GravityModelResult>(evaluable(points[index]));
        expectResultsNear(batched[index], single, 0.0, "batched point " + std::to_string(index));
    }
}

/**
 * The polyhedron's density, normal orientation, and mesh unit have to be honoured by the device-side
 * prefix just as they are by the host.
 */
TEST_F(GravityModelOpenCLTest, ScalingMatchesCpuBackend) {
    using namespace polyhedralGravity;
    constexpr double EPSILON = 1e-4;

    const std::vector<Array3> vertices = _cube.getVertices();
    const std::vector<IndexArray3> faces = _cube.getFaces();

    for (const MetricUnit unit: {MetricUnit::METER, MetricUnit::KILOMETER, MetricUnit::UNITLESS}) {
        const Polyhedron polyhedron{vertices, faces, 2670.0, NormalOrientation::OUTWARDS,
                                    PolyhedronIntegrity::DISABLE, unit};

        const GravityEvaluable openCLEvaluable{polyhedron, ComputeBackend::OPENCL, ComputePrecision::FLOAT32};
        SKIP_UNLESS_OPENCL(openCLEvaluable);
        const GravityEvaluable cpuEvaluable{polyhedron, ComputeBackend::CPU};

        const Array3 point{-3.0, 1.5, 0.25};
        expectResultsNear(std::get<GravityModelResult>(openCLEvaluable(point)),
                          std::get<GravityModelResult>(cpuEvaluable(point)), EPSILON,
                          "scaling with a non-default mesh unit");
    }
}

/**
 * Requesting a backend which cannot be provided must fall back to the CPU rather than fail, which is
 * what makes ComputeBackend::OPENCL usable as the default backend.
 */
TEST_F(GravityModelOpenCLTest, ReportsTheBackendActuallyInUse) {
    using namespace polyhedralGravity;

    const GravityEvaluable cpuEvaluable{_cube, ComputeBackend::CPU};
    EXPECT_EQ(cpuEvaluable.getComputeBackend(), ComputeBackend::CPU);

    // Whether this ends up on OpenCL or on the CPU depends on the machine, but it must never throw
    // and it must report which of the two it settled on
    const GravityEvaluable defaultEvaluable{_cube};
    EXPECT_THAT(defaultEvaluable.getComputeBackend(),
                ::testing::AnyOf(ComputeBackend::CPU, ComputeBackend::OPENCL));

    // The state has to be readable regardless of the backend, since pickling relies on it
    const auto &[polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals] = defaultEvaluable.getState();
    EXPECT_EQ(segmentVectors.size(), _cube.countFaces());
    EXPECT_EQ(planeUnitNormals.size(), _cube.countFaces());
    EXPECT_EQ(segmentUnitNormals.size(), _cube.countFaces());
}
