#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <array>
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/util/UtilityConstants.h"

class PolyhedronMetricTest : public ::testing::Test {
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

    const std::vector<polyhedralGravity::IndexArray3> _cubeFaces{
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

    const polyhedralGravity::Array3 computationPoint{1.0, 1.0, 0.0};
};


TEST_F(PolyhedronMetricTest, MetricUnitConversion) {
    using namespace testing;
    using namespace polyhedralGravity;
    // All normals are pointing outwards, extensive Eros example
    const Polyhedron meterPolyhedron{_cubeVertices, _cubeFaces, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, MetricUnit::METER};
    const Polyhedron kilometerPolyhedron{_cubeVertices, _cubeFaces, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, MetricUnit::KILOMETER};
    const Polyhedron unitlessPolyhedron{_cubeVertices, _cubeFaces, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, MetricUnit::UNITLESS};

    const auto& [meterPot, meterAcc, meterTensor] = GravityModel::evaluate(meterPolyhedron, computationPoint);
    const auto& [kilometerPot, kilometerAcc, kilometerTensor] = GravityModel::evaluate(kilometerPolyhedron, computationPoint);
    const auto& [unitlessPot, unitlessAcc, unitlessTensor] = GravityModel::evaluate(unitlessPolyhedron, computationPoint);

    ASSERT_DOUBLE_EQ(meterPot * 1e-9, kilometerPot);
    ASSERT_DOUBLE_EQ(meterPot, unitlessPot * util::GRAVITATIONAL_CONSTANT);
    ASSERT_DOUBLE_EQ(kilometerPot, unitlessPot * util::GRAVITATIONAL_CONSTANT * 1e-9);
}

TEST_F(PolyhedronMetricTest, MetricUnit) {
    using namespace testing;
    using namespace polyhedralGravity;
    // All normals are pointing outwards, extensive Eros example
    const Polyhedron meterPolyhedron{_cubeVertices, _cubeFaces, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, MetricUnit::METER};
    ASSERT_EQ(meterPolyhedron.getDensityUnit(), "kg/m^3");
    ASSERT_EQ(meterPolyhedron.getMeshUnitAsString(), "m");

    const Polyhedron kilometerPolyhedron{_cubeVertices, _cubeFaces, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, MetricUnit::KILOMETER};
    ASSERT_EQ(kilometerPolyhedron.getDensityUnit(), "kg/km^3");
    ASSERT_EQ(kilometerPolyhedron.getMeshUnitAsString(), "km");

    const Polyhedron unitlessPolyhedron{_cubeVertices, _cubeFaces, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, MetricUnit::UNITLESS};
    ASSERT_EQ(unitlessPolyhedron.getDensityUnit(), "unitless");
    ASSERT_EQ(unitlessPolyhedron.getMeshUnitAsString(), "unitless");
}
