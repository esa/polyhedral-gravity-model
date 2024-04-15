#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <string>
#include "polyhedralGravity/input/TetgenAdapter.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/model/Polyhedron.h"

#include "GravityModelVectorUtility.h"


/**
 * Contains Tests based on the Eros mesh taken from
 * https://github.com/darioizzo/geodesyNets/tree/master/3dmeshes (last accessed: 07.04.2022)
 *
 * The values which are used to check the C++ implementation are calculated by the
 * Tsoulis reference implementation in FORTRAN and saved in the files in test/resources.
 *
 */
class GravityModelBigTest : public ::testing::Test {

protected:

    /**
     * Relative big epsilon due to deviations between FORTRAN implementation and C++ implementation
     */
    static constexpr double LOCAL_TEST_EPSILON = 1e-6;

    static constexpr size_t LOCAL_TEST_COUNT_FACES = 14744;
    static constexpr size_t LOCAL_TEST_COUNT_NODES_PER_FACE = 3;

    polyhedralGravity::Polyhedron _polyhedron{
            polyhedralGravity::TetgenAdapter{
                    {"resources/GravityModelBigTest.node", "resources/GravityModelBigTest.face"}}.getPolyhedron()};

    std::array<double, 3> _computationPoint{0.0, 0.0, 0.0};

    std::vector<std::array<std::array<double, 3>, 3>> expectedGij;

    std::vector<std::array<double, 3>> expectedPlaneUnitNormals;

    std::vector<std::array<std::array<double, 3>, 3>> expectedSegmentUnitNormals;

    std::vector<double> expectedPlaneNormalOrientations;

    std::vector<polyhedralGravity::HessianPlane> expectedHessianPlanes;

    std::vector<double> expectedPlaneDistances;

    std::vector<std::array<double, 3>> expectedOrthogonalProjectionPointsOnPlane;

    std::vector<std::array<double, 3>> expectedSegmentNormalOrientations;

    std::vector<std::array<std::array<double, 3>, 3>> expectedOrthogonalProjectionPointsOnSegment;

    std::vector<std::array<double, 3>> expectedSegmentDistances;

    std::vector<std::array<polyhedralGravity::Distance, 3>> expectedDistancesPerSegmentEndpoint;

    std::vector<std::array<polyhedralGravity::TranscendentalExpression, 3>> expectedTranscendentalExpressions;


    std::vector<std::pair<double, std::array<double, 3>>> expectedSingularityTerms;

public:

    [[nodiscard]] std::vector<std::array<std::array<double, 3>, 3>>
    readTwoDimensionalCartesian(const std::string &filename) const {
        std::vector<std::array<std::array<double, 3>, 3>> result{LOCAL_TEST_COUNT_FACES};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x, y, z;
            if (!(linestream >> x >> y >> z)) {
                break;
            }
            result[i / 3][i % 3] = std::array<double, 3>{x, y, z};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<double, 3>> readOneDimensionalCartesian(const std::string &filename) const {
        std::vector<std::array<double, 3>> result{LOCAL_TEST_COUNT_FACES};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x, y, z;
            if (!(linestream >> x >> y >> z)) {
                break;
            }
            result[i] = std::array<double, 3>{x, y, z};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<double, 3>> readTwoDimensionalValue(const std::string &filename) const {
        std::vector<std::array<double, 3>> result{LOCAL_TEST_COUNT_FACES};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x;
            if (!(linestream >> x)) {
                break;
            }
            result[i / 3][i % 3] = x;
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<double> readOneDimensionalValue(const std::string &filename) const {
        std::vector<double> result(LOCAL_TEST_COUNT_FACES, 0.0);
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x;
            if (!(linestream >> x)) {
                break;
            }
            result[i] = x;
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<polyhedralGravity::HessianPlane>
    readHessianPlanes(const std::string &filename) const {
        std::vector<polyhedralGravity::HessianPlane> result{LOCAL_TEST_COUNT_FACES};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double a, b, c, d;
            if (!(linestream >> a >> b >> c >> d)) {
                break;
            }
            result[i] = polyhedralGravity::HessianPlane{a, b, c, d};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<polyhedralGravity::Distance, 3>>
    readDistances(const std::string &filename) const {
        std::vector<std::array<polyhedralGravity::Distance, 3>> result{LOCAL_TEST_COUNT_FACES};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double l1, l2, s1, s2;
            if (!(linestream >> l1 >> l2 >> s1 >> s2)) {
                break;
            }
            result[i / 3][i % 3] = polyhedralGravity::Distance{l1, l2, s1, s2};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<polyhedralGravity::TranscendentalExpression, 3>>
    readTranscendentalExpressions(const std::string &filename) const {
        std::vector<std::array<polyhedralGravity::TranscendentalExpression, 3>> result{LOCAL_TEST_COUNT_FACES};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double ln, an;
            if (!(linestream >> ln >> an)) {
                break;
            }
            result[i / 3][i % 3] = polyhedralGravity::TranscendentalExpression{ln, an};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<double, 3>>
    readBetaSingularities(const std::string &filename) const {
        std::vector<std::array<double, 3>> result{LOCAL_TEST_COUNT_FACES};
        std::ifstream infile(filename);
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            int i, j;
            double sng;
            if (!(linestream >> i >> j >> sng)) {
                break;
            }
            result[i - 1][j - 1] = sng;
        }
        return result;
    }

    GravityModelBigTest() : ::testing::Test() {
        using namespace polyhedralGravity;
        using namespace util;

        expectedGij = readTwoDimensionalCartesian("resources/GravityModelBigTestExpectedGij.txt");
        expectedPlaneUnitNormals =
                readOneDimensionalCartesian("resources/GravityModelBigTestExpectedPlaneUnitNormals.txt");
        expectedSegmentUnitNormals =
                readTwoDimensionalCartesian("resources/GravityModelBigTestExpectedSegmentUnitNormals.txt");
        expectedPlaneNormalOrientations =
                readOneDimensionalValue("resources/GravityModelBigTestExpectedPlaneOrientation.txt");
        expectedHessianPlanes =
                readHessianPlanes("resources/GravityModelBigTestExpectedHessianPlanes.txt");
        expectedPlaneDistances =
                readOneDimensionalValue("resources/GravityModelBigTestExpectedPlaneDistances.txt");
        expectedOrthogonalProjectionPointsOnPlane =
                readOneDimensionalCartesian("resources/GravityModelBigTestExpectedOrthogonalPlaneProjectionPoints.txt");
        expectedSegmentNormalOrientations =
                readTwoDimensionalValue("resources/GravityModelBigTestExpectedSegmentOrientation.txt");
        expectedOrthogonalProjectionPointsOnSegment = readTwoDimensionalCartesian(
                "resources/GravityModelBigTestExpectedOrthogonalSegmentProjectionPoints.txt");
        expectedSegmentDistances =
                readTwoDimensionalValue("resources/GravityModelBigTestExpectedSegmentDistances.txt");
        expectedDistancesPerSegmentEndpoint =
                readDistances("resources/GravityModelBigTestExpectedDistances.txt");
        expectedTranscendentalExpressions =
                readTranscendentalExpressions("resources/GravityModelBigTestExpectedTranscendentalExpressions.txt");
        auto expectedAlphaSingularityTerms =
                readOneDimensionalValue("resources/GravityModelBigTestExpectedAlphaSingularities.txt");
        auto expectedBetaSingularityTerms =
                readBetaSingularities("resources/GravityModelBigTestExpectedBetaSingularities.txt");

        expectedSingularityTerms.resize(expectedAlphaSingularityTerms.size());
        for (int i = 0; i < expectedAlphaSingularityTerms.size(); ++i) {
            expectedSingularityTerms[i] =
                    std::make_pair(expectedAlphaSingularityTerms[i], expectedBetaSingularityTerms[i]);
        }
    }

};

TEST_F(GravityModelBigTest, GijVectors) {
    using namespace testing;

    auto actualGij = polyhedralGravity::GravityModel::calculateSegmentVectors(_polyhedron);

    ASSERT_THAT(actualGij, ContainerEq(expectedGij));
}

TEST_F(GravityModelBigTest, PlaneUnitNormals) {
    using namespace testing;

    auto actualPlaneUnitNormals = polyhedralGravity::GravityModel::calculatePlaneUnitNormals(expectedGij);

    ASSERT_THAT(actualPlaneUnitNormals, ContainerEq(expectedPlaneUnitNormals));
}

TEST_F(GravityModelBigTest, SegmentUnitNormals) {
    using namespace testing;

    auto actualSegmentUnitNormals = polyhedralGravity::GravityModel::calculateSegmentUnitNormals(expectedGij,
                                                                                                 expectedPlaneUnitNormals);

    ASSERT_THAT(actualSegmentUnitNormals, ContainerEq(expectedSegmentUnitNormals));
}

TEST_F(GravityModelBigTest, PlaneNormalOrientations) {
    using namespace testing;

    auto actualPlaneNormalOrientations = polyhedralGravity::GravityModel::calculatePlaneNormalOrientations(
            _computationPoint, _polyhedron, expectedPlaneUnitNormals);

    ASSERT_THAT(actualPlaneNormalOrientations, ContainerEq(expectedPlaneNormalOrientations));
}

TEST_F(GravityModelBigTest, HessianPlane) {
    using namespace testing;
    using namespace polyhedralGravity;

    auto actualHessianPlane =
            polyhedralGravity::GravityModel::calculateFacesToHessianPlanes(_computationPoint, _polyhedron);

    ASSERT_EQ(actualHessianPlane, expectedHessianPlanes);
}

TEST_F(GravityModelBigTest, PlaneDistances) {
    using namespace testing;

    auto actualPlaneDistances = polyhedralGravity::GravityModel::calculatePlaneDistances(expectedHessianPlanes);

    ASSERT_THAT(actualPlaneDistances, ContainerEq(expectedPlaneDistances));
}

TEST_F(GravityModelBigTest, OrthogonalProjectionPointsOnPlane) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnPlane = polyhedralGravity::GravityModel::calculateOrthogonalProjectionPointsOnPlane(
            expectedHessianPlanes, expectedPlaneUnitNormals, expectedPlaneDistances);

    for (size_t i = 0; i < actualOrthogonalProjectionPointsOnPlane.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(
                    actualOrthogonalProjectionPointsOnPlane[i][j],
                    expectedOrthogonalProjectionPointsOnPlane[i][j])
                                << "Difference for P' of plane=" << i << " and coordinate-Nr.=" << j;
        }
    }
}

TEST_F(GravityModelBigTest, SegmentNormalOrientations) {
    using namespace testing;

    auto actualSegmentNormalOrientations =
            polyhedralGravity::GravityModel::calculateSegmentNormalOrientations(_computationPoint, _polyhedron,
                                                                                expectedSegmentUnitNormals,
                                                                                expectedOrthogonalProjectionPointsOnPlane);

    ASSERT_THAT(actualSegmentNormalOrientations, ContainerEq(expectedSegmentNormalOrientations));
}

TEST_F(GravityModelBigTest, OrthogonalProjectionPointsOnSegment) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnSegment =
            polyhedralGravity::GravityModel::calculateOrthogonalProjectionPointsOnSegments(
                    _computationPoint, _polyhedron,
                    expectedOrthogonalProjectionPointsOnPlane,
                    expectedSegmentNormalOrientations);

    for (size_t i = 0; i < actualOrthogonalProjectionPointsOnSegment.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                EXPECT_NEAR(
                        actualOrthogonalProjectionPointsOnSegment[i][j][k],
                        expectedOrthogonalProjectionPointsOnSegment[i][j][k], 1e16)
                                    << "Difference for P'' of segment=(" << i << ", " << j << ") and coordinate-Nr."
                                    << k;
            }
        }
    }
}

TEST_F(GravityModelBigTest, SegmentDistances) {
    using namespace testing;

    auto actualSegmentDistances =
            polyhedralGravity::GravityModel::calculateSegmentDistances(expectedOrthogonalProjectionPointsOnPlane,
                                                                       expectedOrthogonalProjectionPointsOnSegment);

    ASSERT_THAT(actualSegmentDistances, ContainerEq(expectedSegmentDistances));
}

TEST_F(GravityModelBigTest, DistancesPerSegmentEndpoint) {
    using namespace testing;

    auto actualDistancesPerSegmentEndpoint =
            polyhedralGravity::GravityModel::calculateDistances(_computationPoint, _polyhedron, expectedGij,
                                                                expectedOrthogonalProjectionPointsOnSegment);

    ASSERT_THAT(actualDistancesPerSegmentEndpoint, ContainerEq(expectedDistancesPerSegmentEndpoint));
}

TEST_F(GravityModelBigTest, TranscendentalExpressions) {
    using namespace testing;

    auto actualTranscendentalExpressions =
            polyhedralGravity::GravityModel::calculateTranscendentalExpressions(_computationPoint, _polyhedron,
                                                                                expectedDistancesPerSegmentEndpoint,
                                                                                expectedPlaneDistances,
                                                                                expectedSegmentDistances,
                                                                                expectedSegmentNormalOrientations,
                                                                                expectedOrthogonalProjectionPointsOnPlane);

    ASSERT_EQ(actualTranscendentalExpressions.size(), expectedTranscendentalExpressions.size());

    // For arrays one could use the something like Pointwise(DoubleEqual(), expected_array) for the matcher
    // here, we have no arrays but a custom data structure (the above is just a hint for the future, to safe time)
    for (size_t i = 0; i < actualTranscendentalExpressions.size(); ++i) {
        for (size_t j = 0; j < actualTranscendentalExpressions[i].size(); ++j) {
            ASSERT_NEAR(actualTranscendentalExpressions[i][j].ln,
                        expectedTranscendentalExpressions[i][j].ln, LOCAL_TEST_EPSILON)
                                        << "The LN value differed for transcendental term (i,j) = (" << i << ',' << j
                                        << ')';
            ASSERT_NEAR(actualTranscendentalExpressions[i][j].an,
                        expectedTranscendentalExpressions[i][j].an, LOCAL_TEST_EPSILON)
                                        << "The AN value differed for transcendental term (i,j) = (" << i << ',' << j
                                        << ')';
        }
    }
}

TEST_F(GravityModelBigTest, SingularityTerms) {
    using namespace testing;
    using namespace polyhedralGravity;

    auto actualSingularityTerms =
            polyhedralGravity::GravityModel::calculateSingularityTerms(_computationPoint, _polyhedron, expectedGij,
                                                                       expectedSegmentNormalOrientations,
                                                                       expectedOrthogonalProjectionPointsOnPlane,
                                                                       expectedPlaneDistances,
                                                                       expectedPlaneNormalOrientations,
                                                                       expectedPlaneUnitNormals);

    ASSERT_EQ(actualSingularityTerms.size(), actualSingularityTerms.size());

    for (size_t i = 0; i < actualSingularityTerms.size(); ++i) {
        EXPECT_NEAR(actualSingularityTerms[i].first,
                    expectedSingularityTerms[i].first, LOCAL_TEST_EPSILON)
                            << "The sing A value differed for singularity term (i) = (" << i << ')';
        EXPECT_THAT(actualSingularityTerms[i].second,
                    Pointwise(DoubleNear(LOCAL_TEST_EPSILON), expectedSingularityTerms[i].second))
                            << "The sing B value differed for singularity term (i) = (" << i << ')';
    }
}
