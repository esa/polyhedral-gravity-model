#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <tuple>
#include <fstream>
#include <sstream>
#include <string>
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/model/Polyhedron.h"


/**
 * Contains Tests how the calculation handles a cubic polyhedron
 */
class GravityModelCubeTest :
        public ::testing::TestWithParam<std::tuple<std::array<double, 3>, double, double, std::array<double, 3>>> {

protected:

    /**
     * Small epsilon since we compare to ground truth
     */
    static constexpr double LOCAL_TEST_EPSILON = 1e-20;
    const polyhedralGravity::Polyhedron _cube{{
                                                      {-1.0, -1.0, -1.0},
                                                      {1.0, -1.0, -1.0},
                                                      {1.0, 1.0, -1.0},
                                                      {-1.0, 1.0, -1.0},
                                                      {-1.0, -1.0, 1.0},
                                                      {1.0, -1.0, 1.0},
                                                      {1.0, 1.0, 1.0},
                                                      {-1.0, 1.0, 1.0}},
                                              {
                                                      {1,    3,    2},
                                                      {0,   3,    1},
                                                      {0,   1,   5},
                                                      {0,    5,   4},
                                                      {0,    7,    3},
                                                      {0,   4,    7},
                                                      {1,   2,   6},
                                                      {1,    6,   5},
                                                      {2, 3, 6},
                                                      {3, 7, 6},
                                                      {4, 5, 6},
                                                      {4, 6, 7}}
    };

public:

    [[nodiscard]] static std::vector<std::tuple<std::array<double, 3>, double, double, std::array<double, 3>>>
    readCubePoints(const std::string &filename) {
        std::vector<std::tuple<std::array<double, 3>, double, double, std::array<double, 3>>> result{};
        std::ifstream infile(filename);
        std::string line;

        // The first line contains only one value: The density
        std::getline(infile, line);
        std::istringstream firstLineStream(line);
        double density{1.0};
        firstLineStream >> density;

        // The other lines contain tuples of point, potential and acceleration
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double p1, p2, p3, potential, acc1, acc2, acc3;
            if (!(linestream >> p1 >> p2 >> p3 >> potential >> acc1 >> acc2 >> acc3)) {
                break;
            }
            result.emplace_back(std::array<double, 3>{p1, p2, p3}, density, potential, std::array<double, 3>{acc1, acc2, acc3});
        }
        return result;
    }




};

TEST_P(GravityModelCubeTest, CubePoints) {
    using namespace testing;
    using namespace polyhedralGravity;

    std::tuple<std::array<double, 3>, double, double, std::array<double, 3>> testData = GetParam();
    std::array<double, 3> computationPoint = std::get<0>(testData);
    double density = std::get<1>(testData);
    double expectedPotential = std::get<2>(testData);
    std::array<double, 3> expectedAcceleration = std::get<3>(testData);

    const auto[actualPotential, actualAcceleration, x] = GravityModel::evaluate(_cube, density, computationPoint);

    ASSERT_NEAR(actualPotential, expectedPotential, LOCAL_TEST_EPSILON);
    ASSERT_THAT(actualAcceleration, Pointwise(DoubleNear(LOCAL_TEST_EPSILON), expectedAcceleration));
}

INSTANTIATE_TEST_SUITE_P(CubeGravityModelTest01, GravityModelCubeTest,
                         ::testing::ValuesIn(
                                 GravityModelCubeTest::readCubePoints("resources/analytic_cube_solution_density1.txt")));

INSTANTIATE_TEST_SUITE_P(CubeGravityModelTest42, GravityModelCubeTest,
                         ::testing::ValuesIn(
                                 GravityModelCubeTest::readCubePoints("resources/analytic_cube_solution_density42.txt")));
