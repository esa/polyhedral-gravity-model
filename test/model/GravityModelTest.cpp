#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <utility>
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/model/Polyhedron.h"

#include "GravityModelVectorUtility.h"

/**
 * Contains Tests based on the example from Tsoulis FORTRAN implementation.
 * Hardcoded values taken from his implementation's results.
 */
class GravityModelTest : public ::testing::Test {

protected:

    /**
     * Relative big epsilon due to deviations between FORTRAN implementation and C++ implementation
     */
    static constexpr double LOCAL_TEST_EPSILON = 1e-6;

    //New polyhedron with given vertices and faces
    //this is the base example from Tsoulis
    polyhedralGravity::Polyhedron _polyhedron{
            {
                    {-20, 0, 25},
                    {0, 0, 25},
                    {0, 10, 25},
                    {-20, 10, 25},
                    {-20, 0, 15},
                    {0, 0, 15},
                    {0, 10, 15},
                    {-20, 10, 15}
            },
            {
                    {0,   1, 3},
                    {1, 2, 3},
                    {0, 4,  5},
                    {0,   5,  1},
                    {0,   7, 4},
                    {0, 3, 7},
                    {1, 5,  6},
                    {1,   6,  2},
                    {3, 6, 7},
                    {2, 6, 3},
                    {4, 6, 5},
                    {4, 7, 6}
            }
    };

    std::array<double, 3> _computationPoint{0.0, 0.0, 0.0};

    std::vector<std::array<std::array<double, 3>, 3>> expectedGij{
            std::array<std::array<double, 3>, 3>{{{20.0, 0.0, 0.0}, {-20.0, 10.0, 0.0}, {0.0, -10.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {-20.0, 0.0, 0.0}, {20.0, -10.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, -10.0}, {20.0, 0.0, 0.0}, {-20.0, 0.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{20.0, 0.0, -10.0}, {0.0, 0.0, 10.0}, {-20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, -10.0}, {0.0, -10.0, 0.0}, {0.0, 0.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {0.0, 0.0, -10.0}, {0.0, -10.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, -10.0}, {0.0, 10.0, 0.0}, {0.0, -10.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, -10.0}, {0.0, 0.0, 10.0}, {0.0, -10.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{20.0, 0.0, -10.0}, {-20.0, 0.0, 0.0}, {0.0, 0.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, -10.0}, {-20.0, 0.0, 10.0}, {20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{20.0, 10.0, 0.0}, {0.0, -10.0, 0.0}, {-20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {20.0, 0.0, 0.0}, {-20.0, -10.0, 0.0}}}
    };

    std::vector<std::array<double, 3>> expectedPlaneUnitNormals{
            {0.0,  -0.0, 1.0},
            {0.0,  -0.0, 1.0},
            {0.0,  -1.0, 0.0},
            {0.0,  -1.0, 0.0},
            {-1.0, -0.0, -0.0},
            {-1.0, 0.0,  0.0},
            {1.0,  -0.0, 0.0},
            {1.0,  -0.0, 0.0},
            {0.0,  1.0,  0.0},
            {0.0,  1.0,  0.0},
            {0.0,  0.0,  -1.0},
            {0.0,  0.0,  -1.0}
    };

    std::vector<std::array<std::array<double, 3>, 3>> expectedSegmentUnitNormals{
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -1.0, -0.0}, {0.4472135954999579, 0.8944271909999159, 0.0}, {-1.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{1.0, 0.0, -0.0}, {0.0, 1.0, 0.0}, {-0.4472135954999579, -0.8944271909999159, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-1.0, -0.0, -0.0}, {0.0, 0.0, -1.0}, {0.4472135954999579, 0.0, 0.8944271909999159}}},
            std::array<std::array<double, 3>, 3>{
                    {{-0.4472135954999579, -0.0, -0.8944271909999159}, {1.0, 0.0, -0.0}, {0.0, 0.0, 1.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-0.0, 0.7071067811865475, 0.7071067811865475}, {0.0, 0.0, -1.0}, {0.0, -1.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -0.0, 1.0}, {0.0, 1.0, 0.0}, {-0.0, -0.7071067811865475, -0.7071067811865475}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -1.0, -0.0}, {0.0, 0.0, -1.0}, {0.0, 0.7071067811865475, 0.7071067811865475}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -0.7071067811865475, -0.7071067811865475}, {0.0, 1.0, -0.0}, {0.0, 0.0, 1.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.4472135954999579, -0.0, 0.8944271909999159}, {0.0, 0.0, -1.0}, {-1.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{1.0, -0.0, 0.0}, {-0.4472135954999579, 0.0, -0.8944271909999159}, {0.0, 0.0, 1.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-0.4472135954999579, 0.8944271909999159, 0.0}, {1.0, 0.0, 0.0}, {-0.0, -1.0, -0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-1.0, 0.0, 0.0}, {-0.0, 1.0, 0.0}, {0.4472135954999579, -0.8944271909999159, 0.0}}}
    };

    std::vector<double> expectedPlaneNormalOrientations{1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0};

    std::vector<polyhedralGravity::HessianPlane> expectedHessianPlanes{
            {0.0,    0.0,    200.0,  -5000.0},
            {0.0,    -0.0,   200.0,  -5000.0},
            {0.0,    -200.0, 0.0,    0.0},
            {0.0,    -200.0, 0.0,    0.0},
            {-100.0, 0.0,    0.0,    -2000.0},
            {-100.0, 0.0,    0.0,    -2000.0},
            {100.0,  0.0,    0.0,    -0.0},
            {100.0,  -0.0,   0.0,    0.0},
            {0.0,    200.0,  0.0,    -2000.0},
            {0.0,    200.0,  0.0,    -2000.0},
            {0.0,    0.0,    -200.0, 3000.0},
            {0.0,    0.0,    -200.0, 3000.0}
    };

    std::vector<double> expectedPlaneDistances{25.0,
                                               25.0,
                                               0.0,
                                               0.0,
                                               20.0,
                                               20.0,
                                               0.0,
                                               0.0,
                                               10.0,
                                               10.0,
                                               15.0,
                                               15.0};

    std::vector<std::array<double, 3>> expectedOrthogonalProjectionPointsOnPlane{
            {0.0,   0.0,  25.0},
            {0.0,   0.0,  25.0},
            {0.0,   0.0,  0.0},
            {0.0,   0.0,  0.0},
            {-20.0, 0.0,  0.0},
            {-20.0, 0.0,  0.0},
            {0.0,   0.0,  0.0},
            {0.0,   0.0,  0.0},
            {0.0,   10.0, 0.0},
            {0.0,   10.0, 0.0},
            {0.0,   0.0,  15.0},
            {0.0,   0.0,  15.0}
    };

    std::vector<std::array<double, 3>> expectedSegmentNormalOrientations{
            {0.0,  0.0,  1.0},
            {0.0,  1.0,  0.0},
            {1.0,  -1.0, 1.0},
            {-1.0, 0.0,  1.0},
            {1.0,  -1.0, 0.0},
            {1.0,  1.0,  -1.0},
            {0.0,  -1.0, 1.0},
            {-1.0, 1.0,  1.0},
            {1.0,  -1.0, 1.0},
            {0.0,  -1.0, 1.0},
            {1.0,  0.0,  0.0},
            {1.0,  1.0,  -1.0}
    };

    std::vector<std::array<std::array<double, 3>, 3>> expectedOrthogonalProjectionPointsOnSegment{
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, 25.0}, {0.0, 0.0, 25.0}, {-20.0, -0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, 25.0}, {-0.0, 10.0, 25.0}, {0.0, 0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, -0.0, -0.0}, {-0.0, -0.0, 15.0}, {6.0, -0.0, 12.0}}},
            std::array<std::array<double, 3>, 3>{{{6.0, -0.0, 12.0}, {0.0, 0.0, 0.0}, {-0.0, -0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, 12.5, 12.5}, {-20.0, -0.0, 15.0}, {-20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, -0.0, 25.0}, {-20.0, 10.0, -0.0}, {-20.0, 12.5, 12.5}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, 0.0}, {-0.0, -0.0, 15.0}, {-0.0, 12.5, 12.5}}},
            std::array<std::array<double, 3>, 3>{{{-0.0, 12.5, 12.5}, {-0.0, 10.0, -0.0}, {-0.0, -0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{6.0, 10.0, 12.0}, {-0.0, 10.0, 15.0}, {-20.0, 10.0, -0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {6.0, 10.0, 12.0}, {-0.0, 10.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{-4.0, 8.0, 15.0}, {0.0, 0.0, 15.0}, {0.0, 0.0, 15.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, -0.0, 15.0}, {-0.0, 10.0, 15.0}, {-4.0, 8.0, 15.0}}}
    };

    std::vector<std::array<double, 3>> expectedSegmentDistances{
            {0.0,                0.0,                20.0},
            {0.0,                10.0,               0.0},
            {20.0,               15.0,               13.416407864998739},
            {13.416407864998739, 0.0,                25.0},
            {17.67766952966369,  15.0,               0.0},
            {25.0,               10.0,               17.67766952966369},
            {0.0,                15.0,               17.67766952966369},
            {17.67766952966369,  10.0,               25.0},
            {13.416407864998739, 15.0,               20.0},
            {0.0,                13.416407864998739, 25.0},
            {8.94427190999916,   0.0,                0.0},
            {20.0,               10.0,               8.94427190999916}
    };

    std::vector<std::array<std::array<double, 2>, 3>> expected3DDistancesPerSegmentEndpoint{
            std::array<std::array<double, 2>, 3>{
                    {{32.01562118716424, 25.0}, {25.0, 33.54101966249684}, {33.54101966249684, 32.01562118716424}}},
            std::array<std::array<double, 2>, 3>{
                    {{25.0, 26.92582403567252}, {26.92582403567252, 33.54101966249684}, {33.54101966249684, 25.0}}},
            std::array<std::array<double, 2>, 3>{{{32.01562118716424, 25.0}, {25.0, 15.0}, {15.0, 32.01562118716424}}},
            std::array<std::array<double, 2>, 3>{{{32.01562118716424, 15.0}, {15.0, 25.0}, {25.0, 32.01562118716424}}},
            std::array<std::array<double, 2>, 3>{
                    {{32.01562118716424, 26.92582403567252}, {26.92582403567252, 25.0}, {25.0, 32.01562118716424}}},
            std::array<std::array<double, 2>, 3>{
                    {{32.01562118716424, 33.54101966249684}, {33.54101966249684, 26.92582403567252},
                     {26.92582403567252, 32.01562118716424}}},
            std::array<std::array<double, 2>, 3>{
                    {{-25.0, -15.0}, {15.0, 18.027756377319946}, {18.027756377319946, 25.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{25.0, 18.027756377319946}, {18.027756377319946, 26.92582403567252}, {26.92582403567252, 25.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{33.54101966249684, 18.027756377319946}, {18.027756377319946, 26.92582403567252},
                     {26.92582403567252, 33.54101966249684}}},
            std::array<std::array<double, 2>, 3>{
                    {{26.92582403567252, 18.027756377319946}, {18.027756377319946, 33.54101966249684},
                     {33.54101966249684, 26.92582403567252}}},
            std::array<std::array<double, 2>, 3>{
                    {{25.0, 18.027756377319946}, {18.027756377319946, 15.0}, {15.0, 25.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{25.0, 26.92582403567252}, {26.92582403567252, 18.027756377319946}, {18.027756377319946, 25.0}}}
    };

    std::vector<std::array<std::array<double, 2>, 3>> expected1DDistancesPerSegmentEndpoint{
            std::array<std::array<double, 2>, 3>{{{-20.0, -0.0}, {0.0, 22.360679774997898}, {-10.0, -0.0}}},
            std::array<std::array<double, 2>, 3>{{{0.0, 10.0}, {0.0, 20.0}, {-22.360679774997898, -0.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{-25.0, -15.0}, {-20.0, -0.0}, {6.708203932499369, 29.068883707497267}}},
            std::array<std::array<double, 2>, 3>{
                    {{-29.068883707497267, -6.708203932499369}, {15.0, 25.0}, {0.0, 20.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{-17.67766952966369, -3.5355339059327378}, {-10.0, -0.0}, {15.0, 25.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{0.0, 10.0}, {-25.0, -15.0}, {3.5355339059327378, 17.67766952966369}}},
            std::array<std::array<double, 2>, 3>{
                    {{-25.0, -15.0}, {0.0, 10.0}, {3.5355339059327378, 17.67766952966369}}},
            std::array<std::array<double, 2>, 3>{
                    {{-17.67766952966369, -3.5355339059327378}, {15.0, 25.0}, {-10.0, -0.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{-29.068883707497267, -6.708203932499369}, {0.0, 20.0}, {15.0, 25.0}}},
            std::array<std::array<double, 2>, 3>{
                    {{-25.0, -15.0}, {6.708203932499369, 29.068883707497267}, {-20.0, -0.0}}},
            std::array<std::array<double, 2>, 3>{{{-17.88854381999832, 4.47213595499958}, {-10.0, -0.0}, {0.0, 20.0}}},
            std::array<std::array<double, 2>, 3>{{{0.0, 10.0}, {-20.0, -0.0}, {-4.47213595499958, 17.88854381999832}}}
    };

    std::vector<std::array<polyhedralGravity::Distance, 3>> expectedDistancesPerSegmentEndpoint;

    std::vector<std::array<double, 3>> expectedTranscendentalLN{
            {0.0,                0.0,                 0.30747952872839945},
            {0.0,                0.687362255356451,   0.0},
            {0.3544458320893136, 1.0986122886681098,  1.0345679811316213},
            {1.034567981131622,  0.5108256237659907,  0.7326682560454109},
            {0.4894110007366263, 0.3900353197707153,  0.3544458320893134},
            {0.3074795287283993, 0.33382573681901684, 0.4894110007366262},
            {-0.510825623765990, 0.6251451172504167,  0.6826834766703017},
            {0.6826834766703017, 0.4524679290839864,  0.3900353197707153},
            {0.9286653985398196, 0.9566555518497877,  0.33382573681901667},
            {0.4524679290839866, 0.928665398539819,   0.6873622553564511},
            {1.1518034938098078, 0.0,                 0.0},
            {0.3900353197707153, 0.9566555518497877,  1.1518034938098078}
    };

    std::vector<std::array<double, 3>> expectedTranscendentalAN{
            {0.0,                 0.0,                 0.3567333885140938},
            {0.0,                 0.9799235766494776,  0.0},
            {0.0,                 0.0,                 0.0},
            {0.0,                 0.0,                 0.0},
            {0.4109023045514107,  0.45979025757734426, 0.0},
            {0.23413936163132537, 0.1405746311094993,  0.4109023045514107},
            {0.0,                 0.0,                 0.0},
            {0.0,                 0.0,                 0.0},
            {0.3029908626228055,  0.45979025757734426, 0.08507626483651975},
            {0.0,                 0.3029908626228055,  0.23413936163132537},
            {1.2703024256629791,  0.0,                 0.0},
            {0.27165712367757405, 0.8393489455399783,  1.2703024256629791}
    };

    std::vector<std::array<polyhedralGravity::TranscendentalExpression, 3>> expectedTranscendentalExpressions;


    std::vector<std::pair<double, std::array<double, 3>>> expectedSingularityTerms{
            std::make_pair(-11.591190225020153,
                           std::array<double, 3>{-0.0, -0.0, -0.46364760900080615}),
            std::make_pair(-27.67871794485226,
                           std::array<double, 3>{-0.0, -0.0, -1.1071487177940904}),
            std::make_pair(0.0, std::array<double, 3>{-0.0, 0.0,  -0.0}),
            std::make_pair(0.0, std::array<double, 3>{-0.0, 0.0,  -0.0}),
            std::make_pair(0.0, std::array<double, 3>{0.0,  0.0,  0.0}),
            std::make_pair(0.0, std::array<double, 3>{0.0,  -0.0, -0.0}),
            std::make_pair(0.0, std::array<double, 3>{-0.0, -0.0, -0.0}),
            std::make_pair(0.0, std::array<double, 3>{-0.0, -0.0, -0.0}),
            std::make_pair(0.0, std::array<double, 3>{-0.0, -0.0, -0.0}),
            std::make_pair(0.0, std::array<double, 3>{-0.0, -0.0, -0.0}),
            std::make_pair(-23.5619455575943,
                           std::array<double, 3>{-0.0, -0.0, -1.5707963705062866}),
            std::make_pair(0.0, std::array<double, 3>{-0.0, -0.0, 0.0})
    };

    std::vector<double> expectedAlphaSingularityTerms{-11.591190225020153, -27.67871794485226, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                      0.0, 0.0, 0.0, -23.5619455575943, 0.0};

    std::vector<std::array<double, 3>> expectedBetaSingularityTerms{
            {-0.0, -0.0, -0.46364760900080615},
            {-0.0, -0.0, -1.1071487177940904},
            {-0.0, 0.0,  -0.0},
            {-0.0, 0.0,  -0.0},
            {0.0,  0.0,  0.0},
            {0.0,  -0.0, -0.0},
            {-0.0, -0.0, -0.0},
            {-0.0, -0.0, -0.0},
            {-0.0, -0.0, -0.0},
            {-0.0, -0.0, -0.0},
            {-0.0, -0.0, -1.5707963705062866},
            {-0.0, -0.0, 0.0}
    };

public:

    GravityModelTest() : ::testing::Test() {
        expectedDistancesPerSegmentEndpoint.resize(expected3DDistancesPerSegmentEndpoint.size());
        expectedTranscendentalExpressions.resize(expectedTranscendentalLN.size());
        for (int i = 0; i < expected3DDistancesPerSegmentEndpoint.size(); ++i) {
            for (int j = 0; j < expected3DDistancesPerSegmentEndpoint[i].size(); ++j) {
                expectedDistancesPerSegmentEndpoint[i][j] = polyhedralGravity::Distance{
                        expected3DDistancesPerSegmentEndpoint[i][j][0],
                        expected3DDistancesPerSegmentEndpoint[i][j][1],
                        expected1DDistancesPerSegmentEndpoint[i][j][0],
                        expected1DDistancesPerSegmentEndpoint[i][j][1]
                };
                expectedTranscendentalExpressions[i][j] = polyhedralGravity::TranscendentalExpression{
                        expectedTranscendentalLN[i][j],
                        expectedTranscendentalAN[i][j],
                };
            }
        }
    }

};

TEST_F(GravityModelTest, GijVectors) {
    using namespace testing;

    auto actualGij = polyhedralGravity::GravityModel::calculateSegmentVectors(_polyhedron);

    ASSERT_THAT(actualGij, ContainerEq(expectedGij));
}

TEST_F(GravityModelTest, PlaneUnitNormals) {
    using namespace testing;

    auto actualPlaneUnitNormals = polyhedralGravity::GravityModel::calculatePlaneUnitNormals(expectedGij);

    ASSERT_THAT(actualPlaneUnitNormals, ContainerEq(expectedPlaneUnitNormals));
}

TEST_F(GravityModelTest, SegmentUnitNormals) {
    using namespace testing;

    auto actualSegmentUnitNormals = polyhedralGravity::GravityModel::calculateSegmentUnitNormals(expectedGij,
                                                                                                 expectedPlaneUnitNormals);

    ASSERT_THAT(actualSegmentUnitNormals, ContainerEq(expectedSegmentUnitNormals));
}

TEST_F(GravityModelTest, PlaneNormalOrientations) {
    using namespace testing;
    using namespace polyhedralGravity;

    auto actualPlaneNormalOrientations = GravityModel::calculatePlaneNormalOrientations(
            _computationPoint, _polyhedron, expectedPlaneUnitNormals);

    ASSERT_THAT(actualPlaneNormalOrientations, ContainerEq(expectedPlaneNormalOrientations));
}

TEST_F(GravityModelTest, SimpleHessianPlane) {
    using namespace testing;
    using namespace polyhedralGravity;

    HessianPlane expectedHessian{2, -8, 5, -18};

    auto actualHessianPlane = GravityModel::detail::computeHessianPlane({1, -2, 0}, {3, 1, 4}, {0, -1, 2});

    ASSERT_DOUBLE_EQ(actualHessianPlane.a, expectedHessian.a);
    ASSERT_DOUBLE_EQ(actualHessianPlane.b, expectedHessian.b);
    ASSERT_DOUBLE_EQ(actualHessianPlane.c, expectedHessian.c);
    ASSERT_DOUBLE_EQ(actualHessianPlane.d, expectedHessian.d);
}

TEST_F(GravityModelTest, HessianPlane) {
    using namespace testing;
    using namespace polyhedralGravity;

    auto actualHessianPlane =
            polyhedralGravity::GravityModel::calculateFacesToHessianPlanes(_computationPoint, _polyhedron);

    ASSERT_EQ(actualHessianPlane, expectedHessianPlanes);
}

TEST_F(GravityModelTest, PlaneDistances) {
    using namespace testing;

    auto actualPlaneDistances = polyhedralGravity::GravityModel::calculatePlaneDistances(expectedHessianPlanes);

    ASSERT_THAT(actualPlaneDistances, ContainerEq(expectedPlaneDistances));
}

TEST_F(GravityModelTest, OrthogonalProjectionPointsOnPlane) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnPlane = polyhedralGravity::GravityModel::calculateOrthogonalProjectionPointsOnPlane(
            expectedHessianPlanes, expectedPlaneUnitNormals, expectedPlaneDistances);

    ASSERT_THAT(actualOrthogonalProjectionPointsOnPlane, ContainerEq(expectedOrthogonalProjectionPointsOnPlane));
}

TEST_F(GravityModelTest, SegmentNormalOrientations) {
    using namespace testing;

    auto actualSegmentNormalOrientations =
            polyhedralGravity::GravityModel::calculateSegmentNormalOrientations(_computationPoint, _polyhedron,
                                                                                expectedSegmentUnitNormals,
                                                                                expectedOrthogonalProjectionPointsOnPlane);

    ASSERT_THAT(actualSegmentNormalOrientations, ContainerEq(expectedSegmentNormalOrientations));
}

TEST_F(GravityModelTest, OrthogonalProjectionPointsOnSegment) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnSegment =
            polyhedralGravity::GravityModel::calculateOrthogonalProjectionPointsOnSegments(
                    _computationPoint, _polyhedron,
                    expectedOrthogonalProjectionPointsOnPlane,
                    expectedSegmentNormalOrientations);

    ASSERT_THAT(actualOrthogonalProjectionPointsOnSegment, ContainerEq(expectedOrthogonalProjectionPointsOnSegment));
}

TEST_F(GravityModelTest, SegmentDistances) {
    using namespace testing;

    auto actualSegmentDistances =
            polyhedralGravity::GravityModel::calculateSegmentDistances(expectedOrthogonalProjectionPointsOnPlane,
                                                                       expectedOrthogonalProjectionPointsOnSegment);

    ASSERT_THAT(actualSegmentDistances, ContainerEq(expectedSegmentDistances));
}

TEST_F(GravityModelTest, DistancesPerSegmentEndpoint) {
    using namespace testing;

    auto actualDistancesPerSegmentEndpoint =
            polyhedralGravity::GravityModel::calculateDistances(_computationPoint, _polyhedron, expectedGij,
                                                                expectedOrthogonalProjectionPointsOnSegment);

    ASSERT_THAT(actualDistancesPerSegmentEndpoint, ContainerEq(expectedDistancesPerSegmentEndpoint));
}

TEST_F(GravityModelTest, TranscendentalExpressions) {
    using namespace testing;
    using namespace polyhedralGravity;

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

TEST_F(GravityModelTest, SingularityTerms) {
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
