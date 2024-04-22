#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <set>
#include <chrono>
#include <vector>
#include <utility>
#include <array>
#include <algorithm>
#include <random>
#include <thrust/generate.h>
#include <iostream>
#include "polyhedralGravity/model/Polyhedron.h"

/**
 * Cheks that the detection of wrong faces works as intended.
 */
class PolyhedronBigTest : public ::testing::TestWithParam<std::set<size_t>> {

protected:

    const static inline polyhedralGravity::PolyhedralFiles FILENAMES{"resources/GravityModelBigTest.node", "resources/GravityModelBigTest.face"};

    const static inline polyhedralGravity::Polyhedron CORRECT_POLYHEDRON{
        FILENAMES, 1.0,
        polyhedralGravity::NormalOrientation::OUTWARDS, polyhedralGravity::PolyhedronIntegrity::DISABLE
    };

    static constexpr size_t FACES_COUNT = 14744;
    static constexpr size_t SET_SIZE = 100;
    static constexpr size_t SET_NUMBER = 10;
    static constexpr size_t SEED = 42;

    /**
     * Creates a polyhedron violating the OUTWARDS pointing normals constraint for the given face indices.
     * @param violatinIndices the indices of the faces where the polyhedron's plane unit normals shall point INWARDS
     * @return invalid polyhedron
     */
    polyhedralGravity::Polyhedron createViolatingPolyhedron(const std::set<size_t> &violatinIndices) {
        using namespace polyhedralGravity;
        std::vector<IndexArray3> violatingFaces = CORRECT_POLYHEDRON.getFaces();
        for (const size_t index: violatinIndices) {
            std::swap(violatingFaces[index][0], violatingFaces[index][1]);
        }
        return {
            CORRECT_POLYHEDRON.getVertices(), violatingFaces,
            1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE
        };
    }

public:

    // Function to generate a half sized set
    static std::vector<std::set<size_t>> generateIndices() {
        std::vector<std::set<size_t>> indexSets{};
        indexSets.reserve(SET_NUMBER);

        std::mt19937 engine{SEED};
        std::uniform_int_distribution<size_t> dist{0, PolyhedronBigTest::FACES_COUNT - 1};

        for (size_t i = 0; i < SET_NUMBER; ++i) {
            std::set<size_t> generatedSet;
            thrust::generate_n(std::inserter(generatedSet, generatedSet.end()), SET_SIZE, [&dist, &engine]() {return dist(engine);});
            while(generatedSet.size() < SET_SIZE) {
                generatedSet.insert(dist(engine));
            }
            indexSets.push_back(generatedSet);
        }
        return indexSets;
    }

};


TEST_P(PolyhedronBigTest, BigPolyhedronFindWrongVertices) {
    using namespace testing;
    using namespace polyhedralGravity;

    const std::set<size_t> expectedViolatingIndices = GetParam();
    ASSERT_EQ(expectedViolatingIndices.size(), SET_SIZE);

    Polyhedron invalidPolyhedron = createViolatingPolyhedron(expectedViolatingIndices);

    auto start = std::chrono::high_resolution_clock::now();

    const auto&[actualOrientation, actualViolatingIndices] = invalidPolyhedron.checkPlaneUnitNormalOrientation();

    auto end = std::chrono::high_resolution_clock::now();
    auto dur = end - start;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur);
    std::cout << "Measured time: " << ms.count() << " microseconds ";


    // The orientation changes if violatingIndcies > len(faces)/2 --> So ensure that the size of the test-parameter is smaller than len(faces)/2
    ASSERT_EQ(actualOrientation, NormalOrientation::OUTWARDS);
    ASSERT_THAT(actualViolatingIndices, ContainerEq(expectedViolatingIndices));
}

INSTANTIATE_TEST_SUITE_P(NormalOrientationViolationTest, PolyhedronBigTest, ::testing::ValuesIn(PolyhedronBigTest::generateIndices()));


