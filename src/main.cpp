#include <chrono>
#include "spdlog/spdlog.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/calculation/GravityModel.h"

int main(int argc, char *argv[]) {
    using namespace polyhedralGravity;
    if (argc != 2) {
        SPDLOG_INFO("Wrong program call! Please use the program like that:\n"
                    "./polyhedralGravity [YAML-Configuration-File]\n");
        return 0;
    }

    std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>(argv[1]);
    auto poly = config->getDataSource()->getPolyhedron();
    auto density = config->getDensity();
    auto point = config->getPointsOfInterest()[0];
    auto points = std::vector<Array3>{};
    constexpr size_t M = 1000;
    points.resize(M, point);

    GravityEvaluable evaluable{poly, density};
    constexpr size_t N = 10;
    std::array<size_t, N> time{};
    for (int i = 0; i < N; ++i) {
        //SPDLOG_INFO("The calculation started.");
        auto start = std::chrono::high_resolution_clock::now();
        auto result = evaluable.evaluate(points);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        time[i] = ms.count();
    }

    double avg = 0;
    for (int i = 0; i < N; ++i) {
        avg += time[i];
    }
    avg /= static_cast<double>(N * M);
    SPDLOG_INFO("Average {} microseconds.", avg);

    return 0;
}