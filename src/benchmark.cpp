/**
 * A minimal driver which evaluates one polyhedron at N random computation points on one backend.
 *
 * It exists so that the multi point kernel can be timed and profiled (nsys, ncu) without the Python
 * interpreter in between. Build it with -DBUILD_POLYHEDRAL_GRAVITY_BENCHMARK=ON, then call it as
 *
 *     ./polyhedralGravity_benchmark <node file> <face file> [points] [repetitions] [backend] [precision]
 *
 * where backend is one of serial|cpu|gpu and precision is one of f32|f64. A profiler wants
 * "1 1 gpu f32", i.e. exactly one launch of the kernel under investigation.
 */
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/util/KokkosSession.h"

using namespace polyhedralGravity;

int main(const int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <node file> <face file> [points] [repetitions] [serial|cpu|gpu] [f32|f64]\n";
        return 1;
    }
    const std::string nodeFile{argv[1]};
    const std::string faceFile{argv[2]};
    const size_t pointCount = argc > 3 ? std::strtoul(argv[3], nullptr, 10) : 1000;
    const size_t repetitions = argc > 4 ? std::strtoul(argv[4], nullptr, 10) : 10;
    const std::string backendName = argc > 5 ? argv[5] : "gpu";
    const std::string precisionName = argc > 6 ? argv[6] : "f32";

    const ComputeBackend backend = backendName == "serial" ? ComputeBackend::CPU_SERIAL
            : backendName == "cpu"                         ? ComputeBackend::CPU_PARALLEL
                                                           : ComputeBackend::GPU_PARALLEL;
    const ComputePrecision precision =
            precisionName == "f64" ? ComputePrecision::FLOAT64 : ComputePrecision::FLOAT32;

    const Polyhedron polyhedron{std::vector<std::string>{nodeFile, faceFile}, 1.0, NormalOrientation::OUTWARDS,
                                PolyhedronIntegrity::DISABLE};

    // The same points the Python benchmark uses, i.e. uniformly drawn from the cube [-2, 2]^3
    std::mt19937 generator{42};
    std::uniform_real_distribution<double> distribution{-2.0, 2.0};
    std::vector<Array3> computationPoints(pointCount);
    for (Array3 &point: computationPoints) {
        point = {distribution(generator), distribution(generator), distribution(generator)};
    }

    const GravityEvaluable evaluable{polyhedron, backend, precision};

    // One untimed evaluation, so that neither the lazily created kernel nor a cold cache is measured
    const auto warmup = evaluable(computationPoints);
    std::cout << "Faces: " << polyhedron.countFaces() << ", vertices: " << polyhedron.countVertices()
              << ", points: " << pointCount << ", backend: " << backendName << " " << precisionName << "\n"
              << "Checksum: " << std::get<0>(std::get<std::vector<GravityModelResult>>(warmup).front()) << "\n";

    // Accuracy against a serial double precision evaluation of the very same points, so that a change to
    // the kernel's arithmetic can be judged on accuracy and not only on runtime
    {
        const GravityEvaluable reference{polyhedron, ComputeBackend::CPU_PARALLEL, ComputePrecision::FLOAT64};
        const auto expected = std::get<std::vector<GravityModelResult>>(reference(computationPoints));
        const auto actual = std::get<std::vector<GravityModelResult>>(warmup);
        double worstPotential = 0.0;
        double worstAcceleration = 0.0;
        for (size_t index = 0; index < expected.size(); ++index) {
            const auto relative = [](const double got, const double want) {
                return want == 0.0 ? 0.0 : std::abs(got - want) / std::abs(want);
            };
            worstPotential = std::max(worstPotential,
                                      relative(std::get<0>(actual[index]), std::get<0>(expected[index])));
            // Relative to the magnitude of the vector, not component wise: single components pass through
            // zero, where any component wise relative error is meaningless
            double errorNorm = 0.0;
            double referenceNorm = 0.0;
            for (size_t component = 0; component < 3; ++component) {
                const double difference =
                        std::get<1>(actual[index])[component] - std::get<1>(expected[index])[component];
                errorNorm += difference * difference;
                referenceNorm += std::get<1>(expected[index])[component] * std::get<1>(expected[index])[component];
            }
            worstAcceleration = std::max(worstAcceleration, std::sqrt(errorNorm / referenceNorm));
        }
        std::cout << "MaxRelErr potential: " << worstPotential << ", acceleration: " << worstAcceleration << "\n";
    }

    double best = std::numeric_limits<double>::max();
    double total = 0.0;
    for (size_t repetition = 0; repetition < repetitions; ++repetition) {
        const auto start = std::chrono::high_resolution_clock::now();
        const auto result = evaluable(computationPoints);
        const auto end = std::chrono::high_resolution_clock::now();
        const double milliseconds = std::chrono::duration<double, std::milli>(end - start).count();
        best = std::min(best, milliseconds);
        total += milliseconds;
        (void) result;
    }
    std::cout << "Best:   " << best << " ms (" << best * 1e3 / static_cast<double>(pointCount) << " us/point)\n"
              << "Mean:   " << total / static_cast<double>(repetitions) << " ms ("
              << total * 1e3 / static_cast<double>(repetitions) / static_cast<double>(pointCount)
              << " us/point)\n";
    return 0;
}
