#include "polyhedralGravity/input/YAMLConfigReader.h"
#include <benchmark/benchmark.h>

namespace polyhedralGravity {
    const std::vector<std::string> erosFilePaths{
        "../example-config/data/Eros_scaled-1000", "../example-config/data/Eros_scaled-1732",
        "../example-config/data/Eros_scaled-3000", "../example-config/data/Eros_scaled-5196",
        "../example-config/data/Eros_scaled-9000", "../example-config/data/Eros_scaled-15588",
        "../example-config/data/Eros_scaled-27000", "../example-config/data/Eros_scaled-46765",
        "../example-config/data/Eros_scaled-81000", "../example-config/data/Eros_scaled-140296"
    };

    const std::vector<std::string> sphereFilePaths{
        "../example-config/data/sphere_scaled-1000", "../example-config/data/sphere_scaled-1732",
        "../example-config/data/sphere_scaled-3000", "../example-config/data/sphere_scaled-5196",
        "../example-config/data/sphere_scaled-9000", "../example-config/data/sphere_scaled-15588",
        "../example-config/data/sphere_scaled-27000", "../example-config/data/sphere_scaled-46765",
        "../example-config/data/sphere_scaled-81000", "../example-config/data/sphere_scaled-140296"
    };

    Polyhedron createPolyhedron(const PolyhedralFiles &polyhedralFiles,
                                   const PlaneSelectionAlgorithm::Algorithm &algorithm) {
        return {polyhedralFiles, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, algorithm};
    }

    void BM_Eros_Polyhedron_Tree(benchmark::State &state, const PlaneSelectionAlgorithm::Algorithm &algorithm) {
        const PolyhedralFiles polyhedralFiles = {erosFilePaths.at(state.range(0)) + ".node", erosFilePaths.at(state.range(0)) + ".face"};
        Polyhedron polyhedron = {polyhedralFiles, 1.0, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE, algorithm};//createPolyhedron(polyhedralFiles, algorithm);
        for (auto _: state) {
            auto result = polyhedron.checkPlaneUnitNormalOrientation();
            benchmark::ClobberMemory();
        }
        state.SetComplexityN(static_cast<benchmark::ComplexityN>(polyhedron.countFaces()));
    }

    void BM_Sphere_Polyhedron_Tree(benchmark::State &state, const PlaneSelectionAlgorithm::Algorithm &algorithm) {
        const PolyhedralFiles polyhedralFiles = {sphereFilePaths.at(state.range(0)) + ".node", sphereFilePaths.at(state.range(0)) + ".face"};
        Polyhedron polyhedron = createPolyhedron(polyhedralFiles, algorithm);
        for (auto _: state) {
            auto result = polyhedron.checkPlaneUnitNormalOrientation();
            benchmark::ClobberMemory();
        }
        state.SetComplexityN(static_cast<benchmark::ComplexityN>(polyhedron.countFaces()));
    }

    void BM_Eros_Polyhedron_Tree_Twice(benchmark::State &state) {
        const PolyhedralFiles polyhedralFiles = {erosFilePaths.at(state.range(0)) + ".node", erosFilePaths.at(state.range(0)) + ".face"};
        Polyhedron polyhedron = createPolyhedron(polyhedralFiles, PlaneSelectionAlgorithm::Algorithm::LOG);
        polyhedron.prebuildKDTree();
        for (auto _: state) {
            auto result = polyhedron.checkPlaneUnitNormalOrientation();
            benchmark::ClobberMemory();
        }
        state.SetComplexityN(static_cast<benchmark::ComplexityN>(polyhedron.countFaces()));
    }

    void BM_Sphere_Polyhedron_Tree_Twice(benchmark::State &state) {
        const PolyhedralFiles polyhedralFiles = {sphereFilePaths.at(state.range(0)) + ".node", sphereFilePaths.at(state.range(0)) + ".face"};
        Polyhedron polyhedron = createPolyhedron(polyhedralFiles, PlaneSelectionAlgorithm::Algorithm::LOG);
        polyhedron.prebuildKDTree();
        for (auto _: state) {
            auto result = polyhedron.checkPlaneUnitNormalOrientation();
            benchmark::ClobberMemory();
        }
        state.SetComplexityN(static_cast<benchmark::ComplexityN>(polyhedron.countFaces()));
    }

    void BM_Eros_Polyhedron_Tree_Build(benchmark::State &state, const PlaneSelectionAlgorithm::Algorithm &algorithm) {
        const PolyhedralFiles polyhedralFiles = {erosFilePaths.at(state.range(0)) + ".node", erosFilePaths.at(state.range(0)) + ".face"};
        const auto polyhedron = createPolyhedron(polyhedralFiles, algorithm);
        auto workload = [polyhedron](){polyhedron.prebuildKDTree(); return "finished";};
        for (auto _: state) {
            benchmark::DoNotOptimize(workload());
            benchmark::ClobberMemory();
        }
        state.SetComplexityN(static_cast<benchmark::ComplexityN>(polyhedron.countFaces()));
    }

    void BM_Sphere_Polyhedron_Tree_Build(benchmark::State &state, const PlaneSelectionAlgorithm::Algorithm &algorithm) {
        const PolyhedralFiles polyhedralFiles = {sphereFilePaths.at(state.range(0)) + ".node", sphereFilePaths.at(state.range(0)) + ".face"};
        const auto polyhedron = createPolyhedron(polyhedralFiles, algorithm);
        auto workload = [polyhedron](){polyhedron.prebuildKDTree(); return "finished";};
        for (auto _: state) {
            benchmark::DoNotOptimize(workload());
            benchmark::ClobberMemory();
        }
        state.SetComplexityN(static_cast<benchmark::ComplexityN>(polyhedron.countFaces()));
    }

    // eros mesh benchmarks
    BENCHMARK_CAPTURE(BM_Eros_Polyhedron_Tree, "ErosPolyhedronNoTree", PlaneSelectionAlgorithm::Algorithm::NOTREE)->DenseRange(
        0, static_cast<long>(erosFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Eros_Polyhedron_Tree, "ErosPolyhedronQuadratic",
                      PlaneSelectionAlgorithm::Algorithm::QUADRATIC)->DenseRange(0, static_cast<long>(erosFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Eros_Polyhedron_Tree, "ErosPolyhedronLogSquared",
                      PlaneSelectionAlgorithm::Algorithm::LOGSQUARED)->DenseRange(0, static_cast<long>(erosFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Eros_Polyhedron_Tree, "ErosPolyhedronLog", PlaneSelectionAlgorithm::Algorithm::LOG)->DenseRange(
        0, static_cast<long>(erosFilePaths.size() - 1), 1);
    BENCHMARK(BM_Eros_Polyhedron_Tree_Twice)->Name("ErosPolyhedronSecondRun")->DenseRange(
        0, static_cast<long>(erosFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Eros_Polyhedron_Tree_Build, "ErosPolyhedronBuildTreeSquared", PlaneSelectionAlgorithm::Algorithm::QUADRATIC)->DenseRange(
        0, static_cast<long>(erosFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Eros_Polyhedron_Tree_Build, "ErosPolyhedronBuildTreeLogSquared", PlaneSelectionAlgorithm::Algorithm::LOGSQUARED)->DenseRange(
        0, static_cast<long>(erosFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Eros_Polyhedron_Tree_Build, "ErosPolyhedronBuildTreeLog", PlaneSelectionAlgorithm::Algorithm::LOG)->DenseRange(
        0, static_cast<long>(erosFilePaths.size() - 1), 1);

    // sphere mesh benchmarks
    BENCHMARK_CAPTURE(BM_Sphere_Polyhedron_Tree, "SpherePolyhedronNoTree", PlaneSelectionAlgorithm::Algorithm::NOTREE)->DenseRange(
        0, static_cast<long>(sphereFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Sphere_Polyhedron_Tree, "SpherePolyhedronQuadratic",
                      PlaneSelectionAlgorithm::Algorithm::QUADRATIC)->DenseRange(0, static_cast<long>(sphereFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Sphere_Polyhedron_Tree, "SpherePolyhedronLogSquared",
                      PlaneSelectionAlgorithm::Algorithm::LOGSQUARED)->DenseRange(0, static_cast<long>(sphereFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Sphere_Polyhedron_Tree, "SpherePolyhedronLog", PlaneSelectionAlgorithm::Algorithm::LOG)->DenseRange(
        0, static_cast<long>(sphereFilePaths.size() - 1), 1);
    BENCHMARK(BM_Sphere_Polyhedron_Tree_Twice)->Name("SpherePolyhedronSecondRun")->DenseRange(
        0, static_cast<long>(sphereFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Sphere_Polyhedron_Tree_Build, "SpherePolyhedronBuildTreeSquared", PlaneSelectionAlgorithm::Algorithm::QUADRATIC)->DenseRange(
        0, static_cast<long>(sphereFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Sphere_Polyhedron_Tree_Build, "SpherePolyhedronBuildTreeLogSquared", PlaneSelectionAlgorithm::Algorithm::LOGSQUARED)->DenseRange(
        0, static_cast<long>(sphereFilePaths.size() - 1), 1);
    BENCHMARK_CAPTURE(BM_Sphere_Polyhedron_Tree_Build, "SpherePolyhedronBuildTreeLog", PlaneSelectionAlgorithm::Algorithm::LOG)->DenseRange(
        0, static_cast<long>(sphereFilePaths.size() - 1), 1);
} // namespace polyhedralGravity

BENCHMARK_MAIN();
