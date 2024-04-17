#include <chrono>
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/output/Logging.h"
#include "polyhedralGravity/output/CSVWriter.h"

int main(int argc, char *argv[]) {
    using namespace polyhedralGravity;
    if (argc != 2) {
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Wrong program call! "
                                                                                "Please use the program like this:\n"
                                                                                "./polyhedralGravity [YAML-Configuration-File]\n");
        return 0;
    }

    try {
        std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>(argv[1]);
        auto polyhedralSource = config->getDataSource()->getPolyhedralSource();
        auto density = config->getDensity();
        auto computationPoints = config->getPointsOfInterest();
        auto outputFileName = config->getOutputFileName();
        PolyhedronIntegrity checkPolyhedralInput = config->getMeshInputCheckStatus() ? PolyhedronIntegrity::HEAL : PolyhedronIntegrity::DISABLE;

        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "The calculation started...");
        auto start = std::chrono::high_resolution_clock::now();
        Polyhedron polyhedron{polyhedralSource, density, NormalOrientation::OUTWARDS, checkPolyhedralInput};
        auto result = GravityModel::evaluate(polyhedron, computationPoints, true);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "The calculation finished. It took {} microseconds.", ms.count());
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Or on average {} microseconds/point",
                           static_cast<double>(ms.count()) / static_cast<double>(computationPoints.size()));

        //The results
        if (!outputFileName.empty()) {
            SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                               "Writing results to specified output file {}", outputFileName);
            CSVWriter csvWriter{outputFileName};
            csvWriter.printResult(computationPoints, result);
            SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                               "Writing finished!");
        } else {
            SPDLOG_LOGGER_WARN(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                               "No output filename was specified!");
        }

        return 0;

    } catch (const std::exception &e) {
        SPDLOG_LOGGER_ERROR(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "{}", e.what());
        return -1;
    }
}