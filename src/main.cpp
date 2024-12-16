#include "polyhedralGravity/Info.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/output/CSVWriter.h"
#include "polyhedralGravity/output/Logging.h"
#include <chrono>

int main(int argc, char *argv[]) {
    using namespace polyhedralGravity;

    SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "####################################################################################");
    SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Polyhedral Gravity Model Version:                 {}", POLYHEDRAL_GRAVITY_VERSION);
    SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Polyhedral Gravity Commit Hash:                   {}", POLYHEDRAL_GRAVITY_COMMIT_HASH);
    SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Polyhedral Gravity Model Parallelization Backend: {}", POLYHEDRAL_GRAVITY_PARALLELIZATION);
    SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Polyhedral Gravity Logging Level:                 {}", POLYHEDRAL_GRAVITY_LOGGING_LEVEL);

    if (argc != 2) {
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Wrong program call! "
                                                                                "Please use the program like this:\n"
                                                                                "./polyhedralGravity [YAML-Configuration-File]\n");
        return 0;
    }

    try {
        const std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>(argv[1]);
        const auto polyhedralSource = config->getDataSource()->getPolyhedralSource();
        const auto density = config->getDensity();
        const auto computationPoints = config->getPointsOfInterest();
        const auto outputFileName = config->getOutputFileName();
        const PolyhedronIntegrity checkPolyhedralInput = config->getMeshInputCheckStatus() ? PolyhedronIntegrity::HEAL : PolyhedronIntegrity::DISABLE;

        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "####################################################################################");
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Number of Vertices:                               {}", std::get<0>(polyhedralSource).size());
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Number of Faces:                                  {}", std::get<1>(polyhedralSource).size());
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Number of Computation Points:                     {}", computationPoints.size());
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Mesh Check Enabled:                               {}", config->getMeshInputCheckStatus());
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Output File:                                      {}", outputFileName);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "####################################################################################");

        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "The calculation started...");
        const auto start = std::chrono::high_resolution_clock::now();
        const Polyhedron polyhedron{polyhedralSource, density, NormalOrientation::OUTWARDS, checkPolyhedralInput};
        const auto result = GravityModel::evaluate(polyhedron, computationPoints, true);
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = end - start;
        const auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "The calculation finished. It took {} microseconds.", ms.count());
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Or on average {} microseconds/point",
                           static_cast<double>(ms.count()) / static_cast<double>(computationPoints.size()));

        //The results
        if (!outputFileName.empty()) {
            SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                               "Writing results to specified output file {}", outputFileName);
            const CSVWriter csvWriter{outputFileName};
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