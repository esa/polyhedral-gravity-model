#include <chrono>
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/calculation/MeshChecking.h"
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
        auto poly = config->getDataSource()->getPolyhedron();
        auto density = config->getDensity();
        auto computationPoints = config->getPointsOfInterest();
        auto outputFileName = config->getOutputFileName();
        bool checkPolyhedralInput = config->getMeshInputCheckStatus();

        // Checking that the vertices are correctly set-up in the input if activated
        if (checkPolyhedralInput) {
            SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Checking mesh...");
            if (!MeshChecking::checkTrianglesNotDegenerated(poly)) {
                throw std::runtime_error{
                        "At least on triangle in the mesh is degenerated and its surface area equals zero!"};
            } else if (!MeshChecking::checkNormalsOutwardPointing(poly)) {
                throw std::runtime_error{
                        "The plane unit normals are not pointing outwards! Please check the order "
                        "of the vertices in the polyhedral input source!"};
            }  else {
                SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "The mesh is fine.");
            }
        }

        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "The calculation started...");
        auto start = std::chrono::high_resolution_clock::now();
        auto result = GravityModel::evaluate(poly, density, computationPoints);
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