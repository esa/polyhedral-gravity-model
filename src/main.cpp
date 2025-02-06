#include "polyhedralGravity/Info.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/output/CSVWriter.h"
#include "polyhedralGravity/output/Logging.h"
#include <chrono>

int main(const int argc, char *argv[]) {
    using namespace polyhedralGravity;

    POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Model Version:                 {}", POLYHEDRAL_GRAVITY_VERSION);
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Commit Hash:                   {}", POLYHEDRAL_GRAVITY_COMMIT_HASH);
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Model Parallelization Backend: {}", POLYHEDRAL_GRAVITY_PARALLELIZATION);
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Logging Level:                 {}", POLYHEDRAL_GRAVITY_LOGGING_LEVEL);

    if (argc != 2) {
        POLYHEDRAL_GRAVITY_LOG_ERROR("Wrong program call! "
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

        POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");
        POLYHEDRAL_GRAVITY_LOG_INFO("Number of Vertices:                               {}", std::get<0>(polyhedralSource).size());
        POLYHEDRAL_GRAVITY_LOG_INFO("Number of Faces:                                  {}", std::get<1>(polyhedralSource).size());
        POLYHEDRAL_GRAVITY_LOG_INFO("Number of Computation Points:                     {}", computationPoints.size());
        POLYHEDRAL_GRAVITY_LOG_INFO("Mesh Check Enabled:                               {}", config->getMeshInputCheckStatus());
        POLYHEDRAL_GRAVITY_LOG_INFO("Output File:                                      {}", outputFileName);
        POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");

        POLYHEDRAL_GRAVITY_LOG_INFO("The calculation started...");
        const auto start = std::chrono::high_resolution_clock::now();
        const Polyhedron polyhedron{polyhedralSource, density, NormalOrientation::OUTWARDS, checkPolyhedralInput};
        const auto result = GravityModel::evaluate(polyhedron, computationPoints, true);
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = end - start;
        const auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        POLYHEDRAL_GRAVITY_LOG_INFO("The calculation finished. It took {} microseconds.", ms.count());
        POLYHEDRAL_GRAVITY_LOG_INFO("Or on average {} microseconds/point",
                                    static_cast<double>(ms.count()) / static_cast<double>(computationPoints.size()));

        //The results
        if (!outputFileName.empty()) {
            POLYHEDRAL_GRAVITY_LOG_INFO("Writing results to specified output file {}", outputFileName);
            const CSVWriter csvWriter{outputFileName};
            csvWriter.printResult(computationPoints, result);
            POLYHEDRAL_GRAVITY_LOG_INFO("Writing finished!");
        } else {
            POLYHEDRAL_GRAVITY_LOG_WARN("No output filename was specified!");
        }

        return 0;

    } catch (const std::exception &e) {
        POLYHEDRAL_GRAVITY_LOG_ERROR("{}", e.what());
        return -1;
    }
}