#include "polyhedralGravity/Info.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/model/GravityModel.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/output/CSVWriter.h"
#include "polyhedralGravity/output/Logging.h"
#include <chrono>
#include <sstream>

int main(const int argc, char *argv[]) {
    using namespace polyhedralGravity;

    POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Model Version:                 {}", POLYHEDRAL_GRAVITY_VERSION);
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Commit Hash:                   {}", POLYHEDRAL_GRAVITY_COMMIT_HASH);
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Model Parallelization Backend: {}", POLYHEDRAL_GRAVITY_PARALLELIZATION);
    POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedral Gravity Logging Level:                 {}", POLYHEDRAL_GRAVITY_LOGGING_LEVEL);
    POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");

    if (argc != 2) {
        POLYHEDRAL_GRAVITY_LOG_ERROR("Wrong program call! Please use the program like this: ./polyhedralGravity [YAML-Configuration-File]");
        return 0;
    }

    try {
        const std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>(argv[1]);
        const auto polyhedralSource = config->getPolyhedralSource();
        const auto density = config->getDensity();
        const auto computationPoints = config->getPointsOfInterest();
        const auto outputFileName = config->getOutputFileName();
        const auto metricUnit = config->getMeshUnit();
        const PolyhedronIntegrity checkPolyhedralInput = config->getMeshInputCheckStatus() ? PolyhedronIntegrity::HEAL : PolyhedronIntegrity::DISABLE;

        POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedron creation and check (if enabled) started.");
        const auto startPolyhedron = std::chrono::high_resolution_clock::now();
        const Polyhedron polyhedron{polyhedralSource, density, NormalOrientation::OUTWARDS, checkPolyhedralInput, metricUnit};
        const auto endPolyhedron = std::chrono::high_resolution_clock::now();
        const auto durationPolyhedron = endPolyhedron - startPolyhedron;
        const auto msPolyhedron = std::chrono::duration_cast<std::chrono::microseconds>(durationPolyhedron).count();
        POLYHEDRAL_GRAVITY_LOG_INFO("Polyhedron instantiated and checked. It took {} microseconds.", msPolyhedron);

        POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");
        POLYHEDRAL_GRAVITY_LOG_INFO("Number of Vertices:                               {}", polyhedron.countVertices());
        POLYHEDRAL_GRAVITY_LOG_INFO("Number of Faces:                                  {}", polyhedron.countFaces());
        POLYHEDRAL_GRAVITY_LOG_INFO("Number of Computation Points:                     {}", computationPoints.size());
        POLYHEDRAL_GRAVITY_LOG_INFO("Mesh Check Enabled:                               {}", config->getMeshInputCheckStatus());
        POLYHEDRAL_GRAVITY_LOG_INFO("Mesh Unit:                                        {}", polyhedron.getMeshUnitAsString());
        POLYHEDRAL_GRAVITY_LOG_INFO("Density:                                          {} {}", polyhedron.getDensity(), polyhedron.getDensityUnit());
        POLYHEDRAL_GRAVITY_LOG_INFO("Output File:                                      {}", outputFileName);
        POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");

        POLYHEDRAL_GRAVITY_LOG_INFO("Gravity Evaluation has started!");
        const auto startCalc = std::chrono::high_resolution_clock::now();

        const auto result = GravityModel::evaluate(polyhedron, computationPoints, true);

        const auto endCalc = std::chrono::high_resolution_clock::now();
        const auto durationCalc = endCalc - startCalc;
        const auto msCalc = std::chrono::duration_cast<std::chrono::microseconds>(durationCalc).count();
        POLYHEDRAL_GRAVITY_LOG_INFO("The calculation of the Gravity Model has finished. It took {} microseconds or on average {} microseconds/point",
            msCalc, static_cast<double>(msCalc) / static_cast<double>(computationPoints.size()));
        POLYHEDRAL_GRAVITY_LOG_INFO("####################################################################################");

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