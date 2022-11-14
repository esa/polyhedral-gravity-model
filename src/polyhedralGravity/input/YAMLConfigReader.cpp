#include "YAMLConfigReader.h"

namespace polyhedralGravity {

    std::string YAMLConfigReader::getOutputFileName() {
        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger() ,
                            "Reading the output filename from the configuration file.");
        if (_file[ROOT][OUTPUT] && _file[ROOT][OUTPUT][OUTPUT_FILENAME]) {
            return _file[ROOT][OUTPUT][OUTPUT_FILENAME].as<std::string>();
        } else {
            return "";
        }
    }

    double YAMLConfigReader::getDensity() {
        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger() ,
                            "Reading the density from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_DENSITY]) {
            return _file[ROOT][INPUT][INPUT_DENSITY].as<double>();
        } else {
            throw std::runtime_error{"There happened an error parsing the density from the YAML config file!"};
        }
    }

    std::vector<std::array<double, 3>> YAMLConfigReader::getPointsOfInterest() {
        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger() ,
                            "Reading the computation points from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_POINTS]) {
            return _file[ROOT][INPUT][INPUT_POINTS].as<std::vector<std::array<double, 3>>>();
        } else {
            throw std::runtime_error{
                    "There happened an error parsing the points of interest from the YAML config file!"};
        }
    }

    bool YAMLConfigReader::getMeshInputCheckStatus() {
        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger() ,
                            "Reading the activation of the input mesh sanity check from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_CHECK]) {
            return _file[ROOT][INPUT][INPUT_CHECK].as<bool>();
        } else {
            return false;
        }
    }

    std::shared_ptr<DataSource> YAMLConfigReader::getDataSource() {
        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger() ,
                            "Reading the data sources (file names) from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_POLYHEDRON]) {
            auto vectorOfFiles = _file[ROOT][INPUT][INPUT_POLYHEDRON].as<std::vector<std::string>>();
            return std::make_shared<TetgenAdapter>(vectorOfFiles);
        } else {
            throw std::runtime_error{
                    "There happened an error parsing the DataSource of the Polyhedron from the config file"};
        }
    }

}