#include "YAMLConfigReader.h"

#include "MeshReader.h"

namespace polyhedralGravity {

    std::string YAMLConfigReader::getOutputFileName() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the output filename from the configuration file.");
        if (_file[ROOT][OUTPUT] && _file[ROOT][OUTPUT][OUTPUT_FILENAME]) {
            return _file[ROOT][OUTPUT][OUTPUT_FILENAME].as<std::string>();
        } else {
            return "";
        }
    }

    double YAMLConfigReader::getDensity() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the density from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_DENSITY]) {
            return _file[ROOT][INPUT][INPUT_DENSITY].as<double>();
        } else {
            throw std::runtime_error{"There happened an error parsing the density from the YAML config file!"};
        }
    }

    std::vector<std::array<double, 3>> YAMLConfigReader::getPointsOfInterest() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the computation points from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_POINTS]) {
            return _file[ROOT][INPUT][INPUT_POINTS].as<std::vector<std::array<double, 3>>>();
        } else {
            throw std::runtime_error{"There happened an error parsing the points of interest from the YAML config file!"};
        }
    }

    bool YAMLConfigReader::getMeshInputCheckStatus() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the activation of the input mesh sanity check from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_CHECK]) {
            return _file[ROOT][INPUT][INPUT_CHECK].as<bool>();
        } else {
            return true;
        }
    }

    std::tuple<std::vector<Array3>, std::vector<IndexArray3>> YAMLConfigReader::getPolyhedralSource() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the data sources (file names) from the configuration file.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_POLYHEDRON]) {
            const auto vectorOfFiles = _file[ROOT][INPUT][INPUT_POLYHEDRON].as<std::vector<std::string>>();
            return MeshReader::getPolyhedralSource(vectorOfFiles);
        } else {
            throw std::runtime_error{"There happened an error parsing the DataSource of the Polyhedron from the config file"};
        }
    }

    MetricUnit YAMLConfigReader::getMeshUnit() {
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Reading the unit of the polyhedral mesh.");
        if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_METRIC_UNIT]) {
            const auto unit = _file[ROOT][INPUT][INPUT_METRIC_UNIT].as<std::string>();
            return readMetricUnit(unit);
        } else {
            return MetricUnit::METER;
        }
    }


}