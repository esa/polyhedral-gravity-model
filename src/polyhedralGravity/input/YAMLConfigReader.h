#pragma once

#include "ConfigSource.h"
#include "TetgenAdapter.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "yaml-cpp/yaml.h"
#include <array>
#include <string>
#include <tuple>
#include <vector>

namespace polyhedralGravity {

    /**
     * A class to read configuration parameters from a YAML file.
     * Provides methods to fetch specific configuration data such as density,
     * computation points, mesh check status, and output file details.
     */
    class YAMLConfigReader final : public ConfigSource {

        /*
         * The following static variables contain the names of the YAML nodes.
         * Note: C++20 would allow constexpr std::string which would be more appropriate instead of char[]
         */
        static constexpr char ROOT[] = "gravityModel";
        static constexpr char INPUT[] = "input";
        static constexpr char INPUT_POLYHEDRON[] = "polyhedron";
        static constexpr char INPUT_DENSITY[] = "density";
        static constexpr char INPUT_POINTS[] = "points";
        static constexpr char INPUT_CHECK[] = "check_mesh";
        static constexpr char INPUT_METRIC_UNIT[] = "metric_unit";
        static constexpr char OUTPUT[] = "output";
        static constexpr char OUTPUT_FILENAME[] = "filename";


        /**
         * The member administering the YAML file/ Connection to yaml-cpp
         */
        const YAML::Node _file;

    public:

        /**
         * Creates a new YAML Config Reader.
         * @param filename a reference to a string
         * @throws an exception if the file is malformed or cannot be loaded or if the ROOT node is not found
         */
        explicit YAMLConfigReader(const std::string &filename)
                : _file{YAML::LoadFile(filename)} {
            if (!_file[ROOT]) {
                throw std::runtime_error{"The YAML file does not contain a specification for the \"gravityModel\"!"};
            }
        }

        /**
         * Creates a new YAML Config Reader.
         * @param filename a movable string
         * @throws an exception if the file is malformed or cannot be loaded or if the ROOT node is not found
         */
        explicit YAMLConfigReader(std::string &&filename)
                : _file{YAML::LoadFile(filename)} {
            if (!_file[ROOT]) {
                throw std::runtime_error{"The YAML file does not contain a specification for the \"gravityModel\"!"};
            }
        }

        /**
         * Returns the specified output filename.
         * If none is specified an empty string will be returned.
         * @return a filename a std::string or an empty string if none is specified
         */
        std::string getOutputFileName() override;

        /**
         * Reads the density from the yaml configuration file.
         * @return density as double
         */
        double getDensity() override;

        /**
         * Reads the computation points from the yaml configuration file.
         * @return vector of computation points
         */
        std::vector<std::array<double, 3>> getPointsOfInterest() override;

        /**
         * Reads the enablement of the input sanity check from the yaml file.
         * @return true or false if explicitly enabled, otherwise per-default true
         */
        bool getMeshInputCheckStatus() override;

        /**
         * Reads the DataSource from the yaml configuration file.
         * @return shared_ptr to the DataSource Object created
         */
        std::tuple<std::vector<Array3>, std::vector<IndexArray3>> getPolyhedralSource() override;


        /**
         * Retrieves the metric unit of the mesh as specified in the configuration.
         * Can return meter, kilometer, or unitless depending on the configuration data.
         * @return the metric unit of the mesh as an instance of MetricUnit (if not present defaults to meter)
         */
        MetricUnit getMeshUnit() override;


    };

}
