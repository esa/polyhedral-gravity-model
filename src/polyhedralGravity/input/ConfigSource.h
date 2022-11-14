#pragma once

#include <vector>
#include <memory>
#include <string>
#include "DataSource.h"

namespace polyhedralGravity {

/**
 * Interface defining methods giving the calculation some input.
 * This includes the points of interest as well as the source of the data.
 */
    class ConfigSource {

    public:

        virtual ~ConfigSource() = default;

        /**
         * Returns the specified output file name.
         * @return a std::string with the filename
         */
        virtual std::string getOutputFileName() = 0;

        /**
         * Returns the constant density rho of the given polyhedron in [kg/m^3].
         * The density is required for the calculation.
         * @return density as double
         */
        virtual double getDensity() = 0;

        /**
         * The vector contains the points for which the polyhedral gravity model should
         * be evaluated.
         * @return vector of points
         */
        virtual std::vector<std::array<double, 3>> getPointsOfInterest() = 0;

        /**
         * Returns the activation status of the input polyhedron mesh sanity check.
         * @return true if enabled
         */
        virtual bool getMeshInputCheckStatus() = 0;

        /**
         * The DataSource of the given Polyhedron.
         * @return data source (e. g. a file reader)
         */
        virtual std::shared_ptr<DataSource> getDataSource() = 0;

    };

}
