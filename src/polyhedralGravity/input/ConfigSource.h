#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <array>
#include "polyhedralGravity/model/PolyhedronDefinitions.h"

namespace polyhedralGravity {

/**
 * Interface defining methods giving the calculation some input.
 * This includes the points of interest as well as the source of the data.
 */
    class ConfigSource {

    public:

        /** Default Virtual Destructor */
        virtual ~ConfigSource() = default;

        /**
         * Returns the specified output file name.
         * @return a std::string with the filename
         */
        virtual std::string getOutputFileName() = 0;

        /**
         * Returns the constant density rho of the given polyhedron.
         * The unit must match to the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$.
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
         * The source structure of a Polyhedron consisting of vertices and faces
         * @return a polyhedral source consisting of vertices and faces
         */
        virtual PolyhedralSource getPolyhedralSource() = 0;

        /**
         * Returns the unit of measurement used for the polyhedron mesh.
         * This specifies the physical unit in which the mesh dimensions are defined.
         * Could be metric units such as meters, kilometers, or unitless.
         * @return the unit of the polyhedral mesh as a MetricUnit enum
         */
        virtual MetricUnit getMeshUnit() = 0;

    };

}
