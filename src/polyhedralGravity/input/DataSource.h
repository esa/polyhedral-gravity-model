#pragma once

#include "polyhedralGravity/model/Polyhedron.h"

namespace polyhedralGravity {

/**
 * Interface consisting of a method which returns a polyhedron.
 */
    class DataSource {

    public:

        /**
         * Default destructor of DataSource
         */
        virtual ~DataSource() = default;

        /**
         * Returns a Polyhedron from the underlying source.
         * @return a polyhedron
         */
        virtual Polyhedron getPolyhedron() = 0;

    };

}
