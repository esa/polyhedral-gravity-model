#pragma once

#include <utility>

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "thrust/iterator/transform_iterator.h"

namespace polyhedralGravity {

    /**
     * An iterator transforming the polyhedron's coordinates on demand by a given offset.
     * This function returns a pair of transform iterators (first --> begin(), second --> end()).
     * @param polyhedron - reference to the polyhedron
     * @param offset - the offset to apply
     * @return pair of transform iterators
     */
    inline auto transformPolyhedron(const Polyhedron &polyhedron, const Array3 &offset = {0.0, 0.0, 0.0}) {
        //The offset is captured by value to ensure its lifetime!
        const auto lambdaOffsetApplication = [&polyhedron, offset](const IndexArray3 &face) -> Array3Triplet {
            using namespace util;
            return {polyhedron.getVertex(face[0]) - offset, polyhedron.getVertex(face[1]) - offset,
                    polyhedron.getVertex(face[2]) - offset};
        };
        auto first = thrust::make_transform_iterator(polyhedron.getFaces().begin(), lambdaOffsetApplication);
        auto last = thrust::make_transform_iterator(polyhedron.getFaces().end(), lambdaOffsetApplication);
        return std::make_pair(first, last);
    }
}