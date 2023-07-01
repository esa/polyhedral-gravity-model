#pragma once

#include "thrust/transform.h"
#include "polyhedralGravity/calculation/GravityModel.h"

#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/Polyhedron.h"

// TODO Beautify, At the moment just quickly prototyped

namespace polyhedralGravity {

class GravityEvaluable {


    const Polyhedron _polyhedron;

    const double _density;

    mutable std::vector<Array3Triplet> _segmentVectors;

    mutable std::vector<Array3> _planeUnitNormal;

    mutable std::vector<Array3Triplet> _segmentUnitNormals;



public:

    GravityEvaluable(const Polyhedron &polyhedron, double density) : _polyhedron{polyhedron}, _density{density} {
        auto polyhedronIterator = GravityModel::transformPolyhedron(polyhedron, {0.0, 0.0, 0.0});
        _segmentVectors.resize(polyhedron.countFaces());
        thrust::transform(
                polyhedronIterator.first,
                polyhedronIterator.second,
                _segmentVectors.begin(),
                [](const Array3Triplet &face) {
                    return GravityModel::detail::buildVectorsOfSegments(face[0], face[1], face[2]);
                }
        );
        // ...
    }

    inline GravityModelResult operator()(const Array3 &position) const {
        return GravityModel::evaluate(_polyhedron, _density, position);
    }

};


inline GravityEvaluable prepare(const Polyhedron &polyhedron, double density) {
    return {polyhedron, density};
}

} // namespace polyhedralGravity