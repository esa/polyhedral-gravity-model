#pragma once

#include <tuple>
#include <variant>
#include <string>
#include <sstream>

#include "thrust/transform.h"
#include "polyhedralGravity/calculation/GravityModel.h"

#include "polyhedralGravity/input/TetgenAdapter.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "thrust/execution_policy.h"


namespace polyhedralGravity {

    /**
     * Class for evaluating the polyhedrale gravity model for a given constant density polyhedron.
     * Caches the polyhedron and data which is independent of the computation point P.
     * Provides an operator() for evaluating the polyhedrale gravity model for a given constant density polyhedron
     * at computation point P and choosing between parallel and serial evaluation.
     */
    class GravityEvaluable {

        /** The constant density polyhedron consisting of vertices and triangular faces */
        const Polyhedron _polyhedron;

        /** The constant density of the polyhedron in [kg/m^3] */
        const double _density;

        /** Cache for the segment vectors (segments between vertices of a polyhedral face) */
        mutable std::vector<Array3Triplet> _segmentVectors;

        /** Cache for the plane unit normals (unit normals of the polyhedral faces) */
        mutable std::vector<Array3> _planeUnitNormals;

        /** Cache for the segment unit normals (unit normals of each the polyhedral faces' segments) */
        mutable std::vector<Array3Triplet> _segmentUnitNormals;


    public:

        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron.
         * @param polyhedron - the polyhedron
         * @param density - the constant density in [kg/m^3]
         *
         * @note This is a separate constructor since the polyhedron as a class it not exposed to the user via
         * the Python Interface. Thus, this constructor is only available via the C++ interface.
         */
        GravityEvaluable(const Polyhedron &polyhedron, double density) : _polyhedron{polyhedron}, _density{density} {
            this->prepare();
        }

        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron.
         * @param polyhedralSource - the vertices & faces of the polyhedron as tuple or the filenames of the files
         * @param density - the constant density in [kg/m^3]
         */
        GravityEvaluable(
                const PolyhedralSource &polyhedralSource,
                double density) : _polyhedron{
                std::holds_alternative<std::tuple<std::vector<Array3>, std::vector<IndexArray3>>>(polyhedralSource)
                ? Polyhedron{std::get<std::tuple<std::vector<Array3>, std::vector<IndexArray3>>>(polyhedralSource)}
                : TetgenAdapter(std::get<std::vector<std::string>>(polyhedralSource)).getPolyhedron()
        }, _density{density} {
            this->prepare();
        }

        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
         * point P. Wrapper for evaluate<parallelization>.
         * @param computationPoints - the computation point P or multiple computation points in a vector
         * @param parallelization - if true, the calculation is parallelized
         * @return the GravityModelResult containing the potential, acceleration and second derivative
         */
        inline std::variant<GravityModelResult, std::vector<GravityModelResult>>
        operator()(const std::variant<Array3, std::vector<Array3>> &computationPoints,
                   bool parallelization = true) const {
            if (parallelization) {
                if (std::holds_alternative<Array3>(computationPoints)) {
                    return this->evaluate<true>(std::get<Array3>(computationPoints));
                } else {
                    return this->evaluate<true>(std::get<std::vector<Array3>>(computationPoints));
                }
            } else {
                if (std::holds_alternative<Array3>(computationPoints)) {
                    return this->evaluate<false>(std::get<Array3>(computationPoints));
                } else {
                    return this->evaluate<false>(std::get<std::vector<Array3>>(computationPoints));
                }
            }
        }

        /**
         * Returns a string representation of the GravityEvaluable.
         * @return string representation of the GravityEvaluable
         */
        std::string toString() const;

    private:

        /**
         * Prepares the polyhedron for the evaluation by calculating the segment vectors, the plane unit normals
         * and the segment unit normals.
         * Called by the constructor once.
         */
        void prepare() const;

        /**
        * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
        * point P.
        * @tparam Parallelization - if true, the calculation is parallelized
        * @param computationPoint - the computation Point P
        * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration
        * at computation Point P
        */
        template<bool Parallelization = true>
        GravityModelResult evaluate(const Array3 &computationPoint) const;

        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
         * at multiple computation points.
         * @tparam Parallelization - if true, the calculation is parallelized
         * @param computationPoints - the computation Points
         * @return vector of GravityModelResults containing the potential, the acceleration and the change of acceleration
         */
        template<bool Parallelization = true>
        std::vector<GravityModelResult> evaluate(const std::vector<Array3> &computationPoints) const;


        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation a certain face.
         * @param tuple - consisting of face, segmentVectors, planeUnitNormal and segmentUnitNormals
         * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration which
         * this face contributes to the computation point
         */
        static GravityModelResult
        evaluateFace(const thrust::tuple<Array3Triplet, Array3Triplet, Array3, Array3Triplet> &tuple);

    };

} // namespace polyhedralGravity