#pragma once

#include <tuple>
#include <variant>
#include <string>
#include <optional>
#include <sstream>

#include "thrust/transform.h"
#include "thrust/execution_policy.h"

#include "GravityModelDetail.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/input/TetgenAdapter.h"
#include "GravityModelData.h"
#include "Polyhedron.h"


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

        /** Cache for the segment vectors (segments between vertices of a polyhedral face) */
        mutable std::vector<Array3Triplet> _segmentVectors{};

        /** Cache for the plane unit normals (unit normals of the polyhedral faces) */
        mutable std::vector<Array3> _planeUnitNormals{};

        /** Cache for the segment unit normals (unit normals of each the polyhedral faces' segments) */
        mutable std::vector<Array3Triplet> _segmentUnitNormals{};

    public:
        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron.
         * In contrast to the {@link GravityModel::evaluate}, this evaluate method on the {@link GravityEvaluable}
         * caches intermediate results and input data and subsequent evaluations will be faster.
         *
         * @param polyhedron the constant density polyhedron
         */
        explicit GravityEvaluable(const Polyhedron &polyhedron) :
            _polyhedron{polyhedron} {
            this->prepare();
        }

        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron and caches.
         * This is for restoring a GravityEvaluable from a previous state.
         * @param polyhedron the polyhedron
         * @param segmentVectors the segment vectors
         * @param planeUnitNormals the plane unit normals
         * @param segmentUnitNormals the segment unit normals
         */
        GravityEvaluable(const Polyhedron &polyhedron,
                         const std::vector<Array3Triplet> &segmentVectors,
                         const std::vector<Array3> &planeUnitNormals,
                         const std::vector<Array3Triplet> &segmentUnitNormals) :
            _polyhedron{polyhedron},
            _segmentVectors{segmentVectors},
            _planeUnitNormals{planeUnitNormals},
            _segmentUnitNormals{
                    segmentUnitNormals} {
        }

        /**
         * Evaluates the polyhedral gravity model for a given constant density polyhedron at computation
         * point P. Wrapper for evaluate<parallelization>.
         *
         * The results' units depend on the polyhedron's input units.
         * For example, if the polyhedral mesh is in @f$[m]@f$ and the density in @f$[kg/m^3]@f$, then the potential is in @f$[m^2/s^2]@f$.
         * In case the polyhedron is unitless, the results are not multiplied with the Gravitational Constant @f$G@f$, but returned raw.
         *
         * @param computationPoints the computation point P or multiple computation points in a vector
         * @param parallelization if true, the calculation is parallelized
         * @return the GravityModelResult containing the potential, acceleration, and second derivative
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
        [[nodiscard]] std::string toString() const;

        /**
         * Returns the output units of the GravityEvaluable in order potential, acceleration, second derivative tensor.
         * This depends on the units chosen for the polyhedron.
         * @return human-readable output units of the GravityEvaluable
         */
        [[nodiscard]] std::array<std::string, 3> getOutputMetricUnit() const;


        /**
         * Returns the polyhedron, the density, and the internal caches.
         *
         * @return tuple of polyhedron, density, segmentVectors, planeUnitNormals, and segmentUnitNormals
         */
        std::tuple<Polyhedron, std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>> getState() const;

    private:

        /**
         * Prepares the polyhedron for the evaluation by calculating the segment vectors, the plane unit normals,
         * and the segment unit normals.
         * Called by the constructor once.
         */
        void prepare() const;

        /**
        * Evaluates the polyhedral gravity model for a given constant density polyhedron at computation
        * point P.
        * @tparam Parallelization if true, the calculation is parallelized
        * @param computationPoint the computation Point P
        * @return the GravityModelResult containing the potential, the acceleration, and the change of acceleration
        * at computation Point P
        */
        template<bool Parallelization = true>
        [[nodiscard]] GravityModelResult evaluate(const Array3 &computationPoint) const;


        /**
         * Evaluates the polyhedral gravity model for a given constant density polyhedron at computation
         * at multiple computation points.
         * @tparam Parallelization if true, the calculation is parallelized
         * @param computationPoints the computation Points
         * @return vector of GravityModelResults containing the potential, the acceleration, and the change of acceleration
         */
        template<bool Parallelization = true>
        [[nodiscard]] std::vector<GravityModelResult> evaluate(const std::vector<Array3> &computationPoints) const;

        /**
         * Evaluates the polyhedral gravity model for a given constant density polyhedron at computation a certain face.
         * @param tuple consisting of face, segmentVectors, planeUnitNormal, and segmentUnitNormals
         * @return the GravityModelResult containing the potential, the acceleration, and the change of acceleration which
         * this face contributes to the computation point
         */
        static GravityModelResult
        evaluateFace(const thrust::tuple<Array3Triplet, Array3Triplet, Array3, Array3Triplet> &tuple);

    };

}// namespace polyhedralGravity