#pragma once

#include <array>
#include <memory>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include "GravityModelData.h"
#include "Polyhedron.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"

namespace polyhedralGravity {

    namespace kokkos {
        /* Forward declaration so that the Kokkos headers do not leak into this public header */
        class KokkosEvaluationBase;
    }// namespace kokkos

    /**
     * Class for evaluating the polyhedrale gravity model for a given constant density polyhedron.
     * Caches the polyhedron and data which is independent of the computation point P.
     * Provides an operator() for evaluating the polyhedrale gravity model for a given constant density polyhedron
     * at computation point P and choosing the compute backend the evaluation runs on.
     */
    class GravityEvaluable {

        /** The constant density polyhedron consisting of vertices and triangular faces */
        const Polyhedron _polyhedron;

        /** The floating point precision the kernels compute in */
        ComputePrecision _precision;

        /**
         * The polyhedron as it lives in the memory of the compute devices, together with the caches which only
         * depend on it (the segment vectors, the plane unit normals, and the segment unit normals).
         * Held by shared pointer so that a GravityEvaluable stays copyable and so that copies share the
         * uploaded polyhedron instead of uploading it again.
         */
        std::shared_ptr<kokkos::KokkosEvaluationBase> _evaluation;

    public:
        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron.
         * In contrast to the {@link GravityModel::evaluate}, this evaluate method on the {@link GravityEvaluable}
         * caches intermediate results and input data and subsequent evaluations will be faster.
         *
         * @param polyhedron the constant density polyhedron
         * @param precision the floating point precision the evaluation computes in (default: FLOAT64)
         *
         * @throws std::invalid_argument if the polyhedron has no faces
         *
         * @note FLOAT32 roughly halves the memory traffic of the evaluation, but the polyhedral gravity model
         * cancels large terms against each other, so expect a relative accuracy of only about @f$10^{-4}@f$.
         */
        explicit GravityEvaluable(const Polyhedron &polyhedron,
                                  ComputePrecision precision = ComputePrecision::FLOAT64);

        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron and caches.
         * This is for restoring a GravityEvaluable from a previous state.
         * @param polyhedron the polyhedron
         * @param segmentVectors the segment vectors
         * @param planeUnitNormals the plane unit normals
         * @param segmentUnitNormals the segment unit normals
         * @param precision the floating point precision the evaluation computes in (default: FLOAT64)
         *
         * @throws std::invalid_argument if the polyhedron has no faces or a cache's size does not match
         */
        GravityEvaluable(const Polyhedron &polyhedron,
                         const std::vector<Array3Triplet> &segmentVectors,
                         const std::vector<Array3> &planeUnitNormals,
                         const std::vector<Array3Triplet> &segmentUnitNormals,
                         ComputePrecision precision = ComputePrecision::FLOAT64);

        /**
         * Evaluates the polyhedral gravity model for a given constant density polyhedron at computation
         * point P.
         *
         * The results' units depend on the polyhedron's input units.
         * For example, if the polyhedral mesh is in @f$[m]@f$ and the density in @f$[kg/m^3]@f$, then the potential is in @f$[m^2/s^2]@f$.
         * In case the polyhedron is unitless, the results are not multiplied with the Gravitational Constant @f$G@f$, but returned raw.
         *
         * @param computationPoints the computation point P or multiple computation points in a vector
         * @param backend the compute backend the evaluation runs on (default: CPU_PARALLEL)
         * @return the GravityModelResult containing the potential, acceleration, and second derivative
         *
         * @throws std::runtime_error if GPU_PARALLEL is requested, but this build has no GPU backend
         */
        [[nodiscard]] std::variant<GravityModelResult, std::vector<GravityModelResult>>
        operator()(const std::variant<Array3, std::vector<Array3>> &computationPoints,
                   ComputeBackend backend = ComputeBackend::CPU_PARALLEL) const;

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
         * @return tuple of polyhedron, segmentVectors, planeUnitNormals, and segmentUnitNormals
         */
        [[nodiscard]] std::tuple<Polyhedron, std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
        getState() const;

        /**
         * Returns the floating point precision the evaluation computes in.
         * @return FLOAT32 or FLOAT64
         */
        [[nodiscard]] ComputePrecision getComputePrecision() const;
    };

}// namespace polyhedralGravity
