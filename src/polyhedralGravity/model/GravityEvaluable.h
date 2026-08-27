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

    namespace detail {
        /**
         * The engine behind a {@link GravityEvaluable}, i.e. the mesh with its caches and the kernels.
         * Forward declared because its implementation is templated over the floating point precision and
         * the execution space, both of which the user selects at runtime.
         */
        class EvaluationEngine;
    }// namespace detail

    /**
     * Class for evaluating the polyhedrale gravity model for a given constant density polyhedron.
     * Caches the polyhedron and data which is independent of the computation point P.
     * Provides an operator() for evaluating the polyhedrale gravity model for a given constant density polyhedron
     * at computation point P.
     *
     * The compute backend and the floating point precision are fixed when the GravityEvaluable is created,
     * not when it is called. That is what makes the caching worthwhile: the mesh and the caches of Tsoulis'
     * algorithm are allocated exactly once, in the memory space the chosen backend computes in, and stay
     * there for the whole lifetime of the object. A GravityEvaluable created for
     * {@link ComputeBackend::GPU_PARALLEL} therefore never puts its caches into host memory at all. Use one
     * GravityEvaluable per backend if the same polyhedron is to be evaluated on several of them.
     */
    class GravityEvaluable {

        /** The constant density polyhedron consisting of vertices and triangular faces */
        const Polyhedron _polyhedron;

        /** The compute backend every evaluation of this GravityEvaluable runs on */
        ComputeBackend _backend;

        /** The floating point precision the kernels compute in */
        ComputePrecision _precision;

        /**
         * The mesh and the caches which only depend on the polyhedron (the segment vectors, the plane unit
         * normals, and the segment unit normals) in the memory space of {@link _backend}, together with the
         * kernels evaluating them there.
         *
         * The mesh itself belongs to the polyhedron above and is shared with it rather than uploaded again
         * whenever the memory spaces and the precisions match; only the caches are always allocated here.
         * Held by shared pointer so that a GravityEvaluable stays copyable and so that copies share the
         * caches instead of recomputing them.
         */
        std::shared_ptr<detail::EvaluationEngine> _engine;

    public:
        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron.
         * In contrast to the {@link GravityModel::evaluate}, this evaluate method on the {@link GravityEvaluable}
         * caches intermediate results and input data and subsequent evaluations will be faster.
         *
         * @param polyhedron the constant density polyhedron
         * @param backend the compute backend every evaluation runs on (default: CPU_PARALLEL)
         * @param precision the floating point precision the evaluation computes in (default: FLOAT64)
         *
         * @throws std::invalid_argument if the polyhedron has no faces
         * @throws std::runtime_error if GPU_PARALLEL is requested, but this build has no GPU backend
         *
         * @note FLOAT32 roughly halves the memory traffic of the evaluation, but the polyhedral gravity model
         * cancels large terms against each other, so expect a relative accuracy of only about @f$10^{-4}@f$.
         */
        explicit GravityEvaluable(const Polyhedron &polyhedron,
                                  ComputeBackend backend = ComputeBackend::CPU_PARALLEL,
                                  ComputePrecision precision = ComputePrecision::FLOAT64);

        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron and caches.
         * This is for restoring a GravityEvaluable from a previous state.
         * @param polyhedron the polyhedron
         * @param segmentVectors the segment vectors
         * @param planeUnitNormals the plane unit normals
         * @param segmentUnitNormals the segment unit normals
         * @param backend the compute backend every evaluation runs on (default: CPU_PARALLEL)
         * @param precision the floating point precision the evaluation computes in (default: FLOAT64)
         *
         * @throws std::invalid_argument if the polyhedron has no faces or a cache's size does not match
         * @throws std::runtime_error if GPU_PARALLEL is requested, but this build has no GPU backend
         */
        GravityEvaluable(const Polyhedron &polyhedron,
                         const std::vector<Array3Triplet> &segmentVectors,
                         const std::vector<Array3> &planeUnitNormals,
                         const std::vector<Array3Triplet> &segmentUnitNormals,
                         ComputeBackend backend = ComputeBackend::CPU_PARALLEL,
                         ComputePrecision precision = ComputePrecision::FLOAT64);

        /**
         * Evaluates the polyhedral gravity model for a given constant density polyhedron at computation
         * point P, on the compute backend this GravityEvaluable was created for.
         *
         * The results' units depend on the polyhedron's input units.
         * For example, if the polyhedral mesh is in @f$[m]@f$ and the density in @f$[kg/m^3]@f$, then the potential is in @f$[m^2/s^2]@f$.
         * In case the polyhedron is unitless, the results are not multiplied with the Gravitational Constant @f$G@f$, but returned raw.
         *
         * @param computationPoints the computation point P or multiple computation points in a vector
         * @return the GravityModelResult containing the potential, acceleration, and second derivative
         */
        [[nodiscard]] std::variant<GravityModelResult, std::vector<GravityModelResult>>
        operator()(const std::variant<Array3, std::vector<Array3>> &computationPoints) const;

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
         * Returns the compute backend every evaluation of this GravityEvaluable runs on.
         * @return CPU_SERIAL, CPU_PARALLEL, or GPU_PARALLEL
         */
        [[nodiscard]] ComputeBackend getComputeBackend() const;

        /**
         * Returns the floating point precision the evaluation computes in.
         * @return FLOAT32 or FLOAT64
         */
        [[nodiscard]] ComputePrecision getComputePrecision() const;
    };

}// namespace polyhedralGravity
