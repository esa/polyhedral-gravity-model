#pragma once

#include <memory>
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

    namespace opencl {
        /* Forward declaration so that the OpenCL headers do not leak into this public header */
        class OpenCLEvaluation;
    }// namespace opencl

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

        /**
         * The backend actually in use. This is the backend requested in the constructor unless
         * OpenCL was requested but is unavailable, in which case it is {@link ComputeBackend::CPU}.
         */
        ComputeBackend _backend;

        /** The floating point precision requested for the device-side evaluation */
        ComputePrecision _precision;

        /**
         * The OpenCL engine, or nullptr if the evaluation runs on the host.
         * Held by shared pointer so that a GravityEvaluable stays copyable and copies share the
         * device-resident polyhedron instead of uploading it again.
         */
        std::shared_ptr<opencl::OpenCLEvaluation> _openCLEvaluation{};

        /**
         * Whether the host-side caches above have been filled.
         * On the OpenCL backend they are only needed when the state is read out, so filling them is
         * deferred rather than duplicating work the device already does.
         */
        mutable bool _prepared{false};

    public:
        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron.
         * In contrast to the {@link GravityModel::evaluate}, this evaluate method on the {@link GravityEvaluable}
         * caches intermediate results and input data and subsequent evaluations will be faster.
         *
         * @param polyhedron the constant density polyhedron
         * @param backend the compute backend to evaluate on (default: OPENCL)
         * @param precision the floating point precision the backend computes in (default: FLOAT64)
         *
         * @note If OpenCL is requested but the library was compiled without it, or no OpenCL device
         * supporting the requested precision exists, the evaluation silently falls back to
         * {@link ComputeBackend::CPU}. Query {@link getComputeBackend} to learn what is actually used.
         */
        explicit GravityEvaluable(const Polyhedron &polyhedron,
                                  ComputeBackend backend = ComputeBackend::OPENCL,
                                  ComputePrecision precision = ComputePrecision::FLOAT64);

        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron and caches.
         * This is for restoring a GravityEvaluable from a previous state.
         * @param polyhedron the polyhedron
         * @param segmentVectors the segment vectors
         * @param planeUnitNormals the plane unit normals
         * @param segmentUnitNormals the segment unit normals
         * @param backend the compute backend to evaluate on (default: OPENCL)
         * @param precision the floating point precision the backend computes in (default: FLOAT64)
         */
        GravityEvaluable(const Polyhedron &polyhedron,
                         const std::vector<Array3Triplet> &segmentVectors,
                         const std::vector<Array3> &planeUnitNormals,
                         const std::vector<Array3Triplet> &segmentUnitNormals,
                         ComputeBackend backend = ComputeBackend::OPENCL,
                         ComputePrecision precision = ComputePrecision::FLOAT64);

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

        /**
         * Returns the backend this GravityEvaluable actually evaluates on.
         * This differs from the backend requested in the constructor if OpenCL was requested but is
         * unavailable, in which case {@link ComputeBackend::CPU} is returned.
         * @return the effective compute backend
         */
        [[nodiscard]] ComputeBackend getComputeBackend() const;

        /**
         * Returns the floating point precision the device-side evaluation uses.
         * It is meaningless for {@link ComputeBackend::CPU}, which always computes in double precision.
         * @return the compute precision
         */
        [[nodiscard]] ComputePrecision getComputePrecision() const;

    private:

        /**
         * Prepares the polyhedron for the evaluation by calculating the segment vectors, the plane unit normals,
         * and the segment unit normals.
         * Called once, either by the constructor or, on the OpenCL backend, by {@link ensurePrepared}.
         */
        void prepare() const;

        /**
         * Fills the host-side caches if they are not filled yet.
         * On the OpenCL backend the device derives these values itself, so the host-side caches are
         * only computed once something actually reads them.
         */
        void ensurePrepared() const;

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