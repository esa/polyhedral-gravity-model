#pragma once

#include <memory>
#include <tuple>
#include <vector>

#include "polyhedralGravity/model/PolyhedronDefinitions.h"

namespace polyhedralGravity {

    /* Forward declaration so that this header stays free of the Polyhedron's includes */
    class Polyhedron;

}// namespace polyhedralGravity

namespace polyhedralGravity::kokkos {

    /**
     * The precision-agnostic interface to the device-resident polyhedron and its kernels.
     *
     * A {@link GravityEvaluable} holds one of these. The concrete implementation behind it is templated over
     * the floating point precision, which is why the interface is virtual: the precision is a runtime choice
     * of the user (see {@link ComputePrecision}), while the kernels have to be compiled for both precisions.
     *
     * @note This header deliberately does not include any Kokkos header, so that the public headers of this
     * library do not drag Kokkos into every translation unit of a user of this library.
     */
    class KokkosEvaluationBase {
    public:
        virtual ~KokkosEvaluationBase() = default;

        /**
         * Evaluates the polyhedral gravity model at a single computation point.
         * @param computationPoint the computation point P
         * @param backend the compute backend to run the kernel on
         * @return the potential, the acceleration, and the gradiometric tensor at P
         * @throws std::runtime_error if the backend is not available in this build
         */
        [[nodiscard]] virtual GravityModelResult evaluate(const Array3 &computationPoint,
                                                          ComputeBackend backend) const = 0;

        /**
         * Evaluates the polyhedral gravity model at multiple computation points.
         * The points are evaluated in one kernel launch, with one team of threads per computation point.
         * @param computationPoints the computation points
         * @param backend the compute backend to run the kernel on
         * @return the results foreach computation point, in the order of the input
         * @throws std::runtime_error if the backend is not available in this build
         */
        [[nodiscard]] virtual std::vector<GravityModelResult> evaluate(const std::vector<Array3> &computationPoints,
                                                                       ComputeBackend backend) const = 0;

        /**
         * Copies the polyhedron-dependent caches back from the device into host memory.
         * These are the segment vectors G_pq, the plane unit normals N_p, and the segment unit normals n_pq.
         * @return the three caches, always in double precision
         */
        [[nodiscard]] virtual std::tuple<std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
        getCaches() const = 0;

        /**
         * Returns the floating point precision the kernels of this evaluation compute in.
         * @return FLOAT32 or FLOAT64
         */
        [[nodiscard]] virtual ComputePrecision getPrecision() const = 0;
    };

    /**
     * Uploads a polyhedron to the compute device and computes the caches which only depend on it.
     *
     * @param polyhedron the constant density polyhedron
     * @param precision the floating point precision the kernels should compute in
     * @return the evaluation engine for this polyhedron
     * @throws std::invalid_argument if the polyhedron has no faces
     */
    std::shared_ptr<KokkosEvaluationBase> createKokkosEvaluation(const Polyhedron &polyhedron,
                                                                 ComputePrecision precision);

    /**
     * Uploads a polyhedron together with already known caches, i.e. without recomputing them.
     * This restores a {@link GravityEvaluable} from a previous state, for example after unpickling.
     *
     * @param polyhedron the constant density polyhedron
     * @param precision the floating point precision the kernels should compute in
     * @param segmentVectors the segment vectors G_pq foreach face
     * @param planeUnitNormals the plane unit normals N_p foreach face
     * @param segmentUnitNormals the segment unit normals n_pq foreach face
     * @return the evaluation engine for this polyhedron
     * @throws std::invalid_argument if the polyhedron has no faces or a cache's size does not match
     */
    std::shared_ptr<KokkosEvaluationBase> createKokkosEvaluation(const Polyhedron &polyhedron,
                                                                 ComputePrecision precision,
                                                                 const std::vector<Array3Triplet> &segmentVectors,
                                                                 const std::vector<Array3> &planeUnitNormals,
                                                                 const std::vector<Array3Triplet> &segmentUnitNormals);

}// namespace polyhedralGravity::kokkos
