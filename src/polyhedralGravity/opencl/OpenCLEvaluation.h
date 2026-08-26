#pragma once

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/opencl/OpenCLContext.h"
#include "polyhedralGravity/opencl/OpenCLKernel.h"
#include <cstddef>
#include <mutex>
#include <string>
#include <vector>

namespace polyhedralGravity::opencl {

    /**
     * Evaluates Tsoulis et al.'s polyhedral gravity model on an OpenCL device.
     *
     * The polyhedron is uploaded once and the face properties which do not depend on the computation
     * point (segment vectors, plane unit normals, segment unit normals) are derived on the device and
     * stay resident there, so that evaluating a computation point only costs one kernel launch plus the
     * transfer of the reduced result.
     *
     * The summation over the faces happens on the device as well: the evaluation kernel reduces within
     * each work-group, and for meshes large enough to produce many work-groups a second kernel reduces
     * those partial results once more before they are read back.
     *
     * @note This class is the engine behind {@link ComputeBackend::OPENCL}. Users interact with it
     * through {@link GravityEvaluable} rather than directly.
     */
    class OpenCLEvaluation {

        /**
         * From this number of partial results on, the second reduction stage on the device is cheaper
         * than summing the partial results on the host.
         */
        static constexpr size_t REDUCTION_THRESHOLD = 128;

        /** The device, its context, and its command queue */
        OpenCLContext _context;

        /** Derives the point-independent face properties, run once in the constructor */
        OpenCLKernel _initializationKernel;

        /** Evaluates all faces for one computation point and reduces them per work-group */
        OpenCLKernel _evaluationKernel;

        /** Reduces the evaluation kernel's per-work-group results further */
        OpenCLKernel _reductionKernel;

        /** The number of faces of the polyhedron */
        size_t _numberOfFaces;

        /** The polyhedron's vertices in the device's floating point format */
        cl::Buffer _vertices;

        /** The polyhedron's faces as vertex index triplets */
        cl::Buffer _faces;

        /** The plane unit normals N_p, derived on the device */
        cl::Buffer _planeUnitNormals;

        /** The segment vectors G_pq, derived on the device */
        cl::Buffer _segmentVectors;

        /** The segment unit normals n_pq, derived on the device */
        cl::Buffer _segmentUnitNormals;

        /** One result per work-group of the evaluation kernel */
        cl::Buffer _partialResults;

        /** One result per work-group of the reduction kernel, unused below {@link REDUCTION_THRESHOLD} */
        cl::Buffer _reducedResults;

        /** Whether the second reduction stage is worthwhile for this polyhedron */
        bool _usesReductionKernel;

        /** The number of results transferred back per evaluation */
        size_t _resultCount;

        /** Staging buffer for the results, reused across evaluations to avoid reallocating */
        mutable std::vector<std::byte> _resultStaging;

        /**
         * Serializes evaluations. Setting kernel arguments is not thread-safe in OpenCL, and
         * {@link GravityEvaluable}'s evaluate is const and thus callable from several threads.
         */
        mutable std::mutex _mutex;

        /** GRAVITATIONAL_CONSTANT * density, corrected for normal orientation and mesh unit */
        double _prefix;

    public:
        /**
         * Uploads the polyhedron to a suitable OpenCL device and prepares it for evaluation.
         *
         * @param polyhedron the constant density polyhedron
         * @param precision the floating point precision the device computes in
         *
         * @throws std::runtime_error if no device supporting the precision exists, if a kernel does
         * not compile, or if the polyhedron is too large for the device's addressing
         */
        OpenCLEvaluation(const Polyhedron &polyhedron, ComputePrecision precision);

        /**
         * Evaluates the polyhedral gravity model at the given computation point.
         *
         * @param computationPoint the computation point P
         * @return the potential, the acceleration, and the second derivative tensor at P
         */
        [[nodiscard]] GravityModelResult evaluate(const Array3 &computationPoint) const;

        /**
         * Evaluates the polyhedral gravity model at multiple computation points.
         * The points are processed one after another since the parallelism is already exhausted by
         * spreading a single point's faces across the device.
         *
         * @param computationPoints the computation points
         * @return one result per computation point, in the same order
         */
        [[nodiscard]] std::vector<GravityModelResult> evaluate(const std::vector<Array3> &computationPoints) const;

        /**
         * The name of the device the evaluation runs on, e.g. "Apple M1 Pro".
         * @return the device name
         */
        [[nodiscard]] std::string getDeviceName() const;

        /**
         * The floating point precision the device computes in.
         * @return the precision
         */
        [[nodiscard]] ComputePrecision getPrecision() const;
    };

}// namespace polyhedralGravity::opencl
