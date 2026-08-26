#include "KokkosEvaluation.h"

#include <stdexcept>
#include <string>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "KokkosGravityKernel.h"
#include "KokkosSession.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/output/Logging.h"

namespace polyhedralGravity::kokkos {

    /** The memory space of the GPU (or of the host, if this build has no GPU backend) */
    using DeviceMemory = DeviceSpace::memory_space;
    /** The memory space of the CPU */
    using HostMemory = Kokkos::HostSpace;

    /**
     * Casts a host-side array of doubles into the evaluation's precision.
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam N the array's size
     * @param array the array in double precision
     * @return the same array in the evaluation's precision
     */
    template<typename FloatType, size_t N>
    std::array<FloatType, N> narrow(const std::array<double, N> &array) {
        std::array<FloatType, N> result{};
        for (size_t index = 0; index < N; ++index) {
            result[index] = static_cast<FloatType>(array[index]);
        }
        return result;
    }

    /**
     * Casts an array of the evaluation's precision back into double precision.
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam N the array's size
     * @param array the array in the evaluation's precision
     * @return the same array in double precision
     */
    template<typename FloatType, size_t N>
    std::array<double, N> widen(const std::array<FloatType, N> &array) {
        std::array<double, N> result{};
        for (size_t index = 0; index < N; ++index) {
            result[index] = static_cast<double>(array[index]);
        }
        return result;
    }

    /*
     * The three functions below are the only ones which launch a Kokkos kernel. They are free functions with
     * external linkage on purpose: nvcc does not accept an extended device lambda inside a private member
     * function, which is what they would otherwise be.
     */

    /**
     * Runs Tsoulis' steps 1-01 to 1-03 for every face, i.e. everything which only depends on the polyhedron
     * and can therefore be cached and reused for every computation point.
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the memory space the mesh lives in
     * @param mesh the polyhedron, whose cache views this kernel fills
     * @param faceCount the number of faces of the polyhedron
     */
    template<typename FloatType, typename MemorySpace>
    void runInitializationKernel(const MeshViews<FloatType, MemorySpace> &mesh, const size_t faceCount) {
        using namespace GravityModel::detail;
        Kokkos::parallel_for(
                "polyhedralGravity::initializeFaceProperties",
                Kokkos::RangePolicy<DeviceSpace>(0, faceCount), KOKKOS_LAMBDA(const size_t faceIndex) {
                    const Vector3Triplet<FloatType> face = mesh.resolveFace(faceIndex, Vector3<FloatType>{});
                    //1-01 Step: Compute Segment Vectors G_pq which describe each one the edge between two vertices
                    const Vector3Triplet<FloatType> segmentVectors = buildVectorsOfSegments(face[0], face[1], face[2]);
                    //1-02 Step: Compute the Plane Unit Normals N_p (pointing outside the polyhedron)
                    const Vector3<FloatType> planeUnitNormal = buildUnitNormalOfPlane(segmentVectors[0], segmentVectors[1]);
                    //1-03 Step: Compute Segment Unit Normals n_pq (normal pointing away from each segment)
                    mesh.segmentVectors(faceIndex) = segmentVectors;
                    mesh.planeUnitNormals(faceIndex) = planeUnitNormal;
                    mesh.segmentUnitNormals(faceIndex) = buildUnitNormalOfSegments(segmentVectors, planeUnitNormal);
                });
        Kokkos::fence();
    }

    /**
     * Reduces the contributions of all faces of the polyhedron at one computation point.
     * @tparam ExecutionSpace the Kokkos execution space to run on
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the memory space the mesh lives in, must be accessible from ExecutionSpace
     * @param mesh the polyhedron in the matching memory space
     * @param faceCount the number of faces of the polyhedron
     * @param point the computation point P, already relative to nothing, i.e. in the mesh's coordinates
     * @return the summed contribution of every face
     */
    template<typename ExecutionSpace, typename FloatType, typename MemorySpace>
    FaceContribution<FloatType> runSinglePointKernel(const MeshViews<FloatType, MemorySpace> &mesh,
                                                     const size_t faceCount, const Vector3<FloatType> &point) {
        FaceContribution<FloatType> contribution{};
        Kokkos::parallel_reduce(
                "polyhedralGravity::evaluate", Kokkos::RangePolicy<ExecutionSpace>(0, faceCount),
                KOKKOS_LAMBDA(const size_t faceIndex, FaceContribution<FloatType> &accumulator) {
                    accumulator += evaluateFace(mesh, faceIndex, point);
                },
                contribution);
        return contribution;
    }

    /**
     * Reduces the contributions of all faces of the polyhedron at every given computation point in a single
     * kernel launch. One team of threads handles one computation point, so the kernel saturates a GPU even
     * for a comparatively small polyhedron as long as enough points are asked for at once.
     * @tparam ExecutionSpace the Kokkos execution space to run on
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the memory space the mesh lives in, must be accessible from ExecutionSpace
     * @param mesh the polyhedron in the matching memory space
     * @param faceCount the number of faces of the polyhedron
     * @param points the computation points, in the same memory space as the mesh
     * @return the summed contribution of every face, foreach computation point
     */
    template<typename ExecutionSpace, typename FloatType, typename MemorySpace>
    Kokkos::View<FaceContribution<FloatType> *, MemorySpace> runMultiPointKernel(
            const MeshViews<FloatType, MemorySpace> &mesh, const size_t faceCount,
            const Kokkos::View<Vector3<FloatType> *, MemorySpace> &points) {
        using Policy = Kokkos::TeamPolicy<ExecutionSpace>;
        const size_t pointCount = points.extent(0);
        Kokkos::View<FaceContribution<FloatType> *, MemorySpace> contributions{
                "polyhedralGravity::contributions", pointCount};
        Kokkos::parallel_for(
                "polyhedralGravity::evaluateMultiPoint", Policy(pointCount, Kokkos::AUTO),
                KOKKOS_LAMBDA(const typename Policy::member_type &team) {
                    const auto pointIndex = static_cast<size_t>(team.league_rank());
                    FaceContribution<FloatType> contribution{};
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadRange(team, faceCount),
                            [&](const size_t faceIndex, FaceContribution<FloatType> &accumulator) {
                                accumulator += evaluateFace(mesh, faceIndex, points(pointIndex));
                            },
                            contribution);
                    Kokkos::single(Kokkos::PerTeam(team), [&]() { contributions(pointIndex) = contribution; });
                });
        Kokkos::fence();
        return contributions;
    }

    /**
     * The device-resident polyhedron and the kernels evaluating it, in one fixed floating point precision.
     *
     * The mesh lives twice: once in the memory space of the GPU and once in the memory space of the CPU. On a
     * build without a GPU backend both memory spaces are the same one, so the views alias each other and
     * nothing is duplicated. Holding both lets the very same GravityEvaluable serve CPU_SERIAL, CPU_PARALLEL,
     * and GPU_PARALLEL without re-uploading the polyhedron.
     *
     * @tparam FloatType the floating point precision of the evaluation
     */
    template<typename FloatType>
    class KokkosEvaluation final : public KokkosEvaluationBase {

        /** The polyhedron in the memory space of the device */
        MeshViews<FloatType, DeviceMemory> _deviceMesh{};

        /** The polyhedron in the memory space of the host */
        MeshViews<FloatType, HostMemory> _hostMesh{};

        /** The number of triangular faces of the polyhedron */
        size_t _faceCount{0};

        /** The prefix of Tsoulis' equations, i.e. gravitational constant, density, orientation, and mesh unit */
        double _prefix{0.0};

    public:
        /**
         * Uploads the polyhedron and computes the caches which only depend on it.
         * @param polyhedron the constant density polyhedron
         */
        explicit KokkosEvaluation(const Polyhedron &polyhedron)
            : _faceCount{polyhedron.countFaces()},
              _prefix{polyhedron.getGravityModelScaling()} {
            uploadMesh(polyhedron);
            runInitializationKernel(_deviceMesh, _faceCount);
            copyCachesToHost();
        }

        /**
         * Uploads the polyhedron together with already known caches.
         * @param polyhedron the constant density polyhedron
         * @param segmentVectors the segment vectors G_pq foreach face
         * @param planeUnitNormals the plane unit normals N_p foreach face
         * @param segmentUnitNormals the segment unit normals n_pq foreach face
         */
        KokkosEvaluation(const Polyhedron &polyhedron, const std::vector<Array3Triplet> &segmentVectors,
                         const std::vector<Array3> &planeUnitNormals,
                         const std::vector<Array3Triplet> &segmentUnitNormals)
            : _faceCount{polyhedron.countFaces()},
              _prefix{polyhedron.getGravityModelScaling()} {
            if (segmentVectors.size() != _faceCount || planeUnitNormals.size() != _faceCount ||
                segmentUnitNormals.size() != _faceCount) {
                throw std::invalid_argument{"The given caches do not have one entry per face of the polyhedron!"};
            }
            uploadMesh(polyhedron);
            uploadCaches(segmentVectors, planeUnitNormals, segmentUnitNormals);
        }

        [[nodiscard]] GravityModelResult evaluate(const Array3 &computationPoint,
                                                  const ComputeBackend backend) const override {
            checkBackendAvailable(backend);
            POLYHEDRAL_GRAVITY_LOG_DEBUG("Evaluation for computation point P = [{}, {}, {}] started on {}",
                                         computationPoint[0], computationPoint[1], computationPoint[2],
                                         getExecutionSpaceName(backend));
            const Vector3<FloatType> point = narrow<FloatType>(computationPoint);
            switch (backend) {
                case ComputeBackend::CPU_SERIAL:
                    return finalize(runSinglePointKernel<SerialSpace>(_hostMesh, _faceCount, point));
                case ComputeBackend::CPU_PARALLEL:
                    return finalize(runSinglePointKernel<HostParallelSpace>(_hostMesh, _faceCount, point));
                default:
                    return finalize(runSinglePointKernel<DeviceSpace>(_deviceMesh, _faceCount, point));
            }
        }

        [[nodiscard]] std::vector<GravityModelResult> evaluate(const std::vector<Array3> &computationPoints,
                                                               const ComputeBackend backend) const override {
            checkBackendAvailable(backend);
            POLYHEDRAL_GRAVITY_LOG_DEBUG("Evaluation for {} computation points started on {}",
                                         computationPoints.size(), getExecutionSpaceName(backend));
            if (computationPoints.empty()) {
                return {};
            }
            switch (backend) {
                case ComputeBackend::CPU_SERIAL:
                    return finalizeAll(runMultiPointKernel<SerialSpace>(
                            _hostMesh, _faceCount, uploadPoints<HostMemory>(computationPoints)));
                case ComputeBackend::CPU_PARALLEL:
                    return finalizeAll(runMultiPointKernel<HostParallelSpace>(
                            _hostMesh, _faceCount, uploadPoints<HostMemory>(computationPoints)));
                default:
                    return finalizeAll(runMultiPointKernel<DeviceSpace>(
                            _deviceMesh, _faceCount, uploadPoints<DeviceMemory>(computationPoints)));
            }
        }

        [[nodiscard]] std::tuple<std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
        getCaches() const override {
            std::vector<Array3Triplet> segmentVectors(_faceCount);
            std::vector<Array3> planeUnitNormals(_faceCount);
            std::vector<Array3Triplet> segmentUnitNormals(_faceCount);
            for (size_t index = 0; index < _faceCount; ++index) {
                for (size_t j = 0; j < 3; ++j) {
                    segmentVectors[index][j] = widen(_hostMesh.segmentVectors(index)[j]);
                    segmentUnitNormals[index][j] = widen(_hostMesh.segmentUnitNormals(index)[j]);
                }
                planeUnitNormals[index] = widen(_hostMesh.planeUnitNormals(index));
            }
            return std::make_tuple(segmentVectors, planeUnitNormals, segmentUnitNormals);
        }

        [[nodiscard]] ComputePrecision getPrecision() const override {
            return std::is_same_v<FloatType, float> ? ComputePrecision::FLOAT32 : ComputePrecision::FLOAT64;
        }

    private:
        /**
         * Allocates the mesh in both memory spaces and copies the polyhedron's vertices and faces into them.
         * @param polyhedron the constant density polyhedron
         */
        void uploadMesh(const Polyhedron &polyhedron) {
            const size_t vertexCount = polyhedron.countVertices();

            Kokkos::View<Vector3<FloatType> *, HostMemory> vertices{"polyhedralGravity::stagedVertices", vertexCount};
            Kokkos::View<IndexArray3 *, HostMemory> faces{"polyhedralGravity::stagedFaces", _faceCount};
            for (size_t index = 0; index < vertexCount; ++index) {
                vertices(index) = narrow<FloatType>(polyhedron.getVertex(index));
            }
            for (size_t index = 0; index < _faceCount; ++index) {
                faces(index) = polyhedron.getFace(index);
            }

            _deviceMesh.vertices = Kokkos::View<Vector3<FloatType> *, DeviceMemory>{"polyhedralGravity::vertices", vertexCount};
            _deviceMesh.faces = Kokkos::View<IndexArray3 *, DeviceMemory>{"polyhedralGravity::faces", _faceCount};
            _deviceMesh.segmentVectors = Kokkos::View<Vector3Triplet<FloatType> *, DeviceMemory>{"polyhedralGravity::segmentVectors", _faceCount};
            _deviceMesh.planeUnitNormals = Kokkos::View<Vector3<FloatType> *, DeviceMemory>{"polyhedralGravity::planeUnitNormals", _faceCount};
            _deviceMesh.segmentUnitNormals = Kokkos::View<Vector3Triplet<FloatType> *, DeviceMemory>{"polyhedralGravity::segmentUnitNormals", _faceCount};

            Kokkos::deep_copy(_deviceMesh.vertices, vertices);
            Kokkos::deep_copy(_deviceMesh.faces, faces);

            // On a build without a GPU backend these mirrors are the device views themselves
            _hostMesh.vertices = Kokkos::create_mirror_view(HostMemory{}, _deviceMesh.vertices);
            _hostMesh.faces = Kokkos::create_mirror_view(HostMemory{}, _deviceMesh.faces);
            _hostMesh.segmentVectors = Kokkos::create_mirror_view(HostMemory{}, _deviceMesh.segmentVectors);
            _hostMesh.planeUnitNormals = Kokkos::create_mirror_view(HostMemory{}, _deviceMesh.planeUnitNormals);
            _hostMesh.segmentUnitNormals = Kokkos::create_mirror_view(HostMemory{}, _deviceMesh.segmentUnitNormals);
            Kokkos::deep_copy(_hostMesh.vertices, _deviceMesh.vertices);
            Kokkos::deep_copy(_hostMesh.faces, _deviceMesh.faces);
        }

        /**
         * Fills the caches with values which were computed before instead of running the initialization kernel.
         * @param segmentVectors the segment vectors G_pq foreach face
         * @param planeUnitNormals the plane unit normals N_p foreach face
         * @param segmentUnitNormals the segment unit normals n_pq foreach face
         */
        void uploadCaches(const std::vector<Array3Triplet> &segmentVectors,
                          const std::vector<Array3> &planeUnitNormals,
                          const std::vector<Array3Triplet> &segmentUnitNormals) {
            for (size_t index = 0; index < _faceCount; ++index) {
                for (size_t j = 0; j < 3; ++j) {
                    _hostMesh.segmentVectors(index)[j] = narrow<FloatType>(segmentVectors[index][j]);
                    _hostMesh.segmentUnitNormals(index)[j] = narrow<FloatType>(segmentUnitNormals[index][j]);
                }
                _hostMesh.planeUnitNormals(index) = narrow<FloatType>(planeUnitNormals[index]);
            }
            Kokkos::deep_copy(_deviceMesh.segmentVectors, _hostMesh.segmentVectors);
            Kokkos::deep_copy(_deviceMesh.planeUnitNormals, _hostMesh.planeUnitNormals);
            Kokkos::deep_copy(_deviceMesh.segmentUnitNormals, _hostMesh.segmentUnitNormals);
        }

        /** Mirrors the caches from the device down to the host so that both meshes hold the identical values */
        void copyCachesToHost() {
            Kokkos::deep_copy(_hostMesh.segmentVectors, _deviceMesh.segmentVectors);
            Kokkos::deep_copy(_hostMesh.planeUnitNormals, _deviceMesh.planeUnitNormals);
            Kokkos::deep_copy(_hostMesh.segmentUnitNormals, _deviceMesh.segmentUnitNormals);
        }

        /**
         * Copies the computation points into the given memory space, converting them to the evaluation's precision.
         * @tparam MemorySpace the memory space to allocate in
         * @param computationPoints the computation points
         * @return the view holding the points
         */
        template<typename MemorySpace>
        [[nodiscard]] Kokkos::View<Vector3<FloatType> *, MemorySpace> uploadPoints(
                const std::vector<Array3> &computationPoints) const {
            const size_t pointCount = computationPoints.size();
            Kokkos::View<Vector3<FloatType> *, HostMemory> hostPoints{"polyhedralGravity::stagedPoints", pointCount};
            for (size_t index = 0; index < pointCount; ++index) {
                hostPoints(index) = narrow<FloatType>(computationPoints[index]);
            }
            Kokkos::View<Vector3<FloatType> *, MemorySpace> points{"polyhedralGravity::points", pointCount};
            Kokkos::deep_copy(points, hostPoints);
            return points;
        }

        /**
         * Applies Tsoulis' prefix to the raw sums of one computation point and widens them to double precision.
         * @param contribution the summed contributions of all faces
         * @return the result as the user sees it
         */
        [[nodiscard]] GravityModelResult finalize(const FaceContribution<FloatType> &contribution) const {
            using namespace util;
            if (contribution.numericallyCriticalFaces > 0) {
                POLYHEDRAL_GRAVITY_LOG_WARN("While evaluating {} of the polyhedron's {} planes a significant "
                                            "difference of magnitudes occurred. "
                                            "This may lead to numerically unstable results!",
                                            contribution.numericallyCriticalFaces, _faceCount);
            }
            //9. Step: The prefix consists of the gravitational constant, the density, and the correction
            // factors depending on the alignment of the normals and the polyhedral mesh's unit
            //10. Step: Final expressions after application of the prefix (and a division by 2 for the potential)
            const double potential = (static_cast<double>(contribution.potential) * _prefix) / 2.0;
            const Array3 acceleration = widen(contribution.acceleration) * (-1.0 * _prefix);
            const Array6 gradiometricTensor = widen(contribution.gradiometricTensor) * _prefix;
            return std::make_tuple(potential, acceleration, gradiometricTensor);
        }

        /**
         * Copies the contributions of all computation points back to the host and applies the prefix to each.
         * @tparam MemorySpace the memory space the contributions live in
         * @param contributions one summed contribution per computation point
         * @return the results as the user sees them
         */
        template<typename MemorySpace>
        [[nodiscard]] std::vector<GravityModelResult> finalizeAll(
                const Kokkos::View<FaceContribution<FloatType> *, MemorySpace> &contributions) const {
            const auto hostContributions = Kokkos::create_mirror_view_and_copy(HostMemory{}, contributions);
            std::vector<GravityModelResult> results{};
            results.reserve(hostContributions.extent(0));
            for (size_t index = 0; index < hostContributions.extent(0); ++index) {
                results.push_back(finalize(hostContributions(index)));
            }
            return results;
        }
    };

    /* The kernels are compiled for both precisions the user can choose between */
    template class KokkosEvaluation<float>;
    template class KokkosEvaluation<double>;

    namespace {

        /**
         * Rejects polyhedra which the kernels cannot handle.
         * @param polyhedron the polyhedron to check
         * @throws std::invalid_argument if the polyhedron has no faces
         */
        void checkPolyhedron(const Polyhedron &polyhedron) {
            if (polyhedron.countFaces() == 0) {
                throw std::invalid_argument{"The polyhedron does not have any faces, so there is nothing to evaluate!"};
            }
        }

    }// namespace

    std::shared_ptr<KokkosEvaluationBase> createKokkosEvaluation(const Polyhedron &polyhedron,
                                                                 const ComputePrecision precision) {
        ensureInitialized();
        checkPolyhedron(polyhedron);
        if (precision == ComputePrecision::FLOAT32) {
            return std::make_shared<KokkosEvaluation<float>>(polyhedron);
        }
        return std::make_shared<KokkosEvaluation<double>>(polyhedron);
    }

    std::shared_ptr<KokkosEvaluationBase> createKokkosEvaluation(const Polyhedron &polyhedron,
                                                                 const ComputePrecision precision,
                                                                 const std::vector<Array3Triplet> &segmentVectors,
                                                                 const std::vector<Array3> &planeUnitNormals,
                                                                 const std::vector<Array3Triplet> &segmentUnitNormals) {
        ensureInitialized();
        checkPolyhedron(polyhedron);
        if (precision == ComputePrecision::FLOAT32) {
            return std::make_shared<KokkosEvaluation<float>>(polyhedron, segmentVectors, planeUnitNormals,
                                                             segmentUnitNormals);
        }
        return std::make_shared<KokkosEvaluation<double>>(polyhedron, segmentVectors, planeUnitNormals,
                                                          segmentUnitNormals);
    }

}// namespace polyhedralGravity::kokkos
