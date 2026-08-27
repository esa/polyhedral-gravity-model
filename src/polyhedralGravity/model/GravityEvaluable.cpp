#include "GravityEvaluable.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "polyhedralGravity/model/GravityModelDetail.h"
#include "polyhedralGravity/model/PolyhedralMeshView.h"
#include "polyhedralGravity/output/Logging.h"
#include "polyhedralGravity/util/KokkosSession.h"

namespace polyhedralGravity::detail {

    using namespace polyhedralGravity::kokkos;
    using GravityModel::detail::FaceContribution;

    /**
     * The raw output of the multi point kernel, i.e. one result per computation point before Tsoulis' prefix
     * has been applied.
     *
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the Kokkos memory space the views are allocated in
     */
    template<typename FloatType, typename MemorySpace>
    struct EvaluationResultView {
        /** The gravitational potential V foreach computation point, of the extents @f$(Q)@f$ */
        ScalarView<FloatType, MemorySpace> potential;
        /** The first order derivatives Vx, Vy, Vz foreach computation point, of the extents @f$(Q, 3)@f$ */
        Vector3View<FloatType, MemorySpace> acceleration;
        /** The second order derivatives foreach computation point, of the extents @f$(Q, 6)@f$ */
        Vector6View<FloatType, MemorySpace> gradiometricTensor;
        /** How many faces were numerically critical foreach computation point, of the extents @f$(Q)@f$ */
        ScalarView<int, MemorySpace> numericallyCriticalFaces;

        /**
         * Allocates the four views for a given number of computation points.
         * @param pointCount the number of computation points
         * @return the uninitialized result views
         */
        [[nodiscard]] static EvaluationResultView allocate(const size_t pointCount) {
            return {ScalarView<FloatType, MemorySpace>{
                            Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::potential"), pointCount},
                    Vector3View<FloatType, MemorySpace>{
                            Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::acceleration"), pointCount},
                    Vector6View<FloatType, MemorySpace>{
                            Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::gradiometricTensor"), pointCount},
                    ScalarView<int, MemorySpace>{
                            Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::criticalFaces"), pointCount}};
        }
    };

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

    /*
     * The functions below are the only ones which launch a Kokkos kernel. They are free functions with
     * external linkage on purpose: nvcc does not accept an extended device lambda inside a private member
     * function, which is what they would otherwise be.
     */

    /**
     * Narrows a polyhedral mesh into the evaluation's precision, leaving it in its memory space.
     *
     * If the evaluation computes in the precision the polyhedron is stored in, i.e. in double precision,
     * this hands out the polyhedron's own views and nothing is copied at all.
     *
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam ExecutionSpace the execution space to run the conversion in
     * @param mesh the polyhedron's mesh in double precision
     * @return the same mesh in the evaluation's precision
     */
    template<typename FloatType, typename ExecutionSpace>
    PolyhedralMeshView<FloatType, typename ExecutionSpace::memory_space> narrowMesh(
            const PolyhedralMeshView<double, typename ExecutionSpace::memory_space> &mesh) {
        using MemorySpace = typename ExecutionSpace::memory_space;
        if constexpr (std::is_same_v<FloatType, double>) {
            return mesh;
        } else {
            const size_t vertexCount = mesh.countVertices();
            PolyhedralMeshView<FloatType, MemorySpace> narrowed{};
            narrowed.faces = mesh.faces;
            narrowed.vertices = Vector3View<FloatType, MemorySpace>{
                    Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::vertices"), vertexCount};
            const auto source = mesh.vertices;
            const auto target = narrowed.vertices;
            Kokkos::parallel_for(
                    "polyhedralGravity::narrowVertices", Kokkos::RangePolicy<ExecutionSpace>(0, vertexCount),
                    KOKKOS_LAMBDA(const size_t index) {
                        for (size_t component = 0; component < 3; ++component) {
                            target(index, component) = static_cast<FloatType>(source(index, component));
                        }
                    });
            Kokkos::fence();
            return narrowed;
        }
    }

    /**
     * Runs Tsoulis' steps 1-01 to 1-03 for every face, i.e. everything which only depends on the polyhedron
     * and can therefore be cached and reused for every computation point.
     * @tparam ExecutionSpace the Kokkos execution space to run on
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the memory space the mesh lives in, must be accessible from ExecutionSpace
     * @param mesh the polyhedron, whose cache views this kernel fills
     */
    template<typename ExecutionSpace, typename FloatType, typename MemorySpace>
    void runInitializationKernel(const GravitationalMeshView<FloatType, MemorySpace> &mesh) {
        using namespace GravityModel::detail;
        Kokkos::parallel_for(
                "polyhedralGravity::initializeFaceProperties",
                Kokkos::RangePolicy<ExecutionSpace>(0, mesh.countFaces()), KOKKOS_LAMBDA(const size_t faceIndex) {
                    const Vector3Triplet<FloatType> face = mesh.resolveFace(faceIndex, Vector3<FloatType>{});
                    //1-01 Step: Compute Segment Vectors G_pq which describe each one the edge between two vertices
                    const Vector3Triplet<FloatType> segmentVectors = buildVectorsOfSegments(face[0], face[1], face[2]);
                    //1-02 Step: Compute the Plane Unit Normals N_p (pointing outside the polyhedron)
                    const Vector3<FloatType> planeUnitNormal = buildUnitNormalOfPlane(segmentVectors[0], segmentVectors[1]);
                    //1-03 Step: Compute Segment Unit Normals n_pq (normal pointing away from each segment)
                    mesh.setCaches(faceIndex, segmentVectors, planeUnitNormal,
                                   buildUnitNormalOfSegments(segmentVectors, planeUnitNormal));
                });
        Kokkos::fence();
    }

    /**
     * Reduces the contributions of all faces of the polyhedron at one computation point.
     * @tparam ExecutionSpace the Kokkos execution space to run on
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the memory space the mesh lives in, must be accessible from ExecutionSpace
     * @param mesh the polyhedron in the matching memory space
     * @param point the computation point P, already relative to nothing, i.e. in the mesh's coordinates
     * @return the summed contribution of every face
     */
    template<typename ExecutionSpace, typename FloatType, typename MemorySpace>
    FaceContribution<FloatType> runSinglePointKernel(const GravitationalMeshView<FloatType, MemorySpace> &mesh,
                                                     const Vector3<FloatType> &point) {
        using GravityModel::detail::evaluateFace;
        FaceContribution<FloatType> contribution{};
        Kokkos::parallel_reduce(
                "polyhedralGravity::evaluate", Kokkos::RangePolicy<ExecutionSpace>(0, mesh.countFaces()),
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
     * @param points the computation points, in the same memory space as the mesh
     * @return the raw results, i.e. one entry per computation point in each of the result views
     */
    template<typename ExecutionSpace, typename FloatType, typename MemorySpace>
    EvaluationResultView<FloatType, MemorySpace> runMultiPointKernel(
            const GravitationalMeshView<FloatType, MemorySpace> &mesh,
            const Vector3View<FloatType, MemorySpace> &points) {
        using GravityModel::detail::evaluateFace;
        using Policy = Kokkos::TeamPolicy<ExecutionSpace>;
        const size_t pointCount = points.extent(0);
        const size_t faceCount = mesh.countFaces();
        const EvaluationResultView<FloatType, MemorySpace> results =
                EvaluationResultView<FloatType, MemorySpace>::allocate(pointCount);
        Kokkos::parallel_for(
                "polyhedralGravity::evaluateMultiPoint", Policy(pointCount, Kokkos::AUTO),
                KOKKOS_LAMBDA(const typename Policy::member_type &team) {
                    const auto pointIndex = static_cast<size_t>(team.league_rank());
                    const Vector3<FloatType> point{points(pointIndex, 0), points(pointIndex, 1),
                                                   points(pointIndex, 2)};
                    FaceContribution<FloatType> contribution{};
                    Kokkos::parallel_reduce(
                            Kokkos::TeamThreadRange(team, faceCount),
                            [&](const size_t faceIndex, FaceContribution<FloatType> &accumulator) {
                                accumulator += evaluateFace(mesh, faceIndex, point);
                            },
                            contribution);
                    Kokkos::single(Kokkos::PerTeam(team), [&]() {
                        results.potential(pointIndex) = contribution.potential;
                        results.numericallyCriticalFaces(pointIndex) = contribution.numericallyCriticalFaces;
                        for (size_t component = 0; component < 3; ++component) {
                            results.acceleration(pointIndex, component) = contribution.acceleration[component];
                        }
                        for (size_t component = 0; component < 6; ++component) {
                            results.gradiometricTensor(pointIndex, component) =
                                    contribution.gradiometricTensor[component];
                        }
                    });
                });
        Kokkos::fence();
        return results;
    }

    /**
     * Copies the computation points into the given memory space, converting them to the evaluation's precision.
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the memory space to allocate in
     * @param computationPoints the computation points
     * @return the view holding the points, of the extents @f$(Q, 3)@f$
     */
    template<typename FloatType, typename MemorySpace>
    Vector3View<FloatType, MemorySpace> uploadPoints(const std::vector<Array3> &computationPoints) {
        const size_t pointCount = computationPoints.size();
        Vector3View<FloatType, HostMemory> hostPoints{
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::stagedPoints"), pointCount};
        for (size_t index = 0; index < pointCount; ++index) {
            for (size_t component = 0; component < 3; ++component) {
                hostPoints(index, component) = static_cast<FloatType>(computationPoints[index][component]);
            }
        }
        return Kokkos::create_mirror_view_and_copy(MemorySpace{}, hostPoints);
    }

    /**
     * The precision- and backend-agnostic interface to the polyhedron as a {@link GravityEvaluable} holds it.
     *
     * The implementation behind it is templated over the floating point precision and over the Kokkos
     * execution space, both of which are runtime choices of the user (see {@link ComputePrecision} and
     * {@link ComputeBackend}) while the kernels have to be compiled for each of them. That is what this
     * interface is virtual for.
     */
    class EvaluationEngine {
    public:
        virtual ~EvaluationEngine() = default;

        /**
         * Evaluates the polyhedral gravity model at a single computation point.
         * @param computationPoint the computation point P
         * @return the potential, the acceleration, and the gradiometric tensor at P
         */
        [[nodiscard]] virtual GravityModelResult evaluate(const Array3 &computationPoint) const = 0;

        /**
         * Evaluates the polyhedral gravity model at multiple computation points.
         * The points are evaluated in one kernel launch, with one team of threads per computation point.
         * @param computationPoints the computation points
         * @return the results foreach computation point, in the order of the input
         */
        [[nodiscard]] virtual std::vector<GravityModelResult> evaluate(
                const std::vector<Array3> &computationPoints) const = 0;

        /**
         * Copies the polyhedron-dependent caches into host memory.
         * These are the segment vectors G_pq, the plane unit normals N_p, and the segment unit normals n_pq.
         * @return the three caches, always in double precision
         */
        [[nodiscard]] virtual std::tuple<std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
        getCaches() const = 0;
    };

    /**
     * The polyhedron with its caches and the kernels evaluating it, in one fixed floating point precision
     * and in one fixed execution space.
     *
     * Everything lives in exactly one memory space, namely the one the chosen execution space computes in.
     * An engine created for a GPU backend therefore keeps the mesh and the caches on the GPU for its whole
     * lifetime and only ever moves the computation points up and the results back down. Serving another
     * backend means creating another {@link GravityEvaluable}.
     *
     * The mesh itself belongs to the {@link Polyhedron} and is shared with it rather than copied, as long as
     * the polyhedron already holds it in this memory space and in this precision. Only the caches of
     * Tsoulis' algorithm, which extend a {@link kokkos::PolyhedralMeshView} into a
     * {@link kokkos::GravitationalMeshView}, are always allocated here.
     *
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam ExecutionSpace the Kokkos execution space the kernels run in
     */
    template<typename FloatType, typename ExecutionSpace>
    class TypedEvaluationEngine final : public EvaluationEngine {

        /** The memory space the mesh, the caches, and the results live in */
        using MemorySpace = typename ExecutionSpace::memory_space;

        /** The polyhedron and its caches, in the memory space of the execution space above */
        GravitationalMeshView<FloatType, MemorySpace> _mesh{};

        /** The number of triangular faces of the polyhedron */
        size_t _faceCount{0};

        /** The prefix of Tsoulis' equations, i.e. gravitational constant, density, orientation, and mesh unit */
        double _prefix{0.0};

    public:
        /**
         * Shares the polyhedron's mesh and computes the caches which only depend on it.
         * @param polyhedron the constant density polyhedron
         */
        explicit TypedEvaluationEngine(const Polyhedron &polyhedron)
            : _faceCount{polyhedron.countFaces()},
              _prefix{polyhedron.getGravityModelScaling()} {
            attachMesh(polyhedron);
            runInitializationKernel<ExecutionSpace>(_mesh);
        }

        /**
         * Shares the polyhedron's mesh together with already known caches.
         * @param polyhedron the constant density polyhedron
         * @param segmentVectors the segment vectors G_pq foreach face
         * @param planeUnitNormals the plane unit normals N_p foreach face
         * @param segmentUnitNormals the segment unit normals n_pq foreach face
         */
        TypedEvaluationEngine(const Polyhedron &polyhedron, const std::vector<Array3Triplet> &segmentVectors,
                              const std::vector<Array3> &planeUnitNormals,
                              const std::vector<Array3Triplet> &segmentUnitNormals)
            : _faceCount{polyhedron.countFaces()},
              _prefix{polyhedron.getGravityModelScaling()} {
            if (segmentVectors.size() != _faceCount || planeUnitNormals.size() != _faceCount ||
                segmentUnitNormals.size() != _faceCount) {
                throw std::invalid_argument{"The given caches do not have one entry per face of the polyhedron!"};
            }
            attachMesh(polyhedron);
            uploadCaches(segmentVectors, planeUnitNormals, segmentUnitNormals);
        }

        [[nodiscard]] GravityModelResult evaluate(const Array3 &computationPoint) const override {
            POLYHEDRAL_GRAVITY_LOG_DEBUG("Evaluation for computation point P = [{}, {}, {}] started on {}",
                                         computationPoint[0], computationPoint[1], computationPoint[2],
                                         ExecutionSpace::name());
            return finalize(runSinglePointKernel<ExecutionSpace>(_mesh, narrow<FloatType>(computationPoint)));
        }

        [[nodiscard]] std::vector<GravityModelResult> evaluate(
                const std::vector<Array3> &computationPoints) const override {
            POLYHEDRAL_GRAVITY_LOG_DEBUG("Evaluation for {} computation points started on {}",
                                         computationPoints.size(), ExecutionSpace::name());
            if (computationPoints.empty()) {
                return {};
            }
            return finalizeAll(runMultiPointKernel<ExecutionSpace>(
                    _mesh, uploadPoints<FloatType, MemorySpace>(computationPoints)));
        }

        [[nodiscard]] std::tuple<std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
        getCaches() const override {
            // On a build without a GPU backend, and for every host backend, these mirrors are the caches
            // themselves and nothing is copied
            const GravitationalMeshView<FloatType, HostMemory> hostMesh = mirrorCaches(true);
            std::vector<Array3Triplet> segmentVectors(_faceCount);
            std::vector<Array3> planeUnitNormals(_faceCount);
            std::vector<Array3Triplet> segmentUnitNormals(_faceCount);
            for (size_t index = 0; index < _faceCount; ++index) {
                const Vector3Triplet<FloatType> faceSegmentVectors = hostMesh.getSegmentVectors(index);
                const Vector3Triplet<FloatType> faceSegmentUnitNormals = hostMesh.getSegmentUnitNormals(index);
                for (size_t segment = 0; segment < 3; ++segment) {
                    segmentVectors[index][segment] = widen(faceSegmentVectors[segment]);
                    segmentUnitNormals[index][segment] = widen(faceSegmentUnitNormals[segment]);
                }
                planeUnitNormals[index] = widen(hostMesh.getPlaneUnitNormal(index));
            }
            return std::make_tuple(segmentVectors, planeUnitNormals, segmentUnitNormals);
        }

    private:
        /**
         * Shares the polyhedron's mesh and allocates the caches next to it, both in this engine's memory space.
         * @param polyhedron the constant density polyhedron
         */
        void attachMesh(const Polyhedron &polyhedron) {
            _mesh = GravitationalMeshView<FloatType, MemorySpace>::allocateFor(
                    narrowMesh<FloatType, ExecutionSpace>(polyhedron.getMesh().template getMeshIn<MemorySpace>()));
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
            const GravitationalMeshView<FloatType, HostMemory> hostMesh = mirrorCaches(false);
            for (size_t index = 0; index < _faceCount; ++index) {
                Vector3Triplet<FloatType> faceSegmentVectors{};
                Vector3Triplet<FloatType> faceSegmentUnitNormals{};
                for (size_t segment = 0; segment < 3; ++segment) {
                    faceSegmentVectors[segment] = narrow<FloatType>(segmentVectors[index][segment]);
                    faceSegmentUnitNormals[segment] = narrow<FloatType>(segmentUnitNormals[index][segment]);
                }
                hostMesh.setCaches(index, faceSegmentVectors, narrow<FloatType>(planeUnitNormals[index]),
                                   faceSegmentUnitNormals);
            }
            Kokkos::deep_copy(_mesh.segmentVectors, hostMesh.segmentVectors);
            Kokkos::deep_copy(_mesh.planeUnitNormals, hostMesh.planeUnitNormals);
            Kokkos::deep_copy(_mesh.segmentUnitNormals, hostMesh.segmentUnitNormals);
        }

        /**
         * Creates a host mirror of the three caches, which is the caches themselves for a host backend.
         * @param withValues whether the caches' current values are copied down, which the caller does not
         * want if it is about to overwrite them anyway
         * @return the caches as they are reachable from the host
         */
        [[nodiscard]] GravitationalMeshView<FloatType, HostMemory> mirrorCaches(const bool withValues) const {
            GravitationalMeshView<FloatType, HostMemory> hostMesh{};
            hostMesh.segmentVectors = Kokkos::create_mirror_view(HostMemory{}, _mesh.segmentVectors);
            hostMesh.planeUnitNormals = Kokkos::create_mirror_view(HostMemory{}, _mesh.planeUnitNormals);
            hostMesh.segmentUnitNormals = Kokkos::create_mirror_view(HostMemory{}, _mesh.segmentUnitNormals);
            if (withValues) {
                Kokkos::deep_copy(hostMesh.segmentVectors, _mesh.segmentVectors);
                Kokkos::deep_copy(hostMesh.planeUnitNormals, _mesh.planeUnitNormals);
                Kokkos::deep_copy(hostMesh.segmentUnitNormals, _mesh.segmentUnitNormals);
            }
            return hostMesh;
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
         * Copies the results of all computation points back to the host and applies the prefix to each.
         * @param results one raw result per computation point
         * @return the results as the user sees them
         */
        [[nodiscard]] std::vector<GravityModelResult> finalizeAll(
                const EvaluationResultView<FloatType, MemorySpace> &results) const {
            const auto potential = Kokkos::create_mirror_view_and_copy(HostMemory{}, results.potential);
            const auto acceleration = Kokkos::create_mirror_view_and_copy(HostMemory{}, results.acceleration);
            const auto gradiometricTensor =
                    Kokkos::create_mirror_view_and_copy(HostMemory{}, results.gradiometricTensor);
            const auto criticalFaces =
                    Kokkos::create_mirror_view_and_copy(HostMemory{}, results.numericallyCriticalFaces);
            std::vector<GravityModelResult> gravityResults{};
            gravityResults.reserve(potential.extent(0));
            for (size_t index = 0; index < potential.extent(0); ++index) {
                FaceContribution<FloatType> contribution{};
                contribution.potential = potential(index);
                contribution.numericallyCriticalFaces = criticalFaces(index);
                for (size_t component = 0; component < 3; ++component) {
                    contribution.acceleration[component] = acceleration(index, component);
                }
                for (size_t component = 0; component < 6; ++component) {
                    contribution.gradiometricTensor[component] = gradiometricTensor(index, component);
                }
                gravityResults.push_back(finalize(contribution));
            }
            return gravityResults;
        }
    };

    namespace {

        /**
         * Rejects polyhedra and backends which the kernels cannot handle.
         * @param polyhedron the polyhedron to check
         * @param backend the requested compute backend
         * @throws std::invalid_argument if the polyhedron has no faces
         * @throws std::runtime_error if the backend is not available in this build
         */
        void checkEvaluable(const Polyhedron &polyhedron, const ComputeBackend backend) {
            kokkos::ensureInitialized();
            kokkos::checkBackendAvailable(backend);
            if (polyhedron.countFaces() == 0) {
                throw std::invalid_argument{"The polyhedron does not have any faces, so there is nothing to evaluate!"};
            }
        }

        /**
         * Creates the engine for one precision, dispatching over the requested compute backend.
         * @tparam FloatType the floating point precision of the evaluation
         * @tparam Arguments the trailing constructor arguments, i.e. either nothing or the three caches
         * @param polyhedron the constant density polyhedron
         * @param backend the compute backend the engine's kernels run on
         * @param arguments the trailing constructor arguments
         * @return the engine for this polyhedron
         */
        template<typename FloatType, typename... Arguments>
        std::shared_ptr<EvaluationEngine> createTypedEngine(const Polyhedron &polyhedron,
                                                            const ComputeBackend backend,
                                                            const Arguments &...arguments) {
            switch (backend) {
                case ComputeBackend::CPU_SERIAL:
                    return std::make_shared<TypedEvaluationEngine<FloatType, kokkos::SerialSpace>>(
                            polyhedron, arguments...);
                case ComputeBackend::CPU_PARALLEL:
                    return std::make_shared<TypedEvaluationEngine<FloatType, kokkos::HostParallelSpace>>(
                            polyhedron, arguments...);
                default:
                    return std::make_shared<TypedEvaluationEngine<FloatType, kokkos::DeviceSpace>>(
                            polyhedron, arguments...);
            }
        }

        /**
         * Creates the engine of a {@link GravityEvaluable}, dispatching over precision and compute backend.
         * @tparam Arguments the trailing constructor arguments, i.e. either nothing or the three caches
         * @param polyhedron the constant density polyhedron
         * @param backend the compute backend the engine's kernels run on
         * @param precision the floating point precision the engine's kernels compute in
         * @param arguments the trailing constructor arguments
         * @return the engine for this polyhedron
         */
        template<typename... Arguments>
        std::shared_ptr<EvaluationEngine> createEngine(const Polyhedron &polyhedron, const ComputeBackend backend,
                                                       const ComputePrecision precision,
                                                       const Arguments &...arguments) {
            checkEvaluable(polyhedron, backend);
            if (precision == ComputePrecision::FLOAT32) {
                return createTypedEngine<float>(polyhedron, backend, arguments...);
            }
            return createTypedEngine<double>(polyhedron, backend, arguments...);
        }

    }// namespace

}// namespace polyhedralGravity::detail

namespace polyhedralGravity {

    GravityEvaluable::GravityEvaluable(const Polyhedron &polyhedron, const ComputeBackend backend,
                                       const ComputePrecision precision)
        : _polyhedron{polyhedron},
          _backend{backend},
          _precision{precision},
          _engine{detail::createEngine(polyhedron, backend, precision)} {
    }

    GravityEvaluable::GravityEvaluable(const Polyhedron &polyhedron,
                                       const std::vector<Array3Triplet> &segmentVectors,
                                       const std::vector<Array3> &planeUnitNormals,
                                       const std::vector<Array3Triplet> &segmentUnitNormals,
                                       const ComputeBackend backend, const ComputePrecision precision)
        : _polyhedron{polyhedron},
          _backend{backend},
          _precision{precision},
          _engine{detail::createEngine(polyhedron, backend, precision, segmentVectors, planeUnitNormals,
                                       segmentUnitNormals)} {
    }

    std::variant<GravityModelResult, std::vector<GravityModelResult>>
    GravityEvaluable::operator()(const std::variant<Array3, std::vector<Array3>> &computationPoints) const {
        if (std::holds_alternative<Array3>(computationPoints)) {
            return _engine->evaluate(std::get<Array3>(computationPoints));
        }
        return _engine->evaluate(std::get<std::vector<Array3>>(computationPoints));
    }

    std::string GravityEvaluable::toString() const {
        std::stringstream sstream;
        const auto [unitPotential, unitAcceleration, unitGradiometricTensor] = getOutputMetricUnit();
        sstream << "<polyhedral_gravity.GravityEvaluable, polyhedron = " << _polyhedron.toString()
                << ", output_units = " << unitPotential << ", " << unitAcceleration << ", " << unitGradiometricTensor
                << ", backend = " << _backend << ", precision = " << _precision << ">";
        return sstream.str();
    }

    std::array<std::string, 3> GravityEvaluable::getOutputMetricUnit() const {
        const auto metric = _polyhedron.getMeshUnit();
        if (metric != MetricUnit::UNITLESS) {
            const std::string metricString = _polyhedron.getMeshUnitAsString();
            return {metricString + "^2/s^2", metricString + "/s^2", "1/s^2"};
        } else {
            return {"1/s^2", "1/s^2", "1/s^2"};
        }
    }

    std::tuple<Polyhedron, std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
    GravityEvaluable::getState() const {
        const auto [segmentVectors, planeUnitNormals, segmentUnitNormals] = _engine->getCaches();
        return std::make_tuple(_polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals);
    }

    ComputeBackend GravityEvaluable::getComputeBackend() const {
        return _backend;
    }

    ComputePrecision GravityEvaluable::getComputePrecision() const {
        return _precision;
    }

}// namespace polyhedralGravity
