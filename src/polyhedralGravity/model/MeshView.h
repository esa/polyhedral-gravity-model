#pragma once

#include <cstddef>

#include <Kokkos_Core.hpp>

#include "polyhedralGravity/model/KokkosSession.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/util/UtilityContainer.h"

namespace polyhedralGravity::kokkos {

    /*
     * Every view below is declared with its exact extents, i.e. only the number of vertices, faces, or
     * computation points is a runtime dimension ("*") while the trailing dimensions are compile time
     * constants (3 for a cartesian vector, 6 for the gradiometric tensor, 3x3 for the segments of a
     * triangular face). Two things follow from that:
     *
     *  - Kokkos knows the trailing strides at compile time, so indexing needs no runtime multiplication.
     *  - The memory layout is spelled out as LayoutRight, which is exactly the layout of a C-contiguous
     *    (N, 3) array of NumPy, PyTorch, or JAX. A view can therefore be wrapped around such a buffer
     *    without copying it, no matter whether it lives on the host or on the device.
     *
     * LayoutRight is stated explicitly because Kokkos would otherwise pick LayoutLeft for a device memory
     * space, which would silently transpose a foreign buffer.
     */

    /** A view of one cartesian vector per entry, i.e. of the extents @f$(N, 3)@f$ */
    template<typename Scalar, typename MemorySpace>
    using Vector3View = Kokkos::View<Scalar *[3], Kokkos::LayoutRight, MemorySpace>;

    /** A view of one gradiometric tensor per entry, i.e. of the extents @f$(N, 6)@f$ */
    template<typename Scalar, typename MemorySpace>
    using Vector6View = Kokkos::View<Scalar *[6], Kokkos::LayoutRight, MemorySpace>;

    /** A view of three cartesian vectors per entry, i.e. of the extents @f$(N, 3, 3)@f$ */
    template<typename Scalar, typename MemorySpace>
    using Vector3TripletView = Kokkos::View<Scalar *[3][3], Kokkos::LayoutRight, MemorySpace>;

    /** A view of one scalar per entry, i.e. of the extents @f$(N)@f$ */
    template<typename Scalar, typename MemorySpace>
    using ScalarView = Kokkos::View<Scalar *, Kokkos::LayoutRight, MemorySpace>;

    /**
     * The polyhedral mesh as it lives in the memory of one execution space.
     *
     * These are the elementary properties of a {@link Polyhedron}: which vertices it has and which of them
     * form its triangular faces. Nothing in here depends on the gravity model, which is why a
     * {@link Polyhedron} owns exactly this much and not more.
     *
     * @tparam FloatType the floating point precision the vertices are stored in
     * @tparam MemorySpace the Kokkos memory space the views are allocated in
     */
    template<typename FloatType, typename MemorySpace>
    struct PolyhedralMeshView {
        /** The polyhedron's vertices as cartesian coordinates, of the extents @f$(N, 3)@f$ */
        Vector3View<FloatType, MemorySpace> vertices;
        /** The polyhedron's triangular faces, each referencing three vertices by index, of the extents @f$(M, 3)@f$ */
        Vector3View<size_t, MemorySpace> faces;

        /**
         * The number of vertices of the polyhedron.
         * @return the extent of the vertices' runtime dimension
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION size_t countVertices() const {
            return vertices.extent(0);
        }

        /**
         * The number of triangular faces of the polyhedron.
         * @return the extent of the faces' runtime dimension
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION size_t countFaces() const {
            return faces.extent(0);
        }

        /**
         * Returns the cartesian coordinates of the vertex at the given index.
         * @param vertexIndex the index of the vertex
         * @return the vertex' x, y, and z coordinate
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION Vector3<FloatType> getVertex(const size_t vertexIndex) const {
            return {vertices(vertexIndex, 0), vertices(vertexIndex, 1), vertices(vertexIndex, 2)};
        }

        /**
         * Returns the indices of the three vertices forming the face at the given index.
         * @param faceIndex the index of the face
         * @return the triplet of vertex indices
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION IndexArray3 getFace(const size_t faceIndex) const {
            return {faces(faceIndex, 0), faces(faceIndex, 1), faces(faceIndex, 2)};
        }

        /**
         * Resolves the face at the given index into cartesian coordinates relative to a computation point,
         * i.e. it re-locates the computation point into the origin as Tsoulis' equations require.
         * @param faceIndex the index of the face
         * @param computationPoint the computation point P
         * @return the face's three vertices, each shifted by -P
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION Vector3Triplet<FloatType> resolveFace(
                const size_t faceIndex, const Vector3<FloatType> &computationPoint) const {
            using util::operator-;
            const IndexArray3 face = getFace(faceIndex);
            return {getVertex(face[0]) - computationPoint,
                    getVertex(face[1]) - computationPoint,
                    getVertex(face[2]) - computationPoint};
        }
    };

    /**
     * The polyhedral mesh extended by everything Tsoulis' algorithm caches per face.
     *
     * The three caches only depend on the polyhedron and not on the computation point, so they are filled
     * once by {@link KokkosEvaluation}'s initialization kernel and reused for every subsequent evaluation.
     * In contrast to the mesh itself they are not an elementary property of a polyhedron, which is why they
     * only exist inside a {@link GravityEvaluable} and not inside a {@link Polyhedron}.
     *
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the Kokkos memory space the views are allocated in
     */
    template<typename FloatType, typename MemorySpace>
    struct GravitationalMeshView : PolyhedralMeshView<FloatType, MemorySpace> {
        /** The segment vectors G_pq foreach face, of the extents @f$(M, 3, 3)@f$ */
        Vector3TripletView<FloatType, MemorySpace> segmentVectors;
        /** The plane unit normals N_p foreach face, of the extents @f$(M, 3)@f$ */
        Vector3View<FloatType, MemorySpace> planeUnitNormals;
        /** The segment unit normals n_pq foreach face, of the extents @f$(M, 3, 3)@f$ */
        Vector3TripletView<FloatType, MemorySpace> segmentUnitNormals;

        /**
         * Allocates the three caches for a mesh, i.e. turns a plain mesh view into a gravitational one.
         * @param mesh the polyhedral mesh in the same memory space
         * @return the mesh together with three uninitialized caches
         */
        [[nodiscard]] static GravitationalMeshView allocateFor(
                const PolyhedralMeshView<FloatType, MemorySpace> &mesh) {
            const size_t faceCount = mesh.countFaces();
            GravitationalMeshView result{};
            result.vertices = mesh.vertices;
            result.faces = mesh.faces;
            result.segmentVectors = Vector3TripletView<FloatType, MemorySpace>{
                    Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::segmentVectors"), faceCount};
            result.planeUnitNormals = Vector3View<FloatType, MemorySpace>{
                    Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::planeUnitNormals"), faceCount};
            result.segmentUnitNormals = Vector3TripletView<FloatType, MemorySpace>{
                    Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::segmentUnitNormals"), faceCount};
            return result;
        }

        /**
         * Returns the segment vectors G_pq of the face at the given index.
         * @param faceIndex the index of the face
         * @return the three segment vectors
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION Vector3Triplet<FloatType> getSegmentVectors(const size_t faceIndex) const {
            return readTriplet(segmentVectors, faceIndex);
        }

        /**
         * Returns the plane unit normal N_p of the face at the given index.
         * @param faceIndex the index of the face
         * @return the plane unit normal
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION Vector3<FloatType> getPlaneUnitNormal(const size_t faceIndex) const {
            return {planeUnitNormals(faceIndex, 0), planeUnitNormals(faceIndex, 1), planeUnitNormals(faceIndex, 2)};
        }

        /**
         * Returns the segment unit normals n_pq of the face at the given index.
         * @param faceIndex the index of the face
         * @return the three segment unit normals
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION Vector3Triplet<FloatType> getSegmentUnitNormals(const size_t faceIndex) const {
            return readTriplet(segmentUnitNormals, faceIndex);
        }

        /**
         * Writes the caches of the face at the given index.
         * @param faceIndex the index of the face
         * @param faceSegmentVectors the segment vectors G_pq
         * @param planeUnitNormal the plane unit normal N_p
         * @param faceSegmentUnitNormals the segment unit normals n_pq
         */
        KOKKOS_INLINE_FUNCTION void setCaches(const size_t faceIndex,
                                              const Vector3Triplet<FloatType> &faceSegmentVectors,
                                              const Vector3<FloatType> &planeUnitNormal,
                                              const Vector3Triplet<FloatType> &faceSegmentUnitNormals) const {
            writeTriplet(segmentVectors, faceIndex, faceSegmentVectors);
            writeTriplet(segmentUnitNormals, faceIndex, faceSegmentUnitNormals);
            for (size_t component = 0; component < 3; ++component) {
                planeUnitNormals(faceIndex, component) = planeUnitNormal[component];
            }
        }

    private:
        /**
         * Reads one @f$(3, 3)@f$ entry out of a triplet view.
         * @param view the view to read from
         * @param faceIndex the index of the face
         * @return the three vectors of that entry
         */
        [[nodiscard]] KOKKOS_INLINE_FUNCTION static Vector3Triplet<FloatType> readTriplet(
                const Vector3TripletView<FloatType, MemorySpace> &view, const size_t faceIndex) {
            Vector3Triplet<FloatType> triplet{};
            for (size_t segment = 0; segment < 3; ++segment) {
                for (size_t component = 0; component < 3; ++component) {
                    triplet[segment][component] = view(faceIndex, segment, component);
                }
            }
            return triplet;
        }

        /**
         * Writes one @f$(3, 3)@f$ entry into a triplet view.
         * @param view the view to write into
         * @param faceIndex the index of the face
         * @param triplet the three vectors to write
         */
        KOKKOS_INLINE_FUNCTION static void writeTriplet(const Vector3TripletView<FloatType, MemorySpace> &view,
                                                        const size_t faceIndex,
                                                        const Vector3Triplet<FloatType> &triplet) {
            for (size_t segment = 0; segment < 3; ++segment) {
                for (size_t component = 0; component < 3; ++component) {
                    view(faceIndex, segment, component) = triplet[segment][component];
                }
            }
        }
    };

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

}// namespace polyhedralGravity::kokkos
