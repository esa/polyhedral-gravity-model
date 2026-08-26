#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <Kokkos_Core.hpp>

#include "polyhedralGravity/model/KokkosSession.h"
#include "polyhedralGravity/model/MeshView.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"

namespace polyhedralGravity {

    namespace detail {

        /**
         * Wraps a foreign @f$(N, 3)@f$ buffer in an unmanaged Kokkos view, i.e. without copying it.
         *
         * The resulting view carries no reference count, so the buffer has to outlive it. That is what the
         * keep-alive of {@link PolyhedralMesh} is for.
         *
         * @tparam Scalar the element type of the buffer, which must be the view's element type as well
         * @tparam MemorySpace the memory space the buffer lives in
         * @param data the first element of a C-contiguous @f$(N, 3)@f$ buffer
         * @param count the number of rows N
         * @return a view aliasing the buffer
         *
         * @note The const is cast away because Kokkos' pointer constructor takes a mutable pointer. The
         * library never writes through a wrapped buffer; healing the vertex ordering copies first
         * (see {@link PolyhedralMesh::setFaces}).
         */
        template<typename Scalar, typename MemorySpace>
        kokkos::Vector3View<Scalar, MemorySpace> wrapVector3(const Scalar *data, const size_t count) {
            return kokkos::Vector3View<Scalar, MemorySpace>{const_cast<Scalar *>(data), count};
        }

        /**
         * Turns a foreign @f$(N, 3)@f$ buffer into a view of the library's own element type.
         *
         * If the buffer's element type already is the target type, the buffer is only wrapped and no memory
         * is allocated at all. Otherwise a view is allocated in the same memory space and a kernel casts
         * every element into it, which keeps the data where it is, i.e. a device buffer is converted on the
         * device and never travels through the host.
         *
         * @tparam Target the element type of the resulting view, i.e. double or size_t
         * @tparam ExecutionSpace the execution space to run the conversion in
         * @tparam Source the element type of the buffer
         * @param data the first element of a C-contiguous @f$(N, 3)@f$ buffer
         * @param count the number of rows N
         * @param label the label of the view in case one has to be allocated
         * @return a view of the buffer's content in the target type
         */
        template<typename Target, typename ExecutionSpace, typename Source>
        kokkos::Vector3View<Target, typename ExecutionSpace::memory_space> convertVector3(
                const Source *data, const size_t count, const std::string &label) {
            using MemorySpace = typename ExecutionSpace::memory_space;
            if constexpr (std::is_same_v<Target, Source>) {
                return wrapVector3<Target, MemorySpace>(data, count);
            } else {
                const kokkos::Vector3View<Source, MemorySpace> source = wrapVector3<Source, MemorySpace>(data, count);
                const kokkos::Vector3View<Target, MemorySpace> target{
                        Kokkos::view_alloc(Kokkos::WithoutInitializing, label), count};
                Kokkos::parallel_for(
                        "polyhedralGravity::convertMesh", Kokkos::RangePolicy<ExecutionSpace>(0, count),
                        KOKKOS_LAMBDA(const size_t index) {
                            for (size_t component = 0; component < 3; ++component) {
                                target(index, component) = static_cast<Target>(source(index, component));
                            }
                        });
                Kokkos::fence();
                return target;
            }
        }

    }// namespace detail

    /**
     * The vertices and the triangular faces of a polyhedron, as they live in the memory of the compute
     * devices.
     *
     * This owns the elementary half of the {@link kokkos::PolyhedralMeshView} hierarchy: a mesh is
     * meaningful on its own, whereas the caches of {@link kokkos::GravitationalMeshView} only make sense
     * inside a {@link GravityEvaluable}.
     *
     * The mesh is kept in the memory space of the compute device. Its host mirror is created on demand, so
     * a polyhedron which was handed a device pointer and is only ever evaluated on the device never travels
     * through the host. On a build without a GPU backend both memory spaces are the same one and the mirror
     * is the mesh itself.
     */
    class PolyhedralMesh {

        /** The mesh in the memory space of the compute device, always in double precision */
        kokkos::PolyhedralMeshView<double, kokkos::DeviceMemory> _deviceMesh{};

        /** The mesh in the memory space of the host, created on the first host access */
        mutable kokkos::PolyhedralMeshView<double, kokkos::HostMemory> _hostMesh{};

        /**
         * Keeps a foreign buffer alive for as long as an unmanaged view aliases it.
         * This is empty whenever this mesh owns its memory itself.
         */
        std::shared_ptr<void> _keepAlive{};

    public:
        /** Creates an empty mesh without any vertices or faces */
        PolyhedralMesh() = default;

        /**
         * Copies vertices and faces given as host vectors into freshly allocated views.
         * @param vertices the vertices as cartesian coordinates
         * @param faces the triangular faces, each referencing three vertices by index
         */
        PolyhedralMesh(const std::vector<Array3> &vertices, const std::vector<IndexArray3> &faces);

        /**
         * Wraps the buffers of an array library, i.e. of NumPy, PyTorch, or JAX, without copying them.
         *
         * The buffers are only copied if their element type differs from the library's own, i.e. if the
         * vertices are not double precision or the face indices are not of the platform's @c size_t. Such a
         * conversion happens in the memory space the buffers already live in.
         *
         * @tparam VertexType the element type of the vertex buffer, e.g. float or double
         * @tparam IndexType the element type of the face buffer, e.g. int32_t or size_t
         * @param vertices the first element of a C-contiguous @f$(N, 3)@f$ vertex buffer
         * @param vertexCount the number of vertices N
         * @param faces the first element of a C-contiguous @f$(M, 3)@f$ face buffer
         * @param faceCount the number of faces M
         * @param location whether the buffers live in the host's or in the device's memory
         * @param keepAlive an owner of the buffers which is held for the lifetime of this mesh
         * @return the mesh aliasing the given buffers
         *
         * @throws std::runtime_error if DEVICE is given, but this build has no GPU backend
         *
         * @note If no keepAlive is given, the caller has to guarantee that the buffers outlive every
         * {@link Polyhedron} and {@link GravityEvaluable} built from them.
         */
        template<typename VertexType, typename IndexType>
        [[nodiscard]] static PolyhedralMesh fromBuffers(const VertexType *vertices, size_t vertexCount,
                                                        const IndexType *faces, size_t faceCount,
                                                        MemoryLocation location,
                                                        std::shared_ptr<void> keepAlive = {}) {
            kokkos::ensureInitialized();
            checkBuffers(vertices, faces, location);
            PolyhedralMesh mesh{};
            mesh._keepAlive = std::move(keepAlive);
            if (location == MemoryLocation::DEVICE) {
                mesh._deviceMesh.vertices = detail::convertVector3<double, kokkos::DeviceSpace>(
                        vertices, vertexCount, "polyhedralGravity::vertices");
                mesh._deviceMesh.faces = detail::convertVector3<size_t, kokkos::DeviceSpace>(
                        faces, faceCount, "polyhedralGravity::faces");
                // The host mirror stays empty and is only materialized if somebody asks for host access
                return mesh;
            }
            mesh._hostMesh.vertices = detail::convertVector3<double, Kokkos::DefaultHostExecutionSpace>(
                    vertices, vertexCount, "polyhedralGravity::vertices");
            mesh._hostMesh.faces = detail::convertVector3<size_t, Kokkos::DefaultHostExecutionSpace>(
                    faces, faceCount, "polyhedralGravity::faces");
            mesh.uploadHostMesh();
            return mesh;
        }

        /**
         * The mesh in the memory space of the compute device.
         * @return the device-resident views
         */
        [[nodiscard]] const kokkos::PolyhedralMeshView<double, kokkos::DeviceMemory> &getDeviceMesh() const;

        /**
         * The mesh in the memory space of the host, downloading it from the device if that has not happened
         * yet.
         * @return the host-resident views
         */
        [[nodiscard]] const kokkos::PolyhedralMeshView<double, kokkos::HostMemory> &getHostMesh() const;

        /**
         * The number of vertices of the polyhedron.
         * @return the number of vertices
         */
        [[nodiscard]] size_t countVertices() const;

        /**
         * The number of triangular faces of the polyhedron.
         * @return the number of faces
         */
        [[nodiscard]] size_t countFaces() const;

        /**
         * Returns all vertices as a host vector, i.e. as a copy in the library's public representation.
         * @return the vertices as cartesian coordinates
         */
        [[nodiscard]] std::vector<Array3> getVertices() const;

        /**
         * Returns all faces as a host vector, i.e. as a copy in the library's public representation.
         * @return the faces, each referencing three vertices by index
         */
        [[nodiscard]] std::vector<IndexArray3> getFaces() const;

        /**
         * Replaces all faces and propagates them to the device.
         *
         * The faces are always written into a fresh allocation, so that the buffer of a calling array
         * library is never modified. This is what healing an inconsistent vertex ordering uses.
         *
         * @param faces the new faces, each referencing three vertices by index
         * @throws std::invalid_argument if the number of faces differs from the current one
         */
        void setFaces(const std::vector<IndexArray3> &faces);

        /**
         * Shifts every face index by -1 if the mesh's indexing starts at one instead of at zero.
         *
         * The polyhedral file formats disagree on whether the first vertex is vertex zero or vertex one.
         * Tsoulis' algorithm needs zero based indices, and a mesh whose faces never reference vertex zero
         * is taken to be one based. The shift runs where the faces already are, i.e. on the device for a
         * mesh built from a device pointer, and always writes into a fresh allocation so that the buffer of
         * a calling array library is never modified.
         *
         * @return true if the indexing was shifted, false if it already started at zero
         */
        bool shiftFaceIndicesToZeroBased();

    private:
        /** Copies the device mesh back into the host mirror */
        void downloadDeviceMesh() const;

        /**
         * Copies the host mesh into the device mesh, allocating the latter if necessary.
         * On a build without a GPU backend both are the same view and nothing is copied.
         */
        void uploadHostMesh();

        /**
         * Rejects buffers which cannot be turned into a mesh.
         * @param vertices the vertex buffer
         * @param faces the face buffer
         * @param location the memory space the buffers live in
         * @throws std::invalid_argument if a buffer is a null pointer
         * @throws std::runtime_error if DEVICE is given, but this build has no GPU backend
         */
        static void checkBuffers(const void *vertices, const void *faces, MemoryLocation location);
    };

}// namespace polyhedralGravity
