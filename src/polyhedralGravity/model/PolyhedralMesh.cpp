#include "PolyhedralMesh.h"

#include <limits>
#include <stdexcept>

namespace polyhedralGravity {

    PolyhedralMesh::PolyhedralMesh(const std::vector<Array3> &vertices, const std::vector<IndexArray3> &faces) {
        kokkos::ensureInitialized();
        _hostMesh.vertices = kokkos::Vector3View<double, kokkos::HostMemory>{
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::vertices"), vertices.size()};
        _hostMesh.faces = kokkos::Vector3View<size_t, kokkos::HostMemory>{
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::faces"), faces.size()};
        for (size_t index = 0; index < vertices.size(); ++index) {
            for (size_t component = 0; component < 3; ++component) {
                _hostMesh.vertices(index, component) = vertices[index][component];
            }
        }
        for (size_t index = 0; index < faces.size(); ++index) {
            for (size_t component = 0; component < 3; ++component) {
                _hostMesh.faces(index, component) = faces[index][component];
            }
        }
        uploadHostMesh();
    }

    const kokkos::PolyhedralMeshView<double, kokkos::DeviceMemory> &PolyhedralMesh::getDeviceMesh() const {
        return _deviceMesh;
    }

    const kokkos::PolyhedralMeshView<double, kokkos::HostMemory> &PolyhedralMesh::getHostMesh() const {
        // The mirror is missing exactly if this mesh was built from a device pointer and nobody asked for
        // host access since. On a build without a GPU backend the mirror is the device mesh itself.
        if (_hostMesh.vertices.data() == nullptr && _deviceMesh.vertices.data() != nullptr) {
            downloadDeviceMesh();
        }
        return _hostMesh;
    }

    size_t PolyhedralMesh::countVertices() const {
        return _deviceMesh.countVertices();
    }

    size_t PolyhedralMesh::countFaces() const {
        return _deviceMesh.countFaces();
    }

    std::vector<Array3> PolyhedralMesh::getVertices() const {
        const auto &mesh = getHostMesh();
        std::vector<Array3> vertices(mesh.countVertices());
        for (size_t index = 0; index < vertices.size(); ++index) {
            vertices[index] = mesh.getVertex(index);
        }
        return vertices;
    }

    std::vector<IndexArray3> PolyhedralMesh::getFaces() const {
        const auto &mesh = getHostMesh();
        std::vector<IndexArray3> faces(mesh.countFaces());
        for (size_t index = 0; index < faces.size(); ++index) {
            faces[index] = mesh.getFace(index);
        }
        return faces;
    }

    void PolyhedralMesh::setFaces(const std::vector<IndexArray3> &faces) {
        if (faces.size() != countFaces()) {
            throw std::invalid_argument{"The number of faces of a polyhedral mesh cannot be changed!"};
        }
        // Deliberately a fresh allocation: the previous faces may alias a foreign buffer, which must not be
        // modified by this library
        kokkos::Vector3View<size_t, kokkos::HostMemory> newFaces{
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::faces"), faces.size()};
        for (size_t index = 0; index < faces.size(); ++index) {
            for (size_t component = 0; component < 3; ++component) {
                newFaces(index, component) = faces[index][component];
            }
        }
        _hostMesh.faces = newFaces;
        _deviceMesh.faces = Kokkos::create_mirror_view_and_copy(kokkos::DeviceMemory{}, newFaces);
    }

    bool PolyhedralMesh::shiftFaceIndicesToZeroBased() {
        const kokkos::Vector3View<size_t, kokkos::DeviceMemory> faces = _deviceMesh.faces;
        const size_t faceCount = faces.extent(0);
        if (faceCount == 0) {
            return false;
        }
        // Checking that the vertex with index zero is actually used is the same as asking for the smallest
        // index appearing in any face, since the indices are unsigned
        size_t smallestIndex = std::numeric_limits<size_t>::max();
        Kokkos::parallel_reduce(
                "polyhedralGravity::smallestFaceIndex", Kokkos::RangePolicy<kokkos::DeviceSpace>(0, faceCount),
                KOKKOS_LAMBDA(const size_t index, size_t &accumulator) {
                    for (size_t component = 0; component < 3; ++component) {
                        accumulator = Kokkos::min(accumulator, faces(index, component));
                    }
                },
                Kokkos::Min<size_t>(smallestIndex));
        if (smallestIndex == 0) {
            return false;
        }
        kokkos::Vector3View<size_t, kokkos::DeviceMemory> shiftedFaces{
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "polyhedralGravity::faces"), faceCount};
        Kokkos::parallel_for(
                "polyhedralGravity::shiftFaceIndices", Kokkos::RangePolicy<kokkos::DeviceSpace>(0, faceCount),
                KOKKOS_LAMBDA(const size_t index) {
                    for (size_t component = 0; component < 3; ++component) {
                        shiftedFaces(index, component) = faces(index, component) - 1;
                    }
                });
        Kokkos::fence();
        _deviceMesh.faces = shiftedFaces;
        // Only refresh the host mirror if it exists already, otherwise it stays lazy
        if (_hostMesh.vertices.data() != nullptr) {
            _hostMesh.faces = Kokkos::create_mirror_view_and_copy(kokkos::HostMemory{}, shiftedFaces);
        }
        return true;
    }

    void PolyhedralMesh::downloadDeviceMesh() const {
        _hostMesh.vertices = Kokkos::create_mirror_view_and_copy(kokkos::HostMemory{}, _deviceMesh.vertices);
        _hostMesh.faces = Kokkos::create_mirror_view_and_copy(kokkos::HostMemory{}, _deviceMesh.faces);
    }

    void PolyhedralMesh::uploadHostMesh() {
        _deviceMesh.vertices = Kokkos::create_mirror_view_and_copy(kokkos::DeviceMemory{}, _hostMesh.vertices);
        _deviceMesh.faces = Kokkos::create_mirror_view_and_copy(kokkos::DeviceMemory{}, _hostMesh.faces);
    }

    void PolyhedralMesh::checkBuffers(const void *vertices, const void *faces, const MemoryLocation location) {
        if (vertices == nullptr || faces == nullptr) {
            throw std::invalid_argument{"A polyhedral mesh cannot be built from a null pointer!"};
        }
        if (location == MemoryLocation::DEVICE && !kokkos::GPU_AVAILABLE) {
            throw std::runtime_error{
                    "The polyhedral mesh was given a device pointer, but this build of the polyhedral gravity "
                    "model has no GPU backend. It was compiled with the execution spaces " +
                    kokkos::getEnabledExecutionSpaces() +
                    ". Hand in a host buffer, e.g. a NumPy array or a PyTorch tensor on the CPU, instead."};
        }
    }

}// namespace polyhedralGravity
