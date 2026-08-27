#include "Polyhedron.h"

#include <limits>
#include <stdexcept>

#include <Kokkos_Core.hpp>

#include "polyhedralGravity/util/KokkosSession.h"

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

    Polyhedron::Polyhedron(PolyhedralMesh mesh, const double density, const NormalOrientation &orientation,
                           const PolyhedronIntegrity &integrity, const MetricUnit &metricUnit)
        : _mesh{std::move(mesh)},
          _density{density},
          _orientation{orientation},
          _metricUnit{metricUnit} {
        // Checks that the node with index zero is actually used
        // In case it is not used, the indexing presumably starts mathematically at one
        // In this case, we shift it by -1, so that the indexing start with zero
        if (_mesh.shiftFaceIndicesToZeroBased()) {
            POLYHEDRAL_GRAVITY_LOG_DEBUG("The indexing of the polyhedron's vertices seems to start at 1 instead of 0. The faces array is modfied accordingly!");
        }
        this->runIntegrityMeasures(integrity);
    }

    Polyhedron::Polyhedron(const std::vector<Array3> &vertices,
                           const std::vector<IndexArray3> &faces, const double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity, const MetricUnit& metricUnit)
        : Polyhedron{PolyhedralMesh{vertices, faces}, density, orientation, integrity, metricUnit} {
    }

    Polyhedron::Polyhedron(const PolyhedralSource &polyhedralSource, const double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity, const MetricUnit& metricUnit)
        : Polyhedron{std::get<std::vector<Array3>>(polyhedralSource), std::get<std::vector<IndexArray3>>(polyhedralSource), density, orientation, integrity, metricUnit} {
    }

    Polyhedron::Polyhedron(const PolyhedralFiles &polyhedralFiles, const double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity, const MetricUnit& metricUnit)
        : Polyhedron{MeshReader::getPolyhedralSource(polyhedralFiles), density, orientation, integrity, metricUnit} {
    }

    Polyhedron::Polyhedron(const std::variant<PolyhedralSource, PolyhedralFiles> &polyhedralSource, const double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity, const MetricUnit& metricUnit)
        : Polyhedron{std::holds_alternative<PolyhedralSource>(polyhedralSource) ? std::get<PolyhedralSource>(polyhedralSource) : MeshReader::getPolyhedralSource(std::get<PolyhedralFiles>(polyhedralSource)),
                     density, orientation, integrity, metricUnit} {
    }

    const PolyhedralMesh &Polyhedron::getMesh() const {
        return _mesh;
    }

    std::vector<Array3> Polyhedron::getVertices() const {
        return _mesh.getVertices();
    }

    Array3 Polyhedron::getVertex(size_t index) const {
        return _mesh.getHostMesh().getVertex(index);
    }

    size_t Polyhedron::countVertices() const {
        return _mesh.countVertices();
    }

    std::vector<IndexArray3> Polyhedron::getFaces() const {
        return _mesh.getFaces();
    }

    IndexArray3 Polyhedron::getFace(size_t index) const {
        return _mesh.getHostMesh().getFace(index);
    }

    Array3Triplet Polyhedron::getResolvedFace(size_t index) const {
        return _mesh.getHostMesh().resolveFace(index, Array3{0.0, 0.0, 0.0});
    }

    size_t Polyhedron::countFaces() const {
        return _mesh.countFaces();
    }

    double Polyhedron::getDensity() const {
        return _density;
    }

    void Polyhedron::setDensity(double density) {
        _density = density;
    }

    NormalOrientation Polyhedron::getOrientation() const {
        return _orientation;
    }

    double Polyhedron::getOrientationFactor() const {
        return _orientation == NormalOrientation::OUTWARDS ? 1.0 : -1.0;
    }

    MetricUnit Polyhedron::getMeshUnit() const {
        return _metricUnit;
    }

    std::string Polyhedron::getMeshUnitAsString() const {
        std::stringstream meshUnit{};
        meshUnit << _metricUnit;
        return meshUnit.str();
    }

    std::string Polyhedron::getDensityUnit() const {
        std::stringstream densityUnit{};
        if (_metricUnit != MetricUnit::UNITLESS) {
            densityUnit << "kg/" << _metricUnit << "^3";
        } else {
            densityUnit << _metricUnit;
        }
        return densityUnit.str();
    }


    double Polyhedron::getGravityModelScaling() const {
        switch (_metricUnit) {
            case MetricUnit::UNITLESS:
                return getDensity() * getOrientationFactor();
            case MetricUnit::METER:
                return util::GRAVITATIONAL_CONSTANT * getDensity() * getOrientationFactor();
            case MetricUnit::KILOMETER:
                // Gravitational Constant in km^3/(kg * s^2)
                constexpr double GRAVITATIONAL_CONSTANT_KM = util::GRAVITATIONAL_CONSTANT * 1e-9;
                return GRAVITATIONAL_CONSTANT_KM * getDensity() * getOrientationFactor();
        }
        throw std::invalid_argument{"The metric unit is not supported!"};
    }


    std::string Polyhedron::toString() const {
        std::stringstream sstream{};
        sstream << "<polyhedral_gravity.Polyhedron, density = " << _density << " " << getDensityUnit()
                << ", vertices = " << countVertices()
                << ", faces = " << countFaces()
                << ", orientation = " << _orientation
                << ", mesh_unit = '" << getMeshUnitAsString() << "'"
                << ">";
        return sstream.str();
    }

    std::tuple<std::vector<Array3>, std::vector<IndexArray3>, double, NormalOrientation, MetricUnit> Polyhedron::getState() const {
        return std::make_tuple(_mesh.getVertices(), _mesh.getFaces(), _density, _orientation, _metricUnit);
    }

    std::pair<NormalOrientation, std::set<size_t>> Polyhedron::checkPlaneUnitNormalOrientation() const {
        kokkos::ensureInitialized();
        // 1. Step: Find all indices of normals which vioate the constraint outwards pointing
        const size_t n = this->countFaces();
        // Entry is 1 if the corresponding index VIOLATES the OUTWARDS criteria
        // Entry is 0 if the corresponding index FULFILLS the OUTWARDS criteria
        // This is a char and not a bool since the threads below write into distinct elements concurrently
        std::vector<char> violatesOutwards(n, 0);
        size_t numberOfOutwardsViolations = 0;
        // The ray casting below is host-only code, so it runs on the host execution space no matter which
        // compute backend the user later evaluates the gravity model on
        Kokkos::parallel_reduce(
                "polyhedralGravity::checkPlaneUnitNormalOrientation",
                Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
                [&](const size_t index, size_t &violations) {
                    // If the ray intersects the polyhedron odd number of times the normal points inwards
                    // Hence, violating the OUTWARDS constraint
                    const size_t intersects = this->countRayPolyhedronIntersections(this->getResolvedFace(index));
                    violatesOutwards[index] = static_cast<char>(intersects % 2 != 0);
                    violations += violatesOutwards[index];
                },
                numberOfOutwardsViolations);

        // 2. Step: Create a set with only the indices violating the constraint
        // 3a. Step: If the majority points outwards, return it as the major orientation
        // and the violating faces, i.e. which have inwards pointing normals (and vice versa in 3b.)
        const bool majorityIsOutwards = numberOfOutwardsViolations <= n / 2;
        std::set<size_t> violatingIndices{};
        for (size_t index = 0; index < n; ++index) {
            if (static_cast<bool>(violatesOutwards[index]) == majorityIsOutwards) {
                violatingIndices.insert(index);
            }
        }
        return std::make_pair(majorityIsOutwards ? NormalOrientation::OUTWARDS : NormalOrientation::INWARDS,
                              violatingIndices);
    }

    void Polyhedron::runIntegrityMeasures(const PolyhedronIntegrity &integrity) {
        using util::operator<<;
        switch (integrity) {
            case PolyhedronIntegrity::DISABLE:
                return;
            case PolyhedronIntegrity::AUTOMATIC:
                POLYHEDRAL_GRAVITY_LOG_WARN("The mesh check is enabled and analyzes the polyhedron for degnerated faces & "
                                            "that all plane unit normals point in the specified direction. This checks requires "
                                            "a quadratic runtime cost which is most of the time not desirable. "
                                            "Please explicitly set the integrity_check to either VERIFY, HEAL or DISABLE."
                                            "You can find further details in the documentation!");
            // NO BREAK! AUTOMATIC implies VERIFY, but with a info mesage to explcitly set the option
            case PolyhedronIntegrity::VERIFY:
            // NO BREAK! VERIFY terminates earlier, but does in the beginning the same as HEAL
            case PolyhedronIntegrity::HEAL:
                if (!this->checkTrianglesNotDegenerated()) {
                    throw std::invalid_argument{"At least on triangle in the mesh is degenerated and its surface area equals zero!"};
                }
                const auto &[actualOrientation, violatingIndices] = this->checkPlaneUnitNormalOrientation();
                if (actualOrientation != _orientation || !violatingIndices.empty()) {
                    std::stringstream sstream{};
                    sstream << "The plane unit normals are not all pointing in the specified direction " << _orientation << '\n';
                    if (violatingIndices.empty()) {
                        sstream << "Instead all plane unit normals are pointing "
                                << actualOrientation
                                << ". You can either reconstruct the polyhedron with the orientation set to " << actualOrientation
                                << ". Alternativly, you can reconstruct with the inetgrity_check set to HEAL";
                    } else {
                        sstream << "The actual majority orientation of the polyhedron's normals is " << actualOrientation
                                << ". You can either:\n 1) Fix the ordering of the following faces:\n"
                                << violatingIndices << '\n'
                                << "2) Or you reconstruct the polyhedron using the integrity_check set to HEAL.";
                    }
                    // In case of HEAL, don't throw but repair
                    if (integrity != PolyhedronIntegrity::HEAL) {
                        throw std::invalid_argument(sstream.str());
                    } else {
                        this->healPlaneUnitNormalOrientation(actualOrientation, violatingIndices);
                    }
                }
        }
    }

    bool Polyhedron::checkTrianglesNotDegenerated() const {
        kokkos::ensureInitialized();
        // All triangles surface area needs to be greater than zero
        size_t degenerated = 0;
        Kokkos::parallel_reduce(
                "polyhedralGravity::checkTrianglesNotDegenerated",
                Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, this->countFaces()),
                [&](const size_t index, size_t &accumulator) {
                    accumulator += util::surfaceArea(this->getResolvedFace(index)) > 0.0 ? 0 : 1;
                },
                degenerated);
        return degenerated == 0;
    }

    void Polyhedron::healPlaneUnitNormalOrientation(const NormalOrientation &actualOrientation, const std::set<size_t> &violatingIndices) {
        // Assign the majority plane unit normal orientation
        _orientation = actualOrientation;
        // Fix the vioalting faces by exchaning the vertex ordering (exchaning index 0 with index 1 in the face)
        // The faces are written back as a whole, since they may alias a buffer of the caller which this
        // library must not modify in place
        std::vector<IndexArray3> faces = _mesh.getFaces();
        std::for_each(violatingIndices.cbegin(), violatingIndices.cend(), [&faces](size_t i) {
            std::swap(faces[i][0], faces[i][1]);
        });
        _mesh.setFaces(faces);
    }

    size_t Polyhedron::countRayPolyhedronIntersections(const Array3Triplet &face) const {
        using namespace util;
        // The centroid of the triangular face
        const Array3 centroid = (face[0] + face[1] + face[2]) / 3.0;

        // The normal of the plane calculated with two segments of the triangle
        // The normal is the rayVector starting at the rayOrigin
        const Array3 segmentVector1 = face[1] - face[0];
        const Array3 segmentVector2 = face[2] - face[1];
        const Array3 rayVector = normal(segmentVector1, segmentVector2);

        // The origin of the array has a slight offset in direction of the normal
        const Array3 rayOrigin = centroid + (rayVector * EPSILON_ZERO_OFFSET);

        // Count every triangular face which is intersected by the ray
        std::set<Array3> intersections{};
        for (size_t index = 0; index < this->countFaces(); ++index) {
            const std::unique_ptr<Array3> intersection =
                    rayIntersectsTriangle(rayOrigin, rayVector, this->getResolvedFace(index));
            if (intersection != nullptr) {
                intersections.insert(*intersection);
            }
        }
        return intersections.size();
    }

    std::unique_ptr<Array3> Polyhedron::rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector, const Array3Triplet &triangle) {
        // Adapted Möller–Trumbore intersection algorithm
        // see https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
        using namespace util;
        const Array3 edge1 = triangle[1] - triangle[0];
        const Array3 edge2 = triangle[2] - triangle[0];
        const Array3 h = cross(rayVector, edge2);
        const double a = dot(edge1, h);
        if (a > -EPSILON_ZERO_OFFSET && a < EPSILON_ZERO_OFFSET) {
            return nullptr;
        }

        const double f = 1.0 / a;
        const Array3 s = rayOrigin - triangle[0];
        const double u = f * dot(s, h);
        if (u < 0.0 || u > 1.0) {
            return nullptr;
        }

        const Array3 q = cross(s, edge1);
        const double v = f * dot(rayVector, q);
        if (v < 0.0 || u + v > 1.0) {
            return nullptr;
        }

        const double t = f * dot(edge2, q);
        if (t > EPSILON_ZERO_OFFSET) {
            return std::make_unique<Array3>(rayOrigin + rayVector * t);
        } else {
            return nullptr;
        }
    }

}// namespace polyhedralGravity