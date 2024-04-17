#include "Polyhedron.h"

#include "polyhedralGravity/input/TetgenAdapter.h"


namespace polyhedralGravity {

    Polyhedron::Polyhedron(const std::vector<Array3> &vertices,
            const std::vector<IndexArray3> &faces, double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity)
        : _vertices{vertices},
          _faces{faces},
          _density{density},
          _orientation{orientation} {
        //Checks that the node with index zero is actually used
        if (_faces.end() == std::find_if(_faces.begin(), _faces.end(), [&](auto &face) {
            return face[0] == 0 || face[1] == 0 || face[2] == 0;
        })) {
            throw std::invalid_argument("The node with index zero (0) was never used in any face! This is "
                    "no valid polyhedron. Probable issue: Started numbering the vertices of "
                    "the polyhedron at one (1).");
        }
        this->runIntegrityMeasures(integrity);
    }

    Polyhedron::Polyhedron(const PolyhedralSource &polyhedralSource, double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity)
        : Polyhedron{std::get<std::vector<Array3>>(polyhedralSource), std::get<std::vector<IndexArray3>>(polyhedralSource), density, orientation, integrity} {}

    Polyhedron::Polyhedron(const PolyhedralFiles &polyhedralFiles, double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity)
        : Polyhedron{TetgenAdapter{polyhedralFiles}.getPolyhedralSource(), density, orientation, integrity} {}

    Polyhedron::Polyhedron(const std::variant<PolyhedralSource, PolyhedralFiles> &polyhedralSource, double density, const NormalOrientation &orientation, const PolyhedronIntegrity &integrity)
    : Polyhedron{std::holds_alternative<PolyhedralSource>(polyhedralSource) ? std::get<PolyhedralSource>(polyhedralSource) : TetgenAdapter{std::get<PolyhedralFiles>(polyhedralSource)}.getPolyhedralSource(),
        density, orientation, integrity} {}

    const std::vector<Array3> &Polyhedron::getVertices() const {
        return _vertices;
    }

    const Array3 &Polyhedron::getVertex(size_t index) const {
        return _vertices[index];
    }

    size_t Polyhedron::countVertices() const {
        return _vertices.size();
    }

    const std::vector<IndexArray3> &Polyhedron::getFaces() const {
        return _faces;
    }

    const IndexArray3 &Polyhedron::getFace(size_t index) const {
        return _faces[index];
    }

    Array3Triplet Polyhedron::getResolvedFace(size_t index) const {
        return {_vertices[_faces[index][0]], _vertices[_faces[index][1]], _vertices[_faces[index][2]]};
    }

    size_t Polyhedron::countFaces() const {
        return _faces.size();
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

    std::string Polyhedron::toString() const {
        std::stringstream sstream{};
        sstream << "<polyhedral_gravity.Polyhedron, density=" << _density << ", vertices= "
           << countVertices() << ", faces= " << countFaces() << ", orientation= " << _orientation << ">";
        return sstream.str();
    }

    std::tuple<std::vector<Array3>, std::vector<IndexArray3>, double, NormalOrientation> Polyhedron::getState() const {
        return std::make_tuple(_vertices, _faces, _density, _orientation);
    }

    std::pair<NormalOrientation, std::set<size_t>> Polyhedron::checkPlaneUnitNormalOrientation() const {
        // 1. Step: Find all indices of normals which vioate the constraint outwards pointing
        auto indexIterator = thrust::make_counting_iterator<size_t>(0);
        size_t numberOfFaces = this->countFaces();
        std::set<size_t> violatingOutwards{};
        thrust::copy_if(
                thrust::device,
                indexIterator,
                indexIterator + numberOfFaces,
                std::inserter(violatingOutwards, violatingOutwards.end()),
                [this](const size_t index) {
                    // If the ray intersects the polyhedron odd number if times the normal points inwards
                    // Hence, violating the OUTWARDS constraint
                    const size_t intersects = countRayPolyhedronIntersections(this->getResolvedFace(index));
                    return intersects % 2 != 0;
                });
        // 2. Step: If we have more than n/2 violations, the majority is inwards pointing --> Build complemntary set!
        if (violatingOutwards.size() > numberOfFaces / 2) {
            std::set<size_t> fullSet{};
            std::copy_n(indexIterator, numberOfFaces, std::inserter(fullSet, fullSet.end()));
            std::set<size_t> violatingInwards{};
            std::set_difference(
                    fullSet.begin(), fullSet.end(),
                    violatingOutwards.begin(), violatingOutwards.end(),
                    std::inserter(violatingInwards, violatingInwards.end())
                    );
            // 3a. Step: Return the inwards pointing as major orientation and
            // the violating faces, i.e. which have outwards pointing normals
            return std::make_pair(NormalOrientation::INWARDS, violatingInwards);
        }
        // 3b. Step: Return the outwards pointing as major orientation
        // and the violating faces, i.e. which have inwards pointing normals
        return std::make_pair(NormalOrientation::OUTWARDS, violatingOutwards);
    }

    void Polyhedron::runIntegrityMeasures(const PolyhedronIntegrity &integrity) {
        using util::operator<<;
        switch (integrity) {
            case PolyhedronIntegrity::DISABLE:
                return;
            case PolyhedronIntegrity::AUTOMATIC:
                SPDLOG_LOGGER_WARN(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                                   "The mesh check is enabled and analyzes the polyhedron for degnerated faces & "
                                   "that all plane unit normals point in the specified direction. This checks requires "
                                   "a quadratic runtime cost which is most of the time not desirable. "
                                   "Please explcity enable or disable this by setting 'check' to true orfalse");
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
                                << ". You can either reconstruct the reconstruct the polyhedron with the ortientation set to " << actualOrientation
                                << ". Alternativly, you can reconstruct with the integrity check set to HEAL";
                    } else {
                        sstream << "The actual majority orientation is " << actualOrientation
                                << ". You can either:\n 1) Fix the ordering of the following faces:\n"
                                << violatingIndices << '\n'
                                << "2) Or you reconstruct the polyhedron using the integrity check set to HEAL.";
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
        const auto &[begin, end] = this->transformIterator();
        // All triangles surface area needs to be greater than zero
        return thrust::transform_reduce(
                thrust::device,
                begin, end, [](const Array3Triplet &face) {
                    return util::surfaceArea(face) > 0.0;
                }, true, thrust::logical_and<bool>());
    }

    void Polyhedron::healPlaneUnitNormalOrientation(const NormalOrientation &actualOrientation, const std::set<size_t> &violatingIndices) {
        // Assign the majority plane unit normal orientation
        _orientation = actualOrientation;
        // Fix the vioalting faces by exchaning the vertex ordering (exchaning index 0 with index 1 in the face)
        std::for_each(violatingIndices.cbegin(), violatingIndices.cend(), [this](size_t i) {
            std::swap(this->_faces[i][0], this->_faces[i][1]);
        });
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
        const Array3 rayOrigin = centroid + (rayVector * EPSILON);

        // Count every triangular face which is intersected by the ray
        const auto &[begin, end] = this->transformIterator();
        std::set<Array3> intersections{};
        std::for_each(begin, end, [&rayOrigin, &rayVector, &intersections](const Array3Triplet &otherFace) {
            const std::unique_ptr<Array3> intersection = rayIntersectsTriangle(rayOrigin, rayVector, otherFace);
            if (intersection != nullptr) {
                intersections.insert(*intersection);
            }
        });
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
        if (a > -EPSILON && a < EPSILON) {
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
        if (t > EPSILON) {
            return std::make_unique<Array3>(rayOrigin + rayVector * t);
        } else {
            return nullptr;
        }
    }

}