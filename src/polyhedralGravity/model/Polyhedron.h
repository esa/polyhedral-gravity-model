#pragma once

#include "polyhedralGravity/input/MeshReader.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/PolyhedralMesh.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/output/Logging.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/util/UtilityFloatArithmetic.h"
#include <algorithm>
#include <array>
#include <exception>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace polyhedralGravity {

    /* Forward declaration of Polyhedron */
    class Polyhedron;

    /**
     * Data structure containing the model data of one polyhedron. This includes nodes, edges (faces) and elements.
     * The index always starts with zero!
     *
     * @note Constructing a Polyhedron puts its mesh into the memory of the compute devices and therefore
     * initializes the Kokkos runtime. Kokkos may only be initialized after the program has entered main(),
     * so a Polyhedron must not be a global or a static class member. A function local static, which is
     * constructed on its first use, is fine.
    */
    class Polyhedron {

        /**
         * The vertices and the triangular faces of the polyhedron, as Kokkos views in the memory of the
         * compute devices.
         *
         * Each vertex is an array of size three containing the xyz coordinates. The mesh must be scaled in
         * the same units as the density is given (the unit must match to the mesh, e.g., mesh in @f$[m]@f$
         * requires density in @f$[kg/m^3]@f$).
         *
         * Each face is an array of size three containing the indices of the nodes forming the face.
         * Since every face consists of three nodes, every face consists of three segments. Each segment
         * consists of two nodes. For example, a face consisting of {1, 2, 3} --> segments: {1, 2}, {2, 3},
         * {3, 1}.
         *
         * Copying a Polyhedron shares this mesh instead of duplicating it, since a Kokkos view is a
         * reference counted handle. That is sound because the mesh is only ever modified while a Polyhedron
         * is being constructed, i.e. before anybody else can hold a copy.
         */
        PolyhedralMesh _mesh;

        /** The constant density of the polyhedron (the unit must match to the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$) */
        double _density;

        /**
         * An optional variable representing the orientation of the polyhedron.
         * This variable holds an instance of the NormalOrientation enumeration class,
         * which indicates whether the plane unit normals are pointing outwards or inwards.
         * To set the orientation of the polyhedron, assign a valid `NormalOrientation` value to this variable.
         * If the value is not set, the orientation is considered to be unspecified.
         */
        NormalOrientation _orientation;

        /** Metric Unit of the Vertices Coordinates. One of METER, KILOMETER, or UNITLESS */
        const MetricUnit _metricUnit;

    public:
        /**
         * Generates a polyhedron from nodes and faces.
         * If the indexing of the polyhedron's vertices in the faces' array starts with one, it is shifted so that it starts with zero.
         * @param vertices a vector of nodes
         * @param faces a vector of faces containing the formation of faces off vertices
         * @param density the density of the polyhedron in @f$[kg/X^3]@f$.
         *          It must match the unit of the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         * @param orientation specify if the plane unit normals point outwards or inwards (default: OUTWARDS)
         * @param integrity specify if the mesh input is checked/ healed to fulfill the constraints of Tsoulis' algorithm (see {@link PolyhedronIntegrity})
         * @param metricUnit specify the mesh's coordinate scale's unit. Can be kilometer, meter, or unitless (defaults to meter)
         *
         * @throws std::invalid_argument depending on the {@link integrity} flag
         */
        Polyhedron(
                const std::vector<Array3> &vertices,
                const std::vector<IndexArray3> &faces,
                double density,
                const NormalOrientation &orientation = NormalOrientation::OUTWARDS,
                const PolyhedronIntegrity &integrity = PolyhedronIntegrity::AUTOMATIC,
                const MetricUnit &metricUnit = MetricUnit::METER
                );

        /**
         * Generates a polyhedron from nodes and faces.
         * If the indexing of the polyhedron's vertices in the faces' array starts with one, it is shifted so that it starts with zero.
         * @param polyhedralSource a tuple of vector containing the nodes and triangular faces.
         * @param density the density of the polyhedron in @f$[kg/X^3]@f$.
         *          It must match the unit of the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         * @param orientation specify if the plane unit normals point outwards or inwards (default: OUTWARDS)
         * @param integrity specify if the mesh input is checked/ healed to fulfill the constraints of Tsoulis' algorithm (see {@link PolyhedronIntegrity})
         * @param metricUnit specify the mesh's coordinate scale's unit. Can be kilometer, meter, or unitless (defaults to meter)
         *
         * @throws std::invalid_argument depending on the {@link integrity} flag
         */
        Polyhedron(
                const PolyhedralSource &polyhedralSource,
                double density,
                const NormalOrientation &orientation = NormalOrientation::OUTWARDS,
                const PolyhedronIntegrity &integrity = PolyhedronIntegrity::AUTOMATIC,
                const MetricUnit &metricUnit = MetricUnit::METER
                );

        /**
         * Generates a polyhedron from nodes and faces.
         * If the indexing of the polyhedron's vertices in the faces' array starts with one, it is shifted so that it starts with zero.
         * @param polyhedralFiles a list of files (see {@link TetgenAdapter}
         * @param density the density of the polyhedron in @f$[kg/X^3]@f$.
         *          It must match the unit of the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         * @param orientation specify if the plane unit normals point outwards or inwards (default: OUTWARDS)
         * @param integrity specify if the mesh input is checked/ healed to fulfill the constraints of Tsoulis' algorithm (see {@link PolyhedronIntegrity})
         * @param metricUnit specify the mesh's coordinate scale's unit. Can be kilometer, meter, or unitless (defaults to meter)
         *
         * @throws std::invalid_argument depending on the {@link integrity} flag
         */
        Polyhedron(const PolyhedralFiles &polyhedralFiles, double density,
                   const NormalOrientation &orientation = NormalOrientation::OUTWARDS,
                   const PolyhedronIntegrity &integrity = PolyhedronIntegrity::AUTOMATIC,
                   const MetricUnit &metricUnit = MetricUnit::METER
                                   );

        /**
         * Generates a polyhedron from nodes and faces.
         * If the indexing of the polyhedron's vertices in the faces' array starts with one, it is shifted so that it starts with zero.
         * This constructor using a variant is mainly utilized from the Python Interface.
         * @param polyhedralSource a list of files (see {@link TetgenAdapter} or a tuple of vector containing the nodes and triangular faces.
         * @param density the density of the polyhedron in @f$[kg/X^3]@f$.
         *          It must match the unit of the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         * @param orientation specify if the plane unit normals point outwards or inwards (default: OUTWARDS)
         * @param integrity specify if the mesh input is checked/ healed to fulfill the constraints of Tsoulis' algorithm (see {@link PolyhedronIntegrity})
         * @param metricUnit specify the mesh's coordinate scale's unit. Can be kilometer, meter, or unitless (defaults to meter)
         *
         * @throws std::invalid_argument depending on the {@link integrity} flag
         */
        Polyhedron(const std::variant<PolyhedralSource, PolyhedralFiles> &polyhedralSource, double density,
                   const NormalOrientation &orientation = NormalOrientation::OUTWARDS,
                   const PolyhedronIntegrity &integrity = PolyhedronIntegrity::AUTOMATIC,
                   const MetricUnit &metricUnit = MetricUnit::METER
                   );

        /**
         * Generates a polyhedron from an already built mesh.
         * If the indexing of the polyhedron's vertices in the faces' array starts with one, it is shifted so that it starts with zero.
         * @param mesh the vertices and triangular faces, possibly aliasing a foreign buffer
         * @param density the density of the polyhedron in @f$[kg/X^3]@f$.
         *          It must match the unit of the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         * @param orientation specify if the plane unit normals point outwards or inwards (default: OUTWARDS)
         * @param integrity specify if the mesh input is checked/ healed to fulfill the constraints of Tsoulis' algorithm (see {@link PolyhedronIntegrity})
         * @param metricUnit specify the mesh's coordinate scale's unit. Can be kilometer, meter, or unitless (defaults to meter)
         *
         * @throws std::invalid_argument depending on the {@link integrity} flag
         */
        Polyhedron(PolyhedralMesh mesh, double density,
                   const NormalOrientation &orientation = NormalOrientation::OUTWARDS,
                   const PolyhedronIntegrity &integrity = PolyhedronIntegrity::AUTOMATIC,
                   const MetricUnit &metricUnit = MetricUnit::METER
                   );

        /**
         * Generates a polyhedron directly from the buffers of an array library, i.e. of NumPy, PyTorch, or
         * JAX, without copying them.
         *
         * The buffers are only copied if their element type differs from the library's own, i.e. if the
         * vertices are not double precision or the face indices are not of the platform's @c size_t. Such a
         * conversion happens in the memory space the buffers already live in, so a mesh which is handed in
         * as a device pointer never travels through the host.
         *
         * @tparam VertexType the element type of the vertex buffer, e.g. float or double
         * @tparam IndexType the element type of the face buffer, e.g. int32_t or size_t
         * @param vertices the first element of a C-contiguous @f$(N, 3)@f$ vertex buffer
         * @param vertexCount the number of vertices N
         * @param faces the first element of a C-contiguous @f$(M, 3)@f$ face buffer
         * @param faceCount the number of faces M
         * @param density the density of the polyhedron in @f$[kg/X^3]@f$
         * @param location whether the buffers live in the host's or in the device's memory (default: HOST)
         * @param orientation specify if the plane unit normals point outwards or inwards (default: OUTWARDS)
         * @param integrity specify if the mesh input is checked/ healed to fulfill the constraints of Tsoulis' algorithm (see {@link PolyhedronIntegrity})
         * @param metricUnit specify the mesh's coordinate scale's unit. Can be kilometer, meter, or unitless (defaults to meter)
         * @param keepAlive an owner of the buffers which is held for the lifetime of this polyhedron
         *
         * @throws std::invalid_argument depending on the {@link integrity} flag
         * @throws std::runtime_error if DEVICE is given, but this build has no GPU backend
         *
         * @note If no keepAlive is given, the caller has to guarantee that the buffers outlive this
         * polyhedron and every {@link GravityEvaluable} built from it.
         */
        template<typename VertexType, typename IndexType>
        Polyhedron(const VertexType *vertices, const size_t vertexCount,
                   const IndexType *faces, const size_t faceCount,
                   const double density,
                   const MemoryLocation location = MemoryLocation::HOST,
                   const NormalOrientation &orientation = NormalOrientation::OUTWARDS,
                   const PolyhedronIntegrity &integrity = PolyhedronIntegrity::AUTOMATIC,
                   const MetricUnit &metricUnit = MetricUnit::METER,
                   std::shared_ptr<void> keepAlive = {})
            : Polyhedron{PolyhedralMesh::fromBuffers(vertices, vertexCount, faces, faceCount, location,
                                                     std::move(keepAlive)),
                         density, orientation, integrity, metricUnit} {
        }

        /**
         * Default destructor
         */
        ~Polyhedron() = default;

        /**
         * Returns the vertices of this polyhedron
         * @return vector of cartesian coordinates
         */
        [[nodiscard]] std::vector<Array3> getVertices() const;

        /**
         * Returns the vertex at a specific index
         * @param index size_t
         * @return cartesian coordinates of the vertex at index
         */
        [[nodiscard]] Array3 getVertex(size_t index) const;

        /**
         * The number of points (nodes) that make up the polyhedron.
         * @return a size_t
         */
        [[nodiscard]] size_t countVertices() const;

        /**
         * Returns the triangular faces of this polyhedron
         * @return vector of triangular faces, where each element size_t references a vertex in the vertices vector
         */
        [[nodiscard]] std::vector<IndexArray3> getFaces() const;

        /**
         * Returns the indices of the vertices making up the face at the given index.
         * @param index size_t
         * @return triplet of the vertices indices forming the face
         */
        [[nodiscard]] IndexArray3 getFace(size_t index) const;

        /**
         * Returns the resolved face with its concrete cartesian coordinates at the given index.
         * @param index size_t
         * @return triplet of vertices' cartesian coordinates
         */
        [[nodiscard]] Array3Triplet getResolvedFace(size_t index) const;

        /**
         * Returns the number of faces (triangles) that make up the polyhedral.
         * @return a size_t
         */
        [[nodiscard]] size_t countFaces() const;

        /**
         * Returns the constant density of this polyhedron.
         * Its unit is @f$[kg/X^3]@f$ with X as the metric unit of the mesh.
         * @return the constant density a double
         */
        [[nodiscard]] double getDensity() const;

        /**
         * Sets the density to a new value. The density's unit must match to the scaling of the mesh.
         * In other words, mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$.
         * @param density the new constant density of the polyhedron
         */
        void setDensity(double density);

        /**
         * Returns the orientation of the plane unit normals of this polyhedron
         * @return OUTWARDS or INWARDS
         */
        [[nodiscard]] NormalOrientation getOrientation() const;

        /**
         * Returns the plane unit normal orientation factor.
         * If the unit normals are outwards pointing, it is 1.0 as Tsoulis intended.
         * If the unit normals are inwards pointing, it is -1.0 (reversed).
         * @return 1.0 or -1.0 depending on plane unit orientation
         */
        [[nodiscard]] double getOrientationFactor() const;

        /**
         * Returns the metric unit of the polyhedron's mesh, which is either METER, KILOMETER, or UNITLESS.
         * @return the metric unit of the polyhedron
         */
        [[nodiscard]] MetricUnit getMeshUnit() const;

        /**
         * Returns the metric unit of the polyhedral mesh.
         * It is either in meter, kilometer, or without units.
         * @return unit of mesh (i.e., the unit of the vertices)
         */
        [[nodiscard]] std::string getMeshUnitAsString() const;

        /**
         * Returns the metric unit of the density which is set during construction and depends on the mesh's unit.
         * @return unit of the density
         */
        [[nodiscard]] std::string getDensityUnit() const;

        /**
         * Returns the scaling factor for the gravity model evaluation.
         * If the unit is UNIT_LESS, it returns the product of density and orientation factor.
         * If the unit is METER, it returns the product of density, orientation factor, and util::GRAVITATIONAL_CONSTANT
         * If the unit is KILOMETER, it returns the product of density, orientation factor, and util::GRAVITATIONAL_CONSTANT in @f$[km^3/(kg * s^2)]@f$
         * @return scaling factor in Tsoulis Model (i.e., gravitational constant and density with additional correction depending on alignment and unit)
         */
        [[nodiscard]] double getGravityModelScaling() const;

        /**
         * Returns a string representation of the Polyhedron.
         * Mainly used for the representation method in the Python interface.
         *
         * @return string representation of the Polyhedron
         */
        [[nodiscard]] std::string toString() const;

        /**
         * Returns the internal data structure of Python pickle support.
         * @return tuple of vertices, faces, density, and normal orientation
         */
        [[nodiscard]] std::tuple<std::vector<Array3>, std::vector<IndexArray3>, double, NormalOrientation, MetricUnit> getState() const;

        /**
         * Returns the resolved face at the given index with all of its vertices shifted by an offset.
         * Tsoulis' equations require the computation point P to be relocated into the origin, which is what
         * this offset is used for.
         * @param index size_t
         * @param offset the offset to subtract from every vertex, e.g. the computation point P
         * @return triplet of the face's vertices' cartesian coordinates, each shifted by -offset
         */
        [[nodiscard]] inline Array3Triplet getResolvedFace(const size_t index, const Array3 &offset) const {
            return _mesh.getHostMesh().resolveFace(index, offset);
        }

        /**
         * Returns the polyhedron's mesh, i.e. its vertices and faces as they live in the memory of the
         * compute devices.
         *
         * This is what a {@link GravityEvaluable} extends into a
         * {@link kokkos::GravitationalMeshView} by adding the caches of Tsoulis' algorithm.
         *
         * @return the mesh of this polyhedron
         */
        [[nodiscard]] const PolyhedralMesh &getMesh() const;

        /**
         * This method determines the majority vertex ordering of a polyhedron and the set of faces which
         * violate the majority constraint and need to be adapted.
         * Hence, if the set is empty, all faces obey to the returned ordering/ plane unit normal orientation.
         *
         * @return a pair consisting of majority ordering (OUTWARDS or INWARDS pointing normals)
         *  and a set of face indices which violate the constraint
         */
        [[nodiscard]] std::pair<NormalOrientation, std::set<size_t>> checkPlaneUnitNormalOrientation() const;

    private:
        /**
         * Checks the integrity of the polyhedron depending on the integrity flag.
         *
         * @param integrity the behavior depends on the value, see {@link PolyhedronIntegrity}
         *
         * @throws std::invalid_argument depending on the {@param integrity} flag
         */
        void runIntegrityMeasures(const PolyhedronIntegrity &integrity);


        /**
         * Checks if no triangle is degenerated by checking the surface area being greater than zero.
         * E.g., two points are the same or all three are collinear.
         * @return true if triangles are fine and none of them is degenerate
         */
        [[nodiscard]] bool checkTrianglesNotDegenerated() const;

        /**
         * Fixes the orientation of the plane unit normals for a given set of violating face indices.
         *
         * @param actualOrientation The desired plane unit normal orientation.
         * @param violatingIndices A set of indices representing the faces with violating orientations.
         */
        void healPlaneUnitNormalOrientation(const NormalOrientation &actualOrientation, const std::set<size_t> &violatingIndices);

        /**
         * Calculates how often a vector starting at a specific origin intersects a polyhedron's mesh's triangles.
         * @param face the vector describing the ray
         * @return true if the ray intersects the triangle
         */
        [[nodiscard]] size_t countRayPolyhedronIntersections(const Array3Triplet &face) const;

        /**
         * Calculates how often a vector starting at a specific origin intersects a triangular face.
         * Uses the Möller–Trumbore intersection algorithm.
         * @param rayOrigin the origin of the ray
         * @param rayVector the vector describing the ray
         * @param triangle a triangular face
         * @return intersection point or null
         *
         * @see Adapted from https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
         */
        static std::unique_ptr<Array3> rayIntersectsTriangle(const Array3 &rayOrigin, const Array3 &rayVector, const Array3Triplet &triangle);


    };

}
