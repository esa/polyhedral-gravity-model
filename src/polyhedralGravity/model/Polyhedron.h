#pragma once

#include <utility>
#include <vector>
#include <array>
#include <set>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <tuple>
#include <variant>
#include <string>
#include <sstream>
#include <memory>
#include "thrust/copy.h"
#include "thrust/device_vector.h"
#include "thrust/transform_reduce.h"
#include "thrust/execution_policy.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/Definitions.h"
#include "polyhedralGravity/input/MeshReader.h"
#include "polyhedralGravity/output/Logging.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityFloatArithmetic.h"

namespace polyhedralGravity {

    /* Forward declaration of Polyhedron */
    class Polyhedron;

    /**
     * The orientation of the plane unit normals of the polyhedron.
     * We use this property as the precise definition of the vertices ordering depends on the
     * utilized coordinate system.
     * However, the normal alignment is independent. Tsoulis et al. equations require the
     * normals to point outwards of the polyhedron. If the opposite hold, the result is
     * negated.
     */
    enum class NormalOrientation: char {
        /** Outwards pointing plane unit normals */
        OUTWARDS,
        /** Inwards pointing plane unit normals */
        INWARDS
    };


    /**
     * Overloaded insertion operator to output the string representation of a NormalOrientation enum value.
     *
     * @param os The output stream to write the string representation to.
     * @param orientation The NormalOrientation enum value to output.
     * @return The output stream after writing the string representation.
     */
    inline std::ostream &operator<<(std::ostream &os, const NormalOrientation &orientation) {
        switch (orientation) {
            case NormalOrientation::OUTWARDS:
                os << "OUTWARDS";
                break;
            case NormalOrientation::INWARDS:
                os << "INWARDS";
                break;
            default:
                os << "Unknown";
                break;
        }
        return os;
    }

    /**
     * The three mode the polyhedron class takes in
     * the constructor in order to determine what initialization checks to conduct.
     * This enum is exclusively utilized in the constructor of a {@link Polyhedron} and its private method
     * {@link runIntegrityMeasures}
     */
    enum class PolyhedronIntegrity: char {
        /**
         * All activities regarding MeshChecking are disabled.
         * No runtime overhead!
         */
        DISABLE,
        /**
         * Only verification of the NormalOrientation.
         * A misalignment (e.g., specified OUTWARDS, but is not) leads to a runtime_error.
         * Runtime Cost @f$O(n^2)@f$
         */
        VERIFY,
        /**
         * Like VERIFY, but also informs the user about the option in any case on the runtime costs.
         * This is the implicit default option.
         * Runtime Cost: @f$O(n^2)@f$ and output to stdout in every case!
         */
        AUTOMATIC,
        /**
         * Verification and Automatic Healing of the NormalOrientation.
         * A misalignment does not lead to a runtime_error, but to an internal correction.
         * Runtime Cost: @f$O(n^2)@f$ and a modification of the mesh input!
         */
        HEAL,
    };

    /**
     * Represents the unit of a polyhedron's mesh.
     */
    enum class MetricUnit : char {
        /** The unit meter @f$[m]@f$ */
        METER,
        /** The unit kilometer @f$[km]@f$ */
        KILOMETER,
        /** The mesh is unitless @f$[1]@f$ */
        UNITLESS,
    };


    /**
     * Stream operator for the MetricUnit enum. Prints the enum to a human-readable string.
     * @param os The output stream to write the string representation to.
     * @param metricUnit the metric unit to print
     * @return The output stream after writing the string representation.
     */
    inline std::ostream &operator<<(std::ostream &os, const MetricUnit &metricUnit) {
        switch (metricUnit) {
            case MetricUnit::METER:
                os << "m";
            break;
            case MetricUnit::KILOMETER:
                os << "km";
            break;
            default:
                os << "unitless";
            break;
        }
        return os;
    }

    /**
     * Data structure containing the model data of one polyhedron. This includes nodes, edges (faces) and elements.
     * The index always starts with zero!
    */
    class Polyhedron {

        /**
         * A vector containing the vertices of the polyhedron.
         * Each node is an array of size three containing the xyz coordinates.
         * The mesh must be scaled in the same units as the density is given
         * (the unit must match to the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         */
        const std::vector<Array3> _vertices;

        /**
         * A vector containing the faces (triangles) of the polyhedron.
         * Each face is an array of size three containing the indices of the nodes forming the face.
         * Since every face consists of three nodes, every face consists of three segments. Each segment consists of
         * two nodes.
         * For example, a face consisting of {1, 2, 3} --> segments: {1, 2}, {2, 3}, {3, 1}
         */
        std::vector<IndexArray3> _faces;

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
         * Default destructor
         */
        ~Polyhedron() = default;

        /**
         * Returns the vertices of this polyhedron
         * @return vector of cartesian coordinates
         */
        [[nodiscard]] const std::vector<Array3> &getVertices() const;

        /**
         * Returns the vertex at a specific index
         * @param index size_t
         * @return cartesian coordinates of the vertex at index
         */
        [[nodiscard]] const Array3 &getVertex(size_t index) const;

        /**
         * The number of points (nodes) that make up the polyhedron.
         * @return a size_t
         */
        [[nodiscard]] size_t countVertices() const;

        /**
         * Returns the triangular faces of this polyhedron
         * @return vector of triangular faces, where each element size_t references a vertex in the vertices vector
         */
        [[nodiscard]] const std::vector<IndexArray3> &getFaces() const;

        /**
         * Returns the indices of the vertices making up the face at the given index.
         * @param index size_t
         * @return triplet of the vertices indices forming the face
         */
        [[nodiscard]] const IndexArray3 &getFace(size_t index) const;

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
         * Returns the metric unit of the polyhedron's mesh, which is either METER, KILOMETER, or UNIT_LESS.
         * @return the metric unit of the polyhedron
         */
        [[nodiscard]] MetricUnit getMetricUnit() const;

        /**
         * Returns the metric unit of the polyhedral mesh.
         * It is either in meter, kilometer, or without units.
         * @return unit of mesh (i.e., the unit of the vertices)
         */
        [[nodiscard]] std::string getMeshUnit() const;

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
         * An iterator transforming the polyhedron's coordinates on demand by a given offset.
         * This function returns a pair of transform iterators (first = begin(), second = end()).
         * @param offset the offset to apply
         * @return pair of transform iterators
         */
        [[nodiscard]] inline auto transformIterator(const Array3 &offset = {0.0, 0.0, 0.0}) const {
            //The offset must be captured by value to ensure its lifetime!
            const auto lambdaOffsetApplication = [this, offset](const IndexArray3 &face) -> Array3Triplet {
                using namespace util;
                return {this->getVertex(face[0]) - offset, this->getVertex(face[1]) - offset,
                        this->getVertex(face[2]) - offset};
            };
            auto first = thrust::make_transform_iterator(this->getFaces().begin(), lambdaOffsetApplication);
            auto last = thrust::make_transform_iterator(this->getFaces().end(), lambdaOffsetApplication);
            return std::make_pair(first, last);
        }

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
