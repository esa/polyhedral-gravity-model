#pragma once

#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "thrust/detail/execution_policy.h"
#include "thrust/execution_policy.h"
#include "thrust/system/detail/sequential/for_each.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/iterator/transform_iterator.h>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace polyhedralGravity {

    /**
     * Assigns an integer index to the coordinate axes
     *
     * Used to specify coordinates. E.g. CoordinateArray[Direction]
     */
    enum class Direction {
        X = 0,
        Y = 1,
        Z = 2
    };

    static const std::vector ALL_DIRECTIONS{Direction::X, Direction::Y, Direction::Z};


    /**
     * Returns the normal vector for a direction.
     * @param direction The direction to return the normal vector for.
     * @return The normal vector.
     */
    static Array3 normal(Direction direction);

    /**
     * Number of dimensions for the polyhedron. Also corresponds to the number of elements of the {@link Direction} enum.
     */
    constexpr int DIMENSIONS = 3;

    /**
     * Defines a plane that is parallel to one of the coordinate planes, by taking the fixed axis coordinate value for the plane and the coordinate index ({@link Direction}) that is fixed for every \
     * point on the plane.
     *
     * E.g. Specifying 0.0 and Direction::X would describe the YZ plane that goes through the origin. The direction is equivalent to the coordinate that is 1 in the normal vector of the plane, the others are 0.
     */
    struct Plane {
        /**
         * Each point lying on the plane has to have this value in the dimension specified in the orientation parameter.
         */
        double axisCoordinate;
        /**
         * Specifies which coordinate dimension is fixed for every point on the plane.
         */
        Direction orientation;

        /**
        * Returns the normal vector for this plane.
        * @param returnFlipped Whether to return the normal pointing in the opposite direction.
        * @return The normal vector.
        */
        [[nodiscard]] Array3 normal(bool returnFlipped = false) const;

        /**
         * Returns the origin point of a plane, meaning a point that lies on the plane.
         * @return The origin point.
         */
        [[nodiscard]] Array3 originPoint() const;

        /**
        * Intersects a ray with the splitPlane.
        * @param origin The point where the ray originates from.
        * @param inverseRay The inverse ray direction vector of the ray to be intersected (used for faster calculations).
        * @return Returns the t parameter for the intersection point, with t being from the equation $intersection_point = orig + t * ray$. t is +-infinity if no intersection point is present.
        */
        [[nodiscard]] double rayPlaneIntersection(const Array3 &origin, const Array3 &inverseRay) const;

        /**
        * Equality operator used for testing purposes
        */
        bool operator==(const Plane &other) const;

        /**
        * Inequality operator used for testing purposes
        */
        bool operator!=(const Plane &other) const;

        friend std::ostream& operator<<(std::ostream &os, const Plane &plane);

        Plane() = default;
        Plane(const Array3 &point, Direction direction);
        Plane(double point, Direction direction);
    };

    /**
     * Defines a rectangular box by taking two opposite corner points. First is the point closest to the origin and second is the point farthest away.
     */
    struct Box {
        /**
         * The point closer to the origin, thus minimal
         */
        Array3 minPoint;
        /**
         * The point further away from the origin, thus maximal.
         */
        Array3 maxPoint;

        /**
         * Calculates the intersection points of a ray and a box.
         * @param origin The origin of the ray.
         * @param inverseRay The inverse ray direction vector of the ray to be intersected (used for faster calculations).
         * @return Parameters t of the equation $ intersection_point = origin + t * ray $ for the entry and exit intersection points.
         */
        [[nodiscard]] std::pair<double, double> rayBoxIntersection(const Array3 &origin, const Array3 &inverseRay) const;

        /**
        * Calculates the surface area of a box.
        * @return the surface area
        */
        [[nodiscard]] double surfaceArea() const;

        /**
       * Splits this box into two new boxes.
       * @param plane the plane by which to split the original box.
       * @return a pair of boxes that result by splitting this box.
       */
        [[nodiscard]] std::pair<Box, Box> splitBox(const Plane &plane) const;

        /**
        * Finds the minimal bounding box for a set of vertices.
        * @param vertices the set of vertex coordinates for which to find the box
        * @return the bounding box {@link Box}
        */
        template<typename Container>
        static Box getBoundingBox(const Container &vertices) {
            using namespace util;
            return Box(findMinMaxCoordinates<Container, Array3>(vertices));
        }

        /**
        * Takes points of a face of a polyhedron and clips them to this box. If all the points lie in the box no changes are made but if points lie outside of the box they are linearly interpolated onto the box.
        * Uses the Sutherland-Hodgman-Algorithm.
        * @param points The corner points of the face to be clipped.
        * @return The new corner points of the clipped face.
        */
        [[nodiscard]] std::vector<Array3> clipToVoxel(const std::array<Array3, 3> &points) const;

        explicit Box(const std::pair<Array3, Array3> &pair);
        Box();

    private:
        /**
         * Takes a plane and a set of vertices and clips them accordingly.
         * Used as a sub procedure by the Sutherland-Hodgman-Algorithm.
         * @param plane The plane to split the vertices by
         * @param flipPlaneNormal Specifies which side of the plane is inside (In the direction or opposite of the plane normal).
         * @param source The vertices to be transformed to lie on the inside of the plane.
         * @param dest The transformed vertices.
         */
        static void clipToVoxelPlane(const Plane &plane, bool flipPlaneNormal, const std::vector<Array3> &source, std::vector<Array3> &dest);
    };

    /**
     * A set that stores indices of the faces vector in the KDTree. This effectively corresponds to a set of triangles. For performance purposes a std::vector is used instead of a std::set.
     */
    using TriangleIndexVector = std::vector<size_t>;

    /**
    * Triangle sets contained in an array. Used by the KDTree to divide a bounding boxes included triangles into smaller subsets. For the semantic purpose of the contained sets please refer to the comments in the usage context.
     */
    template<size_t Number>
    using TriangleIndexVectors = std::array<std::unique_ptr<TriangleIndexVector>, Number>;

    /**
    * Used by {@link PlaneEvent} to position the face that generated the event relative to the generated plane.
    */
    enum class PlaneEventType {
        ending = 0,
        planar = 1,
        starting = 2,
    };

    /**
     * Generated when traversing the vector of faces and building their candidate planes.
     */
    struct PlaneEvent {
        PlaneEventType type;
        /**
         * The candidate plane suggested by the face included in this struct.
         */
        Plane plane;
        /**
         * The index of the face that generated this candidate plane.
         */
        unsigned int faceIndex;

        PlaneEvent(PlaneEventType type, Plane plane, unsigned faceIndex);

        PlaneEvent() = default;

        /**
         * Less operator used for sorting an PlaneEvent vector.
         * @param other the PlaneEvent to compare this to.
         * @return true if this should precede the other argument.
         */
        bool operator<(const PlaneEvent &other) const;

        /**
         *Equality operator used for testing purposes
         */
        bool operator==(const PlaneEvent &other) const;
    };

    /**
     * A list of PlaneEvents.
    */
    using PlaneEventVector = std::vector<PlaneEvent>;

    /**
    * An array of PlaneEventLists.
    */
    template<size_t Number>
    using PlaneEventVectors = std::array<std::unique_ptr<PlaneEventVector>, Number>;

    /**
     * Extracts the indices of the triangle faces that are referenced by the given PlaneEvents.
     * @param events The PlaneEvents containing information about the faces.
     * @return A list of face indices.
     */
    TriangleIndexVector convertEventsToFaces(const std::variant<TriangleIndexVector, PlaneEventVector> &events);

    size_t countFaces(const std::variant<TriangleIndexVector, PlaneEventVector> &triangles);

    /**
        * An iterator transforming face indices to vertices and returning both.
        * This function returns a pair of transform iterators (first = begin(), second = end()).
        * @param begin begin iterator of the face indice vector to transform.
        * @param end end iterator of the face indice vector to transform.
        * @param vertices the vector of vertices to look up the indices obtained from the faces vector.
        * @param faces the faces vector to lookup face indices.
        * @return pair of transform iterators.
        */
    [[nodiscard]] static auto transformIterator(const std::vector<size_t>::const_iterator begin, const std::vector<size_t>::const_iterator end, const std::vector<Array3> &vertices, const std::vector<IndexArray3> &faces) {
        //The offset must be captured by value to ensure its lifetime!
        const auto lambdaApplication = [&vertices, &faces](size_t faceIndex) {
            const auto &face = faces[faceIndex];
            Array3Triplet vertexTriplet = {
                    vertices[face[0]],
                    vertices[face[1]],
                    vertices[face[2]]};
            return std::make_pair(faceIndex, vertexTriplet);
        };

        auto first = thrust::make_transform_iterator(begin, lambdaApplication);
        auto last = thrust::make_transform_iterator(end, lambdaApplication);
        return std::make_pair(first, last);
    }

    /**
    * The distance to the tree's root from this node. Used to limit the depth and the size of the tree.
    * @param nodeId The id of the node to determine the depth (distance to root node) of.
    */
    size_t recursionDepth(size_t nodeId);

}// namespace polyhedralGravity