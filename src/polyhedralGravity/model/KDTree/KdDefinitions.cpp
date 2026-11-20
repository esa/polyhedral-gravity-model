#include "polyhedralGravity/model/KDTree/KdDefinitions.h"

namespace polyhedralGravity {
    Array3 normal(const Direction direction) {
        switch (direction) {
            case Direction::X:
                return Array3{1, 0, 0};
            case Direction::Y:
                return Array3{0, 1, 0};
            case Direction::Z:
                return Array3{0, 0, 1};
            default:
                throw std::invalid_argument{"Unknown Direction enum value used during normal fetching."};
        }
    }

    Plane::Plane(const Array3 &point, Direction direction)
        : axisCoordinate(point[static_cast<int>(direction)]), orientation(direction) {
    }

    Plane::Plane(const double point, const Direction direction)
        : axisCoordinate(point), orientation(direction) {
    }

    std::ostream &operator<<(std::ostream &os, const Direction &direction) {
        switch (direction) {
            case Direction::X: os << "X";
                break;
            case Direction::Y: os << "Y";
                break;
            case Direction::Z: os << "Z";
                break;
        }
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Plane &plane) {
        os << "(" << plane.axisCoordinate << ", " << plane.orientation << ")";
        return os;
    }

    Array3 Plane::normal(const bool returnFlipped) const {
        using namespace util;
        return polyhedralGravity::normal(orientation) * (returnFlipped ? -1 : 1);
    }

    Array3 Plane::originPoint() const {
        Array3 point{0.0, 0.0, 0.0};
        point[static_cast<int>(orientation)] = axisCoordinate;
        return point;
    }

    double Plane::rayPlaneIntersection(const Array3 &origin, const Array3 &inverseRay) const {
        const auto &origin_coord = origin[static_cast<int>(orientation)];
        const auto &inverseRay_coord = inverseRay[static_cast<int>(orientation)];

        const auto t = (this->axisCoordinate - origin_coord) * inverseRay_coord;
        // NaN possible through 0./+-inf -> point lies on plane and ray is parallel to plane, t=0 should be returned
        // inverseRay_coord ray is +inf (happens during inverse calculation by 1./0.) if ray is parallel to plane -> If origin additionally not on the plane then algorithm returns +-inf . The sign gives information in which half-space defined by the plane the origin lies.
        return std::isnan(t) ? 0.0 : t;
    }


    bool Plane::operator==(const Plane &other) const {
        return std::fabs(axisCoordinate - other.axisCoordinate) < 1e-15 && orientation == other.orientation;
    }

    bool Plane::operator!=(const Plane &other) const {
        return !(*this == other);
    }

    Box::Box(const std::pair<Array3, Array3> &pair)
        : minPoint{pair.first}, maxPoint{pair.second} {
    }

    Box::Box()
        : minPoint{0.0, 0.0, 0.0}, maxPoint{0.0, 0.0, 0.0} {
    }

    std::pair<double, double> Box::rayBoxIntersection(const Array3 &origin, const Array3 &inverseRay) const {
        //calculate the parameter t in $ origin + t * ray = point $
        const auto lambdaIntersectSlabPoint = [&origin, &inverseRay](const Array3 &point) {
            using namespace util;
            std::array<double, 3> result{};
            for (const auto &direction: ALL_DIRECTIONS) {
                result[static_cast<size_t>(direction)] = Plane(point, direction).rayPlaneIntersection(
                    origin, inverseRay);
            }
            return result;
        };
        //intersections with slabs defined through minPoint
        auto [tx_min, ty_min, tz_min] = lambdaIntersectSlabPoint(minPoint);
        //intersections with slabs defined through maxPoint
        auto [tx_max, ty_max, tz_max] = lambdaIntersectSlabPoint(maxPoint);

        const auto assignEnterExitValue = [](const auto t_min, const auto t_max,
                                             const auto rayDir) -> std::pair<double, double> {
            // if ray is shot in positive ray direction then minPoint slab is hit before max point slab -> min point slab holds the entry point and max point slab the exit point
            // otherwise max point slab is hit first and then min point slab after
            return rayDir < 0 ? std::make_pair(t_max, t_min) : std::make_pair(t_min, t_max);
        };
        //return the parameters in ordered by '<'
        auto [tx_enter, tx_exit] = assignEnterExitValue(tx_min, tx_max, inverseRay[0]);
        auto [ty_enter, ty_exit] = assignEnterExitValue(ty_min, ty_max, inverseRay[1]);
        auto [tz_enter, tz_exit] = assignEnterExitValue(tz_min, tz_max, inverseRay[2]);

        //calculate the point where all slabs have been entered: t_enter
        const double t_enter{std::max(tx_enter, std::max(ty_enter, tz_enter))};
        //calculates the point where the first slab has been exited: t_exit
        const double t_exit{std::min(tx_exit, std::min(ty_exit, tz_exit))};

        return {t_enter, t_exit};
    }

    double Box::surfaceArea() const {
        const double width = std::abs(maxPoint[0] - minPoint[0]);
        const double length = std::abs(maxPoint[1] - minPoint[1]);
        const double height = std::abs(maxPoint[2] - minPoint[2]);
        return 2 * (width * length + width * height + length * height);
    }

    std::pair<Box, Box> Box::splitBox(const Plane &plane) const {
        //clone the original box two times -> modify clones to become child boxes defined by the splitting plane
        Box box1{*this};
        Box box2{*this};
        const Direction &axis{plane.orientation};
        //Shift edges of the boxes to match the plane
        box1.maxPoint[static_cast<int>(axis)] = plane.axisCoordinate;
        box2.minPoint[static_cast<int>(axis)] = plane.axisCoordinate;
        return std::make_pair(box1, box2);
    }

    std::vector<Array3> Box::clipToVoxel(const std::array<Array3, 3> &points) const {
        using namespace util;
        //use clipped as the input vector because the inner for loop swaps input and clipped each iteration,
        //since each iteration needs the output of the previous iteration as input.
        std::vector clipped(points.cbegin(), points.cend());
        std::vector<Array3> input{};
        input.reserve(points.size());
        //every plane defined by the maxPoint has to flip its normal because the normals have to point inside the bounding box.
        bool flipPlane = false;
        for (const Direction direction: ALL_DIRECTIONS) {
            const auto directionPlanes = {Plane(minPoint, direction), Plane(maxPoint, direction)};
            for (const auto &plane: directionPlanes) {
                std::swap(input, clipped);
                clipToVoxelPlane(plane, flipPlane, input, clipped);
                input.clear();
                flipPlane = !flipPlane;
            }
        }
        return clipped;
    }

    void Box::clipToVoxelPlane(const Plane &plane, const bool flipPlaneNormal, const std::vector<Array3> &source,
                               std::vector<Array3> &dest) {
        using namespace util;
        //the distance is interpreted in the normal direction, negative values are in opposite direction of the normal.
        auto distanceMeasures = [&plane, &flipPlaneNormal](
            const Array3 &point) -> double {
            // works for arbitrary plane orientations:
            // return dot(point - plane.originPoint(), plane.normal(flipPlaneNormal));
            // only works for axis aligned planes:
            return (point[static_cast<int>(plane.orientation)] - plane.axisCoordinate) * (flipPlaneNormal ? -1. : 1.);
        };
        static constexpr auto isInside = [](const double distance) { return distance >= 0.0; };
        static constexpr auto intersectionPoint = [](const Array3 &from, const Array3 &to, const double distanceFrom,
                                                     const double distanceTo) {
            // solve for t in $ [(t * from + (1-t) * to ) - origin] * normal = 0 $
            // equation explained: search for a point on the plane defined by $ (point - origin) * normal $, where point is linearly interpolated using vectors from and to.
            const double t{distanceTo / (distanceTo - distanceFrom)};
            return from * t + to * (1.0 - t);
        };
        for (size_t i{0}; i < source.size(); i++) {
            const Array3 &from{source[i]};
            const Array3 &to{source[(i + 1) % source.size()]};
            // $ (from - origin) * normal = cos alpha * |from - origin| * |normal| = cos alpha * |from - origin| * 1 $ ^= distance of from to the plane in the direction of the normal.
            const auto distanceFrom = distanceMeasures(from);
            const auto distanceTo = distanceMeasures(to);
            if (isInside(distanceFrom) && isInside(distanceTo)) {
                dest.push_back(to);
            } else if (isInside(distanceFrom) && !isInside(distanceTo)) {
                dest.emplace_back(intersectionPoint(from, to, distanceFrom, distanceTo));
            } else if (!isInside(distanceFrom) && isInside(distanceTo)) {
                dest.emplace_back(intersectionPoint(from, to, distanceFrom, distanceTo));
                dest.push_back(to);
            } else if (!isInside(distanceFrom) && !isInside(distanceTo)) {
                //do nothing
            }
        }
    }


    PlaneEvent::PlaneEvent(const PlaneEventType type, const Plane plane, const unsigned faceIndex)
        : type{type}, plane{plane}, faceIndex{faceIndex} {
    }

    bool PlaneEvent::operator<(const PlaneEvent &other) const {
        if (this->plane.axisCoordinate != other.plane.axisCoordinate) {
            return this->plane.axisCoordinate < other.plane.axisCoordinate;
        }
        if (this->plane.orientation == other.plane.orientation) {
            return this->type < other.type;
        }
        return static_cast<unsigned>(this->plane.orientation) < static_cast<unsigned>(other.plane.orientation);
    }

    bool PlaneEvent::operator==(const PlaneEvent &other) const {
        return type == other.type && plane == other.plane && faceIndex == other.faceIndex;
    }

    TriangleIndexVector convertEventsToFaces(const std::variant<TriangleIndexVector, PlaneEventVector> &events) {
        if (std::holds_alternative<TriangleIndexVector>(events)) {
            return std::get<TriangleIndexVector>(events);
        }
        const auto &eventList{std::get<PlaneEventVector>(events)};
        TriangleIndexVector triangles{};
        triangles.reserve(eventList.size());
        //used to avoid duplication
        std::unordered_set<size_t> processedFaces{};
        auto insertIfAbsent = [&triangles, &processedFaces](const auto &planeEvent) {
            const auto faceIndex{planeEvent.faceIndex};
            if (processedFaces.find(faceIndex) == processedFaces.end()) {
                processedFaces.insert(faceIndex);
                triangles.push_back(faceIndex);
            }
        };
        std::for_each(eventList.cbegin(), eventList.cend(), insertIfAbsent);
        triangles.shrink_to_fit();
        return triangles;
    }

    size_t countFaces(const std::variant<TriangleIndexVector, PlaneEventVector> &triangles) {
        return std::visit(util::overloaded{
                              [](const TriangleIndexVector &indexList) {
                                  return indexList.size();
                              },
                              [](const PlaneEventVector &eventList) {
                                  std::mutex writeLock{};
                                  size_t count{0};
                                  std::unordered_set<size_t> processedFaces{};
                                  std::for_each(eventList.cbegin(), eventList.cend(),
                                                [&processedFaces, &count, &writeLock](const auto &planeEvent) {
                                                    if (processedFaces.find(planeEvent.faceIndex) == processedFaces.
                                                        end()) {
                                                        processedFaces.insert(planeEvent.faceIndex);
                                                        count++;
                                                    }
                                                });
                                  return count;
                              }
                          },
                          triangles);
    }

    size_t recursionDepth(const size_t nodeId) {
        //t: depth of the current node, N: amount of nodes in the tree if tree is complete
        // $ N = 2^(t+1)-1 $ (geometric series formula for partial sums), assume $ N >= idx + 1 $
        return static_cast<size_t>(std::ceil(std::log2(nodeId + 2))) - 1;
    }
} // namespace polyhedralGravity
