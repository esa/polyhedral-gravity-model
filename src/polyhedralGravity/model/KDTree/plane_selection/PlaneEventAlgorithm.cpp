#include "polyhedralGravity/model/KDTree/plane_selection/PlaneEventAlgorithm.h"

namespace polyhedralGravity {
    TriangleCounter::TriangleCounter(const size_t dimensionCount, const std::array<size_t, 3> &initialValues)
        : dimensionTriangleValues(dimensionCount, initialValues) {
        if (dimensionCount == 0) {
            throw std::invalid_argument("Dimension count must be greater than zero");
        }
    }

    void TriangleCounter::updateMax(Direction direction, const size_t p_planar, const size_t p_end) {
        dimensionTriangleValues.at(static_cast<size_t>(direction) % dimensionTriangleValues.size()).at(1) -= p_planar +
                p_end;
    }

    void TriangleCounter::updateMin(Direction direction, const size_t p_planar, const size_t p_start) {
        dimensionTriangleValues.at(static_cast<size_t>(direction) % dimensionTriangleValues.size()).at(0) += p_planar +
                p_start;
    }

    void TriangleCounter::setPlanar(Direction direction, const size_t p_planar) {
        dimensionTriangleValues.at(static_cast<size_t>(direction) % dimensionTriangleValues.size()).at(2) = p_planar;
    }

    size_t TriangleCounter::getMin(Direction direction) const {
        return dimensionTriangleValues.at(static_cast<size_t>(direction) % dimensionTriangleValues.size()).at(0);
    }

    size_t TriangleCounter::getMax(Direction direction) const {
        return dimensionTriangleValues.at(static_cast<size_t>(direction) % dimensionTriangleValues.size()).at(1);
    }

    size_t TriangleCounter::getPlanar(Direction direction) const {
        return dimensionTriangleValues.at(static_cast<size_t>(direction) % dimensionTriangleValues.size()).at(2);
    }

    PlaneEventVector PlaneEventAlgorithm::generatePlaneEventsFromFaces(const SplitParam &splitParam,
                                                                       std::vector<Direction> directions) {
        // each face has min and max point and each proposes a plane in each of the directions
        PlaneEventVector events{};
        events.reserve(countFaces(splitParam.boundFaces) * 2 * directions.size());
        //mutex used for synchronizing insertions through threads
        std::mutex eventsMutex{};
        if (std::holds_alternative<PlaneEventVector>(splitParam.boundFaces)) {
            return std::get<PlaneEventVector>(splitParam.boundFaces);
        }
        const auto &boundTriangles{std::get<TriangleIndexVector>(splitParam.boundFaces)};
        //transform the faces into vertices
        auto [vertex3_begin, vertex3_end] = transformIterator(boundTriangles.cbegin(), boundTriangles.cend(),
                                                              splitParam.vertices, splitParam.faces);
        thrust::for_each(thrust::device, vertex3_begin, vertex3_end,
                         [&splitParam, &events, &directions, &eventsMutex](const auto &indexAndTriplet) {
                             const auto [index, triplet] = indexAndTriplet;
                             //first clip the triangles vertices to the current bounding box and then get the bounding box of the clipped triangle -> use the box edges as split plane candidates
                             const auto [minPoint, maxPoint] = Box::getBoundingBox<std::vector<Array3> >(
                                 splitParam.boundingBox.clipToVoxel(triplet));
                             std::lock_guard lock(eventsMutex);
                             for (const auto &direction: directions) {
                                 // if the triangle is perpendicular to the split direction, generate a planar event with the candidate plane in which the triangle lies
                                 if (minPoint[static_cast<int>(direction)] == maxPoint[static_cast<int>(direction)]) {
                                     events.emplace_back(
                                         PlaneEventType::planar,
                                         Plane(minPoint, direction),
                                         index);
                                     return;
                                 }
                                 //else create a starting and ending event consisting of the planes defined by the min and max points of the face's bounding box.
                                 events.emplace_back(
                                     PlaneEventType::starting,
                                     Plane(minPoint, direction),
                                     index);
                                 events.emplace_back(
                                     PlaneEventType::ending,
                                     Plane(maxPoint, direction),
                                     index);
                             }
                         });
        //reduce size
        events.shrink_to_fit();
        //sort the events by plane position and then by PlaneEventType. Refer to {@link PlaneEventType} for the specific order
        std::sort(events.begin(), events.end());
        return events;
    }

    std::tuple<Plane, double, bool> PlaneEventAlgorithm::traversePlaneEvents(
        const PlaneEventVector &events, TriangleCounter &triangleCounter, const Box &boundingBox) {
        //initialize the default plane and make it costly
        double cost{std::numeric_limits<double>::infinity()};
        Plane optPlane{};
        bool minSide{true};
        //traverse all the events
        int i{0};
        while (i < events.size()) {
            //poll a plane to test
            const Plane &candidatePlane = events[i].plane;
            //for each plane calculate the faces whose vertices lie in the plane. Differentiate between the face starting in the plane, ending in the plane or all vertices lying in the plane
            size_t p_start{0}, p_end{0}, p_planar{0};
            //count all faces that end in the plane, this works because the PlaneEvents are sorted by position and then by PlaneEventType
            while (i < events.size() && events[i].plane.orientation == candidatePlane.orientation && events[i].plane.
                   axisCoordinate == candidatePlane.axisCoordinate && events[i].type == PlaneEventType::ending) {
                p_end++;
                i++;
            }
            //count all the faces that lie in the plane
            while (i < events.size() && events[i].plane.orientation == candidatePlane.orientation && events[i].plane.
                   axisCoordinate == candidatePlane.axisCoordinate && events[i].type == PlaneEventType::planar) {
                p_planar++;
                i++;
            }
            //count all the faces that start in the plane
            while (i < events.size() && events[i].plane.orientation == candidatePlane.orientation && events[i].plane.
                   axisCoordinate == candidatePlane.axisCoordinate && events[i].type == PlaneEventType::starting) {
                p_start++;
                i++;
            }
            //reference to the triangle counter of the current dimension -> better readability
            triangleCounter.setPlanar(candidatePlane.orientation, p_planar);
            triangleCounter.updateMax(candidatePlane.orientation, p_planar, p_end);
            //evaluate plane and update should the new plane be more efficient
            auto [candidateCost, minSideChosen] = costForPlane(boundingBox, candidatePlane,
                                                               triangleCounter.getMin(candidatePlane.orientation),
                                                               triangleCounter.getMax(candidatePlane.orientation),
                                                               triangleCounter.getPlanar(candidatePlane.orientation));
            // this condition exists to consistently build the same KDTree (choose plane with lower coordinate) by eliminating indeterministic behavior should the cost be equal.
            // this is not important for functionality but for testing purposes
            bool skipEvaluation = candidateCost == cost && optPlane.axisCoordinate < candidatePlane.axisCoordinate;
            if (candidateCost < cost && !skipEvaluation) {
                cost = candidateCost;
                optPlane = candidatePlane;
                minSide = minSideChosen;
            }
            //shift the plane to the next candidate and prepare next iteration
            triangleCounter.updateMin(candidatePlane.orientation, p_planar, p_start);
            triangleCounter.setPlanar(candidatePlane.orientation, 0);
        }
        //return the optimal plane and the plane's cost along with the halfspace to which planar faces are included.
        return {optPlane, cost, minSide};
    }
} // namespace polyhedralGravity
