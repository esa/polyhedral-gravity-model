#include "polyhedralGravity/model/KDTree/plane_selection/LogNSquaredPlane.h"

namespace polyhedralGravity {
    // O(N*log^2(N)) implementation
    std::tuple<Plane, double, std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2> > >
    LogNSquaredPlane::findPlane(const SplitParam &splitParam) {
        Plane optPlane{};
        double cost{std::numeric_limits<double>::infinity()};
        PlaneEventVector optimalEvents{};
        bool minSide{true};
        for (const auto dimension: ALL_DIRECTIONS) {
            splitParam.splitDirection = dimension;
            auto [candidatePlane, candidateCost, events, minSideChosen] = findPlaneForSingleDimension(splitParam);
            // this if clause exists to consistently build the same KDTree (choose plane with lower coordinate) by eliminating indeterministic behavior should the cost be equal.
            // this is not important for functionality but for testing purposes
            if (candidateCost == cost && optPlane.axisCoordinate < candidatePlane.axisCoordinate) {
                continue;
            }
            if (candidateCost < cost) {
                optPlane = candidatePlane;
                cost = candidateCost;
                optimalEvents = events;
                minSide = minSideChosen;
            }
        }
        //generate the triangle index lists for the child bounding boxes and return them along with the optimal plane and the plane's cost.
        return {optPlane, cost, generateTriangleSubsets(optimalEvents, optPlane, minSide)};
    }

    std::tuple<Plane, double, PlaneEventVector, bool> LogNSquaredPlane::findPlaneForSingleDimension(
        const SplitParam &splitParam) {
        const PlaneEventVector events{std::move(generatePlaneEventsFromFaces(splitParam, {splitParam.splitDirection}))};
        TriangleCounter triangleCounter{1, {0, countFaces(splitParam.boundFaces), 0}};
        auto [optPlane, cost, minSide] = traversePlaneEvents(events, triangleCounter, splitParam.boundingBox);
        return {optPlane, cost, events, minSide};
    }


    TriangleIndexVectors<2> LogNSquaredPlane::generateTriangleSubsets(const PlaneEventVector &planeEvents,
                                                                      const Plane &plane, const bool minSide) {
        auto facesMin = std::make_unique<TriangleIndexVector>();
        auto facesMax = std::make_unique<TriangleIndexVector>();
        //set data structure to avoid processing faces twice -> introduces O(1) lookup instead of O(n) lookup using the vectors directly
        std::unordered_set<size_t> facesMinLookup{};
        std::unordered_set<size_t> facesMaxLookup{};
        //each face will most of the time generate two events, the split plane will try to distribute the faces evenly
        //Thus reserving 0.5 * 0.5 * planeEvents.size() for each vector
        facesMin->reserve(planeEvents.size() / 4);
        facesMax->reserve(planeEvents.size() / 4);
        facesMinLookup.reserve(planeEvents.size() / 4);
        facesMaxLookup.reserve(planeEvents.size() / 4);
        std::array<std::mutex, 2> facesMutex{};
        thrust::for_each(thrust::device, planeEvents.cbegin(), planeEvents.cend(),
                         [&facesMin, &facesMax, &plane, minSide, &facesMinLookup, &facesMaxLookup, &facesMutex](
                     const auto &event) {
                             //lambda function to combine lookup and insertion into one place
                             auto insertIfAbsent = [&facesMin, &facesMinLookup, &facesMax, &facesMaxLookup, &facesMutex
                                     ](const size_t faceIndex, const uint8_t index) {
                                 const auto &vector = index == 0 ? facesMin : facesMax;
                                 auto &lookup = index == 1 ? facesMinLookup : facesMaxLookup;
                                 std::lock_guard lock(facesMutex[index]);
                                 if (lookup.find(faceIndex) == lookup.end()) {
                                     lookup.insert(faceIndex);
                                     vector->push_back(faceIndex);
                                 }
                                 // Since each face can only be referenced by max two events (since there are only two planes encasing it),
                                 // after the face has been already been processed once, it can be removed from the lookup buffer after the second time to save space
                                 else {
                                     lookup.erase(lookup.find(faceIndex));
                                 }
                             };
                             //sort the triangles by inferring their position from the event's candidate split plane
                             if (event.plane.axisCoordinate != plane.axisCoordinate) {
                                 insertIfAbsent(event.faceIndex,
                                                event.plane.axisCoordinate < plane.axisCoordinate ? 0 : 1);
                             }
                             //the triangle is in, starting or ending in the plane to split by -> the PlanarEventType signals its position then
                             else if (event.type == PlaneEventType::planar) {
                                 //minSide specifies where to include planar faces
                                 insertIfAbsent(event.faceIndex, minSide ? 0 : 1);
                             }
                             //the face starts in the plane, thus its area overlaps with the bounding box further away from the origin.
                             else if (event.type == PlaneEventType::starting) {
                                 insertIfAbsent(event.faceIndex, 1);
                             }
                             //the face ends in the plane, thus its area overlaps with the bounding box closer to the origin.
                             else {
                                 insertIfAbsent(event.faceIndex, 0);
                             }
                         });
        return {std::move(facesMin), std::move(facesMax)};
    }
} // namespace polyhedralGravity
