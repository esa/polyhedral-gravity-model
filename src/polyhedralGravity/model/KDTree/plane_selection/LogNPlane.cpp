#include "polyhedralGravity/model/KDTree/plane_selection/LogNPlane.h"

namespace polyhedralGravity {
    // O(N*log^2(N)) implementation
    std::tuple<Plane, double, std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2> > > LogNPlane::findPlane(
        const SplitParam &splitParam) {
        const PlaneEventVector events{std::move(generatePlaneEvents(splitParam))};
        TriangleCounter triangleCounter{3, {0, countFaces(splitParam.boundFaces), 0}};
        auto [optPlane, cost, minSide] = traversePlaneEvents(events, triangleCounter, splitParam.boundingBox);
        //generate the triangle index lists for the child bounding boxes and return them along with the optimal plane and the plane's cost.
        return {optPlane, cost, generatePlaneEventSubsets(splitParam, events, optPlane, minSide)};
    }

    PlaneEventVector LogNPlane::generatePlaneEvents(const SplitParam &splitParam) {
        if (std::holds_alternative<TriangleIndexVector>(splitParam.boundFaces)) {
            return generatePlaneEventsFromFaces(splitParam, ALL_DIRECTIONS);
        }
        return std::get<PlaneEventVector>(splitParam.boundFaces);
    }

    PlaneEventVectors<2> LogNPlane::generatePlaneEventSubsets(const SplitParam &splitParam,
                                                              const PlaneEventVector &planeEvents, const Plane &plane,
                                                              const bool minSide) {
        const auto faceClassification{classifyTrianglesRelativeToPlane(planeEvents, plane, minSide)};
        PlaneEventVector planeEventsMin{};
        PlaneEventVector planeEventsMax{};
        TriangleIndexVector facesIndexBoth{};
        planeEventsMin.reserve(planeEvents.size() / 2);
        planeEventsMax.reserve(planeEvents.size() / 2);
        //value estimation taken from source paper
        facesIndexBoth.reserve(std::ceil(std::sqrt(planeEvents.size())));
        std::unordered_set<size_t> processedIndices{};

        auto insertToBothIfAbsent = [&facesIndexBoth, &processedIndices](const auto faceIndex) {
            if (processedIndices.find(faceIndex) == processedIndices.end()) {
                processedIndices.insert(faceIndex);
                facesIndexBoth.push_back(faceIndex);
            }
        };

        //Step 2
        std::for_each(planeEvents.cbegin(), planeEvents.cend(),
                      [&faceClassification, &planeEventsMin, &planeEventsMax, &insertToBothIfAbsent
                      ](const auto &event) {
                          switch (faceClassification.at(event.faceIndex)) {
                              //face of event only contributes to min side event can be added to side without clipping because no overlap with split plane
                              case Locale::MIN_ONLY:
                                  planeEventsMin.push_back(event);
                                  break;
                              //face of event only contributes to max side event can be added to side without clipping because no overlap with split plane
                              case Locale::MAX_ONLY:
                                  planeEventsMax.push_back(event);
                                  break;
                              //face has area on both sides -> event has to be discarded and scheduled for separate event generation
                              case Locale::BOTH:
                              default:
                                  insertToBothIfAbsent(event.faceIndex);
                          }
                      });

        //generate new plane events for straddling faces that were discarded previously in Step 2
        auto [newMinEvents, newMaxEvents] = generatePlaneEventsForClippedFaces(splitParam, facesIndexBoth, plane);
        //merge the new events into the existing sorted lists and return
        return {
            std::move(mergePlaneEventLists(planeEventsMin, newMinEvents)),
            std::move(mergePlaneEventLists(planeEventsMax, newMaxEvents))
        };
    }

    //Step 1
    std::unordered_map<size_t, LogNPlane::Locale> LogNPlane::classifyTrianglesRelativeToPlane(
        const PlaneEventVector &events, const Plane &plane, const bool minSide) {
        std::unordered_map<size_t, Locale> result{};
        //each face generates 6 plane events on average, thus the amount of faces can be roughly estimated.
        result.reserve(events.size() / 6);
        //preparing the map by initializing all faces with them having area in both sub bounding boxes
        std::for_each(events.cbegin(), events.cend(), [&result](const auto &event) {
            result[event.faceIndex] = Locale::BOTH;
        });
        //now search for conditions proving that the faces DO NOT have area in both boxes
        std::for_each(events.begin(), events.end(), [minSide, &result, &plane](const auto &event) {
            if (event.type == PlaneEventType::ending && event.plane.orientation == plane.orientation && event.plane.
                axisCoordinate <= plane.axisCoordinate) {
                result[event.faceIndex] = Locale::MIN_ONLY;
            } else if (event.type == PlaneEventType::starting && event.plane.orientation == plane.orientation && event.
                       plane.axisCoordinate >= plane.axisCoordinate) {
                result[event.faceIndex] = Locale::MAX_ONLY;
            } else if (event.type == PlaneEventType::planar && event.plane.orientation == plane.orientation) {
                if (event.plane.axisCoordinate < plane.axisCoordinate || (
                        event.plane.axisCoordinate == plane.axisCoordinate && minSide)) {
                    result[event.faceIndex] = Locale::MIN_ONLY;
                }
                if (event.plane.axisCoordinate > plane.axisCoordinate || (
                        event.plane.axisCoordinate == plane.axisCoordinate && !minSide)) {
                    result[event.faceIndex] = Locale::MAX_ONLY;
                }
            }
        });
        return result;
    }

    //Step 3
    std::array<PlaneEventVector, 2> LogNPlane::generatePlaneEventsForClippedFaces(
        const SplitParam &splitParam, const TriangleIndexVector &faceIndices, const Plane &plane) {
        auto [minBox, maxBox] = splitParam.boundingBox.splitBox(plane);
        PlaneEventVector minEvents{};
        PlaneEventVector maxEvents{};
        //each face generates six new PlaneEvents and each face has area in both boxes
        minEvents.resize(faceIndices.size() * 6);
        maxEvents.resize(faceIndices.size() * 6);

        //lambda for creating PlaneEvents from a vertex triplet (face) in one of the two sub boxes
        const auto createPlaneEvents = [](const auto &vertices, const auto &boundingBox, const size_t faceIndex,
                                          auto destIt) {
            //clip to the voxel
            auto clipped = boundingBox.clipToVoxel(vertices);
            //create split plane anchor points using the bounding box
            const auto [minPoint, maxPoint] = Box::getBoundingBox(clipped);
            //associate parameters for PlaneEvent creation
            std::array<std::pair<const Array3, PlaneEventType>, 2> planeEventParam{
                std::make_pair(minPoint, PlaneEventType::starting),
                std::make_pair(maxPoint, PlaneEventType::ending)
            };
            //create planes in each dimension, be careful to cluster similar anchor points together.
            size_t planeIndex = 0;
            for (const auto &[point, eventType]: planeEventParam) {
                for (const auto &direction: ALL_DIRECTIONS) {
                    //insert directly for paralellization
                    *(destIt + planeIndex++) = PlaneEvent(eventType, Plane(point, direction), faceIndex);
                }
            }
        };

        //transform faces to vertices
        auto [begin_it, end_it] = transformIterator(faceIndices.cbegin(), faceIndices.cend(), splitParam.vertices,
                                                    splitParam.faces);
        std::atomic_long minIndex{0};
        std::atomic_long maxIndex{0};
        //create new events for each face in both sub boxes
        thrust::for_each(thrust::device, begin_it, end_it,
                         [&minBox, maxBox, &minEvents, &maxEvents, &createPlaneEvents, &minIndex, &maxIndex](
                     const auto &indexAndTriplet) {
                             const auto &[index, vertexTriplet] = indexAndTriplet;
                             //reserve slots of 6 for the threads using the atomic counters. Size fits because of earlier resize
                             createPlaneEvents(vertexTriplet, minBox, index, minEvents.begin() + (minIndex++ * 6));
                             createPlaneEvents(vertexTriplet, maxBox, index, maxEvents.begin() + (maxIndex++ * 6));
                         });

        //sort the lists for later merge sort integration
        std::sort(minEvents.begin(), minEvents.end());
        std::sort(maxEvents.begin(), maxEvents.end());

        return {minEvents, maxEvents};
    }

    //Step 4
    std::unique_ptr<PlaneEventVector> LogNPlane::mergePlaneEventLists(const PlaneEventVector &first,
                                                                      const PlaneEventVector &second) {
        auto first_it{first.cbegin()};
        auto second_it{second.cbegin()};
        auto result{std::make_unique<PlaneEventVector>()};
        result->reserve(first.size() + second.size());

        while (first_it != first.cend() || second_it != second.cend()) {
            if (first_it == first.cend() || second_it == second.cend()) {
                result->push_back(first_it == first.cend() ? *second_it++ : *first_it++);
            } else if (*first_it < *second_it) {
                result->push_back(*first_it++);
            } else {
                result->push_back(*second_it++);
            }
        }
        return result;
    }
} // namespace polyhedralGravity
