#include "polyhedralGravity/model/KDTree/plane_selection/SquaredPlane.h"

namespace polyhedralGravity {
    // O(N^2) implementation
    std::tuple<Plane, double, std::variant<TriangleIndexVectors<2>, PlaneEventVectors<2> > > SquaredPlane::findPlane(
        const SplitParam &splitParam) {
        if (std::holds_alternative<PlaneEventVector>(splitParam.boundFaces)) {
            throw std::invalid_argument("SquaredPlane does not support PlaneEventLists in SplitParam argument");
        }
        const auto &boundFaces = std::get<TriangleIndexVector>(splitParam.boundFaces);
        //initialize the default plane and make it costly
        double cost = std::numeric_limits<double>::infinity();
        Plane optPlane{0, splitParam.splitDirection};
        //store the triangleSets that are implicitly generated during plane testing for later use.
        TriangleIndexVectors<2> optTriangleIndexLists{};
        //each vertex proposes a split plane candidate: test for each of them, store them in buffer set to avoid duplicate testing
        std::unordered_set<double> testedPlaneCoordinates{};
        auto [vertex3_begin, vertex3_end] = transformIterator(boundFaces.cbegin(), boundFaces.cend(),
                                                              splitParam.vertices, splitParam.faces);
        std::mutex optMutex{}, testedPlaneMutex{};
        thrust::for_each(thrust::device, vertex3_begin, vertex3_end,
                         [&splitParam, &optPlane, &cost, &optTriangleIndexLists, &testedPlaneCoordinates, &optMutex, &
                             testedPlaneMutex](
                     const auto &indexAndTriplet) {
                             const auto [index, triplet] = indexAndTriplet;
                             //first clip the triangles vertices to the current bounding box and then get the bounding box of the clipped triangle -> use the box edges as split plane candidates
                             const auto clippedVertices = splitParam.boundingBox.clipToVoxel(triplet);
                             const auto [minPoint, maxPoint] = Box::getBoundingBox<std::vector<Array3> >(
                                 clippedVertices);
                             for (const auto planeSurfacePoint: {minPoint, maxPoint}) {
                                 //constructs the plane that goes through a vertex lying on the bounding box of the face to be checked and spans in a specified direction.
                                 Plane candidatePlane{
                                     planeSurfacePoint[static_cast<int>(splitParam.splitDirection)],
                                     splitParam.splitDirection
                                 }; {
                                     //continue if plane has already been tested
                                     std::lock_guard lock{testedPlaneMutex};
                                     if (testedPlaneCoordinates.find(candidatePlane.axisCoordinate) !=
                                         testedPlaneCoordinates.cend()) {
                                         continue;
                                     }
                                     testedPlaneCoordinates.emplace(candidatePlane.axisCoordinate);
                                 }

                                 auto triangleIndexLists = containedTriangles(splitParam, candidatePlane);

                                 //evaluate the candidate plane and store if it is better than the currently stored result
                                 auto [candidateCost, minSideChosen] = costForPlane(
                                     splitParam.boundingBox, candidatePlane, triangleIndexLists[0]->size(),
                                     triangleIndexLists[1]->size(), triangleIndexLists[2]->size()); {
                                     std::lock_guard lock(optMutex);
                                     // this if clause exists to consistently build the same KDTree (choose plane with lower coordinate) by eliminating indeterministic behavior should the cost be equal.
                                     // this is not important for functionality but for testing purposes
                                     if (candidateCost == cost && optPlane.axisCoordinate < candidatePlane.
                                         axisCoordinate) {
                                         continue;
                                     }
                                     if (candidateCost <= cost) {
                                         cost = candidateCost;
                                         optPlane = candidatePlane;
                                         //planar faces have to be included in one of the two sub boxes.
                                         const auto &includePlanarTo = triangleIndexLists[minSideChosen ? 0 : 1];
                                         includePlanarTo->insert(includePlanarTo->cend(),
                                                                 triangleIndexLists[2]->cbegin(),
                                                                 triangleIndexLists[2]->cend());
                                         optTriangleIndexLists = {
                                             std::move(triangleIndexLists[0]), std::move(triangleIndexLists[1])
                                         };
                                     }
                                 }
                             }
                         });
        return std::make_tuple(optPlane, cost, std::move(optTriangleIndexLists));
    }

    TriangleIndexVectors<3> SquaredPlane::containedTriangles(const SplitParam &splitParam, const Plane &split) {
        using namespace polyhedralGravity;
        if (std::holds_alternative<PlaneEventVector>(splitParam.boundFaces)) {
            throw std::invalid_argument("SquaredPlane does not support PlaneEventLists in SplitParam argument");
        }
        const auto &boundFaces = std::get<TriangleIndexVector>(splitParam.boundFaces);
        //define three sets of triangles: closer to the origin, further away, in the plane
        auto index_less = std::make_unique<TriangleIndexVector>();
        auto index_greater = std::make_unique<TriangleIndexVector>();
        auto index_equal = std::make_unique<TriangleIndexVector>();
        index_less->reserve(boundFaces.size() / 2);
        index_greater->reserve(boundFaces.size() / 2);


        //perform check for every triangle contained in this node's bounding box.
        //transform faceIndices into the vertices
        auto [begin, end] = transformIterator(boundFaces.cbegin(), boundFaces.cend(), splitParam.vertices,
                                              splitParam.faces);
        std::for_each(
            begin, end,
            [&splitParam, &split, &index_greater, &index_less, &index_equal](std::pair<size_t, Array3Triplet> pair) {
                auto [faceIndex, vertices] = pair;
                bool less{false}, greater{false};
                auto clippedVertices = splitParam.boundingBox.clipToVoxel(vertices);
                for (const Array3 vertex: clippedVertices) {
                    //vertex is closer to the origin than the plane
                    if (vertex[static_cast<int>(split.orientation)] < split.axisCoordinate && !less) {
                        less = true;
                        //triangle has area in the closer bounding box and needs to be checked there for intersections
                        index_less->push_back(faceIndex);
                    }
                    //vertex is farther away of the origin than the plane
                    else if (vertex[static_cast<int>(split.orientation)] > split.axisCoordinate && !greater) {
                        greater = true;
                        //triangle has area in the greater bounding box and needs to be checked there for intersections
                        index_greater->push_back(faceIndex);
                    }
                }
                //all vertices of the triangle lie in the plane -> triangle lies in the plane
                if (!less && !greater) {
                    index_equal->push_back(faceIndex);
                }
            });
        return std::array{std::move(index_less), std::move(index_greater), std::move(index_equal)};
    }
} // namespace polyhedralGravity
