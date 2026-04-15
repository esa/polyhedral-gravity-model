#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <variant>
#include <vector>

#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/model/KDTree/KdDefinitions.h"
#include "polyhedralGravity/model/KDTree/SplitParam.h"
#include "polyhedralGravity/model/KDTree/plane_selection/PlaneSelectionAlgorithm.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "thrust/detail/execution_policy.h"
#include "thrust/execution_policy.h"
#include "thrust/iterator/iterator_facade.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/system/detail/sequential/for_each.h"

namespace polyhedralGravity {
struct SplitParam;

    /**
     * Helper class to keep count of triangles that straddle split planes, while PlaneEvents are iterated over.
     */
    class TriangleCounter {
    public:
        /**
         * Update the triangle count for the max side of the plane.
         * @param direction For which dimension to update values.
         * @param p_planar Amount of faces lying in the plane.
         * @param p_end Amount of faces ending in the plane.
         */
        void updateMax(Direction direction, size_t p_planar, size_t p_end);

        /**
        * Update the triangle count for the min side of the plane.
        * @param direction For which dimension to update values.
        * @param p_planar Amount of faces lying in the plane.
        * @param p_start Amount of faces starting in the plane.
        */
        void updateMin(Direction direction, size_t p_planar, size_t p_start);

        /**
         * Sets the amount of planar faces for a specific dimension.
         * @param direction For which dimension to update values.
         * @param p_planar Amount of faces lying in the plane.
         */
        void setPlanar(Direction direction, size_t p_planar);

        /**
         * Returns min value for a specific dimension.
         * @param direction The dimension.
         * @return The min value.
         */
        [[nodiscard]] size_t getMin(Direction direction) const;

        /**
         * Returns max value for a specific dimension.
         * @param direction The dimension.
         * @return The max value.
         */
        [[nodiscard]] size_t getMax(Direction direction) const;

        /**
         * Returns the amount of planar faces for a specific dimension.
         * @param direction The dimension.
         * @return The planar face amount.
         */
        [[nodiscard]] size_t getPlanar(Direction direction) const;

        TriangleCounter(size_t dimensionCount, const std::array<size_t, 3> &initialValues);

    private:
        /**
         * Keeps track of triangle counts for every dimension in order MIN, MAX, PLANAR
         */
        std::vector<std::array<size_t, 3> > dimensionTriangleValues;
    };

    class PlaneEventAlgorithm : public PlaneSelectionAlgorithm {
    protected:
        /**
        * Generates the vector of PlaneEvents comprising all the possible candidate planes using an index list of faces. {@link PlaneEvent}
        * @param splitParam Contains the parameters of the scene to find candidate planes for. {@link SplitParam}
        * @param directions For which directions to generate the events for.
        * @return The vector of PlaneEvents.
        */
        static PlaneEventVector generatePlaneEventsFromFaces(const SplitParam &splitParam,
                                                             std::vector<Direction> directions);

        /**
         * Iterates over PlaneEvents and determines the optimal split plane.
         * @param events The events to base calculations on.
         * @param triangleCounter Used to track face count during iteration.
         * @param boundingBox The current node's bounding box
         * @return Tuple of optimal plane, its cost and where to include planar faces.
         */
        static std::tuple<Plane, double, bool> traversePlaneEvents(const PlaneEventVector &events,
                                                                   TriangleCounter &triangleCounter,
                                                                   const Box &boundingBox);
    };
} // namespace polyhedralGravity
