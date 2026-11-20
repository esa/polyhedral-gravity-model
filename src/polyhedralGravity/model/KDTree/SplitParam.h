#pragma once

#include "polyhedralGravity/model/KDTree/KdDefinitions.h"

namespace polyhedralGravity {
    //forward declaration
    class PlaneSelectionAlgorithm;

    /**
     * Helper struct to bundle important parameters required for splitting a Polyhedron for better readability.
     */
    struct SplitParam {
        /**
         * The vertices that compose the Polyhedron.
         */
        const std::vector<Array3> &vertices;
        /**
         * The faces that connect the vertices to render the Polyhedron.
         */
        const std::vector<IndexArray3> &faces;
        /**
         * Either an index list of faces that are included in the current bounding box of the KDTree or a list of PlaneEvents containing the information about thr bound faces. Important when building deeper levels of a KDTree.
         */
        std::variant<TriangleIndexVector, PlaneEventVector> boundFaces;
        /**
         * The current bounding box that should be divided further by the KDTree.
         */
        Box boundingBox;
        /**
         * The direction in which the current bounding box should be divided by further.
         * Refer to {@link Plane} on how to interpret the Direction.
         */
        mutable Direction splitDirection;
        /**
         * The algorithm used to create new child TreeNodes after splitting the parent.
         */
        const std::shared_ptr<PlaneSelectionAlgorithm> planeSelectionStrategy;

        /**
         * Constructor that initializes all fields. Intended for the use with std::make_unique. See {@link SplitParam} fields for further information.
         *
         */
        SplitParam(const std::vector<Array3> &vertices, const std::vector<IndexArray3> &faces, const Box &boundingBox,
                   const Direction splitDirection,
                   const std::shared_ptr<PlaneSelectionAlgorithm> &planeSelectionStrategy)
            : vertices{vertices}, faces{faces}, boundFaces{TriangleIndexVector(faces.size())}, boundingBox{boundingBox},
              splitDirection{splitDirection}, planeSelectionStrategy{planeSelectionStrategy} {
            auto &indexList = std::get<TriangleIndexVector>(boundFaces);
            std::iota(indexList.begin(), indexList.end(), 0);
        }

        /**
         * Constructor manually initializing boundFaces, used for testing.
         */
        SplitParam(const std::vector<Array3> &vertices, const std::vector<IndexArray3> &faces,
                   const std::variant<TriangleIndexVector, PlaneEventVector> &boundFaces, const Box &boundingBox,
                   const Direction splitDirection,
                   const std::shared_ptr<PlaneSelectionAlgorithm> &planeSelectionStrategy)
            : vertices{vertices}, faces{faces}, boundFaces{boundFaces}, boundingBox{boundingBox},
              splitDirection{splitDirection}, planeSelectionStrategy{planeSelectionStrategy} {
        }
    };
} // namespace polyhedralGravity
