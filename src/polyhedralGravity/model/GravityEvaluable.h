#pragma once

#include <tuple>
#include <variant>
#include <string>
#include <optional>
#include <sstream>

#include "thrust/transform.h"
#include "thrust/execution_policy.h"

#include "GravityModelDetail.h"
#include "MeshChecking.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/input/TetgenAdapter.h"
#include "GravityModelData.h"
#include "Polyhedron.h"


namespace polyhedralGravity {

    /**
     * Class for evaluating the polyhedrale gravity model for a given constant density polyhedron.
     * Caches the polyhedron and data which is independent of the computation point P.
     * Provides an operator() for evaluating the polyhedrale gravity model for a given constant density polyhedron
     * at computation point P and choosing between parallel and serial evaluation.
     */
    class GravityEvaluable {

        /** The constant density polyhedron consisting of vertices and triangular faces */
        const Polyhedron _polyhedron;

        /** The constant density of the polyhedron (the unit must match to the mesh, e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$) */
        const double _density;

        /** The orientation of the normals of the polyhedron, either outwards or inwards pointing the polyhedral source */
        const NormalOrientation _orientation;

        /** Cache for the segment vectors (segments between vertices of a polyhedral face) */
        mutable std::vector<Array3Triplet> _segmentVectors;

        /** Cache for the plane unit normals (unit normals of the polyhedral faces) */
        mutable std::vector<Array3> _planeUnitNormals;

        /** Cache for the segment unit normals (unit normals of each the polyhedral faces' segments) */
        mutable std::vector<Array3Triplet> _segmentUnitNormals;

    public:
        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron.
         * @param polyhedralSource the polyhedron as polyhedron, vertices and faces, or a list of files
         * @param density the constant density (the unit must match to the mesh,
         *                e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         * @param orientation the orientation of the plane unit normals of the polyhedron (either pointing outwards or inwards the polyhedron)
         *                (default: OUTWARDS as this is the original assumption of Tsoulis et al.'s formulation)
         * @param check if true, the polyhedral mesh will be checked for degenerated triangular surfaces and
         *                that the orientation of the normals are outwards pointing/ pointing to as specified by the parameter orientation
         *                The latter check requires @f$O(n^2)$@f operations!
         *                By default the check is enabled and will inform the user about inconstencies and the quadratic runtime.
         *                Explcitly enableing silence the @f$O(n^2)$@f warning, but produces an exception in case of inconstency
         *                Explcity disableing silences everythign and avoids the @f$O(n^2)$@f runtime cost.
         *
         * @note This is a separate constructor since the polyhedron as a class it not exposed to the user via
         * the Python Interface. Thus, this constructor is only available via the C++ interface.
         *
         * @warning The input check requires @f$O(n^2)$@f operations and is enabled by default. Disable it when "you know your mesh"!
         *
         * @throws If inputCheck is true (or not set) and the orientation does not match with the polyhedral match a argument exception
         */
        GravityEvaluable(
                const PolyhedralSource &polyhedralSource,
                double density,
                const NormalOrientation &orientation = NormalOrientation::OUTWARDS,
                const std::optional<bool> check = std::nullopt
                ) :
            _polyhedron{
                    std::visit(util::overloaded{
                                       [](const Polyhedron &arg) { return arg; },
                                       [](const PolyhedralFileSource &arg) { return TetgenAdapter{arg}.getPolyhedron(); },
                                       [](const PolyhedralDataSource &arg) { return Polyhedron{arg}; }
                               }, polyhedralSource)
            },
            _density{density},
            _orientation{orientation} {
            this->runMeshCeck(check);
            this->prepare();
        }

        /**
         * Instantiates a GravityEvaluable with a given constant density polyhedron and caches.
         * This is for restoring a GravityEvaluable from a previous state.
         * @param polyhedron the polyhedron
         * @param density the constant density (the unit must match to the mesh,
         *                e.g., mesh in @f$[m]@f$ requires density in @f$[kg/m^3]@f$)
         * @param orientation the orientation of the plane unit normals of the polyhedron (either pointing outwards or inwards the polyhedron)
         *                (default: OUTWARDS as this is the original assumption of Tsoulis et al.'s formulation)
         * @param segmentVectors the segment vectors
         * @param planeUnitNormals the plane unit normals
         * @param segmentUnitNormals the segment unit normals
         */
        GravityEvaluable(const Polyhedron &polyhedron,
                         double density,
                         const NormalOrientation &orientation,
                         const std::vector<Array3Triplet> &segmentVectors,
                         const std::vector<Array3> &planeUnitNormals,
                         const std::vector<Array3Triplet> &segmentUnitNormals) :
            _polyhedron{polyhedron},
            _density{density},
            _orientation{orientation},
            _segmentVectors{segmentVectors},
            _planeUnitNormals{planeUnitNormals},
            _segmentUnitNormals{
                    segmentUnitNormals} {
        }

        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
         * point P. Wrapper for evaluate<parallelization>.
         * @param computationPoints the computation point P or multiple computation points in a vector
         * @param parallelization if true, the calculation is parallelized
         * @return the GravityModelResult containing the potential, acceleration and second derivative
         */
        inline std::variant<GravityModelResult, std::vector<GravityModelResult>>
        operator()(const std::variant<Array3, std::vector<Array3>> &computationPoints,
                   bool parallelization = true) const {
            if (parallelization) {
                if (std::holds_alternative<Array3>(computationPoints)) {
                    return this->evaluate<true>(std::get<Array3>(computationPoints));
                } else {
                    return this->evaluate<true>(std::get<std::vector<Array3>>(computationPoints));
                }
            } else {
                if (std::holds_alternative<Array3>(computationPoints)) {
                    return this->evaluate<false>(std::get<Array3>(computationPoints));
                } else {
                    return this->evaluate<false>(std::get<std::vector<Array3>>(computationPoints));
                }
            }
        }

        /**
         * Returns a string representation of the GravityEvaluable.
         * @return string representation of the GravityEvaluable
         */
        std::string toString() const;

        /**
         * Returns the polyhedron, the density and the internal caches.
         *
         * @return tuple of polyhedron, density, segmentVectors, planeUnitNormals and segmentUnitNormals
         */
        std::tuple<Polyhedron, double, NormalOrientation, std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
        getState() const;

    private:
        /**
         * Runs the Mesh Sanity Check depending on the flag.
         * If nullopt, the check is run and the user is informed about a potential unwanted @f$O(n^2)$@f runtime cost
         * If true, the check is run, but no runtime cost warning is printed
         * if false, the check is disabled and nothing at all is done
         * @param flag to determine if what is done
         *
         * @throws std::runtime_error in case of wrong mesh properties
         */
        void runMeshCeck(const std::optional<bool> &flag) const;

        /**
         * Prepares the polyhedron for the evaluation by calculating the segment vectors, the plane unit normals
         * and the segment unit normals.
         * Called by the constructor once.
         */
        void prepare() const;

        /**
        * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
        * point P.
        * @tparam Parallelization if true, the calculation is parallelized
        * @param computationPoint the computation Point P
        * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration
        * at computation Point P
        */
        template<bool Parallelization = true>
        GravityModelResult evaluate(const Array3 &computationPoint) const;


        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation
         * at multiple computation points.
         * @tparam Parallelization if true, the calculation is parallelized
         * @param computationPoints the computation Points
         * @return vector of GravityModelResults containing the potential, the acceleration and the change of acceleration
         */
        template<bool Parallelization = true>
        std::vector<GravityModelResult> evaluate(const std::vector<Array3> &computationPoints) const;

        /**
         * Evaluates the polyhedrale gravity model for a given constant density polyhedron at computation a certain face.
         * @param tuple consisting of face, segmentVectors, planeUnitNormal and segmentUnitNormals
         * @return the GravityModelResult containing the potential, the acceleration and the change of acceleration which
         * this face contributes to the computation point
         */
        static GravityModelResult
        evaluateFace(const thrust::tuple<Array3Triplet, Array3Triplet, Array3, Array3Triplet> &tuple);

    };

}// namespace polyhedralGravity