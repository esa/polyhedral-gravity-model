#pragma once

#include <array>
#include <limits>
#include <utility>

#include <Kokkos_Core.hpp>

#include "GravityModelData.h"
#include "polyhedralGravity/model/PolyhedralMeshView.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "polyhedralGravity/util/UtilityFloatArithmetic.h"

/**
 * The single-source implementation of Tsoulis' polyhedral gravity model for one polyhedral face.
 *
 * Every function in here is a KOKKOS_INLINE_FUNCTION and templated over the floating point precision,
 * so the very same code runs inside a Kokkos kernel on a GPU, inside a Kokkos kernel on the CPU, and
 * directly on the host (which is what the unit tests do, always in double precision).
 * Consequently, the functions must not allocate, must not throw, and must only call math functions
 * from the Kokkos:: namespace, which are available in both worlds.
 *
 * The individual steps of Tsoulis' algorithm come first, {@link evaluateFace} assembles them into the
 * body of the kernel which a {@link GravityEvaluable} launches over the faces of a polyhedron.
 */
namespace polyhedralGravity::GravityModel::detail {

    /**
     * The radius around zero which is treated as zero, in the evaluation's own precision.
     *
     * This is {@link util::EPSILON_ZERO_OFFSET} scaled by how much coarser the given precision is than double
     * precision, i.e. it stays at @f$10^{-14}@f$ for double and becomes roughly @f$5 \cdot 10^{-6}@f$ for float.
     * Scaling it is not cosmetic: the singularity cases of Tsoulis' algorithm are detected by comparing
     * distances against this radius, and a threshold below the precision's own resolution makes the
     * evaluation miss them, which produces infinities and NaNs instead of the singularity terms.
     *
     * @tparam FloatType the floating point precision of the evaluation
     */
    template<typename FloatType>
    constexpr FloatType EPSILON_ZERO = static_cast<FloatType>(
            util::EPSILON_ZERO_OFFSET *
            (static_cast<double>(std::numeric_limits<FloatType>::epsilon()) / std::numeric_limits<double>::epsilon()));

    /**
     * The singularity terms sing A and sing B of one plane p of the polyhedron.
     * @tparam FloatType the floating point precision of the evaluation
     */
    template<typename FloatType>
    struct SingularityTerms {
        /** The scalar singularity term sing A used for the potential and the acceleration */
        FloatType alpha;
        /** The vectorial singularity term sing B used for the gradiometric tensor */
        Vector3<FloatType> beta;
    };

    /**
     * Computes the segment vectors G_ij for one plane of the polyhedron according to Tsoulis (18).
     * The segment vectors G_ij represent the vector from one vertex of the face to the neighboring vertex and
     * depict every line segment of the triangular face (A-B-C)
     * @param vertex0 the first vertex A
     * @param vertex1 the second vertex B
     * @param vertex2 the third vertex C
     * @return the segment vectors for a plane
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3Triplet<FloatType> buildVectorsOfSegments(
            const Vector3<FloatType> &vertex0, const Vector3<FloatType> &vertex1, const Vector3<FloatType> &vertex2) {
        using util::operator-;
        //Calculate G_ij
        return {vertex1 - vertex0, vertex2 - vertex1, vertex0 - vertex2};
    }

    /**
     * Computes the plane unit normal N_p for one plane p of the polyhedron according to Tsoulis (19).
     * The plane unit normal is the outward pointing normal of the face from the polyhedron.
     * @param segmentVector1 first edge
     * @param segmentVector2 second edge
     * @return plane unit normal
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3<FloatType> buildUnitNormalOfPlane(
            const Vector3<FloatType> &segmentVector1, const Vector3<FloatType> &segmentVector2) {
        //Calculate N_i as (G_i1 * G_i2) / |G_i1 * G_i2| with * being the cross product
        return util::normal(segmentVector1, segmentVector2);
    }

    /**
     * Computes the segment unit normals n_pq for one plane p of the polyhedron according to Tsoulis (20).
     * The segment unit normal n_pq represent the normal of one line segment of a polyhedrale face.
     * @param segmentVectors the segment vectors of the face G_p(0-2)
     * @param planeUnitNormal the plane unit normal N_p
     * @return segment unit normals n_pq for plane p with q = {0, 1, 2}
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3Triplet<FloatType> buildUnitNormalOfSegments(
            const Vector3Triplet<FloatType> &segmentVectors, const Vector3<FloatType> &planeUnitNormal) {
        Vector3Triplet<FloatType> segmentUnitNormal{};
        //Calculate n_ij as (G_ij * N_i) / |G_ig * N_i| with * being the cross product
        for (size_t j = 0; j < 3; ++j) {
            segmentUnitNormal[j] = util::normal(segmentVectors[j], planeUnitNormal);
        }
        return segmentUnitNormal;
    }

    /**
     * Computes the plane unit normal orientation/ direction sigma_p for one plane p of the polyhedron
     * according to Tsoulis (21).
     * The plane unit normal orientation values represents the relative position of computation point P
     * with respect to the pointing direction of N_p. E. g. if N_p points to the half-space containing P, the
     * inner product of N_p and -G_i1 will be positive, leading to a negative sigma_p.
     * If sigma_p is zero than P and P' lie geometrically in the same plane --> P == P'.
     * @param planeUnitNormal the plane unit normal N_p
     * @param vertex0 the first vertex of the plane
     * @return plane normal orientation
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION FloatType computeUnitNormalOfPlaneDirection(
            const Vector3<FloatType> &planeUnitNormal, const Vector3<FloatType> &vertex0) {
        //Calculate N_i * -G_i1 where * is the dot product and then use the inverted sgn
        //We abstain on the double multiplication with -1 in the line above and beyond since two
        //times multiplying with -1 equals no change
        return static_cast<FloatType>(util::sgn(util::dot(planeUnitNormal, vertex0), EPSILON_ZERO<FloatType>));
    }

    /**
     * Calculates the Hessian Plane form spanned by three given points p, q, and r.
     * @param p first point on the plane
     * @param q second point on the plane
     * @param r third point on the plane
     * @return HessianPlane
     * @related Cross-Product method https://tutorial.math.lamar.edu/classes/calciii/eqnsofplanes.aspx
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION HessianPlaneTemplate<FloatType> computeHessianPlane(
            const Vector3<FloatType> &p, const Vector3<FloatType> &q, const Vector3<FloatType> &r) {
        using util::operator-;
        using util::operator*;
        constexpr Vector3<FloatType> origin{0.0, 0.0, 0.0};
        const Vector3<FloatType> crossProduct = util::cross(p - q, p - r);
        const Vector3<FloatType> res = (origin - p) * crossProduct;
        const FloatType d = res[0] + res[1] + res[2];

        return {crossProduct[0], crossProduct[1], crossProduct[2], d};
    }

    /**
     * Calculates the (plane) distances h_p of computation point P to the plane S_p given in Hessian Form
     * according to the following equation:
     * h_p = D / sqrt(A^2+B^2+C^2)
     * @param hessianPlane Hessian Plane Form of S_p
     * @return plane distance h_p
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION FloatType distanceBetweenOriginAndPlane(const HessianPlaneTemplate<FloatType> &hessianPlane) {
        //Compute h_p as D/sqrt(A^2 + B^2 + C^2)
        return Kokkos::abs(hessianPlane.d / Kokkos::sqrt(hessianPlane.a * hessianPlane.a +
                                                         hessianPlane.b * hessianPlane.b +
                                                         hessianPlane.c * hessianPlane.c));
    }

    /**
     * Computes P' for a given plane p according to equation (22) of Tsoulis paper.
     * P' is the orthogonal projection of the computation point P onto the plane S_p.
     * @param planeUnitNormal the plane unit normal N_p
     * @param planeDistance the distance from P to the plane h_p
     * @param hessianPlane the Hessian Plane Form
     * @return P' for this plane
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3<FloatType> projectPointOrthogonallyOntoPlane(
            const Vector3<FloatType> &planeUnitNormal, FloatType planeDistance,
            const HessianPlaneTemplate<FloatType> &hessianPlane) {
        using util::operator*;
        //Calculate the projection point by (22) P'_ = N_i / norm(N_i) * h_i
        // norm(N_i) is always 1 since N_i is a "normed" vector --> we do not need this division
        Vector3<FloatType> orthogonalProjectionPoint = planeUnitNormal * planeDistance;

        //Calculate alpha, beta and gamma as D/A, D/B and D/C (Notice that we "forget" the minus before those
        // divisions. In consequence, the conditions for signs are reversed below!!!)
        // These values represent the intersections of each polygonal plane with the axes
        // Comparison x == 0.0 is ok, since we only want to avoid nan values
        const Vector3<FloatType> intersections = {
                hessianPlane.a == 0.0 ? static_cast<FloatType>(0.0) : hessianPlane.d / hessianPlane.a,
                hessianPlane.b == 0.0 ? static_cast<FloatType>(0.0) : hessianPlane.d / hessianPlane.b,
                hessianPlane.c == 0.0 ? static_cast<FloatType>(0.0) : hessianPlane.d / hessianPlane.c};

        //Determine the signs of the coordinates of P' according to the intersection values alpha, beta, gamma
        // denoted as __ below, i.e. -alpha, -beta, -gamma denoted -__
        for (size_t index = 0; index < 3; ++index) {
            if (intersections[index] < 0) {
                //If -__ >= 0 --> __ < 0 then coordinates are positive, we calculate abs(orthogonalProjectionPoint[..])
                orthogonalProjectionPoint[index] = Kokkos::abs(orthogonalProjectionPoint[index]);
            } else if (planeUnitNormal[index] > 0) {
                //If -__ < 0 --> __ >= 0 then the coordinate is negative -orthogonalProjectionPoint[..]
                orthogonalProjectionPoint[index] = static_cast<FloatType>(-1.0) * orthogonalProjectionPoint[index];
            }
            //Else the coordinate is positive and stays as it is
        }
        return orthogonalProjectionPoint;
    }

    /**
     * Computes P' for a given plane p according to equation (22) of Tsoulis paper, from the plane's unit
     * normal alone.
     *
     * This is the same projection as the overload above, but it does not need the Hessian form of the plane.
     * The Hessian normal @f$(A, B, C)@f$ is @f$|G_{p1} \times G_{p2}|@f$ times the plane unit normal @f$N_p@f$,
     * and @f$D@f$ is the same positive multiple of the signed plane distance, so the quotients @f$D/A@f$,
     * @f$D/B@f$, and @f$D/C@f$ which decide the signs of P' below are exactly
     * @f$ \tilde{h}_p / N_{p,x} @f$ and its siblings. Since @f$N_p@f$ is cached per face, this turns Tsoulis'
     * steps 1-05 to 1-07 into a single dot product per face and computation point instead of a cross product,
     * a square root, and a division.
     *
     * @param planeUnitNormal the plane unit normal N_p
     * @param planeDistance the distance from P to the plane h_p, i.e. the magnitude of the signed distance
     * @param signedPlaneDistance the signed distance of P to the plane, i.e. -(N_p * v_0) with v_0 relative to P
     * @return P' for this plane
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3<FloatType> projectPointOrthogonallyOntoPlane(
            const Vector3<FloatType> &planeUnitNormal, const FloatType planeDistance,
            const FloatType signedPlaneDistance) {
        using util::operator*;
        //Calculate the projection point by (22) P'_ = N_i / norm(N_i) * h_i
        Vector3<FloatType> orthogonalProjectionPoint = planeUnitNormal * planeDistance;

        //The intersections of the plane with the axes, again without the minus (see the overload above), so
        //the conditions for the signs are the reversed ones. Only the sign of each intersection is ever
        //looked at, and the sign of a quotient is decided by the signs of its two operands, so the three
        //divisions of the overload above are not carried out at all.
        for (size_t index = 0; index < 3; ++index) {
            // Comparison x == 0.0 is ok, since we only want to avoid nan values. A zero numerator or
            // denominator gives an intersection of zero, which is not negative either way.
            const bool intersectionIsNegative = signedPlaneDistance != 0.0 && planeUnitNormal[index] != 0.0 &&
                                                (signedPlaneDistance < 0.0) != (planeUnitNormal[index] < 0.0);
            if (intersectionIsNegative) {
                orthogonalProjectionPoint[index] = Kokkos::abs(orthogonalProjectionPoint[index]);
            } else if (planeUnitNormal[index] > 0) {
                orthogonalProjectionPoint[index] = static_cast<FloatType>(-1.0) * orthogonalProjectionPoint[index];
            }
        }
        return orthogonalProjectionPoint;
    }

    /**
     * Computes the segment normal orientations/ directions sigma_pq for a given plane p.
     * If sigma_pq is negative, this denotes that n_pq points to the half-plane containing P'. Nn case
     * sigma_pq is positive, P' resides in the other half-plane and if sigma_pq is zero, then P' lies directly
     * on the segment pq.
     * @param vertices the vertices of this plane
     * @param projectionPointOnPlane the projection point P' for this plane
     * @param segmentUnitNormalsForPlane the segment unit normals sigma_pq for this plane
     * @return the segment normal orientations for the plane p
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3<FloatType> computeUnitNormalOfSegmentsDirections(
            const Vector3Triplet<FloatType> &vertices, const Vector3<FloatType> &projectionPointOnPlane,
            const Vector3Triplet<FloatType> &segmentUnitNormalsForPlane) {
        using util::operator-;
        Vector3<FloatType> segmentNormalOrientations{};
        //Equation (23)
        //Calculate x_P' - x_ij^1 (x_P' is the projectionPoint and x_ij^1 is the first vertices of one segment,
        //i.e. the coordinates of the training-planes' nodes --> projectionPointOnPlane - vertex
        //Calculate n_ij * x_ij with * being the dot product and use the inverted sgn to determine the value of sigma_pq
        for (size_t j = 0; j < 3; ++j) {
            const FloatType projection = util::dot(segmentUnitNormalsForPlane[j], projectionPointOnPlane - vertices[j]);
            segmentNormalOrientations[j] =
                    static_cast<FloatType>(util::sgn(projection, EPSILON_ZERO<FloatType>)) * static_cast<FloatType>(-1.0);
        }
        return segmentNormalOrientations;
    }

    /**
     * Calculates the point P'' for a given Segment consisting of vertices v1 and v2 and the orthogonal projection
     * point P' for the plane consisting of those vertices. Solves the three equations given in (24), (25) and (26).
     * @param vertex1 first endpoint of segment
     * @param vertex2 second endpoint of segment
     * @param orthogonalProjectionPointOnPlane the orthogonal projection P' of P on this plane
     * @return P'' for this segment
     * @note If sigma_pq is zero then P'' == P', this is not checked by this method, but has to be assured first
     *
     * @note The three equations of (24), (25) and (26) state that P'' lies on the line through the two
     * endpoints and that P'' - P' is perpendicular to that line, i.e. they say nothing more than that P'' is
     * the orthogonal projection of P' onto the line. Solving them by Cramer's rule costs four @f$3 \times 3@f$
     * determinants and three divisions per segment, whereas the foot of the perpendicular
     * @f$ P'' = v_1 + \frac{(P' - v_1) \cdot G}{G \cdot G} G @f$ with @f$ G = v_2 - v_1 @f$ costs two dot
     * products and one division. This is the hottest line of the whole model -- it runs three times per face
     * and computation point -- so the closed form is used instead. It is also the better conditioned of the
     * two, since it never forms the determinant of a matrix built from two nested cross products.
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3<FloatType> projectPointOrthogonallyOntoSegment(
            const Vector3<FloatType> &vertex1, const Vector3<FloatType> &vertex2,
            const Vector3<FloatType> &orthogonalProjectionPointOnPlane) {
        using util::operator-;
        using util::operator+;
        using util::operator*;
        //The segment G_pq spanned by the two endpoints, and P' relative to the first of them
        const Vector3<FloatType> segmentVector = vertex2 - vertex1;
        const Vector3<FloatType> relativeProjectionPoint = orthogonalProjectionPointOnPlane - vertex1;
        //How far along the segment the foot of the perpendicular lies, as a fraction of the segment
        const FloatType relativePosition =
                util::dot(relativeProjectionPoint, segmentVector) / util::dot(segmentVector, segmentVector);
        return vertex1 + segmentVector * relativePosition;
    }

    /**
     * Computes the orthogonal projection Points P'' foreach segment q of a given plane p.
     * @param projectionPointOnPlane the projection Point P'
     * @param segmentNormalOrientations the segment normal orientations sigma_pq for this plane p
     * @param face the vertices of the plane p
     * @return the orthogonal projection points of P on the segment P'' foreach segment q of p
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3Triplet<FloatType> projectPointOrthogonallyOntoSegments(
            const Vector3<FloatType> &projectionPointOnPlane, const Vector3<FloatType> &segmentNormalOrientations,
            const Vector3Triplet<FloatType> &face) {
        Vector3Triplet<FloatType> orthogonalProjectionPointOnSegmentPerPlane{};
        //Running over the segments of this plane
        for (size_t j = 0; j < 3; ++j) {
            //We actually only accept +0.0 or -0.0, so the equal comparison is ok
            if (segmentNormalOrientations[j] == 0.0) {
                //Geometrically trivial case, in neither of the half space --> already on segment
                orthogonalProjectionPointOnSegmentPerPlane[j] = projectionPointOnPlane;
            } else {
                //In one of the half space, evaluate the projection point P'' for the segment
                //with the endpoints v1 and v2
                orthogonalProjectionPointOnSegmentPerPlane[j] =
                        projectPointOrthogonallyOntoSegment(face[j], face[(j + 1) % 3], projectionPointOnPlane);
            }
        }
        return orthogonalProjectionPointOnSegmentPerPlane;
    }

    /**
     * Computes the (segment) distances h_pq between P' for a given plane p and P'' for a given segment q of plane p.
     * @param orthogonalProjectionPointOnPlane the orthogonal projection point P' for p
     * @param orthogonalProjectionPointOnSegments the orthogonal projection points P'' for each segment q of p
     * @return distances h_pq for plane p
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3<FloatType> distancesBetweenProjectionPoints(
            const Vector3<FloatType> &orthogonalProjectionPointOnPlane,
            const Vector3Triplet<FloatType> &orthogonalProjectionPointOnSegments) {
        using util::operator-;
        Vector3<FloatType> segmentDistances{};
        //The inner loop with the running j --> iterating over the segments
        //Using the values P'_i and P''_ij for the calculation of the distance
        for (size_t j = 0; j < 3; ++j) {
            segmentDistances[j] = util::euclideanNorm(
                    orthogonalProjectionPointOnSegments[j] - orthogonalProjectionPointOnPlane);
        }
        return segmentDistances;
    }

    /**
     * Computes the 3D distances l1_pq and l2_pq between the computation point P and the line
     * segment endpoints of each polyhedral segment for one plane.
     * Computes the 1D distances s1_pq and s2_pq between orthogonal projection of P on the line
     * segment P''_pq and the line segment endpoints for each polyhedral segment for one plane
     * @param segmentVectorsForPlane the segment vectors G_pq for plane p
     * @param orthogonalProjectionPointsOnSegmentForPlane the orthogonal projection Points P'' for plane p
     * @param face the vertices of plane p
     * @return distances l1_pq and l2_pq and s1_pq and s2_pq foreach segment q of plane p
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION std::array<DistanceTemplate<FloatType>, 3> distancesToSegmentEndpoints(
            const Vector3Triplet<FloatType> &segmentVectorsForPlane,
            const Vector3Triplet<FloatType> &orthogonalProjectionPointsOnSegmentForPlane,
            const Vector3Triplet<FloatType> &face) {
        using util::operator-;
        std::array<DistanceTemplate<FloatType>, 3> distancesForPlane{};

        //The 3D distance l_pq from P to a vertex is the same for the segment which ends at that vertex and
        //the one which starts at it, so the three norms are taken once instead of six times
        const Vector3<FloatType> vertexNorms{util::euclideanNorm(face[0]), util::euclideanNorm(face[1]),
                                             util::euclideanNorm(face[2])};

        for (size_t j = 0; j < 3; ++j) {
            DistanceTemplate<FloatType> distance{};
            //orthogonal projection point on segment P'' for plane p and segment q
            const Vector3<FloatType> &orthogonalProjectionPointsOnSegment = orthogonalProjectionPointsOnSegmentForPlane[j];

            //Calculate the 3D distances between P (0, 0, 0) and
            // the segment endpoints face[j] and face[(j + 1) % 3])
            distance.l1 = vertexNorms[j];
            distance.l2 = vertexNorms[(j + 1) % 3];
            //Calculate the 1D distances between P'' (every segment has its own) and
            // the segment endpoints face[j] and face[(j + 1) % 3])
            distance.s1 = util::euclideanNorm(orthogonalProjectionPointsOnSegment - face[j]);
            distance.s2 = util::euclideanNorm(orthogonalProjectionPointsOnSegment - face[(j + 1) % 3]);

            /*
             * Additional remark:
             * Details on these conditions are in the second paper referenced in the README.md (Tsoulis, 2021)
             * The numbering of these conditions is equal to the numbering scheme of the paper
             * Assign a sign to those magnitudes depending on the relative position of P'' to the two
             * segment endpoints
             */

            //4. Option: |s1 - l1| == 0 && |s2 - l2| == 0 Computation point P is located from the beginning on
            // the direction of a specific segment (P coincides with P' and P'')
            if (Kokkos::abs(distance.s1 - distance.l1) < EPSILON_ZERO<FloatType> &&
                Kokkos::abs(distance.s2 - distance.l2) < EPSILON_ZERO<FloatType>) {
                //4. Option - Case 2: P is located on the segment from its right side
                // s1 = -|s1|, s2 = -|s2|, l1 = -|l1|, l2 = -|l2|
                if (distance.s2 < distance.s1) {
                    distance.s1 *= -1.0;
                    distance.s2 *= -1.0;
                    distance.l1 *= -1.0;
                    distance.l2 *= -1.0;
                } else if (Kokkos::abs(distance.s2 - distance.s1) < EPSILON_ZERO<FloatType>) {
                    //4. Option - Case 1: P is located inside the segment (s2 == s1)
                    // s1 = -|s1|, s2 = |s2|, l1 = -|l1|, l2 = |l2|
                    distance.s1 *= -1.0;
                    distance.l1 *= -1.0;
                }
                //4. Option - Case 3: P is located on the segment from its left side
                // s1 = |s1|, s2 = |s2|, l1 = |l1|, l2 = |l2| --> Nothing to do!
            } else {
                const FloatType norm = util::euclideanNorm(segmentVectorsForPlane[j]);
                if (distance.s1 < norm && distance.s2 < norm) {
                    //1. Option: |s1| < |G_ij| && |s2| < |G_ij| Point P'' is situated inside the segment
                    // s1 = -|s1|, s2 = |s2|, l1 = |l1|, l2 = |l2|
                    distance.s1 *= -1.0;
                } else if (distance.s2 < distance.s1) {
                    //2. Option: |s2| < |s1| Point P'' is on the right side of the segment
                    // s1 = -|s1|, s2 = -|s2|, l1 = |l1|, l2 = |l2|
                    distance.s1 *= -1.0;
                    distance.s2 *= -1.0;
                }
                //3. Option: |s1| < |s2| Point P'' is on the left side of the segment
                // s1 = |s1|, s2 = |s2|, l1 = |l1|, l2 = |l2| --> Nothing to do!
            }
            distancesForPlane[j] = distance;
        }
        return distancesForPlane;
    }

    /**
     * Calculates the Transcendental Expressions LN_pq and AN_pq for every line segment of the polyhedron for
     * a given plane p.
     * LN_pq is calculated according to (14) using the natural logarithm and AN_pq is calculated according
     * to (15) using the arctan.
     * @param distancesForPlane the distances l1, l2, s1, s2 foreach segment q of plane p
     * @param planeDistance the plane distance h_p for plane p
     * @param segmentDistancesForPlane the segment distance h_pq for segment q of plane p
     * @param segmentNormalOrientationsForPlane the segment normal orientations n_pq for a plane p
     * @param projectionPointVertexNorms the norms of P' and each vertex of plane p
     * @return LN_pq and AN_pq foreach segment q of plane p
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION std::array<TranscendentalExpressionTemplate<FloatType>, 3> computeTranscendentalExpressions(
            const std::array<DistanceTemplate<FloatType>, 3> &distancesForPlane, FloatType planeDistance,
            const Vector3<FloatType> &segmentDistancesForPlane,
            const Vector3<FloatType> &segmentNormalOrientationsForPlane,
            const Vector3<FloatType> &projectionPointVertexNorms) {
        std::array<TranscendentalExpressionTemplate<FloatType>, 3> transcendentalExpressionsForPlane{};

        for (size_t j = 0; j < 3; ++j) {
            //distances l1, l2, s1, s1 for this segment q of plane p
            const DistanceTemplate<FloatType> &distance = distancesForPlane[j];
            //segment distance h_pq for this segment q of plane p
            const FloatType segmentDistance = segmentDistancesForPlane[j];
            //segment normal orientation sigma_pq for this segment q of plane p
            const FloatType segmentNormalOrientation = segmentNormalOrientationsForPlane[j];

            //Result for this segment
            TranscendentalExpressionTemplate<FloatType> transcendentalExpressionPerSegment{};

            //Computation of the norm of P' and segment endpoints
            // If the one of the norms == 0 then P' lies on the corresponding vertex and coincides with P''
            const FloatType r1Norm = projectionPointVertexNorms[(j + 1) % 3];
            const FloatType r2Norm = projectionPointVertexNorms[j];

            //Compute LN_pq according to (14)
            // If sigma_pq == 0 && either of the distances of P' to the two segment endpoints == 0 OR
            // the 1D and 3D distances are smaller than some EPSILON
            // then LN_pq can be set to zero
            if ((segmentNormalOrientation == 0.0 &&
                 (r1Norm < EPSILON_ZERO<FloatType> || r2Norm < EPSILON_ZERO<FloatType>)) ||
                (Kokkos::abs(distance.s1 + distance.s2) < EPSILON_ZERO<FloatType> &&
                 Kokkos::abs(distance.l1 + distance.l2) < EPSILON_ZERO<FloatType>)) {
                transcendentalExpressionPerSegment.ln = 0.0;
            } else {
                //Implementation of
                // log((s2_pq + l2_pq) / (s1_pq + l1_pq))
                transcendentalExpressionPerSegment.ln =
                        Kokkos::log((distance.s2 + distance.l2) / (distance.s1 + distance.l1));
            }

            //Compute AN_pq according to (15)
            // If h_p == 0 or h_pq == 0 then AN_pq is zero, too (distances are always positive!)
            if (planeDistance < EPSILON_ZERO<FloatType> || segmentDistance < EPSILON_ZERO<FloatType>) {
                transcendentalExpressionPerSegment.an = 0.0;
            } else {
                //Implementation of:
                // atan(h_p * s2_pq / h_pq * l2_pq) - atan(h_p * s1_pq / h_pq * l1_pq)
                transcendentalExpressionPerSegment.an =
                        Kokkos::atan((planeDistance * distance.s2) / (segmentDistance * distance.l2)) -
                        Kokkos::atan((planeDistance * distance.s1) / (segmentDistance * distance.l1));
            }

            transcendentalExpressionsForPlane[j] = transcendentalExpressionPerSegment;
        }
        return transcendentalExpressionsForPlane;
    }

    /**
     * Calculates the singularities (correction) terms according to the Flow text for a given plane p.
     * @param segmentVectorsForPlane the segment vectors for a given plane
     * @param segmentNormalOrientationForPlane the segment orientation sigma_pq
     * @param projectionPointVertexNorms the norms of the projection point P' and the plane's vertices
     * @param planeUnitNormal the plane unit normal N_p
     * @param planeDistance the plane distance h_p
     * @param planeNormalOrientation the plane normal orientation sigma_p
     * @return the singularities for a plane p
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION SingularityTerms<FloatType> computeSingularityTerms(
            const Vector3Triplet<FloatType> &segmentVectorsForPlane,
            const Vector3<FloatType> &segmentNormalOrientationForPlane,
            const Vector3<FloatType> &projectionPointVertexNorms, const Vector3<FloatType> &planeUnitNormal,
            FloatType planeDistance, FloatType planeNormalOrientation) {
        using util::operator*;
        //1. Case: If all sigma_pq for a given plane p are 1.0 then P' lies inside the plane S_p
        bool allInside = true;
        for (size_t j = 0; j < 3; ++j) {
            allInside &= segmentNormalOrientationForPlane[j] == 1.0;
        }
        if (allInside) {
            const FloatType factor = static_cast<FloatType>(-1.0 * util::PI2);
            //sing alpha = -2pi*h_p and sing beta = -2pi*sigma_p*N_p
            return {factor * planeDistance, planeUnitNormal * (factor * planeNormalOrientation)};
        }

        //2. Case: If sigma_pq == 0 AND norm(P' - v1) < norm(G_ij) && norm(P' - v2) < norm(G_ij) with G_ij
        // as the vector of v1 and v2
        // then P' is located on one line segment G_p of plane p, but not on any of its vertices
        bool anyOnLine = false;
        for (size_t j = 0; j < 3; ++j) {
            //segmentNormalOrientation != 0.0
            if (Kokkos::abs(segmentNormalOrientationForPlane[j]) > EPSILON_ZERO<FloatType>) {
                continue;
            }
            const FloatType segmentVectorNorm = util::euclideanNorm(segmentVectorsForPlane[j]);
            anyOnLine |= projectionPointVertexNorms[(j + 1) % 3] < segmentVectorNorm &&
                         projectionPointVertexNorms[j] < segmentVectorNorm &&
                         projectionPointVertexNorms[(j + 1) % 3] >= EPSILON_ZERO<FloatType> &&
                         projectionPointVertexNorms[j] >= EPSILON_ZERO<FloatType>;
        }
        if (anyOnLine) {
            const FloatType factor = static_cast<FloatType>(-1.0 * util::PI);
            //sing alpha = -pi*h_p and sing beta = -pi*sigma_p*N_p
            return {factor * planeDistance, planeUnitNormal * (factor * planeNormalOrientation)};
        }

        //3. Case If sigma_pq == 0 AND norm(P' - v1) < 0 || norm(P' - v2) < 0
        // then P' is located at one of G_p's vertices
        for (size_t j = 0; j < 3; ++j) {
            //segmentNormalOrientation != 0.0
            if (Kokkos::abs(segmentNormalOrientationForPlane[j]) > EPSILON_ZERO<FloatType>) {
                continue;
            }
            const FloatType r1Norm = projectionPointVertexNorms[(j + 1) % 3];
            const FloatType r2Norm = projectionPointVertexNorms[j];
            //r1Norm == 0.0 || r2Norm == 0.0
            if (!(r1Norm < EPSILON_ZERO<FloatType> || r2Norm < EPSILON_ZERO<FloatType>)) {
                continue;
            }
            //Two segment vectors G_1 and G_2 of this plane
            const Vector3<FloatType> &g1 = r1Norm < EPSILON_ZERO<FloatType>
                                                   ? segmentVectorsForPlane[j]
                                                   : segmentVectorsForPlane[(j - 1 + 3) % 3];
            const Vector3<FloatType> &g2 = r1Norm < EPSILON_ZERO<FloatType>
                                                   ? segmentVectorsForPlane[(j + 1) % 3]
                                                   : segmentVectorsForPlane[j];
            // theta = arcos((G_2 * -G_1) / (|G_2| * |G_1|))
            const FloatType gdot = util::dot(g1 * static_cast<FloatType>(-1.0), g2);
            const FloatType theta = gdot == 0.0
                                            ? static_cast<FloatType>(util::PI_2)
                                            : Kokkos::acos(gdot / (util::euclideanNorm(g1) * util::euclideanNorm(g2)));
            //sing alpha = -theta*h_p and sing beta = -theta*sigma_p*N_p
            const FloatType factor = static_cast<FloatType>(-1.0) * theta;
            return {factor * planeDistance, planeUnitNormal * (factor * planeNormalOrientation)};
        }

        //4. Case Otherwise P' is located outside the plane S_p and then the singularity equals zero
        return {static_cast<FloatType>(0.0), Vector3<FloatType>{0.0, 0.0, 0.0}};
    }

    /**
     * Computes the L2 norms of the orthogonal projection point P' on a plane p with each vertex of that plane p.
     * The values are later used to determine if P' is situated at a vertex.
     * @param orthogonalProjectionPointOnPlane the orthogonal projection point P'
     * @param face the vertices of plane p
     * @return the norms of p and each vertex
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION Vector3<FloatType> computeNormsOfProjectionPointAndVertices(
            const Vector3<FloatType> &orthogonalProjectionPointOnPlane, const Vector3Triplet<FloatType> &face) {
        using util::operator-;
        return {util::euclideanNorm(orthogonalProjectionPointOnPlane - face[0]),
                util::euclideanNorm(orthogonalProjectionPointOnPlane - face[1]),
                util::euclideanNorm(orthogonalProjectionPointOnPlane - face[2])};
    }

    /**
     * The contribution of one polyhedral face to the gravity model's result at one computation point.
     *
     * This is the value type reduced over all faces of the polyhedron. It intentionally holds the raw sums,
     * i.e. the prefix (gravitational constant, density, normal orientation, and mesh unit) is only applied
     * once on the host after the reduction has finished.
     *
     * @tparam FloatType the floating point precision of the evaluation
     */
    template<typename FloatType>
    struct FaceContribution {
        /** The gravitational potential V, Tsoulis' Equation (11) without the prefix */
        FloatType potential{};
        /** The first order derivatives Vx, Vy, Vz, Tsoulis' Equation (12) without the prefix */
        Vector3<FloatType> acceleration{};
        /** The second order derivatives Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, Tsoulis' Equation (13) without the prefix */
        Vector6<FloatType> gradiometricTensor{};
        /**
         * How many faces were evaluated with magnitudes so far apart that floating point absorption is
         * expected. The host turns a non-zero count into a warning; the face itself cannot log from a kernel.
         */
        int numericallyCriticalFaces{};

        /**
         * Adds another face's contribution to this one. This is the reduction's join operation.
         * @param rhs the other contribution
         * @return this contribution
         */
        KOKKOS_INLINE_FUNCTION FaceContribution &operator+=(const FaceContribution &rhs) {
            potential += rhs.potential;
            for (size_t index = 0; index < 3; ++index) {
                acceleration[index] += rhs.acceleration[index];
            }
            for (size_t index = 0; index < 6; ++index) {
                gradiometricTensor[index] += rhs.gradiometricTensor[index];
            }
            numericallyCriticalFaces += rhs.numericallyCriticalFaces;
            return *this;
        }
    };

    /**
     * Determines whether two values are so far apart in magnitude that adding them absorbs the smaller one.
     *
     * This is the device-capable counterpart of {@link util::isCriticalDifference}: it compares the binary
     * exponents via {@code logb} instead of {@code std::frexp}, which does not exist inside a kernel.
     *
     * @tparam FloatType the floating point precision of the evaluation
     * @param first the first number
     * @param second the second number
     * @return true if the difference in magnitude is too big
     */
    template<typename FloatType>
    KOKKOS_INLINE_FUNCTION bool isCriticalDifference(const FloatType first, const FloatType second) {
        // 50 is the (log2) exponent of the floating point (1 / 1e-15)
        constexpr int maxExponentDifference = 50;
        // logb(0) is negative infinity, but a zero can never absorb anything
        if (first == 0.0 || second == 0.0) {
            return false;
        }
        const int firstExponent = static_cast<int>(Kokkos::logb(first));
        const int secondExponent = static_cast<int>(Kokkos::logb(second));
        return Kokkos::abs(firstExponent - secondExponent) > maxExponentDifference;
    }

    /**
     * Evaluates Tsoulis' polyhedral gravity model for a single face of the polyhedron at a single
     * computation point P.
     *
     * This is the kernel body: it is called once per face by the range policy of the single point evaluation
     * and once per (point, face) pair by the team policy of the multi point evaluation. It runs unchanged on
     * every backend, i.e. serially on the host, on all CPU threads, and on the GPU.
     *
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the memory space of the mesh
     * @param mesh the device-resident polyhedron together with its caches
     * @param faceIndex the index of the face to evaluate
     * @param computationPoint the computation point P
     * @return this face's contribution to the potential, the acceleration, and the gradiometric tensor
     */
    template<typename FloatType, typename MemorySpace>
    KOKKOS_INLINE_FUNCTION FaceContribution<FloatType> evaluateFace(
            const kokkos::GravitationalMeshView<FloatType, MemorySpace> &mesh, const size_t faceIndex,
            const Vector3<FloatType> &computationPoint) {
        using util::operator+;
        using util::operator*;

        //1. Step: Compute the ingredients for the current plane
        //1-01 to 1-03 Step: The segment vectors, plane unit normals, and segment unit normals only depend on
        // the polyhedron and were therefore already computed once by the initialization kernel
        const Vector3Triplet<FloatType> face = mesh.resolveFace(faceIndex, computationPoint);
        const Vector3Triplet<FloatType> segmentVectors = mesh.getSegmentVectors(faceIndex);
        const Vector3<FloatType> planeUnitNormal = mesh.getPlaneUnitNormal(faceIndex);
        const Vector3Triplet<FloatType> segmentUnitNormals = mesh.getSegmentUnitNormals(faceIndex);

        //1-04 to 1-07 Step: The plane normal orientation sigma_p, the plane distance h_p, and the position of
        // P' all follow from projecting the face onto its own plane unit normal, which is already cached. The
        // Hessian form of the plane is therefore never built here: it would recompute, once per computation
        // point, a cross product which only depends on the polyhedron.
        const FloatType planeProjection = util::dot(planeUnitNormal, face[0]);
        //1-04 Step: Compute Plane Normal Orientation sigma_p (direction of N_p in relation to P)
        const auto planeNormalOrientation =
                static_cast<FloatType>(util::sgn(planeProjection, EPSILON_ZERO<FloatType>));
        //1-05 to 1-06 Step: Compute distance h_p between P and P'
        const FloatType planeDistance = Kokkos::abs(planeProjection);
        //1-07 Step: Compute the actual position of P' (projection of P on the plane)
        const Vector3<FloatType> orthogonalProjectionPointOnPlane =
                projectPointOrthogonallyOntoPlane(planeUnitNormal, planeDistance, -planeProjection);
        //1-08 Step: Compute the segment normal orientation sigma_pq (direction of n_pq in relation to P')
        const Vector3<FloatType> segmentNormalOrientations =
                computeUnitNormalOfSegmentsDirections(face, orthogonalProjectionPointOnPlane, segmentUnitNormals);
        //1-09 Step: Compute the orthogonal projection point P'' of P' on each segment
        const Vector3Triplet<FloatType> orthogonalProjectionPointsOnSegmentsForPlane =
                projectPointOrthogonallyOntoSegments(orthogonalProjectionPointOnPlane, segmentNormalOrientations, face);
        //1-10 Step: Compute the segment distances h_pq between P'' and P'
        const Vector3<FloatType> segmentDistances = distancesBetweenProjectionPoints(
                orthogonalProjectionPointOnPlane, orthogonalProjectionPointsOnSegmentsForPlane);
        //1-11 Step: Compute the 3D distances l1, l2 (between P and vertices)
        // and 1D distances s1, s2 (between P'' and vertices)
        const std::array<DistanceTemplate<FloatType>, 3> distances = distancesToSegmentEndpoints(
                segmentVectors, orthogonalProjectionPointsOnSegmentsForPlane, face);
        //1-12 Step: Compute the euclidian Norms of the vectors consisting of P and the vertices
        // they are later used for determining the position of P in relation to the plane
        const Vector3<FloatType> projectionPointVertexNorms =
                computeNormsOfProjectionPointAndVertices(orthogonalProjectionPointOnPlane, face);
        //1-13 Step: Compute the transcendental Expressions LN_pq and AN_pq
        const std::array<TranscendentalExpressionTemplate<FloatType>, 3> transcendentalExpressions =
                computeTranscendentalExpressions(distances, planeDistance, segmentDistances,
                                                 segmentNormalOrientations, projectionPointVertexNorms);
        //1-14 Step: Compute the singularities sing A and sing B if P' is located in the plane,
        // on any vertex, or on one segment (G_pq)
        const SingularityTerms<FloatType> singularities =
                computeSingularityTerms(segmentVectors, segmentNormalOrientations, projectionPointVertexNorms,
                                        planeUnitNormal, planeDistance, planeNormalOrientation);

        //2. Step: Compute Sum 1 used for potential and acceleration (first derivative)
        // sum over: sigma_pq * h_pq * LN_pq
        // --> Equation 11/12 the first summation in the brackets
        FloatType sum1PotentialAcceleration{0.0};
        for (size_t j = 0; j < 3; ++j) {
            sum1PotentialAcceleration += segmentNormalOrientations[j] * segmentDistances[j] * transcendentalExpressions[j].ln;
        }

        //3. Step: Compute Sum 1 used for the gradiometric tensor (second derivative)
        // sum over: n_pq * LN_pq
        // --> Equation 13 the first summation in the brackets
        Vector3<FloatType> sum1Tensor{0.0, 0.0, 0.0};
        for (size_t j = 0; j < 3; ++j) {
            sum1Tensor = sum1Tensor + (segmentUnitNormals[j] * transcendentalExpressions[j].ln);
        }

        //4. Step: Compute Sum 2 which is the same for every result parameter
        // sum over: sigma_pq * AN_pq
        // --> Equation 11/12/13 the second summation in the brackets
        FloatType sum2{0.0};
        for (size_t j = 0; j < 3; ++j) {
            sum2 += segmentNormalOrientations[j] * transcendentalExpressions[j].an;
        }

        //5. Step: Sum for potential and acceleration
        // consisting of: sum1 + h_p * sum2 + sing A
        // --> Equation 11/12 the total sum of the brackets
        const FloatType planeSumPotentialAcceleration =
                sum1PotentialAcceleration + planeDistance * sum2 + singularities.alpha;

        //6. Step: Sum for tensor
        // consisting of: sum1 + sigma_p * N_p * sum2 + sing B
        // --> Equation 13 the total sum of the brackets
        const Vector3<FloatType> subSum =
                (sum1Tensor + (planeUnitNormal * (planeNormalOrientation * sum2))) + singularities.beta;
        // first component: trivial case Vxx, Vyy, Vzz --> just N_p * subSum
        // 00, 11, 22 --> xx, yy, zz with x as 0, y as 1, z as 2
        const Vector3<FloatType> first = planeUnitNormal * subSum;
        // second component: reordering required to build Vxy, Vxz, Vyz
        // 01, 02, 12 --> xy, xz, yz with x as 0, y as 1, z as 2
        const Vector3<FloatType> reorderedNp = {planeUnitNormal[0], planeUnitNormal[0], planeUnitNormal[1]};
        const Vector3<FloatType> reorderedSubSum = {subSum[1], subSum[2], subSum[2]};
        const Vector3<FloatType> second = reorderedNp * reorderedSubSum;

        //7. Step: Assemble this face's contribution
        // Equation (11): sigma_p * h_p * sum
        // Equation (12): N_p * sum
        // Equation (13): already done above, just concat the two components for later summation
        // The multiplication planeDistance * sum2 is not the root cause of a numerical problem, but both
        // numbers are good indicators for the magnitudes appearing during the calculation: planeDistance gets
        // very big when far away, sum2 remains independently very small.
        return {planeNormalOrientation * planeDistance * planeSumPotentialAcceleration,
                planeUnitNormal * planeSumPotentialAcceleration,
                util::concat(first, second),
                isCriticalDifference(planeDistance, sum2) ? 1 : 0};
    }

}// namespace polyhedralGravity::GravityModel::detail

namespace Kokkos {

    /**
     * Teaches Kokkos the neutral element of the FaceContribution reduction, which is an all-zero contribution.
     * @tparam FloatType the floating point precision of the evaluation
     */
    template<typename FloatType>
    struct reduction_identity<polyhedralGravity::GravityModel::detail::FaceContribution<FloatType>> {
        KOKKOS_FORCEINLINE_FUNCTION static polyhedralGravity::GravityModel::detail::FaceContribution<FloatType> sum() {
            return polyhedralGravity::GravityModel::detail::FaceContribution<FloatType>{};
        }
    };

}// namespace Kokkos
