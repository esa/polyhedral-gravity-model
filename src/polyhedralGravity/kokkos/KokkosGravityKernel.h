#pragma once

#include <Kokkos_Core.hpp>

#include "polyhedralGravity/model/GravityModelDetail.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/util/UtilityContainer.h"

namespace polyhedralGravity::kokkos {

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
     * The polyhedron as it lives in the memory of one execution space.
     *
     * Besides the mesh itself, this holds the three caches which only depend on the polyhedron and not on the
     * computation point. They are filled once by {@link KokkosEvaluation}'s initialization kernel and reused
     * for every subsequent evaluation.
     *
     * @tparam FloatType the floating point precision of the evaluation
     * @tparam MemorySpace the Kokkos memory space the views are allocated in
     */
    template<typename FloatType, typename MemorySpace>
    struct MeshViews {
        /** The polyhedron's vertices as cartesian coordinates */
        Kokkos::View<Vector3<FloatType> *, MemorySpace> vertices;
        /** The polyhedron's triangular faces, each referencing three vertices by index */
        Kokkos::View<IndexArray3 *, MemorySpace> faces;
        /** The segment vectors G_pq foreach face */
        Kokkos::View<Vector3Triplet<FloatType> *, MemorySpace> segmentVectors;
        /** The plane unit normals N_p foreach face */
        Kokkos::View<Vector3<FloatType> *, MemorySpace> planeUnitNormals;
        /** The segment unit normals n_pq foreach face */
        Kokkos::View<Vector3Triplet<FloatType> *, MemorySpace> segmentUnitNormals;

        /**
         * Resolves the face at the given index into cartesian coordinates relative to a computation point,
         * i.e. it re-locates the computation point into the origin as Tsoulis' equations require.
         * @param faceIndex the index of the face
         * @param computationPoint the computation point P
         * @return the face's three vertices, each shifted by -P
         */
        KOKKOS_INLINE_FUNCTION Vector3Triplet<FloatType> resolveFace(
                const size_t faceIndex, const Vector3<FloatType> &computationPoint) const {
            using util::operator-;
            const IndexArray3 &face = faces(faceIndex);
            return {vertices(face[0]) - computationPoint,
                    vertices(face[1]) - computationPoint,
                    vertices(face[2]) - computationPoint};
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
     * @param mesh the device-resident polyhedron
     * @param faceIndex the index of the face to evaluate
     * @param computationPoint the computation point P
     * @return this face's contribution to the potential, the acceleration, and the gradiometric tensor
     */
    template<typename FloatType, typename MemorySpace>
    KOKKOS_INLINE_FUNCTION FaceContribution<FloatType> evaluateFace(
            const MeshViews<FloatType, MemorySpace> &mesh, const size_t faceIndex,
            const Vector3<FloatType> &computationPoint) {
        using namespace GravityModel::detail;
        using util::operator+;
        using util::operator*;

        //1. Step: Compute the ingredients for the current plane
        //1-01 to 1-03 Step: The segment vectors, plane unit normals, and segment unit normals only depend on
        // the polyhedron and were therefore already computed once by the initialization kernel
        const Vector3Triplet<FloatType> face = mesh.resolveFace(faceIndex, computationPoint);
        const Vector3Triplet<FloatType> &segmentVectors = mesh.segmentVectors(faceIndex);
        const Vector3<FloatType> &planeUnitNormal = mesh.planeUnitNormals(faceIndex);
        const Vector3Triplet<FloatType> &segmentUnitNormals = mesh.segmentUnitNormals(faceIndex);

        //1-04 Step: Compute Plane Normal Orientation sigma_p (direction of N_p in relation to P)
        const FloatType planeNormalOrientation = computeUnitNormalOfPlaneDirection(planeUnitNormal, face[0]);
        //1-05 Step: Compute Hessian Normal Plane Representation
        const HessianPlaneTemplate<FloatType> hessianPlane = computeHessianPlane(face[0], face[1], face[2]);
        //1-06 Step: Compute distance h_p between P and P'
        const FloatType planeDistance = distanceBetweenOriginAndPlane(hessianPlane);
        //1-07 Step: Compute the actual position of P' (projection of P on the plane)
        const Vector3<FloatType> orthogonalProjectionPointOnPlane =
                projectPointOrthogonallyOntoPlane(planeUnitNormal, planeDistance, hessianPlane);
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

}// namespace polyhedralGravity::kokkos

namespace Kokkos {

    /**
     * Teaches Kokkos the neutral element of the FaceContribution reduction, which is an all-zero contribution.
     * @tparam FloatType the floating point precision of the evaluation
     */
    template<typename FloatType>
    struct reduction_identity<polyhedralGravity::kokkos::FaceContribution<FloatType>> {
        KOKKOS_FORCEINLINE_FUNCTION static polyhedralGravity::kokkos::FaceContribution<FloatType> sum() {
            return polyhedralGravity::kokkos::FaceContribution<FloatType>{};
        }
    };

}// namespace Kokkos
