#include "GravityModel.h"

namespace polyhedralGravity {

    GravityModelResult GravityModel::evaluate(
            const Polyhedron &polyhedron, double density, const Array3 &computationPoint) {
        using namespace util;
        using namespace detail;
        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                            "Evaluation for computation point P = [{}, {}, {}] started, given density = {} kg/m^3",
                            computationPoint[0], computationPoint[1], computationPoint[2], density);
        /*
         * Calculate V and Vx, Vy, Vz and Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
         */
        auto polyhedronIterator = transformPolyhedron(polyhedron, computationPoint);

        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                            "Starting to iterate over the planes...");
        GravityModelResult result{};
        result = thrust::transform_reduce(
                thrust::device,
                polyhedronIterator.first,
                polyhedronIterator.second,
                [](const Array3Triplet &face) {
                    using namespace util;
                    using namespace detail;
                    SPDLOG_LOGGER_TRACE(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                                        "Evaluating the plane with vertices: v1 = [{}, {}, {}], v2 = [{}, {}, {}], "
                                        "v3 = [{}, {}, {}]",
                                        face[0][0], face[0][1], face[0][2],
                                        face[1][0], face[1][1], face[1][2],
                                        face[2][0], face[2][1], face[2][2]);
                    //1. Step: Compute ingredients for current plane
                    //1-01 Step: Compute Segment Vectors G_pq which describe each one the edge between two vertices
                    Array3Triplet segmentVectors = buildVectorsOfSegments(face[0], face[1], face[2]);
                    //1-02 Step: Compute the Plane Unit Normals N_p (pointing outside the polyhedron)
                    Array3 planeUnitNormal = buildUnitNormalOfPlane(segmentVectors[0], segmentVectors[1]);
                    //1-03 Step: Compute Segment Unit Normals n_pq (normal pointing away from each segment)
                    Array3Triplet segmentUnitNormals = buildUnitNormalOfSegments(segmentVectors, planeUnitNormal);
                    //1-04 Step: Compute Plane Normal Orientation sigma_p (direction of N_p in relation to P)
                    double planeNormalOrientation = computeUnitNormalOfPlaneDirection(planeUnitNormal, face[0]);
                    //1-05 Step: Compute Hessian Normal Plane Representation
                    HessianPlane hessianPlane = computeHessianPlane(face[0], face[1], face[2]);
                    //1-06 Step: Compute distance h_p between P and P'
                    double planeDistance = distanceBetweenOriginAndPlane(hessianPlane);
                    //1-07 Step: Compute the actual position of P' (projection of P on the plane)
                    Array3 orthogonalProjectionPointOnPlane = projectPointOrthogonallyOntoPlane(
                            planeUnitNormal,
                            planeDistance,
                            hessianPlane);
                    //1-08 Step: Compute the segment normal orientation sigma_pq (direction of n_pq in relation to P')
                    Array3 segmentNormalOrientations = computeUnitNormalOfSegmentsDirections(
                            face,
                            orthogonalProjectionPointOnPlane,
                            segmentUnitNormals);
                    //1-09 Step: Compute the orthogonal projection point P'' of P' on each segment
                    Array3Triplet orthogonalProjectionPointsOnSegmentsForPlane = projectPointOrthogonallyOntoSegments(
                            orthogonalProjectionPointOnPlane,
                            segmentNormalOrientations,
                            face);
                    //1-10 Step: Compute the segment distances h_pq between P'' and P'
                    Array3 segmentDistances = distancesBetweenProjectionPoints(
                            orthogonalProjectionPointOnPlane,
                            orthogonalProjectionPointsOnSegmentsForPlane);
                    //1-11 Step: Compute the 3D distances l1, l2 (between P and vertices)
                    // and 1D distances s1, s2 (between P'' and vertices)
                    std::array<Distance, 3> distances = distancesToSegmentEndpoints(
                            segmentVectors,
                            orthogonalProjectionPointsOnSegmentsForPlane,
                            face);
                    //1-12 Step: Compute the euclidian Norms of the vectors consisting of P and the vertices
                    // they are later used for determining the position of P in relation to the plane
                    Array3 projectionPointVertexNorms = computeNormsOfProjectionPointAndVertices(
                            orthogonalProjectionPointOnPlane,
                            face);
                    //1-13 Step: Compute the transcendental Expressions LN_pq and AN_pq
                    std::array<TranscendentalExpression, 3> transcendentalExpressions = computeTranscendentalExpressions(
                            distances,
                            planeDistance,
                            segmentDistances,
                            segmentNormalOrientations,
                            projectionPointVertexNorms);
                    //1-14 Step: Compute the singularities sing A and sing B if P' is located in the plane,
                    // on any vertex, or on one segment (G_pq)
                    std::pair<double, Array3> singularities = computeSingularityTerms(
                            segmentVectors,
                            segmentNormalOrientations,
                            projectionPointVertexNorms,
                            planeUnitNormal, planeDistance,
                            planeNormalOrientation);
                    //2. Step: Compute Sum 1 used for potential and acceleration (first derivative)
                    // sum over: sigma_pq * h_pq * LN_pq
                    // --> Equation 11/12 the first summation in the brackets
                    auto zipIteratorSum1PotentialAcceleration =
                            util::zipPair(segmentNormalOrientations, segmentDistances, transcendentalExpressions);
                    const double sum1PotentialAcceleration = std::accumulate(
                            zipIteratorSum1PotentialAcceleration.first,
                            zipIteratorSum1PotentialAcceleration.second,
                            0.0,
                            [](double acc, const auto &tuple) {
                                const double &segmentOrientation = thrust::get<0>(tuple);
                                const double &segmentDistance = thrust::get<1>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<2>(tuple);
                                return acc + segmentOrientation * segmentDistance * transcendentalExpressions.ln;
                            });

                    //3. Step: Compute Sum 1 used for the gradiometric tensor (second derivative)
                    // sum over: n_pq * LN_pq
                    // --> Equation 13 the first summation in the brackets
                    auto zipIteratorSum1Tensor = util::zipPair(segmentUnitNormals, transcendentalExpressions);
                    const Array3 sum1Tensor = std::accumulate(
                            zipIteratorSum1Tensor.first,
                            zipIteratorSum1Tensor.second,
                            Array3{0.0, 0.0, 0.0},
                            [](const Array3 &acc, const auto &tuple) {
                                const Array3 &segmentNormal = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                                return acc + (segmentNormal * transcendentalExpressions.ln);
                            });

                    //4. Step: Compute Sum 2 which is the same for every result parameter
                    // sum over: sigma_pq * AN_pq
                    // --> Equation 11/12/13 the second summation in the brackets
                    auto zipIteratorSum2 = util::zipPair(segmentNormalOrientations, transcendentalExpressions);
                    const double sum2 = std::accumulate(
                            zipIteratorSum2.first,
                            zipIteratorSum2.second,
                            0.0,
                            [](double acc, const auto &tuple) {
                                const double &segmentOrientation = thrust::get<0>(tuple);
                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                                return acc + segmentOrientation * transcendentalExpressions.an;
                            });

                    //5. Step: Sum for potential and acceleration
                    // consisting of: sum1 + h_p * sum2 + sing A
                    // --> Equation 11/12 the total sum of the brackets
                    const double planeSumPotentialAcceleration =
                            sum1PotentialAcceleration + planeDistance * sum2 + singularities.first;

                    //6. Step: Sum for tensor
                    // consisting of: sum1 + sigma_p * N_p * sum2 + sing B
                    // --> Equation 13 the total sum of the brackets
                    const Array3 subSum = (sum1Tensor + (planeUnitNormal * (planeNormalOrientation * sum2))) + singularities.second;
                    // first component: trivial case Vxx, Vyy, Vzz --> just N_p * subSum
                    // 00, 11, 22 --> xx, yy, zz with x as 0, y as 1, z as 2
                    const Array3 first = planeUnitNormal * subSum;
                    // second component: reordering required to build Vxy, Vxz, Vyz
                    // 01, 02, 12 --> xy, xz, yz with x as 0, y as 1, z as 2
                    const Array3 reorderedNp = {planeUnitNormal[0], planeUnitNormal[0], planeUnitNormal[1]};
                    const Array3 reorderedSubSum = {subSum[1], subSum[2], subSum[2]};
                    const Array3 second = reorderedNp * reorderedSubSum;

                    //7. Step: Multiply with prefix
                    // Equation (11): sigma_p * h_p * sum
                    // Equation (12): N_p * sum
                    // Equation (13): already done above, just concat the two components for later summation
                    return GravityModelResult{
                            planeNormalOrientation * planeDistance * planeSumPotentialAcceleration,
                            planeUnitNormal * planeSumPotentialAcceleration,
                            concat(first, second)
                    };
                },
                result,
                [](const GravityModelResult &a, const GravityModelResult &b) {
                    //8. Step: Accumulate the partial sums
                    return GravityModelResult{
                            a.gravitationalPotential + b.gravitationalPotential,
                            a.acceleration + b.acceleration,
                            a.gradiometricTensor + b.gradiometricTensor
                    };
                });

        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                            "Finished the sums. Applying final prefix and eliminating rounding errors.");

        //9. Step: Eliminate rounding errors in the result and set those tiny values to "really" zero
        result.eliminateRoundingErrors();

        //10. Step: Compute prefix consisting of GRAVITATIONAL_CONSTANT * density
        const double prefix = util::GRAVITATIONAL_CONSTANT * density;

        //11. Step: Final expressions after application of the prefix (and a division by 2 for the potential)
        result.gravitationalPotential = (result.gravitationalPotential * prefix) / 2.0;
        result.acceleration = result.acceleration * prefix;
        result.gradiometricTensor = result.gradiometricTensor * prefix;
        return result;
    }

    std::vector<GravityModelResult> GravityModel::evaluate(
            const Polyhedron &polyhedron, double density, const std::vector<Array3> &computationPoints) {
        std::vector<GravityModelResult> result{computationPoints.size()};
        thrust::transform(computationPoints.begin(), computationPoints.end(), result.begin(),
                          [&polyhedron, density](const Array3 &computationPoint) {
                              return evaluate(polyhedron, density, computationPoint);
                          });
        return result;
    }

    Array3Triplet GravityModel::detail::buildVectorsOfSegments(
            const Array3 &vertex0, const Array3 &vertex1, const Array3 &vertex2) {
        using util::operator-;
        //Calculate G_ij
        return {vertex1 - vertex0, vertex2 - vertex1, vertex0 - vertex2};
    }

    Array3 GravityModel::detail::buildUnitNormalOfPlane(const Array3 &segmentVector1, const Array3 &segmentVector2) {
        //Calculate N_i as (G_i1 * G_i2) / |G_i1 * G_i2| with * being the cross product
        return util::normal(segmentVector1, segmentVector2);
    }

    Array3Triplet GravityModel::detail::buildUnitNormalOfSegments(
            const Array3Triplet &segmentVectors, const Array3 &planeUnitNormal) {
        Array3Triplet segmentUnitNormal{};
        //Calculate n_ij as (G_ij * N_i) / |G_ig * N_i| with * being the cross product
        std::transform(segmentVectors.cbegin(), segmentVectors.end(), segmentUnitNormal.begin(),
                       [&planeUnitNormal](const Array3 &segmentVector) -> Array3 {
                           return util::normal(segmentVector, planeUnitNormal);
                       });
        return segmentUnitNormal;
    }

    double
    GravityModel::detail::computeUnitNormalOfPlaneDirection(const Array3 &planeUnitNormal, const Array3 &vertex0) {
        using namespace util;
        //Calculate N_i * -G_i1 where * is the dot product and then use the inverted sgn
        //We abstain on the double multiplication with -1 in the line above and beyond since two
        //times multiplying with -1 equals no change
        return sgn(dot(planeUnitNormal, vertex0), util::EPSILON);
    }

    HessianPlane GravityModel::detail::computeHessianPlane(const Array3 &p, const Array3 &q, const Array3 &r) {
        using namespace util;
        constexpr Array3 origin{0.0, 0.0, 0.0};
        const auto crossProduct = cross(p - q, p - r);
        const auto res = (origin - p) * crossProduct;
        const double d = res[0] + res[1] + res[2];

        return {crossProduct[0], crossProduct[1], crossProduct[2], d};
    }

    double GravityModel::detail::distanceBetweenOriginAndPlane(const HessianPlane &hessianPlane) {
        //Compute h_p as D/sqrt(A^2 + B^2 + C^2)
        return std::abs(hessianPlane.d / std::sqrt(
                hessianPlane.a * hessianPlane.a +
                hessianPlane.b * hessianPlane.b +
                hessianPlane.c * hessianPlane.c));
    }

    Array3 GravityModel::detail::projectPointOrthogonallyOntoPlane(
            const Array3 &planeUnitNormal,
            double planeDistance,
            const HessianPlane &hessianPlane) {
        using namespace util;
        //Calculate the projection point by (22) P'_ = N_i / norm(N_i) * h_i
        // norm(N_i) is always 1 since N_i is a "normed" vector --> we do not need this division
        Array3 orthogonalProjectionPoint = planeUnitNormal * planeDistance;

        //Calculate alpha, beta and gamma as D/A, D/B and D/C (Notice that we "forget" the minus before those
        // divisions. In consequence, the conditions for signs are reversed below!!!)
        // These values represent the intersections of each polygonal plane with the axes
        // Comparison x == 0.0 is ok, since we only want to avoid nan values
        Array3 intersections = {hessianPlane.a == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.a,
                                hessianPlane.b == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.b,
                                hessianPlane.c == 0.0 ? 0.0 : hessianPlane.d / hessianPlane.c};

        //Determine the signs of the coordinates of P' according to the intersection values alpha, beta, gamma
        // denoted as __ below, i.e. -alpha, -beta, -gamma denoted -__
        for (unsigned int index = 0; index < 3; ++index) {
            if (intersections[index] < 0) {
                //If -__ >= 0 --> __ < 0 then coordinates are positive, we calculate abs(orthogonalProjectionPoint[..])
                orthogonalProjectionPoint[index] = std::abs(orthogonalProjectionPoint[index]);
            } else {
                //The coordinates need to be controlled
                if (planeUnitNormal[index] > 0) {
                    //If -__ < 0 --> __ >= 0 then the coordinate is negative -orthogonalProjectionPoint[..]
                    orthogonalProjectionPoint[index] = -1.0 * orthogonalProjectionPoint[index];
                } else {
                    //Else the coordinate is positive orthogonalProjectionPoint[..]
                    orthogonalProjectionPoint[index] = orthogonalProjectionPoint[index];
                }
            }
        }
        return orthogonalProjectionPoint;
    }

    Array3 GravityModel::detail::computeUnitNormalOfSegmentsDirections(
            const Array3Triplet &vertices,
            const Array3 &projectionPointOnPlane,
            const Array3Triplet &segmentUnitNormalsForPlane) {
        using namespace util;
        std::array<double, 3> segmentNormalOrientations{};
        //Equation (23)
        //Calculate x_P' - x_ij^1 (x_P' is the projectionPoint and x_ij^1 is the first vertices of one segment,
        //i.e. the coordinates of the training-planes' nodes --> projectionPointOnPlane - vertex
        //Calculate n_ij * x_ij with * being the dot product and use the inverted sgn to determine the value of sigma_pq
        //--> sgn((dot(segmentNormalOrientation, projectionPointOnPlane - vertex)), util::EPSILON) * -1.0
        std::transform(segmentUnitNormalsForPlane.cbegin(), segmentUnitNormalsForPlane.cend(),
                       vertices.cbegin(), segmentNormalOrientations.begin(),
                       [&projectionPointOnPlane](const Array3 &segmentUnitNormal, const Array3 &vertex) {
                           using namespace util;
                           return sgn((dot(segmentUnitNormal, projectionPointOnPlane - vertex)), util::EPSILON) *
                                  -1.0;
                       });
        return segmentNormalOrientations;
    }

    Array3Triplet GravityModel::detail::projectPointOrthogonallyOntoSegments(
            const Array3 &projectionPointOnPlane,
            const Array3 &segmentNormalOrientations,
            const Array3Triplet &face) {
        auto counterJ = thrust::counting_iterator<unsigned int>(0);
        Array3Triplet orthogonalProjectionPointOnSegmentPerPlane{};

        //Running over the segments of this plane
        thrust::transform(segmentNormalOrientations.begin(), segmentNormalOrientations.end(),
                          counterJ, orthogonalProjectionPointOnSegmentPerPlane.begin(),
                          [&](const double &segmentNormalOrientation, const unsigned int j) {
                              //We actually only accept +0.0 or -0.0, so the equal comparison is ok
                              if (segmentNormalOrientation == 0.0) {
                                  //Geometrically trivial case, in neither of the half space --> already on segment
                                  return projectionPointOnPlane;
                              } else {
                                  //In one of the half space, evaluate the projection point P'' for the segment
                                  //with the endpoints v1 and v2
                                  const auto &vertex1 = face[j];
                                  const auto &vertex2 = face[(j + 1) % 3];
                                  return projectPointOrthogonallyOntoSegment(vertex1, vertex2,
                                                                             projectionPointOnPlane);
                              }
                          });
        return orthogonalProjectionPointOnSegmentPerPlane;
    }

    Array3 GravityModel::detail::projectPointOrthogonallyOntoSegment(const Array3 &vertex1, const Array3 &vertex2,
                                                                     const Array3 &orthogonalProjectionPointOnPlane) {
        using namespace util;
        //Preparing our the planes/ equations in matrix form
        const Array3 matrixRow1 = vertex2 - vertex1;
        const Array3 matrixRow2 = cross(vertex1 - orthogonalProjectionPointOnPlane, matrixRow1);
        const Array3 matrixRow3 = cross(matrixRow2, matrixRow1);
        const Array3 d = {dot(matrixRow1, orthogonalProjectionPointOnPlane),
                          dot(matrixRow2, orthogonalProjectionPointOnPlane),
                          dot(matrixRow3, vertex1)};
        Matrix<double, 3, 3> columnMatrix = transpose(Matrix<double, 3, 3>{matrixRow1, matrixRow2, matrixRow3});
        //Calculation and solving the equations of above
        const double determinant = det(columnMatrix);
        return Array3{
                det(Matrix<double, 3, 3>{d, columnMatrix[1], columnMatrix[2]}),
                det(Matrix<double, 3, 3>{columnMatrix[0], d, columnMatrix[2]}),
                det(Matrix<double, 3, 3>{columnMatrix[0], columnMatrix[1], d})
        } / determinant;
    }

    Array3 GravityModel::detail::distancesBetweenProjectionPoints(const Array3 &orthogonalProjectionPointOnPlane,
                                                                  const Array3Triplet &orthogonalProjectionPointOnSegments) {
        std::array<double, 3> segmentDistances{};
        //The inner loop with the running j --> iterating over the segments
        //Using the values P'_i and P''_ij for the calculation of the distance
        std::transform(orthogonalProjectionPointOnSegments.cbegin(), orthogonalProjectionPointOnSegments.cend(),
                       segmentDistances.begin(), [&](const Array3 &orthogonalProjectionPointOnSegment) {
                    using namespace util;
                    return euclideanNorm(orthogonalProjectionPointOnSegment - orthogonalProjectionPointOnPlane);
                });
        return segmentDistances;
    }

    std::array<Distance, 3> GravityModel::detail::distancesToSegmentEndpoints(
            const Array3Triplet &segmentVectorsForPlane,
            const Array3Triplet &orthogonalProjectionPointsOnSegmentForPlane,
            const Array3Triplet &face) {
        std::array<Distance, 3> distancesForPlane{};
        auto counter = thrust::counting_iterator<unsigned int>(0);
        auto zip = util::zipPair(segmentVectorsForPlane, orthogonalProjectionPointsOnSegmentForPlane);

        thrust::transform(zip.first, zip.second, counter,
                          distancesForPlane.begin(), [&face](const auto &tuple, unsigned int j) {
                    using namespace util;
                    Distance distance{};
                    //segment vector G_pq
                    const Array3 &segmentVector = thrust::get<0>(tuple);
                    //orthogonal projection point on segment P'' for plane p and segment q
                    const Array3 &orthogonalProjectionPointsOnSegment = thrust::get<1>(tuple);

                    //Calculate the 3D distances between P (0, 0, 0) and
                    // the segment endpoints face[j] and face[(j + 1) % 3])
                    distance.l1 = euclideanNorm(face[j]);
                    distance.l2 = euclideanNorm(face[(j + 1) % 3]);
                    //Calculate the 1D distances between P'' (every segment has its own) and
                    // the segment endpoints face[j] and face[(j + 1) % 3])
                    distance.s1 = euclideanNorm(orthogonalProjectionPointsOnSegment - face[j]);
                    distance.s2 = euclideanNorm(orthogonalProjectionPointsOnSegment - face[(j + 1) % 3]);

                    /*
                     * Additional remark:
                     * Details on these conditions are in the second paper referenced in the README.md (Tsoulis, 2021)
                     * The numbering of these conditions is equal to the numbering scheme of the paper
                     * Assign a sign to those magnitudes depending on the relative position of P'' to the two
                     * segment endpoints
                     */

                    //4. Option: |s1 - l1| == 0 && |s2 - l2| == 0 Computation point P is located from the beginning on
                    // the direction of a specific segment (P coincides with P' and P'')
                    if (std::abs(distance.s1 - distance.l1) < EPSILON &&
                        std::abs(distance.s2 - distance.l2) < EPSILON) {
                        //4. Option - Case 2: P is located on the segment from its right side
                        // s1 = -|s1|, s2 = -|s2|, l1 = -|l1|, l2 = -|l2|
                        if (distance.s2 < distance.s1) {
                            distance.s1 *= -1.0;
                            distance.s2 *= -1.0;
                            distance.l1 *= -1.0;
                            distance.l2 *= -1.0;
                            return distance;
                        } else if (std::abs(distance.s2 - distance.s1) < EPSILON) {
                            //4. Option - Case 1: P is located inside the segment (s2 == s1)
                            // s1 = -|s1|, s2 = |s2|, l1 = -|l1|, l2 = |l2|
                            distance.s1 *= -1.0;
                            distance.l1 *= -1.0;
                            return distance;
                        }
                        //4. Option - Case 3: P is located on the segment from its left side
                        // s1 = |s1|, s2 = |s2|, l1 = |l1|, l2 = |l2| --> Nothing to do!
                    } else {
                        const double norm = euclideanNorm(segmentVector);
                        if (distance.s1 < norm && distance.s2 < norm) {
                            //1. Option: |s1| < |G_ij| && |s2| < |G_ij| Point P'' is situated inside the segment
                            // s1 = -|s1|, s2 = |s2|, l1 = |l1|, l2 = |l2|
                            distance.s1 *= -1.0;
                            return distance;
                        } else if (distance.s2 < distance.s1) {
                            //2. Option: |s2| < |s1| Point P'' is on the right side of the segment
                            // s1 = -|s1|, s2 = -|s2|, l1 = |l1|, l2 = |l2|
                            distance.s1 *= -1.0;
                            distance.s2 *= -1.0;
                            return distance;
                        }
                        //3. Option: |s1| < |s2| Point P'' is on the left side of the segment
                        // s1 = |s1|, s2 = |s2|, l1 = |l1|, l2 = |l2| --> Nothing to do!
                    }
                    return distance;
                });
        return distancesForPlane;
    }

    std::array<TranscendentalExpression, 3>
    GravityModel::detail::computeTranscendentalExpressions(
            const std::array<Distance, 3> &distancesForPlane,
            double planeDistance, const Array3 &segmentDistancesForPlane,
            const Array3 &segmentNormalOrientationsForPlane,
            const Array3 &projectionPointVertexNorms) {
        std::array<TranscendentalExpression, 3> transcendentalExpressionsForPlane{};

        //Zip iterator consisting of 3D and 1D distances l1/l2 and s1/2 for this plane | h_pq | sigma_pq for this plane
        auto zip = util::zipPair(distancesForPlane, segmentDistancesForPlane, segmentNormalOrientationsForPlane);
        auto counter = thrust::counting_iterator<unsigned int>(0);

        thrust::transform(zip.first, zip.second, counter,
                          transcendentalExpressionsForPlane.begin(), [&](const auto &tuple, const unsigned int j) {
                    using namespace util;
                    //distances l1, l2, s1, s1 for this segment q of plane p
                    const Distance &distance = thrust::get<0>(tuple);
                    //segment distance h_pq for this segment q of plane p
                    const double segmentDistance = thrust::get<1>(tuple);
                    //segment normal orientation sigma_pq for this segment q of plane p
                    const double segmentNormalOrientation = thrust::get<2>(tuple);

                    //Result for this segment
                    TranscendentalExpression transcendentalExpressionPerSegment{};

                    //Computation of the norm of P' and segment endpoints
                    // If the one of the norms == 0 then P' lies on the corresponding vertex and coincides with P''
                    const double r1Norm = projectionPointVertexNorms[(j + 1) % 3];
                    const double r2Norm = projectionPointVertexNorms[j];

                    //Compute LN_pq according to (14)
                    // If sigma_pq == 0 && either of the distances of P' to the two segment endpoints == 0 OR
                    // the 1D and 3D distances are smaller than some EPSILON
                    // then LN_pq can be set to zero
                    if ((segmentNormalOrientation == 0.0 && (r1Norm < EPSILON || r2Norm < EPSILON)) ||
                        (std::abs(distance.s1 + distance.s2) < EPSILON &&
                         std::abs(distance.l1 + distance.l2) < EPSILON)) {
                        transcendentalExpressionPerSegment.ln = 0.0;
                    } else {
                        //Implementation of
                        // log((s2_pq + l2_pq) / (s1_pq + l1_pq))
                        transcendentalExpressionPerSegment.ln =
                                std::log((distance.s2 + distance.l2) / (distance.s1 + distance.l1));
                    }

                    //Compute AN_pq according to (15)
                    // If h_p == 0 or h_pq == 0 then AN_pq is zero, too (distances are always positive!)
                    if (planeDistance < EPSILON || segmentDistance < EPSILON) {
                        transcendentalExpressionPerSegment.an = 0.0;
                    } else {
                        //Implementation of:
                        // atan(h_p * s2_pq / h_pq * l2_pq) - atan(h_p * s1_pq / h_pq * l1_pq)
                        // in a vectorized manner to reduce the number of atan(..) and divisions(..)

                        //The last subtraction is implemented via -distance.l1 (the minus/ inversion
                        // is sustained through all operations, even the atan(..)
                        xsimd::batch<double> reg1({distance.s2, distance.s1});
                        xsimd::batch<double> reg2({distance.l2, -distance.l1});

                        reg1 = xsimd::mul(reg1, planeDistance);
                        reg2 = xsimd::mul(reg2, segmentDistance);

                        reg1 = xsimd::div(reg1, reg2);
                        reg1 = xsimd::atan(reg1);

                        // "Subtraction" masked as addition (l1 was previously inverted with a minus)
                        transcendentalExpressionPerSegment.an = xsimd::hadd(reg1);
                    }

                    return transcendentalExpressionPerSegment;
                });
        return transcendentalExpressionsForPlane;
    }

    std::pair<double, Array3> GravityModel::detail::computeSingularityTerms(
            const Array3Triplet &segmentVectorsForPlane,
            const Array3 &segmentNormalOrientationForPlane,
            const Array3 &projectionPointVertexNorms,
            const Array3 &planeUnitNormal,
            double planeDistance,
            double planeNormalOrientation) {
        //1. Case: If all sigma_pq for a given plane p are 1.0 then P' lies inside the plane S_p
        if (std::all_of(segmentNormalOrientationForPlane.cbegin(), segmentNormalOrientationForPlane.cend(),
                        [](const double sigma) { return sigma == 1.0; })) {
            using namespace util;
            return std::make_pair(
                    -1.0 * util::PI2 * planeDistance,                               //sing alpha = -2pi*h_p
                    planeUnitNormal * (-1.0 * util::PI2 * planeNormalOrientation)); //sing beta  = -2pi*sigma_p*N_p
        }
        //2. Case: If sigma_pq == 0 AND norm(P' - v1) < norm(G_ij) && norm(P' - v2) < norm(G_ij) with G_ij
        // as the vector of v1 and v2
        // then P' is located on one line segment G_p of plane p, but not on any of its vertices
        auto counter2 = thrust::counting_iterator<unsigned int>(0);
        auto secondCaseBegin = util::zip(segmentVectorsForPlane.begin(), segmentNormalOrientationForPlane.begin(),
                                         counter2);
        auto secondCaseEnd = util::zip(segmentVectorsForPlane.end(), segmentNormalOrientationForPlane.end(),
                                       counter2 + 3);
        if (std::any_of(secondCaseBegin, secondCaseEnd, [&](const auto &tuple) {
            using namespace util;
            const Array3 &segmentVector = thrust::get<0>(tuple);
            const double segmentNormalOrientation = thrust::get<1>(tuple);
            const unsigned int j = thrust::get<2>(tuple);

            //segmentNormalOrientation != 0.0
            if (std::abs(segmentNormalOrientation) > EPSILON) {
                return false;
            }

            const double segmentVectorNorm = euclideanNorm(segmentVector);
            return projectionPointVertexNorms[(j + 1) % 3] < segmentVectorNorm &&
                   projectionPointVertexNorms[j] < segmentVectorNorm;
        })) {
            using namespace util;
            return std::make_pair(
                    -1.0 * util::PI * planeDistance,                                //sing alpha = -pi*h_p
                    planeUnitNormal * (-1.0 * util::PI * planeNormalOrientation));  //sing beta  = -pi*sigma_p*N_p
        }
        //3. Case If sigma_pq == 0 AND norm(P' - v1) < 0 || norm(P' - v2) < 0
        // then P' is located at one of G_p's vertices
        auto counter3 = thrust::counting_iterator<int>(0);
        auto thirdCaseBegin = util::zip(segmentNormalOrientationForPlane.begin(), counter3);
        auto thirdCaseEnd = util::zip(segmentNormalOrientationForPlane.end(), counter3 + 3);
        double r1Norm;
        double r2Norm;
        unsigned int j;
        if (std::any_of(thirdCaseBegin, thirdCaseEnd, [&](const auto &tuple) {
            using namespace util;
            const double segmentNormalOrientation = thrust::get<0>(tuple);
            j = thrust::get<1>(tuple);

            //segmentNormalOrientation != 0.0
            if (std::abs(segmentNormalOrientation) > EPSILON) {
                return false;
            }

            r1Norm = projectionPointVertexNorms[(j + 1) % 3];
            r2Norm = projectionPointVertexNorms[j];
            //r1Norm == 0.0 || r2Norm == 0.0
            return r1Norm < EPSILON || r2Norm < EPSILON;
        })) {
            using namespace util;
            //Two segment vectors G_1 and G_2 of this plane
            const Array3 &g1 = r1Norm == 0.0 ? segmentVectorsForPlane[j] : segmentVectorsForPlane[(j - 1 + 3) % 3];
            const Array3 &g2 = r1Norm == 0.0 ? segmentVectorsForPlane[(j + 1) % 3] : segmentVectorsForPlane[j];
            // theta = arcos((G_2 * -G_1) / (|G_2| * |G_1|))
            const double gdot = dot(g1 * -1.0, g2);
            const double theta =
                    gdot == 0.0 ? util::PI_2 : std::acos(gdot / (euclideanNorm(g1) * euclideanNorm(g2)));
            return std::make_pair(
                    -1.0 * theta * planeDistance,                               //sing alpha = -theta*h_p
                    planeUnitNormal * (-1.0 * theta * planeNormalOrientation)); //sing beta  = -theta*sigma_p*N_p
        }
        //4. Case Otherwise P' is located outside the plane S_p and then the singularity equals zero
        return std::make_pair(
                0.0,                                                            //sing alpha = 0
                Array3{0.0, 0.0, 0.0});                                         //sing beta  = 0
    }

    Array3 GravityModel::detail::computeNormsOfProjectionPointAndVertices(
            const Array3 &orthogonalProjectionPointOnPlane,
            const Array3Triplet &face) {
        using namespace util;
        return {euclideanNorm(orthogonalProjectionPointOnPlane - face[0]),
                euclideanNorm(orthogonalProjectionPointOnPlane - face[1]),
                euclideanNorm(orthogonalProjectionPointOnPlane - face[2])};
    }


}
