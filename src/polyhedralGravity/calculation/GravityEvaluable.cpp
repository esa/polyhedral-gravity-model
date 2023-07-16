#include "GravityEvaluable.h"


namespace polyhedralGravity {

    void GravityEvaluable::prepare() const {
        using namespace GravityModel::detail;
        // Initialize the vectors and allocate the required memory
        const size_t n = _polyhedron.countFaces();
        const auto &vertices = _polyhedron.getVertices();
        const auto &faces = _polyhedron.getFaces();
        _segmentVectors.resize(n);
        _planeUnitNormals.resize(n);
        _segmentUnitNormals.resize(n);

        // Create the iterators for the for_each loop over the polyhedral faces
        thrust::counting_iterator<size_t> begin{0};
        thrust::counting_iterator<size_t> end{n};

        // Compute the segment vectors, the plane unit normals and the segment unit normals
        thrust::for_each(thrust::device, begin, end, [this, &faces, &vertices](size_t index) {
            Array3Triplet face{vertices[faces[index][0]], vertices[faces[index][1]], vertices[faces[index][2]]};
            //1-01 Step: Compute Segment Vectors G_pq which describe each one the edge between two vertices
            _segmentVectors[index] = buildVectorsOfSegments(face[0], face[1], face[2]);
            //1-02 Step: Compute the Plane Unit Normals N_p (pointing outside the polyhedron)
            _planeUnitNormals[index] = buildUnitNormalOfPlane(_segmentVectors[index][0], _segmentVectors[index][1]);
            //1-03 Step: Compute Segment Unit Normals n_pq (normal pointing away from each segment)
            _segmentUnitNormals[index] = buildUnitNormalOfSegments(_segmentVectors[index], _planeUnitNormals[index]);
        });
    }

    template<bool Parallelization>
    GravityModelResult GravityEvaluable::evaluate(const Array3 &computationPoint) const {
        using namespace GravityModel::detail;
        using namespace util;
        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                            "Evaluation for computation point P = [{}, {}, {}] started, given density = {} kg/m^3",
                            computationPoint[0], computationPoint[1], computationPoint[2], _density);
        /*
         * Calculate V and Vx, Vy, Vz and Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
         */
        auto polyhedronIterator = transformPolyhedron(_polyhedron, computationPoint);
        auto zip1 = util::zip(polyhedronIterator.first, _segmentVectors.begin(), _planeUnitNormals.begin(),
                              _segmentUnitNormals.begin());
        auto zip2 = util::zip(polyhedronIterator.second, _segmentVectors.end(), _planeUnitNormals.end(),
                              _segmentUnitNormals.end());

        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                            "Starting to iterate over the planes...");
        GravityModelResult result{};
        auto &[potential, acceleration, gradiometricTensor] = result;

        if constexpr (Parallelization) {
            result = thrust::transform_reduce(thrust::device, zip1, zip2, &GravityEvaluable::evaluateFace, result,
                                              util::operator+<double, Array3, Array6>);
        } else {
            result = thrust::transform_reduce(thrust::host, zip1, zip2, &GravityEvaluable::evaluateFace, result,
                                              util::operator+<double, Array3, Array6>);
        }

        SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                            "Finished the sums. Applying final prefix and eliminating rounding errors.");

        //9. Step: Compute prefix consisting of GRAVITATIONAL_CONSTANT * density
        const double prefix = util::GRAVITATIONAL_CONSTANT * _density;

        //10. Step: Final expressions after application of the prefix (and a division by 2 for the potential)
        potential = (potential * prefix) / 2.0;
        acceleration = acceleration * prefix;
        gradiometricTensor = gradiometricTensor * prefix;
        return result;
    }

    // Explicit template instantiation of the single point evaluate method
    template GravityModelResult GravityEvaluable::evaluate<true>(const Array3 &computationPoints) const;

    template GravityModelResult GravityEvaluable::evaluate<false>(const Array3 &computationPoints) const;

    template<bool Parallelization>
    std::vector<GravityModelResult> GravityEvaluable::evaluate(const std::vector<Array3> &computationPoints) const {
        std::vector<GravityModelResult> result{computationPoints.size()};
        if constexpr (Parallelization) {
            thrust::transform(thrust::device, computationPoints.begin(), computationPoints.end(), result.begin(),
                              [this](const Array3 &computationPoint) {
                                  return this->evaluate<false>(computationPoint);
                              });
        } else {
            thrust::transform(thrust::host, computationPoints.begin(), computationPoints.end(), result.begin(),
                              [this](const Array3 &computationPoint) {
                                  return this->evaluate<false>(computationPoint);
                              });
        }
        return result;
    }

    // Explicit template instantiation of the multipoint evaluate method
    template std::vector<GravityModelResult>
    GravityEvaluable::evaluate<true>(const std::vector<Array3> &computationPoints) const;

    template std::vector<GravityModelResult>
    GravityEvaluable::evaluate<false>(const std::vector<Array3> &computationPoints) const;

    GravityModelResult
    GravityEvaluable::evaluateFace(const thrust::tuple<Array3Triplet, Array3Triplet, Array3, Array3Triplet> &tuple) {
        using namespace util;
        using namespace GravityModel::detail;
        const auto &face = thrust::get<0>(tuple);
        const auto &segmentVectors = thrust::get<1>(tuple);
        const auto &planeUnitNormal = thrust::get<2>(tuple);
        const auto &segmentUnitNormals = thrust::get<3>(tuple);
        SPDLOG_LOGGER_TRACE(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                            "Evaluating the plane with vertices: v1 = [{}, {}, {}], v2 = [{}, {}, {}], "
                            "v3 = [{}, {}, {}]", face[0][0], face[0][1], face[0][2], face[1][0], face[1][1], face[1][2],
                            face[2][0], face[2][1], face[2][2]);
        //1. Step: Compute ingredients for current plane which were not computed before
        //1-04 Step: Compute Plane Normal Orientation sigma_p (direction of N_p in relation to P)
        double planeNormalOrientation = computeUnitNormalOfPlaneDirection(planeUnitNormal, face[0]);
        //1-05 Step: Compute Hessian Normal Plane Representation
        HessianPlane hessianPlane = computeHessianPlane(face[0], face[1], face[2]);
        //1-06 Step: Compute distance h_p between P and P'
        double planeDistance = distanceBetweenOriginAndPlane(hessianPlane);
        //1-07 Step: Compute the actual position of P' (projection of P on the plane)
        Array3 orthogonalProjectionPointOnPlane = projectPointOrthogonallyOntoPlane(planeUnitNormal, planeDistance,
                                                                                    hessianPlane);
        //1-08 Step: Compute the segment normal orientation sigma_pq (direction of n_pq in relation to P')
        Array3 segmentNormalOrientations = computeUnitNormalOfSegmentsDirections(face, orthogonalProjectionPointOnPlane,
                                                                                 segmentUnitNormals);
        //1-09 Step: Compute the orthogonal projection point P'' of P' on each segment
        Array3Triplet orthogonalProjectionPointsOnSegmentsForPlane = projectPointOrthogonallyOntoSegments(
                orthogonalProjectionPointOnPlane, segmentNormalOrientations, face);
        //1-10 Step: Compute the segment distances h_pq between P'' and P'
        Array3 segmentDistances = distancesBetweenProjectionPoints(orthogonalProjectionPointOnPlane,
                                                                   orthogonalProjectionPointsOnSegmentsForPlane);
        //1-11 Step: Compute the 3D distances l1, l2 (between P and vertices)
        // and 1D distances s1, s2 (between P'' and vertices)
        std::array<Distance, 3> distances = distancesToSegmentEndpoints(segmentVectors,
                                                                        orthogonalProjectionPointsOnSegmentsForPlane,
                                                                        face);
        //1-12 Step: Compute the euclidian Norms of the vectors consisting of P and the vertices
        // they are later used for determining the position of P in relation to the plane
        Array3 projectionPointVertexNorms = computeNormsOfProjectionPointAndVertices(orthogonalProjectionPointOnPlane,
                                                                                     face);
        //1-13 Step: Compute the transcendental Expressions LN_pq and AN_pq
        std::array<TranscendentalExpression, 3> transcendentalExpressions = computeTranscendentalExpressions(distances,
                                                                                                             planeDistance,
                                                                                                             segmentDistances,
                                                                                                             segmentNormalOrientations,
                                                                                                             projectionPointVertexNorms);
        //1-14 Step: Compute the singularities sing A and sing B if P' is located in the plane,
        // on any vertex, or on one segment (G_pq)
        std::pair<double, Array3> singularities = computeSingularityTerms(segmentVectors, segmentNormalOrientations,
                                                                          projectionPointVertexNorms, planeUnitNormal,
                                                                          planeDistance, planeNormalOrientation);
        //2. Step: Compute Sum 1 used for potential and acceleration (first derivative)
        // sum over: sigma_pq * h_pq * LN_pq
        // --> Equation 11/12 the first summation in the brackets
        auto zipIteratorSum1PotentialAcceleration = util::zipPair(segmentNormalOrientations, segmentDistances,
                                                                  transcendentalExpressions);
        const double sum1PotentialAcceleration = std::accumulate(zipIteratorSum1PotentialAcceleration.first,
                                                                 zipIteratorSum1PotentialAcceleration.second, 0.0,
                                                                 [](double acc, const auto &tuple) {
                                                                     const double &segmentOrientation = thrust::get<0>(
                                                                             tuple);
                                                                     const double &segmentDistance = thrust::get<1>(
                                                                             tuple);
                                                                     const TranscendentalExpression &transcendentalExpressions = thrust::get<2>(
                                                                             tuple);
                                                                     return acc + segmentOrientation * segmentDistance *
                                                                                  transcendentalExpressions.ln;
                                                                 });

        //3. Step: Compute Sum 1 used for the gradiometric tensor (second derivative)
        // sum over: n_pq * LN_pq
        // --> Equation 13 the first summation in the brackets
        auto zipIteratorSum1Tensor = util::zipPair(segmentUnitNormals, transcendentalExpressions);
        const Array3 sum1Tensor = std::accumulate(zipIteratorSum1Tensor.first, zipIteratorSum1Tensor.second,
                                                  Array3{0.0, 0.0, 0.0}, [](const Array3 &acc, const auto &tuple) {
                    const Array3 &segmentNormal = thrust::get<0>(tuple);
                    const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(tuple);
                    return acc + (segmentNormal * transcendentalExpressions.ln);
                });

        //4. Step: Compute Sum 2 which is the same for every result parameter
        // sum over: sigma_pq * AN_pq
        // --> Equation 11/12/13 the second summation in the brackets
        auto zipIteratorSum2 = util::zipPair(segmentNormalOrientations, transcendentalExpressions);
        const double sum2 = std::accumulate(zipIteratorSum2.first, zipIteratorSum2.second, 0.0,
                                            [](double acc, const auto &tuple) {
                                                const double &segmentOrientation = thrust::get<0>(tuple);
                                                const TranscendentalExpression &transcendentalExpressions = thrust::get<1>(
                                                        tuple);
                                                return acc + segmentOrientation * transcendentalExpressions.an;
                                            });

        //5. Step: Sum for potential and acceleration
        // consisting of: sum1 + h_p * sum2 + sing A
        // --> Equation 11/12 the total sum of the brackets
        const double planeSumPotentialAcceleration =
                sum1PotentialAcceleration + planeDistance * sum2 + singularities.first;

        if (isCriticalDifference(planeDistance, sum2)) {
            // The multiplication planeDistance * sum2 is not the root cause, but both numbers are good
            // indicators for numerical magnitudes appearing during the calculation:
            // planeDistance gets very big when far away, sum2 remains independently very small
            SPDLOG_LOGGER_WARN(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                               "The results of point [{}, {}, {}] may be wrong due to floating point arithmetic");
        }

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
        return std::make_tuple(planeNormalOrientation * planeDistance * planeSumPotentialAcceleration,
                               planeUnitNormal * planeSumPotentialAcceleration, concat(first, second));
    }

    std::string GravityEvaluable::toString() const {
        std::stringstream ss;
        ss  << "<polyhedral_gravity.GravityEvaluable, density=" << _density
            << ", vertices= " << _polyhedron.countVertices()
            << ", faces= " << _polyhedron.countFaces() << ">";
        return ss.str();
    }

}
