#include "GravityEvaluable.h"

#ifdef POLYHEDRAL_GRAVITY_ENABLE_OPENCL
#include "polyhedralGravity/opencl/OpenCLEvaluation.h"
#endif

namespace polyhedralGravity {

    namespace {

        /**
         * Creates the OpenCL engine for the polyhedron, or returns nullptr if the evaluation has to
         * happen on the host.
         *
         * Requesting OpenCL is a preference, not a demand: a library compiled without OpenCL, a machine
         * without a suitable device, and a device lacking @c cl_khr_fp64 for a FLOAT64 request all lead
         * back to the CPU backend rather than to an error, so that the default backend works everywhere.
         *
         * @param polyhedron the polyhedron to upload
         * @param backend the requested backend
         * @param precision the requested precision
         * @return the OpenCL engine, or nullptr to evaluate on the host
         */
        std::shared_ptr<opencl::OpenCLEvaluation> createOpenCLEvaluation(
                [[maybe_unused]] const Polyhedron &polyhedron,
                const ComputeBackend backend,
                [[maybe_unused]] const ComputePrecision precision) {
            if (backend != ComputeBackend::OPENCL) {
                return nullptr;
            }
#ifndef POLYHEDRAL_GRAVITY_ENABLE_OPENCL
            POLYHEDRAL_GRAVITY_LOG_DEBUG("The OpenCL backend was requested, but this library was compiled "
                                         "without OpenCL support. Falling back to the CPU backend.");
            return nullptr;
#else
            try {
                return std::make_shared<opencl::OpenCLEvaluation>(polyhedron, precision);
            } catch (const std::exception &e) {
                POLYHEDRAL_GRAVITY_LOG_WARN("The OpenCL backend is unavailable, falling back to the CPU "
                                            "backend. Reason: {}", e.what());
                return nullptr;
            }
#endif
        }

    }// namespace

    GravityEvaluable::GravityEvaluable(const Polyhedron &polyhedron, const ComputeBackend backend,
                                       const ComputePrecision precision)
        : _polyhedron{polyhedron},
          _backend{backend},
          _precision{precision},
          _openCLEvaluation{createOpenCLEvaluation(polyhedron, backend, precision)} {
        if (_openCLEvaluation == nullptr) {
            _backend = ComputeBackend::CPU;
            // The host-side caches are what the CPU backend evaluates from, so they cannot be deferred
            this->prepare();
            _prepared = true;
        }
    }

    GravityEvaluable::GravityEvaluable(const Polyhedron &polyhedron,
                                       const std::vector<Array3Triplet> &segmentVectors,
                                       const std::vector<Array3> &planeUnitNormals,
                                       const std::vector<Array3Triplet> &segmentUnitNormals,
                                       const ComputeBackend backend, const ComputePrecision precision)
        : _polyhedron{polyhedron},
          _segmentVectors{segmentVectors},
          _planeUnitNormals{planeUnitNormals},
          _segmentUnitNormals{segmentUnitNormals},
          _backend{backend},
          _precision{precision},
          _openCLEvaluation{createOpenCLEvaluation(polyhedron, backend, precision)},
          _prepared{true} {
        if (_openCLEvaluation == nullptr) {
            _backend = ComputeBackend::CPU;
        }
    }

    void GravityEvaluable::ensurePrepared() const {
        if (!_prepared) {
            this->prepare();
            _prepared = true;
        }
    }

    ComputeBackend GravityEvaluable::getComputeBackend() const {
        return _backend;
    }

    ComputePrecision GravityEvaluable::getComputePrecision() const {
        return _precision;
    }

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
        POLYHEDRAL_GRAVITY_LOG_DEBUG("Evaluation for computation point P = [{}, {}, {}] started, given density = {} kg/m^3",
                computationPoint[0], computationPoint[1], computationPoint[2], _polyhedron.getDensity());
#ifdef POLYHEDRAL_GRAVITY_ENABLE_OPENCL
        // The device parallelizes over the polyhedron's faces itself, so the host-side
        // Parallelization flag has no meaning here
        if (_openCLEvaluation != nullptr) {
            return _openCLEvaluation->evaluate(computationPoint);
        }
#endif
        /*
         * Calculate V and Vx, Vy, Vz and Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
         */
        const auto &[polyBegin, polyEnd] = _polyhedron.transformIterator(computationPoint);
        const auto zip1 = zip(polyBegin, _segmentVectors.begin(), _planeUnitNormals.begin(), _segmentUnitNormals.begin());
        const auto zip2 = zip(polyEnd, _segmentVectors.end(), _planeUnitNormals.end(), _segmentUnitNormals.end());

        POLYHEDRAL_GRAVITY_LOG_DEBUG("Starting to iterate over the planes...");
        GravityModelResult result{};
        auto &[potential, acceleration, gradiometricTensor] = result;

        if constexpr (Parallelization) {
            result = thrust::transform_reduce(thrust::device, zip1, zip2, &GravityEvaluable::evaluateFace, result,
                                              util::operator+ <double, Array3, Array6>);
        } else {
            result = thrust::transform_reduce(thrust::host, zip1, zip2, &GravityEvaluable::evaluateFace, result,
                                              util::operator+ <double, Array3, Array6>);
        }

        POLYHEDRAL_GRAVITY_LOG_DEBUG("Finished the sums. Applying final prefix.");

        //9. Step: Compute prefix consisting of GRAVITATIONAL_CONSTANT * density
        //and correction factors depending on alignment of the normals and polyhedral mesh unit
        const double prefix = _polyhedron.getGravityModelScaling();

        //10. Step: Final expressions after application of the prefix (and a division by 2 for the potential)
        potential = (potential * prefix) / 2.0;
        acceleration = acceleration * (-1.0 * prefix);
        gradiometricTensor = gradiometricTensor * prefix;
        return result;
    }

    // Explicit template instantiation of the single point evaluate method
    template GravityModelResult GravityEvaluable::evaluate<true>(const Array3 &computationPoints) const;

    template GravityModelResult GravityEvaluable::evaluate<false>(const Array3 &computationPoints) const;

    template<bool Parallelization>
    std::vector<GravityModelResult> GravityEvaluable::evaluate(const std::vector<Array3> &computationPoints) const {
#ifdef POLYHEDRAL_GRAVITY_ENABLE_OPENCL
        // The points are evaluated one after another: the device is already saturated by a single
        // point's faces, and its command queue and kernel arguments cannot be driven concurrently
        if (_openCLEvaluation != nullptr) {
            return _openCLEvaluation->evaluate(computationPoints);
        }
#endif
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
        POLYHEDRAL_GRAVITY_LOG_TRACE("Evaluating the plane with vertices: v1 = [{}, {}, {}], v2 = [{}, {}, {}], "
                                     "v3 = [{}, {}, {}]",
                                     face[0][0], face[0][1], face[0][2], face[1][0], face[1][1], face[1][2],
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
            POLYHEDRAL_GRAVITY_LOG_WARN("While evaluating the plane with coordinates v1 = [{}, {}, {}], v2 = [{}, {}, {}], "
                                        "v3 = [{}, {}, {}] (with computation point re-located at the origin) a "
                                        "significant difference of magnitudes occurred during the evaluation. "
                                        "This may lead to numerically unstable results!",
                                        face[0][0], face[0][1], face[0][2],
                                        face[1][0], face[1][1], face[1][2], face[2][0], face[2][1], face[2][2]);
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
        std::stringstream sstream;
        const auto[unitPotential, unitAcceleration, unitGradiometricTensor] = getOutputMetricUnit();
        sstream << "<polyhedral_gravity.GravityEvaluable, polyhedron = " << _polyhedron.toString()
                << ", output_units = " << unitPotential << ", " << unitAcceleration << ", " << unitGradiometricTensor
                << ", backend = " << _backend;
#ifdef POLYHEDRAL_GRAVITY_ENABLE_OPENCL
        if (_openCLEvaluation != nullptr) {
            sstream << " (" << _openCLEvaluation->getDeviceName() << ", " << _precision << ")";
        }
#endif
        sstream << ">";
        return sstream.str();
    }

    std::array<std::string, 3> GravityEvaluable::getOutputMetricUnit() const {
        const auto metric = _polyhedron.getMeshUnit();
        if (metric != MetricUnit::UNITLESS) {
            const std::string metricString = _polyhedron.getMeshUnitAsString();
            return {metricString + "^2/s^2", metricString + "/s^2", "1/s^2"};
        } else {
            return {"1/s^2", "1/s^2", "1/s^2"};
        }
    }

    std::tuple<Polyhedron, std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
    GravityEvaluable::getState() const {
        // On the OpenCL backend the caches are filled lazily, and reading the state is what needs them
        this->ensurePrepared();
        return std::make_tuple(_polyhedron, _segmentVectors, _planeUnitNormals, _segmentUnitNormals);
    }

}// namespace polyhedralGravity
