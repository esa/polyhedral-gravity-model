#include "GravityModelVectorUtility.h"

namespace polyhedralGravity {

    std::vector<Array3Triplet> GravityModel::calculateSegmentVectors(const Polyhedron &polyhedron) {
        using namespace detail;
        std::vector<Array3Triplet> segmentVectors{polyhedron.countFaces()};
        //Calculate G_ij for every plane given as input the three vertices of every face
        std::transform(polyhedron.getFaces().cbegin(), polyhedron.getFaces().cend(), segmentVectors.begin(),
                       [&polyhedron](const std::array<size_t, 3> &face) -> Array3Triplet {
                           const Array3 &vertex0 = polyhedron.getVertex(face[0]);
                           const Array3 &vertex1 = polyhedron.getVertex(face[1]);
                           const Array3 &vertex2 = polyhedron.getVertex(face[2]);
                           return buildVectorsOfSegments(vertex0, vertex1, vertex2);
                       });
        return segmentVectors;
    }

    std::vector<Array3> GravityModel::calculatePlaneUnitNormals(const std::vector<Array3Triplet> &segmentVectors) {
        using namespace detail;
        std::vector<Array3> planeUnitNormals{segmentVectors.size()};
        //Calculate N_p for every plane given as input: G_i0 and G_i1 of every plane
        std::transform(segmentVectors.cbegin(), segmentVectors.cend(), planeUnitNormals.begin(),
                       [](const Array3Triplet &segmentVectorsForPlane) -> Array3 {
                           return buildUnitNormalOfPlane(segmentVectorsForPlane[0], segmentVectorsForPlane[1]);
                       });
        return planeUnitNormals;
    }

    std::vector<Array3Triplet> GravityModel::calculateSegmentUnitNormals(
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3> &planeUnitNormals) {
        using namespace detail;
        std::vector<Array3Triplet> segmentUnitNormals{segmentVectors.size()};
        //Loop" over G_i (running i=p) and N_p calculating n_p for every plane
        std::transform(segmentVectors.cbegin(), segmentVectors.cend(), planeUnitNormals.cbegin(),
                       segmentUnitNormals.begin(),
                       [](const Array3Triplet &segmentVectorsForPlane, const Array3 &planeUnitNormal) {
                           return buildUnitNormalOfSegments(segmentVectorsForPlane, planeUnitNormal);
                       });
        return segmentUnitNormals;
    }

    std::vector<double>
    GravityModel::calculatePlaneNormalOrientations(const Array3 &computationPoint, const Polyhedron &polyhedron,
                                                   const std::vector<Array3> &planeUnitNormals) {
        using namespace detail;
        std::vector<double> planeNormalOrientations(planeUnitNormals.size(), 0.0);
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        //Calculate sigma_p for every plane given as input: N_p and vertex0 of every face
        std::transform(planeUnitNormals.cbegin(), planeUnitNormals.cend(), transformedPolyhedronIt.first,
                       planeNormalOrientations.begin(),
                       [](const Array3 &planeUnitNormal, const Array3Triplet &face) {
                           //The first vertices' coordinates of the given face consisting of G_i's
                           return computeUnitNormalOfPlaneDirection(planeUnitNormal, face[0]);
                       });
        return planeNormalOrientations;
    }

    std::vector<HessianPlane>
    GravityModel::calculateFacesToHessianPlanes(const Array3 &computationPoint, const Polyhedron &polyhedron) {
        using namespace detail;
        std::vector<HessianPlane> hessianPlanes{polyhedron.countFaces()};
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        //Calculate for each face/ plane/ triangle (here) the Hessian Plane
        std::transform(transformedPolyhedronIt.first, transformedPolyhedronIt.second, hessianPlanes.begin(),
                       [](const Array3Triplet &face) -> HessianPlane {
                           using namespace util;
                           //The three vertices put up the plane, p is the origin of the reference system default 0,0,0
                           return computeHessianPlane(face[0], face[1], face[2]);
                       });
        return hessianPlanes;
    }

    std::vector<double> GravityModel::calculatePlaneDistances(const std::vector<HessianPlane> &plane) {
        using namespace detail;
        std::vector<double> planeDistances(plane.size(), 0.0);
        //For each plane compute h_p
        std::transform(plane.cbegin(), plane.cend(), planeDistances.begin(),
                       [](const HessianPlane &plane) -> double {
                           return distanceBetweenOriginAndPlane(plane);
                       });
        return planeDistances;
    }

    std::vector<Array3> GravityModel::calculateOrthogonalProjectionPointsOnPlane(
            const std::vector<HessianPlane> &hessianPlanes,
            const std::vector<Array3> &planeUnitNormals,
            const std::vector<double> &planeDistances) {
        using namespace detail;
        std::vector<Array3> orthogonalProjectionPointsOfP{planeUnitNormals.size()};

        //Zip the three required arguments together: Plane normal N_i, Plane Distance h_i and the Hessian Form
        auto zip = util::zipPair(planeUnitNormals, planeDistances, hessianPlanes);

        //Calculates the Projection Point P' for every plane p
        thrust::transform(zip.first, zip.second, orthogonalProjectionPointsOfP.begin(), [](const auto &tuple) {
            using namespace util;
            const Array3 &planeUnitNormal = thrust::get<0>(tuple);
            const double planeDistance = thrust::get<1>(tuple);
            const HessianPlane &hessianPlane = thrust::get<2>(tuple);
            return projectPointOrthogonallyOntoPlane(planeUnitNormal, planeDistance, hessianPlane);
        });
        return orthogonalProjectionPointsOfP;
    }

    std::vector<Array3> GravityModel::calculateSegmentNormalOrientations(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentUnitNormals,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane) {
        using namespace detail;
        std::vector<Array3> segmentNormalOrientations{segmentUnitNormals.size()};

        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(transformedPolyhedronIt.first,
                               orthogonalProjectionPointsOnPlane.begin(), segmentUnitNormals.begin());
        auto last = util::zip(transformedPolyhedronIt.second,
                              orthogonalProjectionPointsOnPlane.end(), segmentUnitNormals.end());

        //Calculates the segment normal orientation sigma_pq for every plane p
        thrust::transform(first, last, segmentNormalOrientations.begin(), [](const auto &tuple) {
            const Array3Triplet &face = thrust::get<0>(tuple);
            const Array3 &projectionPointOnPlaneForPlane = thrust::get<1>(tuple);
            const Array3Triplet &segmentUnitNormalsForPlane = thrust::get<2>(tuple);

            return computeUnitNormalOfSegmentsDirections(face, projectionPointOnPlaneForPlane,
                                                         segmentUnitNormalsForPlane);
        });
        return segmentNormalOrientations;
    }

    std::vector<Array3Triplet> GravityModel::calculateOrthogonalProjectionPointsOnSegments(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3> &segmentNormalOrientation) {
        using namespace detail;
        std::vector<Array3Triplet> orthogonalProjectionPointsOnSegments{orthogonalProjectionPointsOnPlane.size()};

        //Zip the three required arguments together: P' for every plane, sigma_pq for every segment, the faces
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(orthogonalProjectionPointsOnPlane.begin(), segmentNormalOrientation.begin(),
                               transformedPolyhedronIt.first);
        auto last = util::zip(orthogonalProjectionPointsOnPlane.end(), segmentNormalOrientation.end(),
                              transformedPolyhedronIt.second);

        //The outer loop with the running i --> the planes
        thrust::transform(first, last, orthogonalProjectionPointsOnSegments.begin(), [](const auto &tuple) {
            //P' for plane i, sigma_pq[i] with fixed i, the nodes making up plane i
            const Array3 &orthogonalProjectionPointOnPlane = thrust::get<0>(tuple);
            const Array3 &segmentNormalOrientationsForPlane = thrust::get<1>(tuple);
            const Array3Triplet &face = thrust::get<2>(tuple);
            return projectPointOrthogonallyOntoSegments(
                    orthogonalProjectionPointOnPlane, segmentNormalOrientationsForPlane, face);
        });

        return orthogonalProjectionPointsOnSegments;
    }

    std::vector<Array3> GravityModel::calculateSegmentDistances(
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        using namespace detail;
        std::vector<Array3> segmentDistances{orthogonalProjectionPointsOnPlane.size()};
        //Iterating over planes (P'_i and P''_i are the parameters of the lambda)
        std::transform(orthogonalProjectionPointsOnPlane.cbegin(), orthogonalProjectionPointsOnPlane.cend(),
                       orthogonalProjectionPointsOnSegment.cbegin(), segmentDistances.begin(),
                       [](const Array3 projectionPointOnPlane, const Array3Triplet &projectionPointOnSegments) {
                           return distancesBetweenProjectionPoints(projectionPointOnPlane, projectionPointOnSegments);
                       });
        return segmentDistances;
    }

    std::vector<std::array<Distance, 3>> GravityModel::calculateDistances(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        using namespace detail;
        std::vector<std::array<Distance, 3>> distances{segmentVectors.size()};

        //Zip the three required arguments together: G_ij for every segment, P'' for every segment
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(segmentVectors.begin(), orthogonalProjectionPointsOnSegment.begin(),
                               transformedPolyhedronIt.first);
        auto last = util::zip(segmentVectors.end(), orthogonalProjectionPointsOnSegment.end(),
                              transformedPolyhedronIt.second);

        thrust::transform(first, last, distances.begin(), [](const auto &tuple) {
            const Array3Triplet &segmentVectorsForPlane = thrust::get<0>(tuple);
            const Array3Triplet &orthogonalProjectionPointsOnSegmentForPlane = thrust::get<1>(tuple);
            const Array3Triplet &face = thrust::get<2>(tuple);
            return distancesToSegmentEndpoints(segmentVectorsForPlane, orthogonalProjectionPointsOnSegmentForPlane,
                                               face);
        });
        return distances;
    }

    std::vector<std::array<TranscendentalExpression, 3>> GravityModel::calculateTranscendentalExpressions(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<std::array<Distance, 3>> &distances,
            const std::vector<double> &planeDistances,
            const std::vector<Array3> &segmentDistances,
            const std::vector<Array3> &segmentNormalOrientation,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane) {
        using namespace detail;
        std::vector<std::array<TranscendentalExpression, 3>> transcendentalExpressions{distances.size()};

        //Zip iterator consisting of  3D and 1D distances l1/l2 and s1/2 | h_p | h_pq | sigma_pq | P'_p | faces
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(distances.begin(), planeDistances.begin(), segmentDistances.begin(),
                               segmentNormalOrientation.begin(), orthogonalProjectionPointsOnPlane.begin(),
                               transformedPolyhedronIt.first);
        auto last = util::zip(distances.end(), planeDistances.end(), segmentDistances.end(),
                              segmentNormalOrientation.end(), orthogonalProjectionPointsOnPlane.end(),
                              transformedPolyhedronIt.second);

        thrust::transform(first, last, transcendentalExpressions.begin(), [](const auto &tuple) {
            const std::array<Distance, 3> &distancesForPlane = thrust::get<0>(tuple);
            const double planeDistance = thrust::get<1>(tuple);
            const Array3 &segmentDistancesForPlane = thrust::get<2>(tuple);
            const Array3 &segmentNormalOrientationsForPlane = thrust::get<3>(tuple);
            const Array3 &projectionPointOnPlane = thrust::get<4>(tuple);
            const Array3Triplet &face = thrust::get<5>(tuple);

            const Array3 projectionPointVertexNorms = computeNormsOfProjectionPointAndVertices(
                    projectionPointOnPlane, face);

            return computeTranscendentalExpressions(distancesForPlane, planeDistance, segmentDistancesForPlane,
                                                    segmentNormalOrientationsForPlane,
                                                    projectionPointVertexNorms);
        });
        return transcendentalExpressions;
    }

    std::vector<std::pair<double, Array3>> GravityModel::calculateSingularityTerms(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3> &segmentNormalOrientation,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<double> &planeDistances,
            const std::vector<double> &planeNormalOrientations,
            const std::vector<Array3> &planeUnitNormals) {
        using namespace detail;
        //The result
        std::vector<std::pair<double, Array3>> singularities{planeDistances.size()};

        //Zip iterator consisting of G_ij vectors | sigma_pq | faces | P' | h_p | sigma_p | N_i
        auto transformedPolyhedronIt = transformPolyhedron(polyhedron, computationPoint);
        auto first = util::zip(segmentVectors.begin(), segmentNormalOrientation.begin(),
                               orthogonalProjectionPointsOnPlane.begin(), planeUnitNormals.begin(),
                               planeDistances.begin(), planeNormalOrientations.begin(), transformedPolyhedronIt.first);
        auto last = util::zip(segmentVectors.end(), segmentNormalOrientation.end(),
                              orthogonalProjectionPointsOnPlane.end(), planeUnitNormals.end(),
                              planeDistances.end(), planeNormalOrientations.end(), transformedPolyhedronIt.second);

        thrust::transform(first, last, singularities.begin(), [&](const auto &tuple) {
            const Array3Triplet &segmentVectorsForPlane = thrust::get<0>(tuple);
            const Array3 segmentNormalOrientationForPlane = thrust::get<1>(tuple);
            const Array3 &orthogonalProjectionPointOnPlane = thrust::get<2>(tuple);
            const Array3 &planeUnitNormal = thrust::get<3>(tuple);
            const double planeDistance = thrust::get<4>(tuple);
            const double planeNormalOrientation = thrust::get<5>(tuple);
            const Array3Triplet &face = thrust::get<6>(tuple);

            const Array3 projectionPointVertexNorms = computeNormsOfProjectionPointAndVertices(
                    orthogonalProjectionPointOnPlane, face);

            return computeSingularityTerms(segmentVectorsForPlane, segmentNormalOrientationForPlane,
                                           projectionPointVertexNorms, planeUnitNormal, planeDistance,
                                           planeNormalOrientation);
        });
        return singularities;
    }

}