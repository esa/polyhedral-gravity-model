#include "GravityModelVectorUtility.h"

namespace polyhedralGravity {

    using namespace GravityModel::detail;

    std::vector<Array3Triplet> GravityModel::calculateSegmentVectors(const Polyhedron &polyhedron) {
        std::vector<Array3Triplet> segmentVectors{polyhedron.countFaces()};
        //Calculate G_ij for every plane given as input the three vertices of every face
        for (size_t i = 0; i < polyhedron.countFaces(); ++i) {
            const IndexArray3 &face = polyhedron.getFace(i);
            segmentVectors[i] = buildVectorsOfSegments(polyhedron.getVertex(face[0]), polyhedron.getVertex(face[1]),
                                                       polyhedron.getVertex(face[2]));
        }
        return segmentVectors;
    }

    std::vector<Array3> GravityModel::calculatePlaneUnitNormals(const std::vector<Array3Triplet> &segmentVectors) {
        std::vector<Array3> planeUnitNormals{segmentVectors.size()};
        //Calculate N_p for every plane given as input: G_i0 and G_i1 of every plane
        for (size_t i = 0; i < segmentVectors.size(); ++i) {
            planeUnitNormals[i] = buildUnitNormalOfPlane(segmentVectors[i][0], segmentVectors[i][1]);
        }
        return planeUnitNormals;
    }

    std::vector<Array3Triplet> GravityModel::calculateSegmentUnitNormals(
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3> &planeUnitNormals) {
        std::vector<Array3Triplet> segmentUnitNormals{segmentVectors.size()};
        //Loop over G_i (running i=p) and N_p calculating n_p for every plane
        for (size_t i = 0; i < segmentVectors.size(); ++i) {
            segmentUnitNormals[i] = buildUnitNormalOfSegments(segmentVectors[i], planeUnitNormals[i]);
        }
        return segmentUnitNormals;
    }

    std::vector<double>
    GravityModel::calculatePlaneNormalOrientations(const Array3 &computationPoint, const Polyhedron &polyhedron,
                                                   const std::vector<Array3> &planeUnitNormals) {
        std::vector<double> planeNormalOrientations(planeUnitNormals.size(), 0.0);
        //Calculate sigma_p for every plane given as input: N_p and vertex0 of every face
        for (size_t i = 0; i < planeUnitNormals.size(); ++i) {
            //The first vertices' coordinates of the given face consisting of G_i's
            const Array3Triplet face = polyhedron.getResolvedFace(i, computationPoint);
            planeNormalOrientations[i] = computeUnitNormalOfPlaneDirection(planeUnitNormals[i], face[0]);
        }
        return planeNormalOrientations;
    }

    std::vector<HessianPlane>
    GravityModel::calculateFacesToHessianPlanes(const Array3 &computationPoint, const Polyhedron &polyhedron) {
        std::vector<HessianPlane> hessianPlanes{polyhedron.countFaces()};
        //Calculate for each face/ plane/ triangle (here) the Hessian Plane
        for (size_t i = 0; i < polyhedron.countFaces(); ++i) {
            //The three vertices put up the plane, p is the origin of the reference system default 0,0,0
            const Array3Triplet face = polyhedron.getResolvedFace(i, computationPoint);
            hessianPlanes[i] = computeHessianPlane(face[0], face[1], face[2]);
        }
        return hessianPlanes;
    }

    std::vector<double> GravityModel::calculatePlaneDistances(const std::vector<HessianPlane> &plane) {
        std::vector<double> planeDistances(plane.size(), 0.0);
        //For each plane compute h_p
        for (size_t i = 0; i < plane.size(); ++i) {
            planeDistances[i] = distanceBetweenOriginAndPlane(plane[i]);
        }
        return planeDistances;
    }

    std::vector<Array3> GravityModel::calculateOrthogonalProjectionPointsOnPlane(
            const std::vector<HessianPlane> &hessianPlanes,
            const std::vector<Array3> &planeUnitNormals,
            const std::vector<double> &planeDistances) {
        std::vector<Array3> orthogonalProjectionPointsOfP{planeUnitNormals.size()};
        //Calculates the Projection Point P' for every plane p
        for (size_t i = 0; i < planeUnitNormals.size(); ++i) {
            orthogonalProjectionPointsOfP[i] =
                    projectPointOrthogonallyOntoPlane(planeUnitNormals[i], planeDistances[i], hessianPlanes[i]);
        }
        return orthogonalProjectionPointsOfP;
    }

    std::vector<Array3> GravityModel::calculateSegmentNormalOrientations(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentUnitNormals,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane) {
        std::vector<Array3> segmentNormalOrientations{segmentUnitNormals.size()};
        //Calculates the segment normal orientation sigma_pq for every plane p
        for (size_t i = 0; i < segmentUnitNormals.size(); ++i) {
            const Array3Triplet face = polyhedron.getResolvedFace(i, computationPoint);
            segmentNormalOrientations[i] = computeUnitNormalOfSegmentsDirections(
                    face, orthogonalProjectionPointsOnPlane[i], segmentUnitNormals[i]);
        }
        return segmentNormalOrientations;
    }

    std::vector<Array3Triplet> GravityModel::calculateOrthogonalProjectionPointsOnSegments(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3> &segmentNormalOrientation) {
        std::vector<Array3Triplet> orthogonalProjectionPointsOnSegments{orthogonalProjectionPointsOnPlane.size()};
        //The outer loop with the running i --> the planes
        for (size_t i = 0; i < orthogonalProjectionPointsOnPlane.size(); ++i) {
            //P' for plane i, sigma_pq[i] with fixed i, the nodes making up plane i
            const Array3Triplet face = polyhedron.getResolvedFace(i, computationPoint);
            orthogonalProjectionPointsOnSegments[i] = projectPointOrthogonallyOntoSegments(
                    orthogonalProjectionPointsOnPlane[i], segmentNormalOrientation[i], face);
        }
        return orthogonalProjectionPointsOnSegments;
    }

    std::vector<Array3> GravityModel::calculateSegmentDistances(
            const std::vector<Array3> &orthogonalProjectionPointsOnPlane,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        std::vector<Array3> segmentDistances{orthogonalProjectionPointsOnPlane.size()};
        //Iterating over planes (P'_i and P''_i are the parameters)
        for (size_t i = 0; i < orthogonalProjectionPointsOnPlane.size(); ++i) {
            segmentDistances[i] = distancesBetweenProjectionPoints(orthogonalProjectionPointsOnPlane[i],
                                                                   orthogonalProjectionPointsOnSegment[i]);
        }
        return segmentDistances;
    }

    std::vector<std::array<Distance, 3>> GravityModel::calculateDistances(
            const Array3 &computationPoint,
            const Polyhedron &polyhedron,
            const std::vector<Array3Triplet> &segmentVectors,
            const std::vector<Array3Triplet> &orthogonalProjectionPointsOnSegment) {
        std::vector<std::array<Distance, 3>> distances{segmentVectors.size()};
        for (size_t i = 0; i < segmentVectors.size(); ++i) {
            const Array3Triplet face = polyhedron.getResolvedFace(i, computationPoint);
            distances[i] = distancesToSegmentEndpoints(segmentVectors[i], orthogonalProjectionPointsOnSegment[i], face);
        }
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
        std::vector<std::array<TranscendentalExpression, 3>> transcendentalExpressions{distances.size()};
        for (size_t i = 0; i < distances.size(); ++i) {
            const Array3Triplet face = polyhedron.getResolvedFace(i, computationPoint);
            const Array3 projectionPointVertexNorms =
                    computeNormsOfProjectionPointAndVertices(orthogonalProjectionPointsOnPlane[i], face);
            transcendentalExpressions[i] = computeTranscendentalExpressions(
                    distances[i], planeDistances[i], segmentDistances[i], segmentNormalOrientation[i],
                    projectionPointVertexNorms);
        }
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
        //The result
        std::vector<std::pair<double, Array3>> singularities{planeDistances.size()};
        for (size_t i = 0; i < planeDistances.size(); ++i) {
            const Array3Triplet face = polyhedron.getResolvedFace(i, computationPoint);
            const Array3 projectionPointVertexNorms =
                    computeNormsOfProjectionPointAndVertices(orthogonalProjectionPointsOnPlane[i], face);
            const SingularityTerms<double> terms = computeSingularityTerms(
                    segmentVectors[i], segmentNormalOrientation[i], projectionPointVertexNorms, planeUnitNormals[i],
                    planeDistances[i], planeNormalOrientations[i]);
            singularities[i] = std::make_pair(terms.alpha, terms.beta);
        }
        return singularities;
    }

}
