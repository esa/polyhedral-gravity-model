/**
 * Evaluates Tsoulis et al.'s polyhedral gravity model for one computation point P.
 *
 * Every work-item handles one polyhedral face and computes that face's contribution to Equations
 * (11), (12), and (13). The contributions are then summed within the work-group so that only one
 * partial result per work-group leaves the device.
 *
 * This is the device counterpart of GravityEvaluable::evaluateFace() and the step numbering below
 * refers to the very same steps.
 */

/**
 * The signum function with a cutoff around zero.
 * @param value the value to take the sign of
 * @return -1, 0, or 1
 */
int sgn(const FloatType value) {
    if (value < -EPSILON_ZERO_OFFSET) {
        return -1;
    }
    if (value > EPSILON_ZERO_OFFSET) {
        return 1;
    }
    return 0;
}

/**
 * Transposes a 3x3 matrix given as three row vectors in place.
 * @param matrix the matrix to transpose
 */
void transpose(FloatType3 matrix[3]) {
    const FloatType3 copy[3] = {
            matrix[0], matrix[1], matrix[2],
    };

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            matrix[i][j] = copy[j][i];
        }
    }
}

/**
 * The determinant of a 3x3 matrix given as three row vectors.
 * @param matrix the matrix
 * @return the determinant
 */
FloatType det(const FloatType3 matrix[3]) {
    return matrix[0][0] * matrix[1][1] * matrix[2][2] +
           matrix[0][1] * matrix[1][2] * matrix[2][0] +
           matrix[0][2] * matrix[1][0] * matrix[2][1] -
           matrix[0][2] * matrix[1][1] * matrix[2][0] -
           matrix[0][0] * matrix[1][2] * matrix[2][1] -
           matrix[0][1] * matrix[1][0] * matrix[2][2];
}

/**
 * The determinant of the 3x3 matrix built from the three given row vectors.
 * @param a first row
 * @param b second row
 * @param c third row
 * @return the determinant
 */
FloatType detOfRows(const FloatType3 a, const FloatType3 b, const FloatType3 c) {
    const FloatType3 matrix[3] = {a, b, c};
    return det(matrix);
}

/** The 3D distances l1_pq, l2_pq and the 1D distances s1_pq, s2_pq of one segment, see Tsoulis Equation (9)/(10) */
typedef struct {
    FloatType l1;
    FloatType l2;
    FloatType s1;
    FloatType s2;
} Distance;

/** The transcendental expressions LN_pq and AN_pq of one segment, see Tsoulis Equation (14)/(15) */
typedef struct {
    FloatType ln;
    FloatType an;
} TranscendentalExpression;

/**
 * Computes the angle theta of the singularity terms sing A and sing B for one plane, i.e. determines
 * where the orthogonal projection point P' is located in relation to the plane.
 *
 * @param faceIndex the index of the evaluated face
 * @param segmentNormalOrientations the segment normal orientations sigma_pq of this plane
 * @param projectionPointVertexNorms the norms between P' and this plane's vertices
 * @param segmentVectors all faces' segment vectors G_pq
 * @return theta, from which the caller derives sing alpha and sing beta
 */
FloatType computeSingularityAngle(
        const int faceIndex,
        const int3 segmentNormalOrientations,
        const FloatType3 projectionPointVertexNorms,
        global const FloatType3 *segmentVectors) {

    //1. Case: If all sigma_pq for a given plane p are 1.0 then P' lies inside the plane S_p
    bool allInside = true;
    for (uint index = 0; index < 3; ++index) {
        allInside &= segmentNormalOrientations[index] == 1;
    }
    if (allInside) {
        return PI2;
    }

    //2. Case: If sigma_pq == 0 AND norm(P' - v1) < norm(G_ij) && norm(P' - v2) < norm(G_ij) with G_ij
    // as the vector of v1 and v2, then P' is located on one line segment G_p of plane p,
    // but not on any of its vertices -- the latter is what the epsilon comparisons exclude, so that
    // such points are handled by the third case below instead.
    bool anyOnLine = false;
    for (uint index = 0; index < 3; ++index) {
        if (segmentNormalOrientations[index] != 0) {
            continue;
        }
        const FloatType segmentVectorNorm = length(segmentVectors[faceIndex * 3 + index]);
        anyOnLine |= projectionPointVertexNorms[(index + 1) % 3] < segmentVectorNorm &&
                     projectionPointVertexNorms[index] < segmentVectorNorm &&
                     projectionPointVertexNorms[(index + 1) % 3] >= EPSILON_ZERO_OFFSET &&
                     projectionPointVertexNorms[index] >= EPSILON_ZERO_OFFSET;
    }
    if (anyOnLine) {
        return PI;
    }

    //3. Case: If sigma_pq == 0 AND norm(P' - v1) == 0 || norm(P' - v2) == 0
    // then P' is located at one of G_p's vertices
    for (uint index = 0; index < 3; ++index) {
        if (segmentNormalOrientations[index] != 0) {
            continue;
        }

        const FloatType r1Norm = projectionPointVertexNorms[(index + 1) % 3];
        const FloatType r2Norm = projectionPointVertexNorms[index];

        if (!(r1Norm < EPSILON_ZERO_OFFSET || r2Norm < EPSILON_ZERO_OFFSET)) {
            continue;
        }

        //Two segment vectors G_1 and G_2 of this plane
        const FloatType3 g1 = r1Norm < EPSILON_ZERO_OFFSET
                                      ? segmentVectors[faceIndex * 3 + index]
                                      : segmentVectors[faceIndex * 3 + (index - 1 + 3) % 3];
        const FloatType3 g2 = r1Norm < EPSILON_ZERO_OFFSET
                                      ? segmentVectors[faceIndex * 3 + (index + 1) % 3]
                                      : segmentVectors[faceIndex * 3 + index];

        // theta = arccos((G_2 * -G_1) / (|G_2| * |G_1|))
        const FloatType gdot = dot(-g1, g2);
        return gdot == 0.0 ? PI_2 : acos(gdot / (length(g1) * length(g2)));
    }

    //4. Case: Otherwise P' is located outside the plane S_p and then the singularity equals zero
    return 0.0;
}

kernel void evaluateFaces(
        global const FloatType3 *vertices,
        global const int3 *faces,
        global const FloatType3 *planeUnitNormals,
        global const FloatType3 *segmentVectors,
        global const FloatType3 *segmentUnitNormals,
        global FloatType16 *results,
        const int numberOfFaces,
        const FloatType pointX,
        const FloatType pointY,
        const FloatType pointZ,
        local FloatType *scratch) {

    const FloatType3 point = {pointX, pointY, pointZ};

    // Out-of-range work-items must NOT return early: reduceOverWorkGroup() below contains barriers
    // which every work-item of the work-group has to reach, otherwise the behaviour is undefined and
    // implementations deadlock on the work-items that already left. Instead the out-of-range items
    // redundantly evaluate face 0 and contribute zero to the reduction.
    const int globalIndex = get_global_id(0);
    const bool inRange = globalIndex < numberOfFaces;
    const int faceIndex = inRange ? globalIndex : 0;

    const FloatType3 planeUnitNormal = planeUnitNormals[faceIndex];

    // The face's vertices, relocated so that the computation point P sits at the origin
    const FloatType3 face[3] = {
            vertices[faces[faceIndex][0]] - point,
            vertices[faces[faceIndex][1]] - point,
            vertices[faces[faceIndex][2]] - point,
    };

    //1-04 Step: Compute Plane Normal Orientation sigma_p (direction of N_p in relation to P)
    const int planeNormalOrientation = sgn(dot(planeUnitNormal, face[0]));

    //1-05 Step: Compute Hessian Normal Plane Representation
    FloatType4 hessianPlane;
    {
        const FloatType3 origin = {0.0, 0.0, 0.0};
        const FloatType3 crossProduct = cross(face[0] - face[1], face[0] - face[2]);
        const FloatType3 res = (origin - face[0]) * crossProduct;

        hessianPlane.xyz = crossProduct;
        hessianPlane.w = res[0] + res[1] + res[2];
    }

    //1-06 Step: Compute distance h_p between P and P'
    const FloatType planeDistance = fabs(hessianPlane.w / length(hessianPlane.xyz));

    //1-07 Step: Compute the actual position of P' (projection of P on the plane)
    FloatType3 orthogonalProjectionPointOnPlane = planeUnitNormal * planeDistance;
    {
        const FloatType3 intersections = {
                hessianPlane.x == 0.0 ? 0.0 : hessianPlane.w / hessianPlane.x,
                hessianPlane.y == 0.0 ? 0.0 : hessianPlane.w / hessianPlane.y,
                hessianPlane.z == 0.0 ? 0.0 : hessianPlane.w / hessianPlane.z,
        };

        for (uint index = 0; index < 3; ++index) {
            if (intersections[index] < 0) {
                orthogonalProjectionPointOnPlane[index] = fabs(orthogonalProjectionPointOnPlane[index]);
            } else if (orthogonalProjectionPointOnPlane[index] > 0) {
                orthogonalProjectionPointOnPlane[index] = -orthogonalProjectionPointOnPlane[index];
            }
        }
    }

    //1-08 Step: Compute the segment normal orientation sigma_pq (direction of n_pq in relation to P')
    int3 segmentNormalOrientations;
    for (uint index = 0; index < 3; ++index) {
        const FloatType inner = dot(segmentUnitNormals[faceIndex * 3 + index],
                                    orthogonalProjectionPointOnPlane - face[index]);
        segmentNormalOrientations[index] = -sgn(inner);
    }

    //1-09 Step: Compute the orthogonal projection point P'' of P' on each segment
    FloatType3 orthogonalProjectionPointsOnSegmentsForPlane[3];
    for (uint index = 0; index < 3; ++index) {
        if (segmentNormalOrientations[index] == 0) {
            orthogonalProjectionPointsOnSegmentsForPlane[index] = orthogonalProjectionPointOnPlane;
        } else {
            const FloatType3 vertex1 = face[index];
            const FloatType3 vertex2 = face[(index + 1) % 3];

            const FloatType3 matrixRow1 = vertex2 - vertex1;
            const FloatType3 matrixRow2 = cross(vertex1 - orthogonalProjectionPointOnPlane, matrixRow1);
            const FloatType3 matrixRow3 = cross(matrixRow2, matrixRow1);

            const FloatType3 d = {
                    dot(matrixRow1, orthogonalProjectionPointOnPlane),
                    dot(matrixRow2, orthogonalProjectionPointOnPlane),
                    dot(matrixRow3, vertex1),
            };

            FloatType3 columnMatrix[3] = {matrixRow1, matrixRow2, matrixRow3};
            transpose(columnMatrix);

            const FloatType determinant = det(columnMatrix);

            if (determinant != 0.0) {
                const FloatType3 r = {
                        detOfRows(d, columnMatrix[1], columnMatrix[2]),
                        detOfRows(columnMatrix[0], d, columnMatrix[2]),
                        detOfRows(columnMatrix[0], columnMatrix[1], d),
                };
                orthogonalProjectionPointsOnSegmentsForPlane[index] = r / determinant;
            }
        }
    }

    //1-10 Step: Compute the segment distances h_pq between P'' and P'
    FloatType3 segmentDistances;
    for (uint index = 0; index < 3; ++index) {
        segmentDistances[index] = length(orthogonalProjectionPointsOnSegmentsForPlane[index] -
                                         orthogonalProjectionPointOnPlane);
    }

    //1-11 Step: Compute the 3D distances l1, l2 (between P and vertices)
    // and the 1D distances s1, s2 (between P'' and vertices)
    Distance distances[3];
    for (uint index = 0; index < 3; ++index) {
        distances[index].l1 = length(face[index]);
        distances[index].l2 = length(face[(index + 1) % 3]);

        distances[index].s1 = length(orthogonalProjectionPointsOnSegmentsForPlane[index] - face[index]);
        distances[index].s2 = length(orthogonalProjectionPointsOnSegmentsForPlane[index] - face[(index + 1) % 3]);

        if (fabs(distances[index].s1 - distances[index].l1) < EPSILON_ZERO_OFFSET &&
            fabs(distances[index].s2 - distances[index].l2) < EPSILON_ZERO_OFFSET) {
            if (distances[index].s2 < distances[index].s1) {
                distances[index].s1 *= -1.0;
                distances[index].s2 *= -1.0;
                distances[index].l1 *= -1.0;
                distances[index].l2 *= -1.0;
            } else if (fabs(distances[index].s2 - distances[index].s1) < EPSILON_ZERO_OFFSET) {
                distances[index].s1 *= -1.0;
                distances[index].l1 *= -1.0;
            }
        } else {
            const FloatType norm = length(segmentVectors[faceIndex * 3 + index]);
            if (distances[index].s1 < norm && distances[index].s2 < norm) {
                distances[index].s1 *= -1.0;
            } else if (distances[index].s2 < distances[index].s1) {
                distances[index].s1 *= -1.0;
                distances[index].s2 *= -1.0;
            }
        }
    }

    //1-12 Step: Compute the euclidian norms of the vectors consisting of P' and the vertices,
    // they are later used for determining the position of P' in relation to the plane
    const FloatType3 projectionPointVertexNorms = {
            length(orthogonalProjectionPointOnPlane - face[0]),
            length(orthogonalProjectionPointOnPlane - face[1]),
            length(orthogonalProjectionPointOnPlane - face[2]),
    };

    //1-13 Step: Compute the transcendental expressions LN_pq and AN_pq
    TranscendentalExpression transcendentalExpressions[3];
    for (uint index = 0; index < 3; ++index) {
        const FloatType r1Norm = projectionPointVertexNorms[(index + 1) % 3];
        const FloatType r2Norm = projectionPointVertexNorms[index];

        if ((segmentNormalOrientations[index] == 0 &&
             (r1Norm < EPSILON_ZERO_OFFSET || r2Norm < EPSILON_ZERO_OFFSET)) ||
            (fabs(distances[index].s1 + distances[index].s2) < EPSILON_ZERO_OFFSET &&
             fabs(distances[index].l1 + distances[index].l2) < EPSILON_ZERO_OFFSET)) {
            transcendentalExpressions[index].ln = 0.0;
        } else {
            const FloatType numerator = distances[index].s2 + distances[index].l2;
            const FloatType denominator = distances[index].s1 + distances[index].l1;

            transcendentalExpressions[index].ln = numerator <= 0.0 || denominator <= 0.0
                                                          ? 0.0
                                                          : log(numerator / denominator);
        }

        if (planeDistance < EPSILON_ZERO_OFFSET || segmentDistances[index] < EPSILON_ZERO_OFFSET) {
            transcendentalExpressions[index].an = 0.0;
        } else {
            const FloatType fraction1 = (planeDistance * distances[index].s2) /
                                        (segmentDistances[index] * distances[index].l2);
            const FloatType fraction2 = (planeDistance * distances[index].s1) /
                                        (segmentDistances[index] * distances[index].l1);

            transcendentalExpressions[index].an = atan(fraction1) - atan(fraction2);
        }
    }

    //1-14 Step: Compute the singularities sing A and sing B if P' is located in the plane,
    // on any vertex, or on one segment (G_pq)
    const FloatType singularityAngle = computeSingularityAngle(faceIndex, segmentNormalOrientations,
                                                               projectionPointVertexNorms, segmentVectors);
    const FloatType singularityAlpha = -planeDistance * singularityAngle;
    const FloatType3 singularityBeta = planeUnitNormal * ((FloatType) -1.0 * singularityAngle * planeNormalOrientation);

    //2. Step: Compute Sum 1 used for potential and acceleration, sum over: sigma_pq * h_pq * LN_pq
    FloatType sum1PotentialAcceleration = 0.0;
    for (uint index = 0; index < 3; ++index) {
        sum1PotentialAcceleration += segmentNormalOrientations[index] * segmentDistances[index] *
                                     transcendentalExpressions[index].ln;
    }

    //3. Step: Compute Sum 1 used for the gradiometric tensor, sum over: n_pq * LN_pq
    FloatType3 sum1Tensor = {0.0, 0.0, 0.0};
    for (uint index = 0; index < 3; ++index) {
        sum1Tensor += segmentUnitNormals[faceIndex * 3 + index] * transcendentalExpressions[index].ln;
    }

    //4. Step: Compute Sum 2 which is the same for every result parameter, sum over: sigma_pq * AN_pq
    FloatType sum2 = 0.0;
    for (uint index = 0; index < 3; ++index) {
        sum2 += segmentNormalOrientations[index] * transcendentalExpressions[index].an;
    }

    //5. Step: Sum for potential and acceleration, consisting of: sum1 + h_p * sum2 + sing A
    const FloatType planeSumPotentialAcceleration = sum1PotentialAcceleration + planeDistance * sum2 + singularityAlpha;

    //6. Step: Sum for the tensor, consisting of: sum1 + sigma_p * N_p * sum2 + sing B
    const FloatType3 subSum = (sum1Tensor + (planeUnitNormal * (planeNormalOrientation * sum2))) + singularityBeta;
    // first component: trivial case Vxx, Vyy, Vzz --> just N_p * subSum
    const FloatType3 first = planeUnitNormal * subSum;
    // second component: reordering required to build Vxy, Vxz, Vyz
    const FloatType3 reorderedNp = {planeUnitNormal[0], planeUnitNormal[0], planeUnitNormal[1]};
    const FloatType3 reorderedSubSum = {subSum[1], subSum[2], subSum[2]};
    const FloatType3 second = reorderedNp * reorderedSubSum;

    //7. Step: Multiply with the prefix and store this face's contribution.
    // Components 0-2 hold the acceleration, component 3 the potential, and components 4-9 the tensor.
    FloatType16 value = (FloatType16) (0.0);
    if (inRange) {
        value.xyz = planeUnitNormal * planeSumPotentialAcceleration;
        value.w = planeNormalOrientation * planeDistance * planeSumPotentialAcceleration;
        value.s456 = first;
        value.s789 = second;
    }

    //8. Step: Sum this work-group's faces so that only one partial result per work-group leaves the device
    const FloatType16 result = reduceOverWorkGroup(value, scratch);

    if (get_local_id(0) == 0) {
        results[get_group_id(0)] = result;
    }
}
