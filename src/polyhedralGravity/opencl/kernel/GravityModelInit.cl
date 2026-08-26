/**
 * Computes the properties of every polyhedral face which are independent of the computation point P:
 * the segment vectors G_pq, the plane unit normals N_p, and the segment unit normals n_pq.
 *
 * This mirrors GravityEvaluable::prepare() on the host and is run exactly once per polyhedron,
 * which is why the results stay resident on the device across evaluations.
 */
kernel void initializeFaceProperties(
        global const FloatType3 *vertices,
        global const int3 *faces,
        global FloatType3 *planeUnitNormals,
        global FloatType3 *segmentVectors,
        global FloatType3 *segmentUnitNormals,
        const int numberOfFaces) {

    const int index = get_global_id(0);

    // This kernel contains no barrier, so out-of-range work-items may return early
    if (index >= numberOfFaces) {
        return;
    }

    const FloatType3 face[3] = {
            vertices[faces[index].x],
            vertices[faces[index].y],
            vertices[faces[index].z],
    };

    //1-01 Step: Compute Segment Vectors G_pq which describe each one the edge between two vertices
    const FloatType3 segmentVector[3] = {
            face[1] - face[0],
            face[2] - face[1],
            face[0] - face[2],
    };

    //1-02 Step: Compute the Plane Unit Normal N_p (pointing outside the polyhedron)
    const FloatType3 planeUnitNormal = normalize(cross(segmentVector[0], segmentVector[1]));
    planeUnitNormals[index] = planeUnitNormal;

    segmentVectors[index * 3 + 0] = segmentVector[0];
    segmentVectors[index * 3 + 1] = segmentVector[1];
    segmentVectors[index * 3 + 2] = segmentVector[2];

    //1-03 Step: Compute Segment Unit Normals n_pq (normal pointing away from each segment)
    segmentUnitNormals[index * 3 + 0] = normalize(cross(segmentVector[0], planeUnitNormal));
    segmentUnitNormals[index * 3 + 1] = normalize(cross(segmentVector[1], planeUnitNormal));
    segmentUnitNormals[index * 3 + 2] = normalize(cross(segmentVector[2], planeUnitNormal));
}
