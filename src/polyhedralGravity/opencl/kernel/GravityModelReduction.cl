/**
 * Second stage of the summation over the polyhedral faces.
 *
 * The evaluation kernel already reduces within each of its work-groups and emits one partial result
 * per work-group. For meshes with many faces the number of those partial results is still large enough
 * that summing them on the host dominates the runtime, so this kernel reduces them once more.
 */
kernel void reducePartialResults(
        global const FloatType16 *partialResults,
        global FloatType16 *reducedResults,
        const int numberOfPartialResults,
        local FloatType *scratch) {

    const int index = get_global_id(0);

    // Out-of-range work-items must not return early, see reduceOverWorkGroup()
    const FloatType16 value = index < numberOfPartialResults ? partialResults[index] : (FloatType16) (0.0);

    const FloatType16 result = reduceOverWorkGroup(value, scratch);

    if (get_local_id(0) == 0) {
        reducedResults[get_group_id(0)] = result;
    }
}
