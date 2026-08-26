/**
 * Definitions shared by all polyhedral gravity OpenCL kernels.
 *
 * This source is not a standalone program. It is prepended to every kernel source by the host
 * (see OpenCLProgram), which is why it declares no kernel of its own.
 *
 * The floating point type is injected through the build options (-D FloatType=double, ...) so that
 * one and the same kernel source can be compiled in single and in double precision.
 */

#ifdef POLYHEDRAL_GRAVITY_DOUBLE_PRECISION
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

/**
 * The number of meaningful components of the per-face result vector:
 * the potential (1), the acceleration (3), and the gradiometric tensor (6).
 * The result is stored in a 16-component vector since OpenCL has no 10-component vector type;
 * the remaining components are padding and never read.
 */
#define RESULT_COMPONENTS 10

#define EPSILON_ZERO_OFFSET ((FloatType) 1e-14)
#define PI    ((FloatType) 3.1415926535897932384626433832795028841971693993751058209749445923)
#define PI2   ((FloatType) 6.2831853071795864769252867665590057683943387987502116419498891846)
#define PI_2  ((FloatType) 1.5707963267948966192313216916397514420985846996875529104874722961)

/**
 * Sums the RESULT_COMPONENTS leading components of value across the entire work-group.
 * The total is returned to the work-item with local id 0; what the other work-items receive is undefined.
 *
 * @param value the per-work-item contribution
 * @param scratch a local buffer of at least RESULT_COMPONENTS * get_local_size(0) elements
 * @return the work-group's sum in the work-item with local id 0
 *
 * @note This deliberately avoids OpenCL 2.0's work_group_reduce_add(). Apple's implementation -- the
 * only one available on macOS -- is stuck at OpenCL 1.2 and does not provide the work-group collectives.
 * @note This function contains barriers, so every work-item of the work-group has to reach it,
 * the out-of-range ones included. They contribute a zero vector instead of returning early.
 * @note The reduction assumes get_local_size(0) to be a power of two, which the host guarantees.
 */
FloatType16 reduceOverWorkGroup(FloatType16 value, local FloatType *scratch) {
    const uint localId = get_local_id(0);
    const uint localSize = get_local_size(0);

    for (uint component = 0; component < RESULT_COMPONENTS; ++component) {
        scratch[component * localSize + localId] = value[component];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (uint offset = localSize / 2; offset > 0; offset >>= 1) {
        if (localId < offset) {
            for (uint component = 0; component < RESULT_COMPONENTS; ++component) {
                scratch[component * localSize + localId] += scratch[component * localSize + localId + offset];
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    FloatType16 result = (FloatType16) (0.0);
    for (uint component = 0; component < RESULT_COMPONENTS; ++component) {
        result[component] = scratch[component * localSize];
    }
    return result;
}
