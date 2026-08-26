#pragma once

/**
 * We target the OpenCL 1.2 API. Apple's implementation -- the only one available on macOS -- never
 * advanced beyond it, and the polyhedral gravity kernels require no feature introduced afterwards.
 * (Notably, the work-group reduction is written by hand rather than using OpenCL 2.0's
 * work_group_reduce_add(), see GravityModelCommon.cl.)
 */
#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/opencl.hpp>

#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include <cstddef>

namespace polyhedralGravity::opencl {

    /**
     * The number of meaningful components of the per-face result vector:
     * the potential (1), the acceleration (3), and the gradiometric tensor (6).
     * @note Must match RESULT_COMPONENTS in GravityModelCommon.cl
     */
    constexpr size_t RESULT_COMPONENTS = 10;

    /** The index of the acceleration's first component inside the result vector */
    constexpr size_t RESULT_ACCELERATION_OFFSET = 0;

    /** The index of the potential inside the result vector */
    constexpr size_t RESULT_POTENTIAL_INDEX = 3;

    /** The index of the gradiometric tensor's first component inside the result vector */
    constexpr size_t RESULT_TENSOR_OFFSET = 4;

    /**
     * The host-side counterpart of OpenCL's @c float3 / @c double3.
     * OpenCL pads a three-component vector to four components, hence the explicit padding member.
     * @tparam FloatType either @c cl_float or @c cl_double
     */
    template<typename FloatType>
    struct DeviceVector3 {
        FloatType x;
        FloatType y;
        FloatType z;
        FloatType padding;
    };

    /**
     * The host-side counterpart of OpenCL's @c float16 / @c double16, used for the per-work-group
     * results. Only the leading {@link RESULT_COMPONENTS} components carry data.
     * @tparam FloatType either @c cl_float or @c cl_double
     */
    template<typename FloatType>
    struct DeviceVector16 {
        FloatType components[16];
    };

    /** The host-side counterpart of OpenCL's @c int3, holding the vertex indices of one face. */
    struct DeviceIndex3 {
        cl_int x;
        cl_int y;
        cl_int z;
        cl_int padding;
    };

    static_assert(sizeof(DeviceVector3<cl_float>) == sizeof(cl_float3), "DeviceVector3 must match cl_float3");
    static_assert(sizeof(DeviceVector3<cl_double>) == sizeof(cl_double3), "DeviceVector3 must match cl_double3");
    static_assert(sizeof(DeviceVector16<cl_float>) == sizeof(cl_float16), "DeviceVector16 must match cl_float16");
    static_assert(sizeof(DeviceVector16<cl_double>) == sizeof(cl_double16), "DeviceVector16 must match cl_double16");
    static_assert(sizeof(DeviceIndex3) == sizeof(cl_int3), "DeviceIndex3 must match cl_int3");

    /**
     * The size in bytes one scalar occupies on the device for the given precision.
     * @param precision the precision
     * @return 4 for FLOAT32, 8 for FLOAT64
     */
    constexpr size_t scalarSize(const ComputePrecision precision) {
        return precision == ComputePrecision::FLOAT64 ? sizeof(cl_double) : sizeof(cl_float);
    }

    /**
     * The size in bytes one three-component vector occupies on the device for the given precision.
     * @param precision the precision
     * @return the size of @c cl_float3 respectively @c cl_double3
     */
    constexpr size_t vector3Size(const ComputePrecision precision) {
        return precision == ComputePrecision::FLOAT64 ? sizeof(cl_double3) : sizeof(cl_float3);
    }

    /**
     * The size in bytes one sixteen-component vector occupies on the device for the given precision.
     * @param precision the precision
     * @return the size of @c cl_float16 respectively @c cl_double16
     */
    constexpr size_t vector16Size(const ComputePrecision precision) {
        return precision == ComputePrecision::FLOAT64 ? sizeof(cl_double16) : sizeof(cl_float16);
    }

}// namespace polyhedralGravity::opencl
