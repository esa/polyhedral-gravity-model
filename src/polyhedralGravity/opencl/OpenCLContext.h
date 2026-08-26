#pragma once

#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/opencl/OpenCLDefinitions.h"
#include <cstddef>
#include <string>

namespace polyhedralGravity::opencl {

    /**
     * The OpenCL device the polyhedral gravity model is evaluated on, together with its context
     * and command queue.
     *
     * The device is picked automatically: the first GPU capable of the requested
     * {@link ComputePrecision} wins, and only if no such GPU exists any other OpenCL device
     * (i.e. the host CPU, if its vendor ships an OpenCL runtime) is used.
     */
    class OpenCLContext {

        /** The selected device */
        cl::Device _device;

        /** The context owning the device */
        cl::Context _context;

        /** The in-order command queue of the device */
        mutable cl::CommandQueue _queue;

        /** The precision this context was created for */
        ComputePrecision _precision;

    public:
        /**
         * Selects a device supporting the given precision and creates a context and command queue for it.
         * @param precision the floating point precision the kernels are compiled in
         * @throws std::runtime_error if no OpenCL device supporting the precision exists
         */
        explicit OpenCLContext(ComputePrecision precision);

        /**
         * Checks whether an OpenCL device supporting the given precision exists on this machine.
         * In contrast to the constructor, this never throws and never propagates an OpenCL error, which
         * makes it suitable for deciding whether to fall back to {@link ComputeBackend::CPU}.
         * @param precision the floating point precision
         * @return true if a suitable device exists
         */
        static bool isAvailable(ComputePrecision precision) noexcept;

        /**
         * Returns whether a device is able to compute in the given precision.
         * FLOAT64 requires the device to support the @c cl_khr_fp64 extension.
         *
         * @param device the device to check
         * @param precision the floating point precision
         * @return true if the device supports the precision
         *
         * @note We query CL_DEVICE_DOUBLE_FP_CONFIG rather than trying to build a double precision
         * kernel: Apple's OpenCL implementation reports success from @c clBuildProgram for a double
         * precision kernel and only fails later in @c clCreateKernel, so a successful build is no
         * proof that double precision actually works.
         */
        static bool supportsPrecision(const cl::Device &device, ComputePrecision precision);

        /**
         * Returns the device this context runs on.
         * @return the device
         */
        [[nodiscard]] const cl::Device &getDevice() const;

        /**
         * Returns the OpenCL context.
         * @return the context
         */
        [[nodiscard]] const cl::Context &getContext() const;

        /**
         * Returns the command queue of the device.
         * @return the command queue
         */
        [[nodiscard]] cl::CommandQueue &getQueue() const;

        /**
         * Returns the precision the kernels of this context are compiled in.
         * @return the precision
         */
        [[nodiscard]] ComputePrecision getPrecision() const;

        /**
         * Returns the name of the device, e.g. "Apple M1 Pro".
         * @return the device name
         */
        [[nodiscard]] std::string getDeviceName() const;

        /**
         * Returns the amount of local (i.e. work-group shared) memory of the device in bytes.
         * @return the local memory size
         */
        [[nodiscard]] size_t getLocalMemorySize() const;

        /**
         * Returns the build options the kernels have to be compiled with, which inject the
         * floating point type into the kernel sources.
         * @return the build options
         */
        [[nodiscard]] std::string getBuildOptions() const;

        /**
         * Returns a human-readable listing of every OpenCL platform and device found on this machine,
         * used to explain why no suitable device could be selected.
         * @param precision the precision the devices are judged against
         * @return a multi-line listing
         */
        static std::string describeAvailableDevices(ComputePrecision precision);
    };

}// namespace polyhedralGravity::opencl
