#pragma once

#include "polyhedralGravity/opencl/OpenCLContext.h"
#include <cstddef>
#include <string>

namespace polyhedralGravity::opencl {

    /**
     * A compiled OpenCL kernel together with the work-group size it may actually be launched with.
     *
     * Bundling the two is what makes the kernel usable without further ceremony: the nominal
     * {@link DESIRED_WORK_GROUP_SIZE} is only achievable where the kernel is small enough for it, and
     * exceeding either the kernel's or the device's limit makes @c clEnqueueNDRangeKernel fail with
     * @c CL_INVALID_WORK_GROUP_SIZE. The evaluation kernel in particular is register hungry enough
     * that devices report a much smaller CL_KERNEL_WORK_GROUP_SIZE than CL_DEVICE_MAX_WORK_GROUP_SIZE.
     */
    class OpenCLKernel {

        /**
         * The work-group size aimed for. Whether it is reached depends on the device and the kernel,
         * see {@link getWorkGroupSize}.
         */
        static constexpr size_t DESIRED_WORK_GROUP_SIZE = 256;

        /** The context the kernel is compiled for and enqueued on */
        const OpenCLContext &_context;

        /** The built program */
        cl::Program _program;

        /** The kernel itself */
        mutable cl::Kernel _kernel;

        /** The largest legal work-group size for this kernel on this device, always a power of two */
        size_t _workGroupSize;

    public:
        /**
         * Compiles the given kernel source for the context's device and looks up the named kernel.
         *
         * @param context the context to compile for
         * @param source the kernel source, to which the shared definitions are prepended
         * @param name the name of the kernel function inside the source
         *
         * @throws std::runtime_error carrying the build log if the source does not compile
         */
        OpenCLKernel(const OpenCLContext &context, const std::string &source, const std::string &name);

        /**
         * Returns the underlying kernel, e.g. to set its arguments.
         * @return the kernel
         */
        [[nodiscard]] cl::Kernel &get() const;

        /**
         * The largest work-group size this kernel can be launched with on this device.
         * It is a power of two, which the hand-written work-group reduction relies on and which also
         * keeps it a clean multiple of the sub-group width.
         * @return the work-group size
         */
        [[nodiscard]] size_t getWorkGroupSize() const;

        /**
         * The size in bytes of the local scratch buffer the work-group reduction needs.
         * @return the scratch size, to be passed as a @c cl::Local kernel argument
         */
        [[nodiscard]] size_t getLocalScratchSize() const;

        /**
         * The number of work-groups the given number of work items is spread over.
         * @param workItems the number of work items
         * @return the number of work-groups, i.e. the number of partial results the kernel emits
         */
        [[nodiscard]] size_t countWorkGroups(size_t workItems) const;

        /**
         * Enqueues the kernel for the given number of work items on the context's queue.
         * The global size is padded to a multiple of the work-group size, so the kernel has to guard
         * against work items beyond {@link workItems}.
         * @param workItems the number of work items
         */
        void enqueue(size_t workItems) const;

        /**
         * Enqueues the kernel for the given number of work items without specifying a work-group size,
         * leaving the partitioning to the implementation. Only valid for kernels which contain no
         * barrier and no work-group wide reduction.
         * @param workItems the number of work items
         */
        void enqueueWithoutWorkGroups(size_t workItems) const;
    };

}// namespace polyhedralGravity::opencl
