#include "polyhedralGravity/opencl/OpenCLKernel.h"

#include "polyhedralGravity/output/Logging.h"
#include <algorithm>
#include <sstream>
#include <stdexcept>

// Generated at build time from kernel/GravityModelCommon.cl, see compile_opencl()
#include "GravityModelCommonKernel.h"

namespace polyhedralGravity::opencl {

    namespace {

        /**
         * Rounds a size down to the next power of two.
         * @param size the size, may be zero
         * @return the largest power of two less than or equal to size, at least one
         */
        size_t roundDownToPowerOfTwo(const size_t size) {
            size_t powerOfTwo{1};
            while (powerOfTwo * 2 <= size) {
                powerOfTwo *= 2;
            }
            return powerOfTwo;
        }

    }// namespace

    OpenCLKernel::OpenCLKernel(const OpenCLContext &context, const std::string &source, const std::string &name)
        : _context{context},
          // The shared definitions are prepended rather than #include'd since an OpenCL program has no
          // include path unless the host provides one
          _program{context.getContext(), std::string{KERNEL_COMMON} + source},
          _kernel{},
          _workGroupSize{1} {
        const cl::Device &device = context.getDevice();
        try {
            _program.build({device}, context.getBuildOptions().c_str());
        } catch (const cl::Error &e) {
            std::stringstream message{};
            message << "Building the OpenCL kernel '" << name << "' for the device '"
                    << context.getDeviceName() << "' failed (" << e.what() << "):\n"
                    << _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
            throw std::runtime_error{message.str()};
        }
        _kernel = cl::Kernel{_program, name.c_str()};

        // Every constraint below can only shrink the work-group size, so we simply take the minimum:
        // what the kernel permits, what the device permits, and what fits into local memory.
        const size_t bytesPerWorkItem = RESULT_COMPONENTS * scalarSize(context.getPrecision());
        const auto kernelStaticLocalMemory = static_cast<size_t>(_kernel.getWorkGroupInfo<CL_KERNEL_LOCAL_MEM_SIZE>(device));
        const size_t localMemory = context.getLocalMemorySize() > kernelStaticLocalMemory
                                           ? context.getLocalMemorySize() - kernelStaticLocalMemory
                                           : 0;

        const size_t limit = std::min({DESIRED_WORK_GROUP_SIZE,
                                       _kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device),
                                       static_cast<size_t>(device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()),
                                       localMemory / bytesPerWorkItem});
        _workGroupSize = roundDownToPowerOfTwo(limit);

        POLYHEDRAL_GRAVITY_LOG_DEBUG("OpenCL kernel '{}' compiled with a work-group size of {} "
                                     "(kernel limit {}, device limit {}, local memory limit {})",
                                     name, _workGroupSize,
                                     _kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device),
                                     device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>(),
                                     localMemory / bytesPerWorkItem);
    }

    cl::Kernel &OpenCLKernel::get() const {
        return _kernel;
    }

    size_t OpenCLKernel::getWorkGroupSize() const {
        return _workGroupSize;
    }

    size_t OpenCLKernel::getLocalScratchSize() const {
        return RESULT_COMPONENTS * scalarSize(_context.getPrecision()) * _workGroupSize;
    }

    size_t OpenCLKernel::countWorkGroups(const size_t workItems) const {
        return (workItems + _workGroupSize - 1) / _workGroupSize;
    }

    void OpenCLKernel::enqueue(const size_t workItems) const {
        const size_t globalSize = countWorkGroups(workItems) * _workGroupSize;
        _context.getQueue().enqueueNDRangeKernel(_kernel, cl::NullRange,
                                                 cl::NDRange{globalSize}, cl::NDRange{_workGroupSize});
    }

    void OpenCLKernel::enqueueWithoutWorkGroups(const size_t workItems) const {
        _context.getQueue().enqueueNDRangeKernel(_kernel, cl::NullRange, cl::NDRange{workItems}, cl::NullRange);
    }

}// namespace polyhedralGravity::opencl
