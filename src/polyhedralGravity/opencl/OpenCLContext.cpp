#include "polyhedralGravity/opencl/OpenCLContext.h"

#include "polyhedralGravity/output/Logging.h"
#include <sstream>
#include <stdexcept>
#include <vector>

namespace polyhedralGravity::opencl {

    namespace {

        /**
         * Returns every device of every platform, ordered so that GPUs come first.
         * Platforms which expose no device of the requested type are skipped rather than treated as an
         * error: systems with several runtimes installed (e.g. Intel's oneAPI, which advertises a CPU
         * runtime, an FPGA emulation, and a GPU runtime as separate platforms) list them in an
         * arbitrary order, so the GPU is frequently not on the first platform.
         *
         * @param type the device type to collect
         * @return the collected devices
         */
        std::vector<cl::Device> collectDevices(const cl_device_type type) {
            std::vector<cl::Platform> platforms{};
            cl::Platform::get(&platforms);

            std::vector<cl::Device> devices{};
            for (const cl::Platform &platform: platforms) {
                std::vector<cl::Device> platformDevices{};
                try {
                    platform.getDevices(type, &platformDevices);
                } catch (const cl::Error &) {
                    // CL_DEVICE_NOT_FOUND merely means this platform has no device of that type
                    continue;
                }
                devices.insert(devices.end(), platformDevices.begin(), platformDevices.end());
            }
            return devices;
        }

        /**
         * Selects the device to evaluate on, preferring a GPU over any other device type.
         * @param precision the required precision
         * @return the selected device
         * @throws std::runtime_error if no device supports the precision
         */
        cl::Device selectDevice(const ComputePrecision precision) {
            // A GPU is preferred; only if none supports the precision any other device type
            // (i.e. the host CPU, where its vendor ships an OpenCL runtime) is considered
            for (const cl_device_type type: {cl_device_type{CL_DEVICE_TYPE_GPU}, cl_device_type{CL_DEVICE_TYPE_ALL}}) {
                for (const cl::Device &device: collectDevices(type)) {
                    if (OpenCLContext::supportsPrecision(device, precision)) {
                        return device;
                    }
                }
            }
            std::stringstream message{};
            message << "No OpenCL device supporting " << precision << " was found on this machine.\n"
                    << OpenCLContext::describeAvailableDevices(precision);
            throw std::runtime_error{message.str()};
        }

    }// namespace

    OpenCLContext::OpenCLContext(const ComputePrecision precision)
        : _device{selectDevice(precision)},
          _context{_device},
          _queue{_context, _device},
          _precision{precision} {
        POLYHEDRAL_GRAVITY_LOG_INFO("OpenCL backend uses the device '{}' in {} precision",
                                    getDeviceName(), _precision == ComputePrecision::FLOAT64 ? "float64" : "float32");
    }

    bool OpenCLContext::isAvailable(const ComputePrecision precision) noexcept {
        try {
            for (const cl::Device &device: collectDevices(CL_DEVICE_TYPE_ALL)) {
                if (supportsPrecision(device, precision)) {
                    return true;
                }
            }
        } catch (...) {
            // A machine without any OpenCL runtime makes the ICD loader fail rather than report zero
            // platforms, which is not an error here but simply means "no device available".
            return false;
        }
        return false;
    }

    bool OpenCLContext::supportsPrecision(const cl::Device &device, const ComputePrecision precision) {
        if (precision != ComputePrecision::FLOAT64) {
            return true;
        }
        // Both indicators are consulted because neither is reliable on its own: PoCL executes double
        // precision correctly on arm64 macOS while reporting a zeroed CL_DEVICE_DOUBLE_FP_CONFIG, and
        // conversely a device may list the extension without a usable configuration.
        cl_device_fp_config doubleConfig{0};
        if (device.getInfo(CL_DEVICE_DOUBLE_FP_CONFIG, &doubleConfig) == CL_SUCCESS && doubleConfig != 0) {
            return true;
        }
        std::string extensions{};
        if (device.getInfo(CL_DEVICE_EXTENSIONS, &extensions) != CL_SUCCESS) {
            return false;
        }
        return extensions.find("cl_khr_fp64") != std::string::npos;
    }

    const cl::Device &OpenCLContext::getDevice() const {
        return _device;
    }

    const cl::Context &OpenCLContext::getContext() const {
        return _context;
    }

    cl::CommandQueue &OpenCLContext::getQueue() const {
        return _queue;
    }

    ComputePrecision OpenCLContext::getPrecision() const {
        return _precision;
    }

    std::string OpenCLContext::getDeviceName() const {
        return _device.getInfo<CL_DEVICE_NAME>();
    }

    size_t OpenCLContext::getLocalMemorySize() const {
        return static_cast<size_t>(_device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>());
    }

    std::string OpenCLContext::getBuildOptions() const {
        // The floating point type is a build-time parameter of the kernels so that one and the same
        // kernel source serves both precisions, see GravityModelCommon.cl
        if (_precision == ComputePrecision::FLOAT64) {
            return "-cl-std=CL1.2 -D POLYHEDRAL_GRAVITY_DOUBLE_PRECISION"
                   " -D FloatType=double -D FloatType3=double3 -D FloatType4=double4 -D FloatType16=double16";
        }
        return "-cl-std=CL1.2"
               " -D FloatType=float -D FloatType3=float3 -D FloatType4=float4 -D FloatType16=float16";
    }

    std::string OpenCLContext::describeAvailableDevices(const ComputePrecision precision) {
        std::stringstream listing{};
        listing << "Detected OpenCL devices:\n";
        std::vector<cl::Platform> platforms{};
        try {
            cl::Platform::get(&platforms);
        } catch (const cl::Error &e) {
            listing << "  <querying the OpenCL platforms failed: " << e.what() << ">\n";
            return listing.str();
        }
        if (platforms.empty()) {
            listing << "  <none, this machine exposes no OpenCL platform>\n";
        }
        for (const cl::Platform &platform: platforms) {
            listing << "  Platform '" << platform.getInfo<CL_PLATFORM_NAME>() << "' ("
                    << platform.getInfo<CL_PLATFORM_VERSION>() << ")\n";
            std::vector<cl::Device> devices{};
            try {
                platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
            } catch (const cl::Error &) {
                listing << "    <no devices>\n";
                continue;
            }
            for (const cl::Device &device: devices) {
                listing << "    Device '" << device.getInfo<CL_DEVICE_NAME>() << "' ("
                        << device.getInfo<CL_DEVICE_VERSION>() << "): "
                        << (supportsPrecision(device, precision) ? "supports" : "does NOT support")
                        << " " << precision << "\n";
            }
        }
        return listing.str();
    }

}// namespace polyhedralGravity::opencl
