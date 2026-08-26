#include "polyhedralGravity/opencl/OpenCLEvaluation.h"

#include "polyhedralGravity/output/Logging.h"
#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>

// Generated at build time from the sources in kernel/, see compile_opencl()
#include "GravityModelEvaluationKernel.h"
#include "GravityModelInitKernel.h"
#include "GravityModelReductionKernel.h"

namespace polyhedralGravity::opencl {

    namespace {

        /**
         * Converts the host's coordinates into the device's floating point format.
         * @tparam FloatType either cl_float or cl_double
         * @param vertices the vertices
         * @return the vertices in the device's layout
         */
        template<typename FloatType>
        std::vector<DeviceVector3<FloatType>> toDeviceVectors(const std::vector<Array3> &vertices) {
            std::vector<DeviceVector3<FloatType>> deviceVertices{};
            deviceVertices.reserve(vertices.size());
            for (const Array3 &vertex: vertices) {
                deviceVertices.push_back({static_cast<FloatType>(vertex[0]),
                                          static_cast<FloatType>(vertex[1]),
                                          static_cast<FloatType>(vertex[2]),
                                          FloatType{0}});
            }
            return deviceVertices;
        }

        /**
         * Converts the faces' vertex indices into the device's int3 layout.
         * @param faces the faces
         * @return the faces in the device's layout
         */
        std::vector<DeviceIndex3> toDeviceIndices(const std::vector<IndexArray3> &faces) {
            std::vector<DeviceIndex3> deviceFaces{};
            deviceFaces.reserve(faces.size());
            for (const IndexArray3 &face: faces) {
                deviceFaces.push_back({static_cast<cl_int>(face[0]),
                                       static_cast<cl_int>(face[1]),
                                       static_cast<cl_int>(face[2]),
                                       0});
            }
            return deviceFaces;
        }

        /**
         * Creates a device buffer holding a copy of the given host data.
         * @tparam T the element type
         * @param context the context to allocate in
         * @param data the host data
         * @return the buffer
         */
        template<typename T>
        cl::Buffer createBufferFrom(const cl::Context &context, std::vector<T> &data) {
            return cl::Buffer{context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, data.size() * sizeof(T), data.data()};
        }

        /**
         * Sums the partial results read back from the device and applies the model's prefix.
         *
         * @tparam FloatType either cl_float or cl_double
         * @param staging the raw bytes read back from the device
         * @param count the number of partial results the staging buffer holds
         * @param prefix GRAVITATIONAL_CONSTANT * density, corrected for orientation and mesh unit
         * @return the potential, the acceleration, and the second derivative tensor
         */
        template<typename FloatType>
        GravityModelResult accumulate(const std::vector<std::byte> &staging, const size_t count, const double prefix) {
            const auto *partialResults = reinterpret_cast<const DeviceVector16<FloatType> *>(staging.data());

            GravityModelResult result{};
            auto &[potential, acceleration, gradiometricTensor] = result;
            for (size_t index = 0; index < count; ++index) {
                const FloatType *components = partialResults[index].components;
                potential += static_cast<double>(components[RESULT_POTENTIAL_INDEX]);
                for (size_t component = 0; component < acceleration.size(); ++component) {
                    acceleration[component] += static_cast<double>(components[RESULT_ACCELERATION_OFFSET + component]);
                }
                for (size_t component = 0; component < gradiometricTensor.size(); ++component) {
                    gradiometricTensor[component] += static_cast<double>(components[RESULT_TENSOR_OFFSET + component]);
                }
            }

            // Final expressions after application of the prefix (and a division by 2 for the potential),
            // matching step 9 and 10 of GravityEvaluable::evaluate
            potential = (potential * prefix) / 2.0;
            for (double &component: acceleration) {
                component *= -1.0 * prefix;
            }
            for (double &component: gradiometricTensor) {
                component *= prefix;
            }
            return result;
        }

    }// namespace

    OpenCLEvaluation::OpenCLEvaluation(const Polyhedron &polyhedron, const ComputePrecision precision)
        : _context{precision},
          _initializationKernel{_context, KERNEL_INIT, "initializeFaceProperties"},
          _evaluationKernel{_context, KERNEL_EVALUATION, "evaluateFaces"},
          _reductionKernel{_context, KERNEL_REDUCTION, "reducePartialResults"},
          _numberOfFaces{polyhedron.countFaces()},
          _usesReductionKernel{false},
          _resultCount{0},
          _prefix{polyhedron.getGravityModelScaling()} {
        // The kernels address vertices and faces with a signed 32 bit index
        constexpr size_t maximumIndex = static_cast<size_t>(std::numeric_limits<cl_int>::max());
        if (polyhedron.countVertices() > maximumIndex || _numberOfFaces > maximumIndex) {
            throw std::runtime_error{"The polyhedron exceeds the number of vertices or faces the OpenCL "
                                     "backend can address. Use ComputeBackend::CPU for this polyhedron."};
        }
        if (_numberOfFaces == 0) {
            throw std::runtime_error{"The OpenCL backend cannot evaluate a polyhedron without faces."};
        }

        const cl::Context &context = _context.getContext();

        std::vector<DeviceIndex3> faces = toDeviceIndices(polyhedron.getFaces());
        _faces = createBufferFrom(context, faces);
        if (precision == ComputePrecision::FLOAT64) {
            auto vertices = toDeviceVectors<cl_double>(polyhedron.getVertices());
            _vertices = createBufferFrom(context, vertices);
        } else {
            auto vertices = toDeviceVectors<cl_float>(polyhedron.getVertices());
            _vertices = createBufferFrom(context, vertices);
        }

        // The point-independent face properties are derived on the device and stay resident there
        const size_t vectorSize = vector3Size(precision);
        _planeUnitNormals = cl::Buffer{context, CL_MEM_READ_WRITE, _numberOfFaces * vectorSize};
        _segmentVectors = cl::Buffer{context, CL_MEM_READ_WRITE, _numberOfFaces * vectorSize * 3};
        _segmentUnitNormals = cl::Buffer{context, CL_MEM_READ_WRITE, _numberOfFaces * vectorSize * 3};

        const size_t partialResultCount = _evaluationKernel.countWorkGroups(_numberOfFaces);
        _usesReductionKernel = partialResultCount > REDUCTION_THRESHOLD;
        _resultCount = _usesReductionKernel ? _reductionKernel.countWorkGroups(partialResultCount) : partialResultCount;

        _partialResults = cl::Buffer{context, CL_MEM_READ_WRITE, partialResultCount * vector16Size(precision)};
        // A zero-sized buffer is invalid per the OpenCL specification, so allocate at least one element
        // even when the second reduction stage is unused
        _reducedResults = cl::Buffer{context, CL_MEM_READ_WRITE,
                                     std::max(_usesReductionKernel ? _resultCount : size_t{0}, size_t{1}) *
                                             vector16Size(precision)};
        _resultStaging.resize(_resultCount * vector16Size(precision));

        cl::Kernel &initializationKernel = _initializationKernel.get();
        initializationKernel.setArg(0, _vertices);
        initializationKernel.setArg(1, _faces);
        initializationKernel.setArg(2, _planeUnitNormals);
        initializationKernel.setArg(3, _segmentVectors);
        initializationKernel.setArg(4, _segmentUnitNormals);
        initializationKernel.setArg(5, static_cast<cl_int>(_numberOfFaces));

        cl::Kernel &evaluationKernel = _evaluationKernel.get();
        evaluationKernel.setArg(0, _vertices);
        evaluationKernel.setArg(1, _faces);
        evaluationKernel.setArg(2, _planeUnitNormals);
        evaluationKernel.setArg(3, _segmentVectors);
        evaluationKernel.setArg(4, _segmentUnitNormals);
        evaluationKernel.setArg(5, _partialResults);
        evaluationKernel.setArg(6, static_cast<cl_int>(_numberOfFaces));
        evaluationKernel.setArg(10, cl::Local(_evaluationKernel.getLocalScratchSize()));

        cl::Kernel &reductionKernel = _reductionKernel.get();
        reductionKernel.setArg(0, _partialResults);
        reductionKernel.setArg(1, _reducedResults);
        reductionKernel.setArg(2, static_cast<cl_int>(partialResultCount));
        reductionKernel.setArg(3, cl::Local(_reductionKernel.getLocalScratchSize()));

        // The initialization kernel contains no barrier, so the implementation may partition it freely
        _initializationKernel.enqueueWithoutWorkGroups(_numberOfFaces);
        _context.getQueue().finish();

        POLYHEDRAL_GRAVITY_LOG_DEBUG("OpenCL backend prepared {} faces on '{}', producing {} partial "
                                     "result(s) per evaluation (second reduction stage {})",
                                     _numberOfFaces, getDeviceName(), _resultCount,
                                     _usesReductionKernel ? "enabled" : "disabled");
    }

    GravityModelResult OpenCLEvaluation::evaluate(const Array3 &computationPoint) const {
        const std::lock_guard<std::mutex> lock{_mutex};

        cl::Kernel &evaluationKernel = _evaluationKernel.get();
        if (_context.getPrecision() == ComputePrecision::FLOAT64) {
            evaluationKernel.setArg(7, static_cast<cl_double>(computationPoint[0]));
            evaluationKernel.setArg(8, static_cast<cl_double>(computationPoint[1]));
            evaluationKernel.setArg(9, static_cast<cl_double>(computationPoint[2]));
        } else {
            evaluationKernel.setArg(7, static_cast<cl_float>(computationPoint[0]));
            evaluationKernel.setArg(8, static_cast<cl_float>(computationPoint[1]));
            evaluationKernel.setArg(9, static_cast<cl_float>(computationPoint[2]));
        }

        _evaluationKernel.enqueue(_numberOfFaces);

        cl::CommandQueue &queue = _context.getQueue();
        if (_usesReductionKernel) {
            _reductionKernel.enqueue(_evaluationKernel.countWorkGroups(_numberOfFaces));
            queue.enqueueReadBuffer(_reducedResults, CL_TRUE, 0, _resultStaging.size(), _resultStaging.data());
        } else {
            queue.enqueueReadBuffer(_partialResults, CL_TRUE, 0, _resultStaging.size(), _resultStaging.data());
        }
        queue.finish();

        return _context.getPrecision() == ComputePrecision::FLOAT64
                       ? accumulate<cl_double>(_resultStaging, _resultCount, _prefix)
                       : accumulate<cl_float>(_resultStaging, _resultCount, _prefix);
    }

    std::vector<GravityModelResult> OpenCLEvaluation::evaluate(const std::vector<Array3> &computationPoints) const {
        std::vector<GravityModelResult> results{};
        results.reserve(computationPoints.size());
        for (const Array3 &computationPoint: computationPoints) {
            results.push_back(this->evaluate(computationPoint));
        }
        return results;
    }

    std::string OpenCLEvaluation::getDeviceName() const {
        return _context.getDeviceName();
    }

    ComputePrecision OpenCLEvaluation::getPrecision() const {
        return _context.getPrecision();
    }

}// namespace polyhedralGravity::opencl
