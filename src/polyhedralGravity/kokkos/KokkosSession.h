#pragma once

#include <string>
#include <type_traits>

#include <Kokkos_Core.hpp>

#include "polyhedralGravity/model/PolyhedronDefinitions.h"

/**
 * Everything which directly touches Kokkos, i.e. the runtime session, the device-resident polyhedron,
 * and the kernels evaluating Tsoulis' polyhedral gravity model.
 */
namespace polyhedralGravity::kokkos {

    /**
     * The execution space behind {@link ComputeBackend::CPU_SERIAL}.
     * The Serial backend is unconditionally enabled by the CMake configuration.
     */
    using SerialSpace = Kokkos::Serial;

    /**
     * The execution space behind {@link ComputeBackend::CPU_PARALLEL}.
     * This is OpenMP whenever an OpenMP installation was found while configuring, otherwise Kokkos'
     * default host space, which then degrades to Serial (i.e., no parallelization on the host).
     */
#ifdef KOKKOS_ENABLE_OPENMP
    using HostParallelSpace = Kokkos::OpenMP;
#else
    using HostParallelSpace = Kokkos::DefaultHostExecutionSpace;
#endif

    /**
     * The execution space behind {@link ComputeBackend::GPU_PARALLEL}, i.e. CUDA, HIP, or SYCL.
     * If the library was built without a device backend, this is a host space and
     * {@link GPU_AVAILABLE} is false.
     */
    using DeviceSpace = Kokkos::DefaultExecutionSpace;

    /**
     * Whether the library was compiled with a GPU backend.
     * Kokkos only makes its default execution space differ from the default host execution space if
     * one of the device backends (CUDA, HIP, SYCL) is enabled.
     */
    constexpr bool GPU_AVAILABLE = !std::is_same_v<DeviceSpace, Kokkos::DefaultHostExecutionSpace>;

    /**
     * Initializes the Kokkos runtime if it is not initialized yet.
     *
     * This is called by everything that touches a Kokkos::View, so a user of this library never has to
     * think about Kokkos' lifecycle. The runtime is finalized again when the process exits, which happens
     * strictly after all Views of this library have been released.
     *
     * @note If the embedding application initializes Kokkos itself, this function does nothing and, in
     * particular, does not finalize the application's session.
     */
    void ensureInitialized();

    /**
     * Returns the Kokkos execution spaces this library was compiled with as a human-readable list,
     * for example @code "Serial, OpenMP, Cuda" @endcode.
     * @return the comma-separated names of the enabled execution spaces
     */
    std::string getEnabledExecutionSpaces();

    /**
     * Returns the name of the execution space a given compute backend maps to,
     * for example @code "OpenMP" @endcode for {@link ComputeBackend::CPU_PARALLEL}.
     * @param backend the compute backend
     * @return the name of the corresponding Kokkos execution space
     */
    std::string getExecutionSpaceName(ComputeBackend backend);

    /**
     * Throws if the given compute backend cannot be served by this build.
     * The only backend which can be unavailable is {@link ComputeBackend::GPU_PARALLEL}, namely if the
     * library was compiled without a CUDA, HIP, or SYCL backend.
     * @param backend the requested compute backend
     * @throws std::runtime_error if the backend is not available
     */
    void checkBackendAvailable(ComputeBackend backend);

}// namespace polyhedralGravity::kokkos
