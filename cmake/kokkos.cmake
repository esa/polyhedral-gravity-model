message(STATUS "Setting up Kokkos")
set(POLYHEDRAL_GRAVITY_KOKKOS_VERSION 5.1.1)

include(CheckLanguage)

# The Serial backend implements ComputeBackend::CPU_SERIAL and is therefore always required
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable the Kokkos Serial backend" FORCE)

##########################################################
# Host Backend: OpenMP implements ComputeBackend::CPU_PARALLEL
##########################################################
# Apple's Clang needs to be pointed at the libomp of Homebrew/ MacPorts as it does not ship one itself
if (APPLE AND NOT DEFINED OpenMP_ROOT)
    foreach (HINT IN ITEMS /opt/homebrew/opt/libomp /usr/local/opt/libomp /opt/local)
        if (EXISTS "${HINT}/include/omp.h")
            set(OpenMP_ROOT "${HINT}")
            message(STATUS "Kokkos: Found a libomp installation at ${OpenMP_ROOT}")
            break()
        endif ()
    endforeach ()
endif ()

find_package(OpenMP QUIET COMPONENTS CXX)
if (OpenMP_CXX_FOUND)
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable the Kokkos OpenMP backend" FORCE)
else ()
    message(WARNING "Kokkos: No OpenMP installation was found. ComputeBackend::CPU_PARALLEL will fall back to the "
            "Serial backend, i.e., it will not parallelize. Install libomp (macOS) or libomp-dev (Linux) to fix this.")
endif ()

##################################################################
# Device Backend: The native vendor paradigm of the available GPU
##################################################################
set(POLYHEDRAL_GRAVITY_DEVICE_BACKEND "AUTO" CACHE STRING "The Kokkos device backend used for ComputeBackend::GPU_PARALLEL
 (AUTO = detect the vendor paradigm, NONE = CPU-only build, or one of CUDA, HIP, SYCL)")
set_property(CACHE POLYHEDRAL_GRAVITY_DEVICE_BACKEND PROPERTY STRINGS AUTO NONE CUDA HIP SYCL)

if (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "AUTO")
    # NVIDIA GPUs are programmed via CUDA, AMD GPUs via HIP, and Intel GPUs via SYCL.
    # Whichever toolchain is installed determines which GPU is actually attached to this machine.
    check_language(CUDA)
    check_language(HIP)
    if (CMAKE_CUDA_COMPILER)
        set(POLYHEDRAL_GRAVITY_DEVICE_BACKEND "CUDA")
    elseif (CMAKE_HIP_COMPILER)
        set(POLYHEDRAL_GRAVITY_DEVICE_BACKEND "HIP")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        set(POLYHEDRAL_GRAVITY_DEVICE_BACKEND "SYCL")
    else ()
        set(POLYHEDRAL_GRAVITY_DEVICE_BACKEND "NONE")
    endif ()
endif ()

if (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "CUDA")
    enable_language(CUDA)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable the Kokkos CUDA backend" FORCE)
elseif (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "HIP")
    enable_language(HIP)
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable the Kokkos HIP backend" FORCE)
elseif (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "SYCL")
    set(Kokkos_ENABLE_SYCL ON CACHE BOOL "Enable the Kokkos SYCL backend" FORCE)
endif ()

message(STATUS "Kokkos: Host Backend   ${CMAKE_CXX_COMPILER_ID} (Serial + OpenMP: ${OpenMP_CXX_FOUND})")
message(STATUS "Kokkos: Device Backend ${POLYHEDRAL_GRAVITY_DEVICE_BACKEND}")

##########################################
# Get Kokkos itself, preferring a local one
##########################################
find_package(Kokkos 4.6.02...${POLYHEDRAL_GRAVITY_KOKKOS_VERSION} QUIET)

if (Kokkos_FOUND)
    message(STATUS "Found existing Kokkos libraries: ${Kokkos_DIR}")
else ()
    message(STATUS "Using Kokkos from GitHub Release ${POLYHEDRAL_GRAVITY_KOKKOS_VERSION}")
    include(FetchContent)

    # For the CPU code always optimize for the machine being built on (use vectorization, etc.)
    set(Kokkos_ARCH_NATIVE ON CACHE BOOL "Always build for the machine on which is being compiled" FORCE)

    FetchContent_Declare(
            Kokkos
            URL https://github.com/kokkos/kokkos/archive/refs/tags/${POLYHEDRAL_GRAVITY_KOKKOS_VERSION}.tar.gz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(Kokkos)
endif ()
