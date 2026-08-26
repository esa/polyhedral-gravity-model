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

# Kokkos compiles the device code as CXX, using nvcc_wrapper/ hipcc/ a CUDA-capable Clang, and passes the
# architecture flags itself. So we must NOT enable_language(CUDA/HIP) here: that would put the language into
# ENABLED_LANGUAGES with an empty CMAKE_CUDA_ARCHITECTURES, which makes Kokkos' own architecture
# auto-detection fail to even generate its try_run project ("CUDA_ARCHITECTURES is empty for target ...").
if (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "CUDA")
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable the Kokkos CUDA backend" FORCE)
    # The kernels call constexpr host functions (e.g. std::array::operator[]) from device code. Kokkos
    # defaults this to ON for Clang but to OFF for nvcc, where it becomes -expt-relaxed-constexpr.
    set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE BOOL "Allow constexpr host functions in device code" FORCE)
    # Kokkos detects the compute capability by briefly enabling the CUDA language itself and building a
    # probe. That probe needs a non-empty architecture list, so give it one unless the user picked their own.
    # "all-major" is used rather than "native" because it does not require a visible GPU at configure time.
    if (NOT CMAKE_CUDA_ARCHITECTURES AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.23)
        set(CMAKE_CUDA_ARCHITECTURES "all-major")
    endif ()
elseif (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "HIP")
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable the Kokkos HIP backend" FORCE)
elseif (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "SYCL")
    set(Kokkos_ENABLE_SYCL ON CACHE BOOL "Enable the Kokkos SYCL backend" FORCE)
endif ()

##################################################################
# Which compilers actually translate the host and the device code
##################################################################
set(POLYHEDRAL_GRAVITY_HOST_COMPILER "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "NONE")
    set(POLYHEDRAL_GRAVITY_DEVICE_COMPILER "None")
elseif (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "CUDA" AND CMAKE_CUDA_COMPILER)
    # Kokkos compiles the device code as CXX, so the CUDA language is never enabled (see below) and CMake
    # never identifies nvcc itself. Asking the binary is the only way to get a version out of it.
    execute_process(COMMAND "${CMAKE_CUDA_COMPILER}" --version OUTPUT_VARIABLE NVCC_VERSION_OUTPUT ERROR_QUIET)
    if (NVCC_VERSION_OUTPUT MATCHES "release ([0-9]+\\.[0-9]+)")
        set(POLYHEDRAL_GRAVITY_DEVICE_COMPILER "NVIDIA nvcc ${CMAKE_MATCH_1}")
    else ()
        set(POLYHEDRAL_GRAVITY_DEVICE_COMPILER "NVIDIA nvcc")
    endif ()
elseif (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "HIP" AND CMAKE_HIP_COMPILER)
    execute_process(COMMAND "${CMAKE_HIP_COMPILER}" --version OUTPUT_VARIABLE HIPCC_VERSION_OUTPUT ERROR_QUIET)
    if (HIPCC_VERSION_OUTPUT MATCHES "HIP version: ([0-9][^\n]*)")
        set(POLYHEDRAL_GRAVITY_DEVICE_COMPILER "AMD hipcc ${CMAKE_MATCH_1}")
    else ()
        set(POLYHEDRAL_GRAVITY_DEVICE_COMPILER "AMD hipcc")
    endif ()
else ()
    # SYCL, and a CUDA/ HIP build whose device code the C++ compiler understands by itself
    set(POLYHEDRAL_GRAVITY_DEVICE_COMPILER "${POLYHEDRAL_GRAVITY_DEVICE_BACKEND} via ${POLYHEDRAL_GRAVITY_HOST_COMPILER}")
endif ()

message(STATUS "Kokkos: Host Backend   ${CMAKE_CXX_COMPILER_ID} (Serial + OpenMP: ${OpenMP_CXX_FOUND})")
message(STATUS "Kokkos: Device Backend ${POLYHEDRAL_GRAVITY_DEVICE_BACKEND}")

# The device code is compiled by the C++ compiler, so it has to be one that understands the device paradigm
if (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "CUDA" AND
        NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND
        NOT CMAKE_CXX_COMPILER_ID STREQUAL "NVIDIA" AND
        NOT CMAKE_CXX_COMPILER MATCHES "nvcc_wrapper")
    message(WARNING "Kokkos: The CUDA backend compiles the device code with the C++ compiler "
            "(${CMAKE_CXX_COMPILER_ID}), which Kokkos only supports for Clang and NVIDIA's nvcc_wrapper. "
            "If the build fails, configure with CXX pointing to nvcc_wrapper, or set "
            "POLYHEDRAL_GRAVITY_DEVICE_BACKEND=NONE for a CPU-only build.")
endif ()

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
