##################################################################################
# Fast Math: trade the last digits of the FLOAT32 evaluation for a faster kernel
##################################################################################
# The polyhedral gravity model divides and takes square roots several dozen times per face and computation
# point. Compiled for strict IEEE-754 semantics, each of those becomes a hardware approximation followed by a
# correction step and a call into a slow path, which is why the option below exists.
#
# It only ever affects single precision: neither nvcc's -use_fast_math nor Clang's -ffast-math has a fast path
# for double precision arithmetic. ComputePrecision::FLOAT64 therefore produces bit-identical results with and
# without it, while ComputePrecision::FLOAT32 -- which is documented to be good for about 1e-4 relative anyway
# -- becomes noticeably faster. Measured for the Eros mesh (14744 faces) on an RTX 4060, the GPU evaluation of
# 2000 points goes from 2.83 to 1.97 microseconds per point, i.e. it gains about 30%.
option(POLYHEDRAL_GRAVITY_FAST_MATH
        "Compile the FLOAT32 evaluation with fast, less accurate division, square root, and transcendentals (Default: OFF)"
        OFF)

if (POLYHEDRAL_GRAVITY_FAST_MATH)
    # Which flag is right is decided by whoever compiles the *device* code, not by CMAKE_CXX_COMPILER_ID.
    # With the CUDA backend, Kokkos substitutes its own nvcc_wrapper as the compile rule while the cache
    # entry still names the host compiler, so CMAKE_CXX_COMPILER_ID reads "GNU" for a build whose kernels
    # nvcc translates. Keying on it would hand nvcc_wrapper a -ffast-math, which it forwards to the host
    # compiler, leaving the FLOAT32 kernels -- the only thing this option is about -- untouched.
    if (POLYHEDRAL_GRAVITY_DEVICE_BACKEND STREQUAL "CUDA" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # nvcc, reached either through nvcc_wrapper or directly. Applies to the device code only.
        set(POLYHEDRAL_GRAVITY_FAST_MATH_FLAG "-use_fast_math")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|IntelLLVM|GNU")
        # Clang compiling CUDA itself, the HIP and SYCL backends, and a CPU-only build, where this is the
        # only way to reach the FLOAT32 evaluation at all. Unlike nvcc's flag, this one also applies to the
        # host code of the targets it is attached to.
        set(POLYHEDRAL_GRAVITY_FAST_MATH_FLAG "-ffast-math")
    else ()
        # Never downgrade silently: whoever asked for this has to learn that they did not get it
        message(FATAL_ERROR "POLYHEDRAL_GRAVITY_FAST_MATH is ON, but this project does not know a fast math "
                "flag for ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER}) with device backend "
                "${POLYHEDRAL_GRAVITY_DEVICE_BACKEND}. Configure with -DPOLYHEDRAL_GRAVITY_FAST_MATH=OFF, or "
                "add the compiler to cmake/fast_math.cmake.")
    endif ()
    # Deliberately NOT add_compile_options(): that would reach the bundled dependencies too, and TetGen's
    # exact geometric predicates depend on the IEEE semantics this flag relaxes. The flag is attached to the
    # project's own targets only, in src/CMakeLists.txt.
    message(STATUS "Fast Math: compiling the polyhedral gravity targets with ${POLYHEDRAL_GRAVITY_FAST_MATH_FLAG}")
endif ()

# Spelled as a C++ bool literal for Info.h
if (POLYHEDRAL_GRAVITY_FAST_MATH)
    set(POLYHEDRAL_GRAVITY_FAST_MATH_LITERAL "true")
else ()
    set(POLYHEDRAL_GRAVITY_FAST_MATH_LITERAL "false")
endif ()
