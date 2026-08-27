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
    if (CMAKE_CXX_COMPILER MATCHES "nvcc_wrapper" OR CMAKE_CXX_COMPILER_ID STREQUAL "NVIDIA")
        # nvcc applies this to the device code only, the host compiler behind it never sees it
        add_compile_options(-use_fast_math)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
        # Covers a Clang CUDA build as well as the HIP and SYCL backends. Unlike nvcc's flag, this one also
        # applies to the host code compiled from the same translation unit.
        add_compile_options(-ffast-math)
    else ()
        message(WARNING "POLYHEDRAL_GRAVITY_FAST_MATH is ON, but ${CMAKE_CXX_COMPILER_ID} is not one of the "
                "compilers this project knows a fast math flag for. It has no effect.")
        set(POLYHEDRAL_GRAVITY_FAST_MATH OFF)
    endif ()
endif ()

# Spelled as a C++ bool literal for Info.h
if (POLYHEDRAL_GRAVITY_FAST_MATH)
    set(POLYHEDRAL_GRAVITY_FAST_MATH_LITERAL "true")
else ()
    set(POLYHEDRAL_GRAVITY_FAST_MATH_LITERAL "false")
endif ()
