# Locates the OpenCL runtime and the OpenCL C++ bindings (<CL/opencl.hpp>) and provides
# the compile_opencl() function which embeds *.cl kernel sources into generated C++ headers.
find_package(OpenCL REQUIRED)

# The OpenCL C++ bindings are a separate Khronos package and are frequently not shipped
# alongside the runtime: Apple's OpenCL.framework contains the C headers only (and under the
# <OpenCL/...> prefix rather than <CL/...>), and several Linux distributions package the ICD
# loader without the clhpp headers. Hence, we search for the C and C++ headers separately and
# inject whatever we find into the imported target.
find_path(OpenCL_CXX_INCLUDE_DIR
        NAMES CL/opencl.hpp
        HINTS
        ${OpenCL_INCLUDE_DIRS}
        /opt/homebrew/opt/opencl-clhpp-headers/include
        /usr/local/opt/opencl-clhpp-headers/include
        /opt/homebrew/include
        /usr/local/include
        /opt/local/include
        DOC "Directory containing the OpenCL C++ bindings (CL/opencl.hpp)"
)

find_path(OpenCL_C_INCLUDE_DIR
        NAMES CL/cl.h
        HINTS
        ${OpenCL_INCLUDE_DIRS}
        /opt/homebrew/opt/opencl-headers/include
        /usr/local/opt/opencl-headers/include
        /opt/homebrew/include
        /usr/local/include
        /opt/local/include
        DOC "Directory containing the OpenCL C headers (CL/cl.h)"
)

if (NOT OpenCL_CXX_INCLUDE_DIR)
    message(FATAL_ERROR
            "Could not find the OpenCL C++ bindings (CL/opencl.hpp). They are required by "
            "POLYHEDRAL_GRAVITY_ENABLE_OPENCL. Install them via e.g. 'brew install opencl-clhpp-headers' (macOS), "
            "'apt install opencl-clhpp-headers' (Debian/Ubuntu), or point OpenCL_CXX_INCLUDE_DIR at them.")
endif ()
if (NOT OpenCL_C_INCLUDE_DIR)
    message(FATAL_ERROR
            "Could not find the OpenCL C headers (CL/cl.h) which CL/opencl.hpp includes. "
            "Install them via e.g. 'brew install opencl-headers' (macOS), 'apt install opencl-headers' "
            "(Debian/Ubuntu), or point OpenCL_C_INCLUDE_DIR at them.")
endif ()

message(STATUS "Found OpenCL C headers at:      ${OpenCL_C_INCLUDE_DIR}")
message(STATUS "Found OpenCL C++ bindings at:   ${OpenCL_CXX_INCLUDE_DIR}")

target_include_directories(OpenCL::OpenCL INTERFACE "${OpenCL_C_INCLUDE_DIR}" "${OpenCL_CXX_INCLUDE_DIR}")

include(compile_opencl)

#############################################################
# Embedding the OpenCL kernels into generated C++ headers
#############################################################
# The kernels are shipped inside the binary rather than as loose *.cl files so that the library stays
# relocatable and the Python wheel needs no data files.
set(POLYHEDRAL_GRAVITY_OPENCL_KERNEL_DIR "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/opencl/kernel")
set(POLYHEDRAL_GRAVITY_OPENCL_GENERATED_DIR "${PROJECT_BINARY_DIR}/generated/opencl")

set(POLYHEDRAL_GRAVITY_OPENCL_HEADERS "")
foreach (KERNEL IN ITEMS Common Init Evaluation Reduction)
    string(TOUPPER ${KERNEL} KERNEL_UPPERCASE)
    compile_opencl(
            INFILE "${POLYHEDRAL_GRAVITY_OPENCL_KERNEL_DIR}/GravityModel${KERNEL}.cl"
            OUTFILE "GravityModel${KERNEL}Kernel.h"
            PATH "${POLYHEDRAL_GRAVITY_OPENCL_GENERATED_DIR}"
            VAR_NAME "KERNEL_${KERNEL_UPPERCASE}"
            NAMESPACE "polyhedralGravity::opencl"
    )
    list(APPEND POLYHEDRAL_GRAVITY_OPENCL_HEADERS
            "${POLYHEDRAL_GRAVITY_OPENCL_GENERATED_DIR}/GravityModel${KERNEL}Kernel.h")
endforeach ()

# Every target compiling the OpenCL sources depends on this target so that the headers exist first
add_custom_target(polyhedral_gravity_opencl_kernels DEPENDS ${POLYHEDRAL_GRAVITY_OPENCL_HEADERS})

# Carries everything a target needs to compile and link the OpenCL backend
add_library(polyhedral_gravity_opencl INTERFACE)
target_link_libraries(polyhedral_gravity_opencl INTERFACE OpenCL::OpenCL)
target_include_directories(polyhedral_gravity_opencl INTERFACE "${POLYHEDRAL_GRAVITY_OPENCL_GENERATED_DIR}")
target_compile_definitions(polyhedral_gravity_opencl INTERFACE POLYHEDRAL_GRAVITY_ENABLE_OPENCL)
