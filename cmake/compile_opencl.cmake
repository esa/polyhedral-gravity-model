# - Provides a function to embed an OpenCL kernel (*.cl) file into a C++ header file (*.h)
# Usage:
#   compile_opencl(INFILE <input_file>
#                  OUTFILE <header_name.h>
#                  PATH <output_dir>
#                  VAR_NAME <name_of_variable_in_header>
#                  NAMESPACE <namespace_of_variable>
#   )
#  VAR_NAME (defaults to "KERNEL_SOURCE") and NAMESPACE (defaults to "polyhedralGravity::opencl")
#  are optional arguments.
# Produces a file ${PATH}/${OUTFILE} with the content:
#   #pragma once
#   namespace <NAMESPACE> {
#       inline constexpr const char <VAR_NAME>[] = R"( ... file contents ... )";
#   }
function(compile_opencl)
    set(options)
    set(oneValueArgs INFILE OUTFILE PATH VAR_NAME NAMESPACE)
    set(multiValueArgs)
    cmake_parse_arguments(COMPILE_OCL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT COMPILE_OCL_INFILE)
        message(FATAL_ERROR "compile_opencl: INFILE is required")
    endif ()
    if (NOT COMPILE_OCL_OUTFILE)
        message(FATAL_ERROR "compile_opencl: OUTFILE is required")
    endif ()
    if (NOT COMPILE_OCL_PATH)
        message(FATAL_ERROR "compile_opencl: PATH is required")
    endif ()

    if (NOT DEFINED COMPILE_OCL_VAR_NAME OR COMPILE_OCL_VAR_NAME STREQUAL "")
        set(COMPILE_OCL_VAR_NAME "KERNEL_SOURCE")
    endif ()
    if (NOT DEFINED COMPILE_OCL_NAMESPACE OR COMPILE_OCL_NAMESPACE STREQUAL "")
        set(COMPILE_OCL_NAMESPACE "polyhedralGravity::opencl")
    endif ()

    set(OUTPUT_PATH "${COMPILE_OCL_PATH}/${COMPILE_OCL_OUTFILE}")

    # Make the input path absolute
    if (NOT IS_ABSOLUTE "${COMPILE_OCL_INFILE}")
        set(INPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${COMPILE_OCL_INFILE}")
    else ()
        set(INPUT_PATH "${COMPILE_OCL_INFILE}")
    endif ()

    # Make the output directory at configure time
    file(MAKE_DIRECTORY "${COMPILE_OCL_PATH}")

    # Create a custom command that regenerates the header whenever the kernel source changes
    add_custom_command(
            OUTPUT "${OUTPUT_PATH}"
            COMMAND ${CMAKE_COMMAND}
            -D "INPUT_FILE=${INPUT_PATH}"
            -D "OUTPUT_FILE=${OUTPUT_PATH}"
            -D "VAR_NAME=${COMPILE_OCL_VAR_NAME}"
            -D "NAMESPACE=${COMPILE_OCL_NAMESPACE}"
            -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/scripts/opencl_to_header.cmake"
            # The generator script is a dependency too, so that changing how the kernels are embedded
            # regenerates the headers even when the kernel sources themselves are untouched
            DEPENDS "${INPUT_PATH}" "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/scripts/opencl_to_header.cmake"
            COMMENT "Generating OpenCL kernel header ${COMPILE_OCL_OUTFILE}"
            VERBATIM
    )
endfunction()
