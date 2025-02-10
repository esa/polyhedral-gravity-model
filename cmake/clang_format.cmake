file(GLOB_RECURSE CLANG_FORMAT_SRC
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.h"
        "${PROJECT_SOURCE_DIR}/test/*.cpp"
        "${PROJECT_SOURCE_DIR}/test/*.h"
)

# Define a variable for clang-format command
find_program(CLANG_FORMAT clang-format)

# Ensure clang-format was found
if(NOT CLANG_FORMAT)
    message(STATUS "clang-format not found. Please install it to use clang-format via CMake")
else()
    message(STATUS "clang-format found. You can format all source files via `cmake --build . --target format`")
    add_custom_command(
            OUTPUT format_all_files
            COMMAND ${CLANG_FORMAT} -i ${CLANG_FORMAT_SRC}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Formatting all source and test files with clang-format"
            VERBATIM
    )

    add_custom_target(format
            DEPENDS format_all_files
    )
endif()
