find_program(CMAKE_FORMAT_EXECUTABLE NAMES cmake-format)

if (NOT CMAKE_FORMAT_EXECUTABLE)
    message(STATUS "cmake-format not found. Please install it to use cmake-format via CMake")
else ()
    message(STATUS "cmake-format found. You can format all source files via `cmake --build . --target format_cmake`")

    # Create a list of all .cmake files in the project
    file(GLOB_RECURSE CMAKE_FORMAT_FILES
            "${PROJECT_SOURCE_DIR}/*.cmake"
            "${PROJECT_SOURCE_DIR}/CMakeLists.txt"
    )

    # Add a custom target to format all .cmake files
    add_custom_command(
            OUTPUT format_cmake_files
            COMMAND ${CMAKE_FORMAT_EXECUTABLE} --in-place ${CMAKE_FORMAT_FILES}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Formatting all cmake files in the project with cmake-format"
            VERBATIM
    )
    add_custom_target(format_cmake
            DEPENDS format_cmake_files
    )
endif ()