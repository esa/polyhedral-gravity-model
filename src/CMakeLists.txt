###############################
# Building the Polyhedral Libray
################################
if (BUILD_POLYHEDRAL_GRAVITY_LIBRARY)
    message(STATUS "Building the Polyhedral Gravity Library")

    # YAML is only required for the library target, but not for e.g. the Python interface
    include(yaml)

    add_subdirectory(polyhedralGravity)

    # Adds the include Path PROJECT_SOURCE_DIR/src to the target polyhedralGravity
    target_include_directories(${PROJECT_NAME}_lib PUBLIC
            "${CMAKE_CURRENT_SOURCE_DIR}"
            )

    # Link libraries needed in all targets
    target_link_libraries(${PROJECT_NAME}_lib
            spdlog::spdlog
            yaml-cpp
            tetgen
            xsimd
            Thrust
            )
endif ()
#####################################
# Building the Polyhedral Executable
#####################################
if (BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE)
    message(STATUS "Building the Polyhedral Gravity Executable")

    # Building the standalone Executable
    add_executable(${PROJECT_NAME} main.cpp)

    # Link executable with library
    target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib)

    # Place the executable in the top level directory
    set_target_properties(${PROJECT_NAME} PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            )
endif()

##########################################
# Building the Polyhedral Python Interface
##########################################
if(BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE)
    message(STATUS "Building the Polyhedral Gravity Python Interface")

    # Include pybind11
    include(pybind11)

    add_subdirectory(polyhedralGravityPython)

    target_include_directories(polyhedral_gravity PUBLIC
            "${CMAKE_CURRENT_SOURCE_DIR}"
            )

    target_link_libraries(polyhedral_gravity PUBLIC
            spdlog::spdlog
            tetgen
            xsimd
            Thrust
            )
endif()