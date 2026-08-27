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

    # Fast, less accurate FLOAT32 arithmetic, if requested (see cmake/fast_math.cmake). Attached per target so
    # that the bundled dependencies, TetGen's exact predicates above all, keep their IEEE semantics.
    if (POLYHEDRAL_GRAVITY_FAST_MATH)
        target_compile_options(${PROJECT_NAME}_lib PRIVATE ${POLYHEDRAL_GRAVITY_FAST_MATH_FLAG})
    endif ()

    # Link libraries needed in all targets
    target_link_libraries(${PROJECT_NAME}_lib
            spdlog::spdlog
            yaml-cpp::yaml-cpp
            tetgen_lib
            Kokkos::kokkos
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

###########################################
# Building the Polyhedral Benchmark Driver
###########################################
if (BUILD_POLYHEDRAL_GRAVITY_BENCHMARK)
    message(STATUS "Building the Polyhedral Gravity Benchmark Driver")

    add_executable(${PROJECT_NAME}_benchmark ${PROJECT_SOURCE_DIR}/script/benchmark.cpp)
    target_link_libraries(${PROJECT_NAME}_benchmark PUBLIC ${PROJECT_NAME}_lib)
    set_target_properties(${PROJECT_NAME}_benchmark PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
            )
endif ()

##########################################
# Building the Polyhedral Python Interface
##########################################
if(BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE)
    message(STATUS "Building the Polyhedral Gravity Python Interface")

    # Include pybind11
    include(pybind11)

    add_subdirectory(polyhedralGravityPython)

    target_include_directories(_core PUBLIC
            "${CMAKE_CURRENT_SOURCE_DIR}"
            )

    if (POLYHEDRAL_GRAVITY_FAST_MATH)
        target_compile_options(_core PRIVATE ${POLYHEDRAL_GRAVITY_FAST_MATH_FLAG})
    endif ()

    target_link_libraries(_core PUBLIC
            spdlog::spdlog
            tetgen_lib
            Kokkos::kokkos
            )
endif()