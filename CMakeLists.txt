cmake_minimum_required(VERSION 3.16)
project(polyhedralGravity)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Appends the the module path to contain additional CMake modules for this project
# and include everything necessary
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CMakeDependentOption)

#####################################
# PARALLELIZATION OPTIONS AND LOGGING
#####################################
# Parallelization on the HOST
set(POLYHEDRAL_GRAVITY_PARALLELIZATION "CPP" CACHE STRING "Host parallelization chosen by the user
 (CPP= Serial, OMP = OpenMP, TBB = Intel Threading Building Blocks")
set_property(CACHE POLYHEDRAL_GRAVITY_PARALLELIZATION PROPERTY STRINGS CPP, OMP, TBB)

# Set the Logging Level
set(POLYHEDRAL_GRAVITY_LOGGING_LEVEL_LIST "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "CRITICAL" "OFF")
set(POLYHEDRAL_GRAVITY_LOGGING_LEVEL "INFO" CACHE STRING "Set the Logging level, default (INFO), available options: TRACE, DEBUG, INFO, WARN, ERROR, CRITICAL, OFF")
set_property(CACHE POLYHEDRAL_GRAVITY_LOGGING_LEVEL PROPERTY STRINGS ${POLYHEDRAL_GRAVITY_LOGGING_LEVEL_LIST})
# Convert the logging level string to its corresponding number
list(FIND POLYHEDRAL_GRAVITY_LOGGING_LEVEL_LIST ${POLYHEDRAL_GRAVITY_LOGGING_LEVEL} LOGGING_LEVEL_INDEX)
if (${LOGGING_LEVEL_INDEX} EQUAL -1)
    message(FATAL_ERROR "Invalid logging level: ${POLYHEDRAL_GRAVITY_LOGGING_LEVEL}")
endif ()
add_compile_definitions(SPDLOG_ACTIVE_LEVEL=${LOGGING_LEVEL_INDEX})

#########################################################
# What actually to build? - Options, Versions and Output
#########################################################
# Build docs
option(BUILD_POLYHEDRAL_GRAVITY_DOCS "Builds the documentation (Default: OFF)" OFF)
# Build C++ executable
option(BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE "Builds the C++ executable (Default: ON)" ON)
# Build library (default ON), if the executable or tests are built this forced to ON
cmake_dependent_option(BUILD_POLYHEDRAL_GRAVITY_LIBRARY "Builds the library (Default: ON)" ON
        "NOT BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE AND NOT BUILD_POLYHEDRAL_GRAVITY_TESTS" ON)
# Option to build the python interface
option(BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE "Set this to on if the python interface should be built (Default: ON)" ON)
# Option to build tests or not
option(BUILD_POLYHEDRAL_GRAVITY_TESTS "Set to on if the tests should be built (Default: ON)" ON)

IF(_LIBCPP_DISABLE_AVAILABILITY)
    message(STATUS "Disabling availability macros for libc++")
    add_definitions(-D_LIBCPP_DISABLE_AVAILABILITY)
endif ()

# Resolves missing fmt symbols when working with spdlog (bundled via brew/ conda on Arm architecture)
# Refer to https://github.com/gabime/spdlog/issues/660
add_compile_definitions(FMT_HEADER_ONLY)

include(git)
include(version.cmake)

message(STATUS "#################################################################")
message(STATUS "Polyhedral Gravity Version          ${POLYHEDRAL_GRAVITY_VERSION}")
message(STATUS "Polyhedral Gravity Commit Hash      ${POLYHEDRAL_GRAVITY_COMMIT_HASH}")
message(STATUS "Polyhedral Parallelization Backend  ${POLYHEDRAL_GRAVITY_PARALLELIZATION}")
message(STATUS "Polyhedral Gravity Logging Level    ${POLYHEDRAL_GRAVITY_LOGGING_LEVEL}")
message(STATUS "#################################################################")
message(STATUS "Polyhedral Gravity Documentation    ${BUILD_POLYHEDRAL_GRAVITY_DOCS}")
message(STATUS "Polyhedral Gravity Library          ${BUILD_POLYHEDRAL_GRAVITY_LIBRARY}")
message(STATUS "Polyhedral Gravity C++ Executable   ${BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE}")
message(STATUS "Polyhedral Gravity Python Interface ${BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE}")
message(STATUS "Polyhedral Gravity Tests            ${BUILD_POLYHEDRAL_GRAVITY_TESTS}")
message(STATUS "#################################################################")
#######################################################
# Including dependencies needed across multiple targets
#######################################################

include(thrust)
include(spdlog)
include(tetgen)
include(xsimd)
include(clang_format)

###############################
# Thrust Parallelization Set-Up
###############################
# Get a version of tbb from the github repository, simplifies compilation for the user since tbb does not need to be
# preinstalled but rather gets automatically set up via CMake
# Nevertheless, there is still the option to enforce to use a local installation if one exists
if (${POLYHEDRAL_GRAVITY_PARALLELIZATION} STREQUAL "TBB")
    include(tbb)
    thrust_set_TBB_target(TBB::tbb)
    add_compile_definitions(POLYHEDRAL_GRAVITY_TBB)
elseif (${POLYHEDRAL_GRAVITY_PARALLELIZATION} STREQUAL "OMP")
    add_compile_definitions(POLYHEDRAL_GRAVITY_OMP)
endif ()

# Thrust set-up i.e. the parallelization library, create targets according to the users specification
thrust_create_target(Thrust HOST CPP DEVICE ${POLYHEDRAL_GRAVITY_PARALLELIZATION})
message(STATUS "Set Parallelization: HOST CPP DEVICE ${POLYHEDRAL_GRAVITY_PARALLELIZATION}")

############################################################
# Polyhedral Gravity Library & Executable & Python Interface
############################################################
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

##############################
# Building the Polyhedral Docs
##############################
if (BUILD_POLYHEDRAL_GRAVITY_DOCS)
    add_subdirectory(${PROJECT_SOURCE_DIR}/docs)
endif ()

###############################
# Building the Polyhedral Tests
###############################
if (BUILD_POLYHEDRAL_GRAVITY_TESTS)
    message(STATUS "Building the Polyhedral Gravity Tests")
    # Enables CTest, must be in the top-level CMakeList.txt, otherwise it won't work
    enable_testing()

    # Subdirectory where the tests are located
    add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif ()
