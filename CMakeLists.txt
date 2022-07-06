cmake_minimum_required(VERSION 3.16)
project(polyhedralGravity)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CMakeDependentOption)

#####################################
# PARALLELIZATION OPTIONS AND LOGGING
#####################################
# Parallelization on the HOST
set(PARALLELIZATION_HOST "CPP" CACHE STRING "Host parallelization chosen by the user
 (CPP= Serial, OMP = OpenMP, TBB = Intel Threading Building Blocks")
set_property(CACHE PARALLELIZATION_HOST PROPERTY STRINGS CPP, OMP, TBB)

# Parallelization on the DEVICE
set(PARALLELIZATION_DEVICE "CPP" CACHE STRING "Device parallelization chosen by the user
 (CPP= Serial, OMP = OpenMP, TBB = Intel Threading Building Blocks")
set_property(CACHE PARALLELIZATION_DEVICE PROPERTY STRINGS CPP, OMP, TBB)

# Enforce to use an already installed tbb library instead of compiling from source
option(USE_LOCAL_TBB "Uses the local tbb installation rather than on using the automatically fetched version from
GitHub via CMake (Default: OFF)" OFF)

# Set the Logging Level
set(LOGGING_LEVEL "2" CACHE STRING "Set the Logging level, default (INFO=2), available options:
TRACE=0, DEBUG=1, INFO=2, WARN=3, ERROR=4, CRITICAL=5, OFF=6")
set_property(CACHE LOGGING_LEVEL PROPERTY STRINGS 0, 1, 2, 3, 4, 5, 6)
add_compile_definitions(SPDLOG_ACTIVE_LEVEL=${LOGGING_LEVEL})

###################################
# What actually to build? - Options
###################################

# Build docs
option(BUILD_POLYHEDRAL_GRAVITY_DOCS "Builds the documentation (Default: OFF)" OFF)
message(STATUS "BUILD_POLYHEDRAL_GRAVITY_DOCS = ${BUILD_POLYHEDRAL_GRAVITY_DOCS}")
# Build C++ executable
option(BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE "Builds the C++ executable (Default: ON)" ON)
message(STATUS "BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE = ${BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE}")
# Build library (default ON), if the executable or tests are built this forced to ON
cmake_dependent_option(BUILD_POLYHEDRAL_GRAVITY_LIBRARY "Builds the library (Default: ON)" ON
        "NOT BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE AND NOT BUILD_POLYHEDRAL_GRAVITY_TESTS" ON)
message(STATUS "BUILD_POLYHEDRAL_GRAVITY_LIBRARY = ${BUILD_POLYHEDRAL_GRAVITY_LIBRARY}")
# Option to build the python interface
option(BUILD_POLYHEDRAL_PYTHON_INTERFACE "Set this to on if the python interface should be built (Default: ON)" ON)
message(STATUS "BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE = ${BUILD_POLYHEDRAL_PYTHON_INTERFACE}")
# Option to build tests or not
option(BUILD_POLYHEDRAL_GRAVITY_TESTS "Set to on if the tests should be built (Default: ON)" ON)
message(STATUS "BUILD_POLYHEDRAL_GRAVITY_TESTS = ${BUILD_POLYHEDRAL_GRAVITY_TESTS}")


#######################################################
# Including dependencies needed across multiple targets
#######################################################
# Appends the the module path to contain additional CMake modules for this project
# and include everything necessary
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Include dependecies
include(thrust)
include(spdlog)
include(tetgen)
include(xsimd)

###############################
# Thrust Parallelization Set-Up
###############################
# Get a version of tbb from the github repository, simplifies compilation for the user since tbb does not need to be
# preinstalled but rather gets automatically set up via CMake
# Nevertheless, there is still the option to enforce to use a local installation if one exists
if (NOT USE_LOCAL_TBB AND (${PARALLELIZATION_HOST} STREQUAL "TBB" OR ${PARALLELIZATION_DEVICE} STREQUAL "TBB"))
    include(tbb)
    thrust_set_TBB_target(tbb)
endif ()

# Thrust set-up i.e. the parallelization library, create targets according to the users specification
thrust_create_target(Thrust HOST ${PARALLELIZATION_HOST} DEVICE ${PARALLELIZATION_DEVICE})
message(STATUS "Set Parallelization: HOST ${PARALLELIZATION_HOST} DEVICE ${PARALLELIZATION_DEVICE}")

# This definitions will set the main eval(..) routine to be executed on the DEVICE if the users has given
# an option parallelize on the Device
if(NOT PARALLELIZATION_DEVICE STREQUAL "CPP")
    add_compile_definitions(DEVICE)
endif()

################################
# Building the Polyhedral Libray
################################
if (BUILD_POLYHEDRAL_GRAVITY_LIBRARY)
    # YAML is only required for the library target, but not for e.g. the Python interface
    include(yaml)

    # Including the actual Polyhedral source files
    file(GLOB_RECURSE SRC
        "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.h"
        "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.cpp")

    add_library(${PROJECT_NAME}_lib OBJECT ${SRC})

    # Adds the include Path PROJECT_SOURCE_DIR/src to the target polyhedralGravity
    target_include_directories(${PROJECT_NAME}_lib PUBLIC
            "${PROJECT_SOURCE_DIR}/src"
            )

    # Link libaries
    target_link_libraries(${PROJECT_NAME}_lib
            spdlog::spdlog
            tetgen
            yaml-cpp
            xsimd
            )

    # Adds the include Path PROJECT_SOURCE_DIR/src to the target polyhedralGravity
    target_include_directories(${PROJECT_NAME}_lib PUBLIC
            "${PROJECT_SOURCE_DIR}/src"
            )

    # Link libraries needed in all targets
    target_link_libraries(${PROJECT_NAME}_lib
            spdlog::spdlog
            yaml-cpp
            tetgen
            xsimd
            Thrust
            )

endif()

#####################################
# Building the Polyhedral Executable
#####################################
if (BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE)
    # Main file, containing the gravity model's driver
    set(mainFile ${PROJECT_SOURCE_DIR}/src/main.cpp)

    # Building the standalone Executable
    add_executable(${PROJECT_NAME} ${mainFile})

    # Link executable with library
    target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib)
endif()

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

##########################################
# Building the Polyhedral Python Interface
##########################################
if(BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE)
    message(STATUS "Building the Polyhedral Gravity Python Interface")
    add_subdirectory(${PROJECT_SOURCE_DIR}/python-binding)
endif()