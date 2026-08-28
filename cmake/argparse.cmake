include(FetchContent)

message(STATUS "Setting up argparse")
set(ARGPARSE_VERSION 3.2)

# argparse is only used by the benchmark driver (see src/benchmark.cpp), which is why this module is
# included from within the BUILD_POLYHEDRAL_GRAVITY_BENCHMARK branch and not from the top-level file.
find_package(argparse ${ARGPARSE_VERSION} QUIET)

if (${argparse_FOUND})
    message(STATUS "Found existing argparse Library: ${argparse_DIR}")
else ()
    message(STATUS "Using argparse Library from GitHub Release ${ARGPARSE_VERSION}")

    FetchContent_Declare(argparse
            GIT_REPOSITORY https://github.com/p-ranav/argparse.git
            GIT_TAG v${ARGPARSE_VERSION}
    )
    # Disable everything we don't need
    set(ARGPARSE_BUILD_SAMPLES OFF CACHE BOOL "" FORCE)
    set(ARGPARSE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(ARGPARSE_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(argparse)

    # Disable warnings from the headers of this interface library
    get_target_property(propval argparse INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(argparse SYSTEM INTERFACE "${propval}")
endif ()
