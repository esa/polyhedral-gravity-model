include(FetchContent)

message(STATUS "Setting up Pybind11 Library")
set(PYBIND11_VERSION 3.0.1)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 ${PYBIND11_VERSION} QUIET)

if(${pybind11_FOUND})
    message(STATUS "Found existing Pybind11 Library: ${pybind11_DIR}")
else()
    message(STATUS "Using Pybind11 Library from GitHub Release ${PYBIND11_VERSION}")

    FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG v${PYBIND11_VERSION}
    )
    FetchContent_MakeAvailable(pybind11)
endif()