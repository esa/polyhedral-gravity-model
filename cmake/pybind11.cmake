include(FetchContent)

message(STATUS "Setting up pybind11")

find_package(pybind11 2.12.0 QUIET)

if (${pybind11_FOUND})

    message(STATUS "Using local pybind11 installation")

else()

    message(STATUS "Using pybind11 from git repository")

    FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG v2.12.0
    )

    FetchContent_MakeAvailable(pybind11)

endif()
