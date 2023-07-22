include(FetchContent)

message(STATUS "Setting up pybind11")

find_package(pybind11 2.10.4 QUIET)

if (${pybind11_FOUND})

    message(STATUS "Using local pybind11 installation")

else()

    message(STATUS "Using pybind11 from git repository")

    #Fetches the version 2.10.4 from the official github of pybind11
    FetchContent_Declare(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG v2.10.4
            )

    FetchContent_MakeAvailable(pybind11)

endif()
