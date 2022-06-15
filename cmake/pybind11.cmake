include(FetchContent)

message(STATUS "Setting up pybind11")

#Fetches the version 2.9.2 from the official github of pybind11
FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v2.9.2
        )

FetchContent_MakeAvailable(pybind11)
