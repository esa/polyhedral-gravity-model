include(FetchContent)

message(STATUS "Setting up xsimd via CMake")

#Fetches the version 8.1.0 from the official github of tbb
FetchContent_Declare(xsimd
        GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
        GIT_TAG 8.1.0
        )

FetchContent_MakeAvailable(xsimd)
