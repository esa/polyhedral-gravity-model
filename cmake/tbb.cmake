include(FetchContent)

message(STATUS "Setting up tbb via CMake")

#Fetches the version v2021.5.0 from the official github of tbb
FetchContent_Declare(tbb
        GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
        GIT_TAG v2021.5.0
        )

# Disable tests & and do not treat tbb-compile errors as warnings
option(TBB_TEST "Enable testing" OFF)
option(TBB_STRICT "Treat compiler warnings as errors" OFF)

FetchContent_MakeAvailable(tbb)
