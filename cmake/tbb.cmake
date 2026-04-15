include(FetchContent)

message(STATUS "Setting up tbb")
set(TBB_VERSION 2021.12.0)

find_package(TBB QUIET HINTS /opt/homebrew/Cellar/tbb)

if(${TBB_FOUND})
    message(STATUS "Found existing TBB library: ${TBB_DIR}")
else()
    message(STATUS "Using TBB from GitHub Release ${TBB_VERSION}")

    #Fetches the version v2021.12.0 from the official github of tbb
    FetchContent_Declare(tbb
            GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
            GIT_TAG v${TBB_VERSION}
    )

    # Disable tests & and do not treat tbb-compile errors as warnings
    set(TBB_TEST OFF CACHE BOOL "" FORCE)
    set(TBB_STRICT OFF CACHE BOOl "" FORCE)

    FetchContent_MakeAvailable(tbb)
endif()
