include(FetchContent)

message(STATUS "Setting up xsimd via CMake")


find_package(xsimd 9.0.1 QUIET)

if (${xsimd_FOUND})

    message(STATUS "Using existing xsimd installation")

else()

    message(STATUS "Using xsimd from git repository")

    #Fetches the version 9.0.1 from the official github of tbb
    FetchContent_Declare(xsimd
            GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
            GIT_TAG 9.0.1
            )

    FetchContent_MakeAvailable(xsimd)

endif()
