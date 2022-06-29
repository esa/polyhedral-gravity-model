include(FetchContent)

message(STATUS "Setting up xsimd via CMake")


find_package(xsimd 8.1.0 QUIET)

if (${xsimd_FOUND})

    message(STATUS "Using local xsimd installation")

else()

    #Fetches the version 8.1.0 from the official github of tbb
    FetchContent_Declare(xsimd
            GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
            GIT_TAG 8.1.0
            )

    FetchContent_MakeAvailable(xsimd)

endif()
