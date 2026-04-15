include(FetchContent)

message(STATUS "Setting up xsimd Library")
set(XSIMD_VERSION 11.1.0)

find_package(xsimd ${XSIMD_VERSION} QUIET)

if (${xsimd_FOUND})
    message(STATUS "Found existing xsimd Library: ${xsimd_DIR}")
else()
    message(STATUS "Using xsimd Library from GitHub Release ${XSIMD_VERSION}")
    FetchContent_Declare(xsimd
            GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
            GIT_TAG ${XSIMD_VERSION}
            )

    FetchContent_MakeAvailable(xsimd)
endif()
