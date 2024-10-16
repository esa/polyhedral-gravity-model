include(FetchContent)

message(STATUS "Setting up thrust")
set(THRUST_VERSION 1.16.0)

# Set custom variables, policies, etc.
# Disable stuff not needed
set(THRUST_ENABLE_HEADER_TESTING "OFF")
set(THRUST_ENABLE_TESTING "OFF")
set(THRUST_ENABLE_EXAMPLES "OFF")
# Set standard CPP Dialect to 17 (default of thrust would be 14)
set(THRUST_CPP_DIALECT 17)

find_package(Thrust ${THRUST_VERSION} QUIET)


if (${Thrust_FOUND})
    message(STATUS "Found existing thrust installation: ${Thrust_DIR}")
else()
    message(STATUS "Using thrust from GitHub Release ${THRUST_VERSION}")
    FetchContent_Declare(thrust
            GIT_REPOSITORY https://github.com/NVIDIA/thrust.git
            GIT_TAG ${THRUST_VERSION}
            )
    FetchContent_MakeAvailable(thrust)
endif()
