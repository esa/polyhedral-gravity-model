include(FetchContent)

message(STATUS "Setting up thrust")

# Set custom variables, policies, etc.
# Disable stuff not needed
set(THRUST_ENABLE_HEADER_TESTING "OFF")
set(THRUST_ENABLE_TESTING "OFF")
set(THRUST_ENABLE_EXAMPLES "OFF")

# Set standard CPP Dialect to 17 (default of thrust would be 14)
set(THRUST_CPP_DIALECT 17)

find_package(Thrust 1.16.0 QUIET)

if (${Thrust_FOUND})

    message(STATUS "Using existing thrust installation")

else()
    message(STATUS "Using thrust from git repository")
    # Fetches the version 1.16.0 of the official NVIDIA Thrust repository
    FetchContent_Declare(thrust
            GIT_REPOSITORY https://github.com/NVIDIA/thrust.git
            GIT_TAG 1.16.0
            )
    FetchContent_MakeAvailable(thrust)
endif()
