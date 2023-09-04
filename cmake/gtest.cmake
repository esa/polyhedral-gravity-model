include(FetchContent)

message(STATUS "Setting up gtest")

#Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/testing/googletest.html
#Fetches the version 1.13.0 from the official github for googletest
FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.13.0
        )

FetchContent_MakeAvailable(googletest)

# Disable warnings from the library target
target_compile_options(gtest_main PRIVATE -w)
# Disable warnings from included headers
get_target_property(propval gtest_main INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(gtest_main SYSTEM PUBLIC "${propval}")