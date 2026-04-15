# Includes the gtest.cmake
include(gtest)

# Recursively collects all the test files in the test Directory
file(GLOB_RECURSE TEST_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
        )

# Creates the Test Target
add_executable(${PROJECT_NAME}_test ${TEST_SRC})

# Links target breakupModelTest against gtest_main (entry point, so no own main method is required)
target_link_libraries(${PROJECT_NAME}_test
        GTest::gtest_main
        GTest::gtest
        GTest::gmock
        ${PROJECT_NAME}_lib
        )

# Copies the test resources into a directory next to the test executable
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/resources" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

include(GoogleTest)

# Adds the Tests to CTest by querying the test target executable
gtest_discover_tests(${PROJECT_NAME}_test)