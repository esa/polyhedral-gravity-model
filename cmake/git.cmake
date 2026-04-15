find_package(Git QUIET REQUIRED)

function(get_git_commit_hash OUTPUT_VAR)
    # Run a Git command to get the first 8 characters of the current commit hash
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short=8 HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_VARIABLE GIT_ERROR
            ERROR_STRIP_TRAILING_WHITESPACE
    )

    # Check if the Git command was successful
    if (NOT GIT_COMMIT_HASH OR GIT_ERROR)
        message(WARNING "Failed to retrieve Git commit hash: ${GIT_ERROR}")
        set(${OUTPUT_VAR} "UNKNOWN" PARENT_SCOPE)
    else()
        # Pass the short hash back to the calling scope
        set(${OUTPUT_VAR} "${GIT_COMMIT_HASH}" PARENT_SCOPE)
    endif()
endfunction()

function(is_git_working_tree_clean OUTPUT_VAR)
    # Run a Git command to check if the working tree is clean
    execute_process(
            COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD --
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            RESULT_VARIABLE GIT_STATUS_RESULT
            ERROR_VARIABLE GIT_ERROR
            ERROR_STRIP_TRAILING_WHITESPACE
    )

    # Check the result of the Git command
    if (NOT GIT_ERROR AND GIT_STATUS_RESULT EQUAL 0)
        # Working tree is clean
        set(${OUTPUT_VAR} TRUE PARENT_SCOPE)
    else()
        # Working tree has uncommitted changes or an error occurred
        set(${OUTPUT_VAR} FALSE PARENT_SCOPE)

        if (GIT_ERROR)
            message(WARNING "Error while checking Git working tree: ${GIT_ERROR}")
        endif()
    endif()
endfunction()