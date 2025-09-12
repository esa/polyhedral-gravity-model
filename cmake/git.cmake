find_package(Git QUIET)

function(get_git_commit_hash OUTPUT_VAR)
    # Run a Git command to get the first 8 characters of the current commit hash
    if(NOT GIT_FOUND OR NOT GIT_EXECUTABLE)
        set(${OUTPUT_VAR} "Unknown" PARENT_SCOPE)
        return()
    endif()

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
    if(NOT GIT_FOUND OR NOT GIT_EXECUTABLE)
        set(${OUTPUT_VAR} "Unknown" PARENT_SCOPE)
        return()
    endif()

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

function(get_git_version_tag OUTPUT_VAR)
    # Runs a git command to return the latest version number (given the most similar git tag)
    if(NOT GIT_FOUND OR NOT GIT_EXECUTABLE)
        set(${OUTPUT_VAR} "Unknown" PARENT_SCOPE)
        return()
    endif()

    execute_process(
            COMMAND git describe --tags --abbrev=0
            --match v[0-9]*.[0-9]*.[0-9]*        # plain
            --match v[0-9]*.[0-9]*.[0-9]*a*      # alpha
            --match v[0-9]*.[0-9]*.[0-9]*b*      # beta
            --match v[0-9]*.[0-9]*.[0-9]*rc*     # release candidate
            --match v[0-9]*.[0-9]*.[0-9]*.post*  # post
            --match v[0-9]*.[0-9]*.[0-9]*.dev*   # dev
            OUTPUT_VARIABLE GIT_TAG
            ERROR_VARIABLE GIT_TAG_ERR
            RESULT_VARIABLE GIT_TAG_RES
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    string(LENGTH ${GIT_TAG} TAG_LENGTH)
    string(SUBSTRING ${GIT_TAG} 1 ${TAG_LENGTH} VERSION_NUMBER)
    set(${OUTPUT_VAR} "${VERSION_NUMBER}" PARENT_SCOPE)
endfunction()