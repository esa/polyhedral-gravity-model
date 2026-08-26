# opencl_to_header.cmake - Script executed at build time by compile_opencl()
# Embeds an OpenCL kernel source file (*.cl) into a C++ header as a raw string literal.
# Do not call this script by itself.
file(READ "${INPUT_FILE}" FILE_CONTENTS)

# The C++ standard limits a raw string delimiter to 16 characters
set(RAW_DELIM "clKernel")

set(HEADER_TEXT
        "#pragma once

// This file is generated at build time from ${INPUT_FILE}. Do not edit it manually.

namespace ${NAMESPACE} {

    inline constexpr const char ${VAR_NAME}[] = R\"${RAW_DELIM}(
${FILE_CONTENTS}
)${RAW_DELIM}\";

}// namespace ${NAMESPACE}
")

file(WRITE "${OUTPUT_FILE}" "${HEADER_TEXT}")
