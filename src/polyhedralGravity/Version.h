#pragma once

#include <string_view>

namespace polyhedralGravity {

    /**
     * The API version of the polyhedral gravity model's interface.
     * This value is utilized across C++ and Python Interface and the single value
     * to change in case of updates.
     */
    constexpr std::string_view POLYHEDRAL_GRAVITY_VERSION = "3.2.1";

    /**
     * The API parallelization backend of the polyhedral gravity model's interface.
     * This value depends on the CMake configuration.
     */
    constexpr std::string_view POLYHEDRAL_GRAVITY_PARALLELIZATION =
#ifdef POLYHEDRAL_GRAVITY_TBB
            "TBB"
#elifdef POLYHEDRAL_GRAVITY_OMP
            "OMP"
#else
            "CPP"
#endif
            ;


}// namespace polyhedralGravity
