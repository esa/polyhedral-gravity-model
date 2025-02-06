#pragma once

#include <string>
#include <algorithm>

namespace polyhedralGravity::util {

    /**
     * Returns true if a string ends with the given suffix.
     * @param str the string whose suffix to check
     * @param suffix the suffix
     * @return true if the string ends with suffix, otherwise false
     */
    inline bool ends_with(const std::string& str, const std::string& suffix) {
        if (suffix.size() > str.size()) {
            return false;
        }
        return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
    }
}