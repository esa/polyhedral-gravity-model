#pragma once

#include <string>
#include <algorithm>

namespace polyhedralGravity::util {

    /**
     * Returns true if a string ends with the given suffix or any of the given suffices.
     * @param str the string whose suffix to check
     * @param suffix the suffix
     * @param suffixes more suffixes arguments.
     * @return true if the string ends with suffix, otherwise false
     */
    template <typename... Args>
    inline bool ends_with(const std::string& str, const std::string& suffix, const Args&... suffixes) {
        if (suffix.size() <= str.size() && std::equal(suffix.rbegin(), suffix.rend(), str.rbegin())) {
            return true;
        }
        if constexpr (sizeof...(suffixes) > 0) {
            return ends_with(str, suffixes...);
        }
        return false;
    }
}