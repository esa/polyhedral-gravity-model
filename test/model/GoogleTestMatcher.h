#pragma once

#include "gtest/gtest.h"
#include "gmock/gmock.h"


#include "polyhedralGravity/util/UtilityFloatArithmetic.h"

MATCHER_P2(FloatContainter2D, container, ulpDiff, "Comparing 2D Containers") {
    if (container.size() != arg.size()) {
        *result_listener << "The top-level container sizes do not match. Sizes: " << container.size() << " != " << arg.size();
        return false;
    }
    for (size_t idx = 0; idx < std::size(container); ++idx) {
        if (container[idx].size() != arg[idx].size()) {
            *result_listener << "The inner container sizes at idx = " << idx << " do not match. Sizes: " << container[idx].size() << " != " << arg[idx].size();
            return false;
        }
        for (size_t idy = 0; idy < std::size(container[idx]); ++idy) {
            if (!polyhedralGravity::util::floatNear(container[idx][idy], arg[idx][idy], ulpDiff)) {
                *result_listener << "The elements at idx = " << idx << " and idy = " << idy << " do not match. Values: " << container[idx][idy] << " != " << arg[idx][idy];
                return false;
            }
        }
    }
    return true;
}