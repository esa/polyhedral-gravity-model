#pragma once

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "polyhedralGravity/util/UtilityFloatArithmetic.h"

/*
 * This file provides several matchers to be used with ASSERT_TAHT(..) statements for multi-dimensional containers
 * with flotaing points as content.
 * GoogleTest provides matchers like DoubleNear(..), ContainerEq(..) and Pointwise(..). These work great when working
 * with 1D containers. However, there is a lack of these matchers for multi-dimensional containers.
 * This files provides Matcher for comparing nested containers with floating points where the comparison is
 * conducted with an epsilon to counter the effect of potential rounding erros due to floating point arithmetic.
 */


MATCHER_P(FloatContainter1D, container, "Comparing 1D Containers") {
    if (container.size() != arg.size()) {
        *result_listener << "The container sizes do not match. Sizes: " << container.size() << " != " << arg.size();
        return false;
    }
    for (size_t idx = 0; idx < std::size(container); ++idx) {
        if (!polyhedralGravity::util::almostEqualRelative(container[idx], arg[idx])) {
            *result_listener
                    << "The elements at idx = " << idx << " do not match. Values: "
                    << container[idx] << " != " << arg[idx];
            return false;
        }
    }
    return true;
}

// Using the FloatContainter1D would be nice, but this leads to issues with template instatantation
// Hence, this is the easy way and a double-for-loop (but at least better fitting messages)
MATCHER_P(FloatContainter2D, container, "Comparing 2D Containers") {
    if (container.size() != arg.size()) {
        *result_listener << "The container sizes do not match. Sizes: " << container.size() << " != " << arg.size();
        return false;
    }
    for (size_t idx = 0; idx < std::size(container); ++idx) {
        if (container[idx].size() != arg[idx].size()) {
            *result_listener
                    << "The container sizes at idx = " << idx << " do not match. Sizes: "
                    << container[idx].size() << " != " << arg[idx].size();
            return false;
        }
        for (size_t idy = 0; idy < std::size(container[idx]); ++idy) {
            if (!polyhedralGravity::util::almostEqualRelative(container[idx][idy], arg[idx][idy])) {
                *result_listener
                        << "The elements at idx = " << idx << " and idy = " << idy << " do not match. Values: "
                        << container[idx][idy] << " != " << arg[idx][idy];
                return false;
            }
        }
    }
    return true;
}

// Using the FloatContainter2D would be nice, but this leads to issues with template instatantation
// Hence, this is the easy way and a triple-for-loop (but at least better fitting messages)
MATCHER_P(FloatContainter3D, container, "Comparing 3D Containers") {
    if (container.size() != arg.size()) {
        *result_listener << "The container sizes do not match. Sizes: " << container.size() << " != " << arg.size();
        return false;
    }
    for (size_t idx = 0; idx < std::size(container); ++idx) {
        if (container[idx].size() != arg[idx].size()) {
            *result_listener
                    << "The container sizes at idx = " << idx << " do not match. Sizes: "
                    << container[idx].size() << " != " << arg[idx].size();
            return false;
        }
        for (size_t idy = 0; idy < std::size(container[idx]); ++idy) {
            if (container[idx][idy].size() != arg[idx][idy].size()) {
                *result_listener
                        << "The container sizes at idx = " << idx << " do not match. Sizes: "
                        << container[idx].size() << " != " << arg[idx].size();
                return false;
            }
            for (size_t idz = 0; idz < std::size(container[idx][idy]); ++idz) {
                if (!polyhedralGravity::util::almostEqualRelative(container[idx][idy][idz], arg[idx][idy][idz])) {
                    *result_listener
                            << "The elements at idx = " << idx << " and idy = " << idy << " and idz = " << idz << " do not match. Values: "
                            << container[idx][idy][idz] << " != " << arg[idx][idy][idz];
                    return false;
                }
            }
        }
    }
    return true;
}
