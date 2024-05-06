#include "UtilityFloatArithmetic.h"


namespace polyhedralGravity::util {

    bool floatNear(double lhs, double rhs, int ulpDistance) {
        if (lhs == rhs) {
            return true;
        }
        if (lhs < 0.0 && rhs > 0.0 || lhs > 0.0 && rhs < 0.0) {
            return false;
        }
        return reinterpret_cast<std::int64_t&>(rhs) - reinterpret_cast<std::int64_t&>(lhs) <= ulpDistance;
    }

    bool floatNear(float lhs, float rhs, int ulpDistance) {
        if (lhs == rhs) {
            return true;
        }
        if (lhs < 0.0 && rhs > 0.0 || lhs > 0.0 && rhs < 0.0) {
            return false;
        }
        return reinterpret_cast<std::int32_t&>(rhs) - reinterpret_cast<std::int32_t&>(lhs) <= ulpDistance;
    }


}