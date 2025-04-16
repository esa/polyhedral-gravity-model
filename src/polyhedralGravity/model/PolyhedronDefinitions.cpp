#include "PolyhedronDefinitions.h"

namespace polyhedralGravity {

    std::ostream &operator<<(std::ostream &os, const NormalOrientation &orientation) {
        switch (orientation) {
            case NormalOrientation::OUTWARDS:
                os << "OUTWARDS";
            break;
            case NormalOrientation::INWARDS:
                os << "INWARDS";
            break;
            default:
                os << "Unknown";
            break;
        }
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const MetricUnit &metricUnit) {
        switch (metricUnit) {
            case MetricUnit::METER:
                os << "m";
            break;
            case MetricUnit::KILOMETER:
                os << "km";
            break;
            default:
                os << "unitless";
            break;
        }
        return os;
    }

    MetricUnit readMetricUnit(const std::string &unit) {
        if (unit == "m") {
            return MetricUnit::METER;
        } else if (unit == "km") {
            return MetricUnit::KILOMETER;
        } else if (unit == "unitless") {
            return MetricUnit::UNITLESS;
        } else {
            throw std::runtime_error{"The unit of the mesh is not supported! Must be either 'm', 'km' or 'unitless'"};
        }
    }
}