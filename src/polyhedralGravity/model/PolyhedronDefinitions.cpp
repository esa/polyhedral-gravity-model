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

    std::ostream &operator<<(std::ostream &os, const ComputeBackend &backend) {
        switch (backend) {
            case ComputeBackend::CPU:
                os << "CPU";
            break;
            case ComputeBackend::OPENCL:
                os << "OpenCL";
            break;
            default:
                os << "Unknown";
            break;
        }
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const ComputePrecision &precision) {
        switch (precision) {
            case ComputePrecision::FLOAT32:
                os << "float32";
            break;
            case ComputePrecision::FLOAT64:
                os << "float64";
            break;
            default:
                os << "Unknown";
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