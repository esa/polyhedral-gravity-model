#include "GravityEvaluable.h"

#include <sstream>

#include "polyhedralGravity/kokkos/KokkosEvaluation.h"

namespace polyhedralGravity {

    GravityEvaluable::GravityEvaluable(const Polyhedron &polyhedron, const ComputePrecision precision)
        : _polyhedron{polyhedron},
          _precision{precision},
          _evaluation{kokkos::createKokkosEvaluation(polyhedron, precision)} {
    }

    GravityEvaluable::GravityEvaluable(const Polyhedron &polyhedron,
                                       const std::vector<Array3Triplet> &segmentVectors,
                                       const std::vector<Array3> &planeUnitNormals,
                                       const std::vector<Array3Triplet> &segmentUnitNormals,
                                       const ComputePrecision precision)
        : _polyhedron{polyhedron},
          _precision{precision},
          _evaluation{kokkos::createKokkosEvaluation(polyhedron, precision, segmentVectors, planeUnitNormals,
                                                     segmentUnitNormals)} {
    }

    std::variant<GravityModelResult, std::vector<GravityModelResult>>
    GravityEvaluable::operator()(const std::variant<Array3, std::vector<Array3>> &computationPoints,
                                 const ComputeBackend backend) const {
        if (std::holds_alternative<Array3>(computationPoints)) {
            return _evaluation->evaluate(std::get<Array3>(computationPoints), backend);
        }
        return _evaluation->evaluate(std::get<std::vector<Array3>>(computationPoints), backend);
    }

    std::string GravityEvaluable::toString() const {
        std::stringstream sstream;
        const auto [unitPotential, unitAcceleration, unitGradiometricTensor] = getOutputMetricUnit();
        sstream << "<polyhedral_gravity.GravityEvaluable, polyhedron = " << _polyhedron.toString()
                << ", output_units = " << unitPotential << ", " << unitAcceleration << ", " << unitGradiometricTensor
                << ", precision = " << _precision << ">";
        return sstream.str();
    }

    std::array<std::string, 3> GravityEvaluable::getOutputMetricUnit() const {
        const auto metric = _polyhedron.getMeshUnit();
        if (metric != MetricUnit::UNITLESS) {
            const std::string metricString = _polyhedron.getMeshUnitAsString();
            return {metricString + "^2/s^2", metricString + "/s^2", "1/s^2"};
        } else {
            return {"1/s^2", "1/s^2", "1/s^2"};
        }
    }

    std::tuple<Polyhedron, std::vector<Array3Triplet>, std::vector<Array3>, std::vector<Array3Triplet>>
    GravityEvaluable::getState() const {
        const auto [segmentVectors, planeUnitNormals, segmentUnitNormals] = _evaluation->getCaches();
        return std::make_tuple(_polyhedron, segmentVectors, planeUnitNormals, segmentUnitNormals);
    }

    ComputePrecision GravityEvaluable::getComputePrecision() const {
        return _precision;
    }

}// namespace polyhedralGravity
