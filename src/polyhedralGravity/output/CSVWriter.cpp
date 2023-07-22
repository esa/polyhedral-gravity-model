#include "CSVWriter.h"

namespace polyhedralGravity {

    void CSVWriter::printResult(const std::vector<std::array<double, 3>> &computationPoints,
                                const std::vector<GravityModelResult> &gravityResults) const {
        _logger->info("Point P,Potential [m^2/s^2],Acceleration [m/s^2],Second Derivative Gravity Tensor [1/s^2]");
        for (size_t i = 0; i < computationPoints.size() && i < gravityResults.size(); ++i) {
            const std::array<double, 3> &computationPoint = computationPoints[i];
            const auto &[potential, acceleration, secondDerivative] = gravityResults[i];
            //Because of ADL the overload operator<< for arrays (in UtilityContainer) does not work here
            //@related https://en.cppreference.com/w/cpp/language/adl
            _logger->info("[{} {} {}],{},[{} {} {}],[{} {} {} {} {} {}]",
                          computationPoint[0], computationPoint[1], computationPoint[2], potential,
                          acceleration[0], acceleration[1], acceleration[2], secondDerivative[0], secondDerivative[1],
                          secondDerivative[2], secondDerivative[3], secondDerivative[4], secondDerivative[5]);
        }
    }

}

