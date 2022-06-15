#pragma once

#include <iostream>
#include <string>
#include <exception>
#include <utility>
#include <memory>
#include "polyhedralGravity/model/GravityModelData.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/spdlog.h"

namespace polyhedralGravity {

/**
 * By using spdlog, this class provides methods to print the result the polyhedral gravity model into an CSV file.
 * @note The used Logger level is "info", so this should be activated at least if output is wished!
 */
    class CSVWriter {

        /**
        * The logger to write to the file sink
        */
        std::shared_ptr<spdlog::logger> _logger;


    public:

        /**
         * Creates a new CSVWriter.
         * Results are written to "polyhedralGravityModel.csv".
         */
        CSVWriter() : CSVWriter("polyhedralGravityModel.csv") {}

        /**
         * Creates a new CSVWriter.
         * Results are written to file with filename as name.
         * @param filename - a string
         */
        explicit CSVWriter(const std::string &filename)
                : _logger{
                spdlog::basic_logger_mt<spdlog::synchronous_factory>("CSVWriter_" + filename, filename, true)} {
            _logger->set_pattern("%v");
        }

        /**
         * De-Registers the logger of the CSVWriter, to ensure that a similar CSVWriter can be constructed once again.
         */
        ~CSVWriter() {
            spdlog::drop(_logger->name());
        }

        /**
         * Prints the result to the specified file
         * @param computationPoints - vector of computation points
         * @param gravityResults - vector of gravity results
         */
        void printResult(const std::vector<std::array<double, 3>> &computationPoints,
                         const std::vector<GravityModelResult> &gravityResults) const;

    };

}
