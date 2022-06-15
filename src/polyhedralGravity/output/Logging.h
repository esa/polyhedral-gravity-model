#pragma once

#include <memory>
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/stdout_color_sinks.h"


namespace polyhedralGravity {

    class PolyhedralGravityLogger;

    /**
     * Wrapper Class for spdlog logger
     */
    class PolyhedralGravityLogger {

    public:

        /**
         * The default logger which is used in the whole implementation for every logging.
         */
        const static PolyhedralGravityLogger DEFAULT_LOGGER;

    private:

        /**
         * The actual spdlog::logger
         */
        const std::shared_ptr<spdlog::logger> _logger;

    public:
        /**
         * Constructs a new PolyhedralGravityLogger. Further, it registers the new logger in  spdlog's registry with
         * the name POLYHEDRAL_GRAVITY_LOGGER.
         */
        PolyhedralGravityLogger()
                : _logger(spdlog::stdout_color_mt<spdlog::synchronous_factory>("POLYHEDRAL_GRAVITY_LOGGER")) {
            _logger->set_level(spdlog::level::trace);
        }

        [[nodiscard]] inline std::shared_ptr<spdlog::logger> getLogger() const {
            return _logger;
        }


    };

}