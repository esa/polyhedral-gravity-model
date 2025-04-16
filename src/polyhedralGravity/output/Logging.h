#pragma once

#include <memory>
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/stdout_color_sinks.h"


namespace polyhedralGravity {
    class PolyhedralGravityLogger;

    /**
     * Wrapper Class for a spdlog logger.
     * It is used to configure the static DEFAULT_LOGGER of the polyhedral gravity model
     * within the constructor.
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

        /**
         * Returns the underlying pointer to spdlog's logger object.
         * @return a shared pointer to a spdlog's logger object
         */
        [[nodiscard]] inline std::shared_ptr<spdlog::logger> getLogger() const {
            return _logger;
        }


    };

#define POLYHEDRAL_GRAVITY_LOG_TRACE(msg, ...) SPDLOG_LOGGER_TRACE(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), msg, ##__VA_ARGS__)
#define POLYHEDRAL_GRAVITY_LOG_DEBUG(msg, ...) SPDLOG_LOGGER_DEBUG(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), msg, ##__VA_ARGS__)
#define POLYHEDRAL_GRAVITY_LOG_INFO(msg, ...) SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), msg, ##__VA_ARGS__)
#define POLYHEDRAL_GRAVITY_LOG_WARN(msg, ...) SPDLOG_LOGGER_WARN(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), msg, ##__VA_ARGS__)
#define POLYHEDRAL_GRAVITY_LOG_ERROR(msg, ...) SPDLOG_LOGGER_ERROR(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), msg, ##__VA_ARGS__)
#define POLYHEDRAL_GRAVITY_LOG_CRITICAL(msg, ...) SPDLOG_LOGGER_CRITICAL(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), msg, ##__VA_ARGS__)

}