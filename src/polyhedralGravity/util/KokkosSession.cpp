#include "polyhedralGravity/util/KokkosSession.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "polyhedralGravity/output/Logging.h"

namespace polyhedralGravity::kokkos {

    namespace {

        /**
         * Redirects everything written to @c std::cout and @c std::cerr for as long as it exists into the
         * library's log, at the DEBUG level.
         *
         * Kokkos reports its initialization diagnostics by writing to @c std::cerr directly -- most
         * prominently the OMP_PROC_BIND warning, which it emits whenever that environment variable is
         * unset. That is reasonable for an application, but this is a library: its diagnostics belong in
         * its own log and not in the output of whichever program happens to link it, which for the Python
         * interface is an interactive interpreter.
         *
         * The redirection is process-global, so it is kept to the shortest possible scope: the single
         * call to @c Kokkos::initialize() or @c Kokkos::finalize() it wraps.
         *
         * @note The obvious alternative, setting @c OMP_PROC_BIND from within the process so that the
         * warning has nothing to complain about, does not work. libgomp parses @c OMP_PROC_BIND and
         * @c OMP_PLACES in a load-time constructor, i.e. before any code of this library runs, so a
         * @c setenv() here silences the warning without binding a single thread. The variables have to be
         * set in the environment of the process before it starts to have any effect at all.
         */
        class CapturedStandardStreams {
        public:
            CapturedStandardStreams()
                : _outBuffer{std::cout.rdbuf(_capture.rdbuf())}, _errBuffer{std::cerr.rdbuf(_capture.rdbuf())} {}

            CapturedStandardStreams(const CapturedStandardStreams &) = delete;
            CapturedStandardStreams &operator=(const CapturedStandardStreams &) = delete;

            ~CapturedStandardStreams() {
                std::cout.rdbuf(_outBuffer);
                std::cerr.rdbuf(_errBuffer);
                // A diagnostic must never be the reason a process terminates, and this destructor may well
                // run while an exception from Kokkos is propagating
                try {
                    std::istringstream lines{_capture.str()};
                    std::string line{};
                    while (std::getline(lines, line)) {
                        // Kokkos pads its warnings with blank lines, which carry nothing worth a log record
                        if (line.find_first_not_of(" \t\r") != std::string::npos) {
                            POLYHEDRAL_GRAVITY_LOG_DEBUG("Kokkos: {}", line);
                        }
                    }
                } catch (...) {
                }
            }

        private:
            /** Collects both streams, so that their relative order is preserved */
            std::ostringstream _capture{};
            std::streambuf *_outBuffer;
            std::streambuf *_errBuffer;
        };

        /**
         * Owns the Kokkos runtime for the lifetime of the process.
         *
         * Kokkos may only be initialized once and it must be finalized after the last Kokkos::View has been
         * released. Tying it to the destructor of a function-local static gives exactly that ordering: static
         * destructors run at process exit, i.e. after any GravityEvaluable held by C++ code has gone out of
         * scope and after the Python interpreter has torn down its objects.
         */
        struct SessionGuard {
            /** Whether this guard, and not the embedding application, initialized Kokkos */
            bool owned{false};

            SessionGuard() {
                if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
                    {
                        const CapturedStandardStreams capture{};
                        Kokkos::initialize();
                    }// the capture forwards whatever Kokkos printed to the log when it goes out of scope
                    owned = true;
                    POLYHEDRAL_GRAVITY_LOG_DEBUG("Initialized the Kokkos runtime with the execution spaces {}",
                                                 getEnabledExecutionSpaces());
                }
            }

            ~SessionGuard() {
                if (owned && Kokkos::is_initialized() && !Kokkos::is_finalized()) {
                    const CapturedStandardStreams capture{};
                    Kokkos::finalize();
                }
            }
        };

    }// namespace

    void ensureInitialized() {
        static SessionGuard guard{};
        // Reading the guard keeps compilers from optimizing the static away
        (void) guard.owned;
    }

    std::string getEnabledExecutionSpaces() {
        std::vector<std::string> spaces{};
#ifdef KOKKOS_ENABLE_SERIAL
        spaces.emplace_back(Kokkos::Serial::name());
#endif
#ifdef KOKKOS_ENABLE_OPENMP
        spaces.emplace_back(Kokkos::OpenMP::name());
#endif
#ifdef KOKKOS_ENABLE_THREADS
        spaces.emplace_back(Kokkos::Threads::name());
#endif
        if constexpr (GPU_AVAILABLE) {
            spaces.emplace_back(DeviceSpace::name());
        }
        std::stringstream sstream{};
        for (size_t index = 0; index < spaces.size(); ++index) {
            sstream << (index == 0 ? "" : ", ") << spaces[index];
        }
        return sstream.str();
    }

    std::string getExecutionSpaceName(const ComputeBackend backend) {
        switch (backend) {
            case ComputeBackend::CPU_SERIAL:
                return SerialSpace::name();
            case ComputeBackend::CPU_PARALLEL:
                return HostParallelSpace::name();
            case ComputeBackend::GPU_PARALLEL:
                return GPU_AVAILABLE ? DeviceSpace::name() : "None";
            default:
                return "Unknown";
        }
    }

    void checkBackendAvailable(const ComputeBackend backend) {
        if (backend != ComputeBackend::GPU_PARALLEL || GPU_AVAILABLE) {
            return;
        }
        throw std::runtime_error{
                "The compute backend GPU_PARALLEL was requested, but this build of the polyhedral gravity model "
                "has no GPU backend. It was compiled with the execution spaces " + getEnabledExecutionSpaces() +
                ". Re-configure with a CUDA (NVIDIA), ROCm/HIP (AMD), or oneAPI/SYCL (Intel) toolchain available, "
                "or use the compute backend CPU_PARALLEL instead."};
    }

}// namespace polyhedralGravity::kokkos
