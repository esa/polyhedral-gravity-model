#include "KokkosSession.h"

#include <sstream>
#include <stdexcept>
#include <vector>

#include "polyhedralGravity/output/Logging.h"

namespace polyhedralGravity::kokkos {

    namespace {

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
                    Kokkos::initialize();
                    owned = true;
                    POLYHEDRAL_GRAVITY_LOG_DEBUG("Initialized the Kokkos runtime with the execution spaces {}",
                                                 getEnabledExecutionSpaces());
                }
            }

            ~SessionGuard() {
                if (owned && Kokkos::is_initialized() && !Kokkos::is_finalized()) {
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
