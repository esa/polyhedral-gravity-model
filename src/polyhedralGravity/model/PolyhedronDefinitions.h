#pragma once

#include <array>
#include <vector>
#include <string>
#include <tuple>
#include <ostream>

namespace polyhedralGravity {

    /**
     * Alias for an array of size 3 for x, y, z coordinates in the given floating point precision.
     * The evaluation of the gravity model is templated over this precision (see {@link ComputePrecision}),
     * whereas everything the user hands in or gets back is always in double precision.
     */
    template<typename FloatType>
    using Vector3 = std::array<FloatType, 3>;

    /**
     * Alias for an array of size 6 for xx, yy, zz, xy, xz, yz second derivatives
     * in the given floating point precision.
     */
    template<typename FloatType>
    using Vector6 = std::array<FloatType, 6>;

    /**
     * Alias for a triplet of arrays of size 3 for the segments of a triangular face
     * in the given floating point precision.
     */
    template<typename FloatType>
    using Vector3Triplet = std::array<Vector3<FloatType>, 3>;

    /**
     * Alias for an array of size 3 (double) for x, y, z coordinates.
     */
    using Array3 = Vector3<double>;

    /**
     * Alias for an array of size 3 (size_t) for the vertex indices in a triangular face.
     */
    using IndexArray3 = std::array<size_t, 3>;

    /**
     * Alias for an array of size 6 for xx, yy, zz, xy, xz, yz second derivatives.
     */
    using Array6 = Vector6<double>;

    /**
     * Alias for a triplet of arrays of size 3 for the segment of a triangular face
     */
    using Array3Triplet = Vector3Triplet<double>;

    /**
     * Contains in the order of the tuple:
     * The gravitational potential in [m^2/s^2] <--> [J/kg] at point P.
     * Related are Equation (1) and (11) of Tsoulis Paper, here referred as V.
     *
     * The first order derivatives of the gravitational potential in [m/s^2].
     * The array contains the derivatives depending on the coordinates x-y-z in this order.
     * Related are Equation (2) and (12) of Tsoulis Paper, here referred as Vx, Vy, Vz.
     *
     * The second order derivatives or also called gradiometric Tensor in [1/s^2].
     * The array contains the second order derivatives in the following order xx, yy, zz, xy, xz, yz.
     * Related are Equation (3) and (13) of Tsoulis Paper, here referred as Vxx, Vyy, Vzz, Vxy, Vxz, Vyz.
     */
    using GravityModelResult = std::tuple<double, Array3, Array6>;

    /** A polyhedron defined by a set of filenames, as available in {@link TetgenAdapter} */
    using PolyhedralFiles = std::vector<std::string>;

    /** A polyhedron defined by a set of vertices and face indices */
    using PolyhedralSource = std::tuple<std::vector<Array3>, std::vector<IndexArray3>>;


    /**
     * The orientation of the plane unit normals of the polyhedron.
     * We use this property as the precise definition of the vertices ordering depends on the
     * utilized coordinate system.
     * However, the normal alignment is independent. Tsoulis et al. equations require the
     * normals to point outwards of the polyhedron. If the opposite hold, the result is
     * negated.
     */
    enum class NormalOrientation: char {
        /** Outwards pointing plane unit normals */
        OUTWARDS,
        /** Inwards pointing plane unit normals */
        INWARDS
    };


    /**
     * Overloaded insertion operator to output the string representation of a NormalOrientation enum value.
     *
     * @param os The output stream to write the string representation to.
     * @param orientation The NormalOrientation enum value to output.
     * @return The output stream after writing the string representation.
     */
    std::ostream &operator<<(std::ostream &os, const NormalOrientation &orientation);

    /**
     * The three mode the polyhedron class takes in
     * the constructor in order to determine what initialization checks to conduct.
     * This enum is exclusively utilized in the constructor of a {@link Polyhedron} and its private method
     * {@link runIntegrityMeasures}
     */
    enum class PolyhedronIntegrity: char {
        /**
         * All activities regarding MeshChecking are disabled.
         * No runtime overhead!
         */
        DISABLE,
        /**
         * Only verification of the NormalOrientation.
         * A misalignment (e.g., specified OUTWARDS, but is not) leads to a runtime_error.
         * Runtime Cost @f$O(n^2)@f$
         */
        VERIFY,
        /**
         * Like VERIFY, but also informs the user about the option in any case on the runtime costs.
         * This is the implicit default option.
         * Runtime Cost: @f$O(n^2)@f$ and output to stdout in every case!
         */
        AUTOMATIC,
        /**
         * Verification and Automatic Healing of the NormalOrientation.
         * A misalignment does not lead to a runtime_error, but to an internal correction.
         * Runtime Cost: @f$O(n^2)@f$ and a modification of the mesh input!
         */
        HEAL,
    };

    /**
     * Represents the unit of a polyhedron's mesh.
     */
    enum class MetricUnit : char {
        /** The unit meter @f$[m]@f$ */
        METER,
        /** The unit kilometer @f$[km]@f$ */
        KILOMETER,
        /** The mesh is unitless @f$[1]@f$ */
        UNITLESS,
    };

    /**
     * Stream operator for the MetricUnit enum. Prints the enum to a human-readable string.
     * @param os The output stream to write the string representation to.
     * @param metricUnit the metric unit to print
     * @return The output stream after writing the string representation.
     */
    std::ostream &operator<<(std::ostream &os, const MetricUnit &metricUnit);

    /**
     * The compute backend on which the polyhedral gravity model is evaluated.
     * Every backend is a Kokkos execution space; which of them are compiled in is decided
     * by the CMake configuration (see {@code cmake/kokkos.cmake}).
     */
    enum class ComputeBackend : char {
        /** Evaluation on a single CPU thread, i.e. the Kokkos Serial execution space */
        CPU_SERIAL,
        /** Evaluation on all CPU threads, i.e. the Kokkos OpenMP execution space */
        CPU_PARALLEL,
        /**
         * Evaluation on the GPU using the vendor's native paradigm, i.e. the Kokkos CUDA (NVIDIA),
         * HIP (AMD), or SYCL (Intel) execution space.
         * @throws std::runtime_error if the library was not compiled with a GPU backend
         */
        GPU_PARALLEL,
    };

    /**
     * Stream operator for the ComputeBackend enum. Prints the enum to a human-readable string.
     * @param os The output stream to write the string representation to.
     * @param backend the compute backend to print
     * @return The output stream after writing the string representation.
     */
    std::ostream &operator<<(std::ostream &os, const ComputeBackend &backend);

    /**
     * The floating point precision in which the polyhedral gravity model is evaluated.
     * The polyhedron's mesh and the results are always given in double precision;
     * this enum only selects the precision of the evaluation itself.
     */
    enum class ComputePrecision : char {
        /** Single precision @f$[32\ bit]@f$ */
        FLOAT32,
        /** Double precision @f$[64\ bit]@f$ */
        FLOAT64,
    };

    /**
     * Stream operator for the ComputePrecision enum. Prints the enum to a human-readable string.
     * @param os The output stream to write the string representation to.
     * @param precision the compute precision to print
     * @return The output stream after writing the string representation.
     */
    std::ostream &operator<<(std::ostream &os, const ComputePrecision &precision);

    /**
     * Where a raw pointer handed to the library points to.
     *
     * This is what makes a {@link Polyhedron} constructible from the buffer of an array library without
     * copying it: a NumPy array lives in {@link MemoryLocation::HOST}, a PyTorch or JAX array on an
     * accelerator in {@link MemoryLocation::DEVICE}.
     */
    enum class MemoryLocation : char {
        /** The pointer addresses the CPU's memory */
        HOST,
        /**
         * The pointer addresses the compute device's memory, i.e. the GPU's.
         * @throws std::runtime_error if the library was compiled without a GPU backend
         */
        DEVICE,
    };

    /**
     * Stream operator for the MemoryLocation enum. Prints the enum to a human-readable string.
     * @param os The output stream to write the string representation to.
     * @param location the memory location to print
     * @return The output stream after writing the string representation.
     */
    std::ostream &operator<<(std::ostream &os, const MemoryLocation &location);

    /**
     * Converts a given string representation of a metric unit into the corresponding MetricUnit enum value.
     * This function maps input strings to enum values such as "METER", "KILOMETER", or "UNITLESS".
     * If the input does not match a valid metric unit representation, appropriate handling is expected.
     *
     * @param unit The input string representing the metric unit.
     * @return The corresponding MetricUnit enum value for the provided string.
     */
    MetricUnit readMetricUnit(const std::string &unit);

}