#pragma once

#include <array>
#include <vector>
#include <string>
#include <tuple>
#include <ostream>

namespace polyhedralGravity {

    /**
     * Alias for an array of size 3 (double) for x, y, z coordinates.
     */
    using Array3 = std::array<double, 3>;

    /**
     * Alias for an array of size 3 (size_t) for the vertex indices in a triangular face.
     */
    using IndexArray3 = std::array<size_t, 3>;

    /**
     * Alias for an array of size 6 for xx, yy, zz, xy, xz, yz second derivatives.
     */
    using Array6 = std::array<double, 6>;

    /**
     * Alias for a triplet of arrays of size 3 for the segment of a triangular face
     */
    using Array3Triplet = std::array<Array3, 3>;

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
     * Converts a given string representation of a metric unit into the corresponding MetricUnit enum value.
     * This function maps input strings to enum values such as "METER", "KILOMETER", or "UNITLESS".
     * If the input does not match a valid metric unit representation, appropriate handling is expected.
     *
     * @param unit The input string representing the metric unit.
     * @return The corresponding MetricUnit enum value for the provided string.
     */
    MetricUnit readMetricUnit(const std::string &unit);

    /**
     * The compute backend used to evaluate the polyhedral gravity model.
     * The chosen backend only influences how the result is computed, not the result itself
     * (up to the floating point precision selected via {@link ComputePrecision}).
     */
    enum class ComputeBackend : char {
        /**
         * Evaluation on the host using the parallelization technology the library was compiled with
         * (see POLYHEDRAL_GRAVITY_PARALLELIZATION).
         */
        CPU,
        /**
         * Evaluation on an OpenCL device, i.e. typically a GPU.
         * Requires the library to be compiled with POLYHEDRAL_GRAVITY_ENABLE_OPENCL.
         * Falls back to {@link ComputeBackend::CPU} if no OpenCL device supporting the requested
         * {@link ComputePrecision} is available.
         */
        OPENCL,
    };

    /**
     * Stream operator for the ComputeBackend enum. Prints the enum to a human-readable string.
     * @param os The output stream to write the string representation to.
     * @param backend the compute backend to print
     * @return The output stream after writing the string representation.
     */
    std::ostream &operator<<(std::ostream &os, const ComputeBackend &backend);

    /**
     * The floating point precision in which a {@link ComputeBackend} evaluates the gravity model.
     *
     * This exclusively concerns the device-side computation. The host-side interface is always
     * double precision, and the {@link ComputeBackend::CPU} backend always computes in double
     * precision, so this setting is ignored for it.
     */
    enum class ComputePrecision : char {
        /**
         * Single precision (32 bit) evaluation.
         * Considerably faster on most GPUs and the only option on devices without @c cl_khr_fp64,
         * such as every Apple Silicon GPU. Note that Tsoulis' algorithm evaluates differences of
         * transcendental terms which are prone to cancellation, so single precision noticeably
         * reduces the accuracy of the result.
         */
        FLOAT32,
        /**
         * Double precision (64 bit) evaluation, matching the accuracy of the CPU backend.
         * Requires the OpenCL device to support the @c cl_khr_fp64 extension.
         */
        FLOAT64,
    };

    /**
     * Stream operator for the ComputePrecision enum. Prints the enum to a human-readable string.
     * @param os The output stream to write the string representation to.
     * @param precision the compute precision to print
     * @return The output stream after writing the string representation.
     */
    std::ostream &operator<<(std::ostream &os, const ComputePrecision &precision);

}