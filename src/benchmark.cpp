/**
 * A driver which evaluates one polyhedron at N random computation points, for a configurable set of
 * backends and precisions.
 *
 * It exists so that the multi point kernel can be timed and profiled (nsys, ncu) without the Python
 * interpreter in between. Build it with -DBUILD_POLYHEDRAL_GRAVITY_BENCHMARK=ON, then call it as
 *
 *     ./polyhedralGravity_benchmark [options] <mesh file>...
 *
 * The mesh is given in any format the library reads, i.e. a .node/.face pair or a single .obj/.tab/
 * .off/.ply/.stl/.mesh file. Every remaining knob has a default, and the point sample is drawn from a
 * fixed seed, so that two runs of the same command line measure the very same work.
 *
 * A profiler wants "-n 1 -r 1 -b gpu -p f32 --skip-accuracy-check", i.e. exactly one launch of the
 * kernel under investigation and nothing else.
 *
 * All output goes through the library's logger on INFO, so it is visible with the default build and is
 * silenced by configuring with -DPOLYHEDRAL_GRAVITY_LOGGING_LEVEL=WARN or higher.
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <argparse/argparse.hpp>

#include "polyhedralGravity/model/GravityEvaluable.h"
#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/PolyhedronDefinitions.h"
#include "polyhedralGravity/output/Logging.h"
#include "polyhedralGravity/util/KokkosSession.h"

using namespace polyhedralGravity;

namespace {

    /**
     * The command line spelling of every enum the driver exposes. These maps are the single source of
     * both the help text listing the accepted values and the validation of what was actually given.
     */
    const std::map<std::string, ComputeBackend> BACKENDS{
            {"serial", ComputeBackend::CPU_SERIAL},
            {"cpu", ComputeBackend::CPU_PARALLEL},
            {"gpu", ComputeBackend::GPU_PARALLEL},
    };

    const std::map<std::string, ComputePrecision> PRECISIONS{
            {"f32", ComputePrecision::FLOAT32},
            {"f64", ComputePrecision::FLOAT64},
    };

    const std::map<std::string, NormalOrientation> ORIENTATIONS{
            {"outwards", NormalOrientation::OUTWARDS},
            {"inwards", NormalOrientation::INWARDS},
    };

    const std::map<std::string, PolyhedronIntegrity> INTEGRITIES{
            {"disable", PolyhedronIntegrity::DISABLE},
            {"verify", PolyhedronIntegrity::VERIFY},
            {"automatic", PolyhedronIntegrity::AUTOMATIC},
            {"heal", PolyhedronIntegrity::HEAL},
    };

    const std::map<std::string, MetricUnit> METRIC_UNITS{
            {"meter", MetricUnit::METER},
            {"kilometer", MetricUnit::KILOMETER},
            {"unitless", MetricUnit::UNITLESS},
    };

    /** "a, b, c", i.e. the accepted spellings of one of the maps above for a help text. */
    template<typename Enum>
    std::string joinChoices(const std::map<std::string, Enum> &mapping) {
        std::string joined{};
        for (const auto &[name, value]: mapping) {
            joined.append(joined.empty() ? "" : ", ").append(name);
        }
        return joined;
    }

    /**
     * The enum one of the maps above associates with the given spelling.
     *
     * The validation is done here rather than through argparse's choices(), since those are checked
     * against the string representation of the default value, which an option taking several values
     * (--backend, --precision) does not have.
     *
     * @throws std::invalid_argument naming the accepted spellings if the given one is not among them
     */
    template<typename Enum>
    Enum lookup(const std::map<std::string, Enum> &mapping, const std::string &name, const std::string &option) {
        const auto entry = mapping.find(name);
        if (entry == mapping.end()) {
            throw std::invalid_argument{option + ": '" + name + "' is not one of " + joinChoices(mapping)};
        }
        return entry->second;
    }

    /** The axis aligned bounding box of the mesh, as the two corners (minimum, maximum). */
    std::pair<Array3, Array3> boundingBox(const std::vector<Array3> &vertices) {
        Array3 minimum{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::max()};
        Array3 maximum{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                       std::numeric_limits<double>::lowest()};
        for (const Array3 &vertex: vertices) {
            for (size_t axis = 0; axis < 3; ++axis) {
                minimum[axis] = std::min(minimum[axis], vertex[axis]);
                maximum[axis] = std::max(maximum[axis], vertex[axis]);
            }
        }
        return {minimum, maximum};
    }

    /**
     * The computation points, drawn uniformly from a box and reproducible for a fixed seed.
     *
     * Every point is drawn from the same generator in the same order, so that the sample only depends
     * on the seed and on the box, but not on which configurations are being measured.
     */
    std::vector<Array3> samplePoints(const size_t count, const Array3 &minimum, const Array3 &maximum,
                                     const unsigned int seed) {
        std::mt19937 generator{seed};
        std::array<std::uniform_real_distribution<double>, 3> distributions{
                std::uniform_real_distribution<double>{minimum[0], maximum[0]},
                std::uniform_real_distribution<double>{minimum[1], maximum[1]},
                std::uniform_real_distribution<double>{minimum[2], maximum[2]},
        };
        std::vector<Array3> points(count);
        for (Array3 &point: points) {
            point = {distributions[0](generator), distributions[1](generator), distributions[2](generator)};
        }
        return points;
    }

    /** The worst relative error of a result against a reference, for the potential and the acceleration. */
    std::pair<double, double> worstRelativeError(const std::vector<GravityModelResult> &actual,
                                                 const std::vector<GravityModelResult> &expected) {
        double worstPotential = 0.0;
        double worstAcceleration = 0.0;
        for (size_t index = 0; index < expected.size(); ++index) {
            const auto relative = [](const double got, const double want) {
                return want == 0.0 ? 0.0 : std::abs(got - want) / std::abs(want);
            };
            worstPotential =
                    std::max(worstPotential, relative(std::get<0>(actual[index]), std::get<0>(expected[index])));
            // Relative to the magnitude of the vector, not component wise: single components pass through
            // zero, where any component wise relative error is meaningless
            double errorNorm = 0.0;
            double referenceNorm = 0.0;
            for (size_t component = 0; component < 3; ++component) {
                const double difference =
                        std::get<1>(actual[index])[component] - std::get<1>(expected[index])[component];
                errorNorm += difference * difference;
                referenceNorm += std::get<1>(expected[index])[component] * std::get<1>(expected[index])[component];
            }
            if (referenceNorm != 0.0) {
                worstAcceleration = std::max(worstAcceleration, std::sqrt(errorNorm / referenceNorm));
            }
        }
        return {worstPotential, worstAcceleration};
    }

    /** One row of the result table, i.e. what was measured for one backend/ precision combination. */
    struct Measurement {
        std::string backend;
        std::string precision;
        double best;
        double mean;
        double worstPotential;
        double worstAcceleration;
        double checksum;
    };

}// namespace

int main(const int argc, char *argv[]) {
    argparse::ArgumentParser program{"polyhedralGravity_benchmark"};
    program.add_description("Times the polyhedral gravity model for one mesh at N random computation points.");

    program.add_argument("mesh")
            .help("the mesh file(s) describing the polyhedron: a .node and a .face file (in that order), "
                  "or a single .obj/.tab/.off/.ply/.stl/.mesh file")
            .nargs(argparse::nargs_pattern::at_least_one);

    program.add_argument("-n", "--points")
            .help("the number of computation points to evaluate")
            .default_value(std::size_t{1000})
            .scan<'u', std::size_t>();

    program.add_argument("-r", "--repetitions")
            .help("how often each configuration is timed; the best and the mean run are reported")
            .default_value(std::size_t{10})
            .scan<'u', std::size_t>();

    program.add_argument("-b", "--backend")
            .help("the compute backend(s) to time, one measurement per backend (" + joinChoices(BACKENDS) + ")")
            .default_value(std::vector<std::string>{"gpu"})
            .nargs(argparse::nargs_pattern::at_least_one);

    program.add_argument("-p", "--precision")
            .help("the floating point precision(s) to time, one measurement per precision (" +
                  joinChoices(PRECISIONS) + ")")
            .default_value(std::vector<std::string>{"f32"})
            .nargs(argparse::nargs_pattern::at_least_one);

    program.add_argument("-s", "--seed")
            .help("the seed of the point sample, so that a run is reproducible")
            .default_value(42u)
            .scan<'u', unsigned int>();

    program.add_argument("-d", "--density")
            .help("the density of the polyhedron; it must match the unit of the mesh")
            .default_value(1.0)
            .scan<'g', double>();

    program.add_argument("--orientation")
            .help("the direction the plane unit normals point in (" + joinChoices(ORIENTATIONS) + ")")
            .default_value(std::string{"outwards"});

    program.add_argument("--integrity")
            .help("whether the mesh is checked/ healed before it is evaluated (" + joinChoices(INTEGRITIES) + ")")
            .default_value(std::string{"disable"});

    program.add_argument("--metric-unit")
            .help("the unit of the mesh's coordinates (" + joinChoices(METRIC_UNITS) + ")")
            .default_value(std::string{"meter"});

    program.add_argument("--sample-box")
            .help("the box the computation points are drawn from, as MIN MAX applied to every axis "
                  "(default: the mesh's bounding box, scaled by --sample-scale)")
            .nargs(2)
            .scan<'g', double>();

    program.add_argument("--sample-scale")
            .help("the factor the mesh's bounding box is scaled by to obtain the sample box; a factor "
                  "above one puts points outside of the body as well. Ignored with --sample-box")
            .default_value(2.0)
            .scan<'g', double>();

    program.add_argument("--skip-accuracy-check")
            .help("skip the reference evaluation each configuration is compared against; the reference "
                  "costs one additional CPU_PARALLEL/FLOAT64 evaluation per run")
            .flag();

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception &exception) {
        POLYHEDRAL_GRAVITY_LOG_ERROR("{}", exception.what());
        POLYHEDRAL_GRAVITY_LOG_ERROR("{} (see --help)", program.usage());
        return 1;
    }

    const auto meshFiles = program.get<std::vector<std::string>>("mesh");
    const auto pointCount = program.get<std::size_t>("--points");
    const auto repetitions = program.get<std::size_t>("--repetitions");
    const auto backendNames = program.get<std::vector<std::string>>("--backend");
    const auto precisionNames = program.get<std::vector<std::string>>("--precision");
    const auto seed = program.get<unsigned int>("--seed");
    const auto density = program.get<double>("--density");
    const auto skipAccuracyCheck = program.get<bool>("--skip-accuracy-check");

    // Every spelling of an enum is resolved before any work is done, so that a typo fails immediately
    // and not only once the mesh has been read
    NormalOrientation orientation{};
    PolyhedronIntegrity integrity{};
    MetricUnit metricUnit{};
    std::vector<ComputeBackend> backends{};
    std::vector<ComputePrecision> precisions{};
    try {
        orientation = lookup(ORIENTATIONS, program.get<std::string>("--orientation"), "--orientation");
        integrity = lookup(INTEGRITIES, program.get<std::string>("--integrity"), "--integrity");
        metricUnit = lookup(METRIC_UNITS, program.get<std::string>("--metric-unit"), "--metric-unit");
        for (const std::string &name: backendNames) {
            backends.push_back(lookup(BACKENDS, name, "--backend"));
        }
        for (const std::string &name: precisionNames) {
            precisions.push_back(lookup(PRECISIONS, name, "--precision"));
        }
    } catch (const std::exception &exception) {
        POLYHEDRAL_GRAVITY_LOG_ERROR("{}", exception.what());
        POLYHEDRAL_GRAVITY_LOG_ERROR("{} (see --help)", program.usage());
        return 1;
    }

    const Polyhedron polyhedron{meshFiles, density, orientation, integrity, metricUnit};

    // The sample box is either given explicitly, or it is the mesh's bounding box about its centre, so
    // that the same command line produces a sensible sample for a unit cube and for a 500 m asteroid
    Array3 minimum{};
    Array3 maximum{};
    if (auto explicitBox = program.present<std::vector<double>>("--sample-box")) {
        minimum = {(*explicitBox)[0], (*explicitBox)[0], (*explicitBox)[0]};
        maximum = {(*explicitBox)[1], (*explicitBox)[1], (*explicitBox)[1]};
    } else {
        const double scale = program.get<double>("--sample-scale");
        const auto [meshMinimum, meshMaximum] = boundingBox(polyhedron.getVertices());
        for (size_t axis = 0; axis < 3; ++axis) {
            const double centre = 0.5 * (meshMinimum[axis] + meshMaximum[axis]);
            const double halfExtent = 0.5 * (meshMaximum[axis] - meshMinimum[axis]) * scale;
            minimum[axis] = centre - halfExtent;
            maximum[axis] = centre + halfExtent;
        }
    }
    for (size_t axis = 0; axis < 3; ++axis) {
        if (minimum[axis] > maximum[axis]) {
            POLYHEDRAL_GRAVITY_LOG_ERROR("--sample-box: the minimum must not be larger than the maximum");
            POLYHEDRAL_GRAVITY_LOG_ERROR("{} (see --help)", program.usage());
            return 1;
        }
    }
    const std::vector<Array3> computationPoints = samplePoints(pointCount, minimum, maximum, seed);

    POLYHEDRAL_GRAVITY_LOG_INFO("Mesh:        {} faces, {} vertices", polyhedron.countFaces(),
                                polyhedron.countVertices());
    POLYHEDRAL_GRAVITY_LOG_INFO("Points:      {}, seed {}, drawn from [{}, {}] x [{}, {}] x [{}, {}]", pointCount,
                                seed, minimum[0], maximum[0], minimum[1], maximum[1], minimum[2], maximum[2]);
    POLYHEDRAL_GRAVITY_LOG_INFO("Repetitions: {}", repetitions);

    // The reference every configuration is judged against, so that a change to the kernel's arithmetic
    // can be judged on accuracy and not only on runtime. It is evaluated once, not once per row.
    std::vector<GravityModelResult> reference{};
    if (!skipAccuracyCheck) {
        const GravityEvaluable referenceEvaluable{polyhedron, ComputeBackend::CPU_PARALLEL,
                                                  ComputePrecision::FLOAT64};
        reference = std::get<std::vector<GravityModelResult>>(referenceEvaluable(computationPoints));
    }

    std::vector<Measurement> measurements{};
    for (size_t backendIndex = 0; backendIndex < backends.size(); ++backendIndex) {
        for (size_t precisionIndex = 0; precisionIndex < precisions.size(); ++precisionIndex) {
            const std::string &backendName = backendNames[backendIndex];
            const std::string &precisionName = precisionNames[precisionIndex];
            const GravityEvaluable evaluable{polyhedron, backends[backendIndex], precisions[precisionIndex]};

            // One untimed evaluation, so that neither the lazily created kernel nor a cold cache is measured
            const auto warmup = std::get<std::vector<GravityModelResult>>(evaluable(computationPoints));

            const auto [worstPotential, worstAcceleration] =
                    skipAccuracyCheck ? std::pair{0.0, 0.0} : worstRelativeError(warmup, reference);

            double best = std::numeric_limits<double>::max();
            double total = 0.0;
            for (size_t repetition = 0; repetition < repetitions; ++repetition) {
                const auto start = std::chrono::high_resolution_clock::now();
                const auto result = evaluable(computationPoints);
                const auto end = std::chrono::high_resolution_clock::now();
                const double milliseconds = std::chrono::duration<double, std::milli>(end - start).count();
                best = std::min(best, milliseconds);
                total += milliseconds;
                (void) result;
            }

            measurements.push_back({backendName, precisionName, best, total / static_cast<double>(repetitions),
                                    worstPotential, worstAcceleration, std::get<0>(warmup.front())});
        }
    }

    const auto perPoint = [pointCount](const double milliseconds) {
        return milliseconds * 1e3 / static_cast<double>(pointCount);
    };

    // The table is assembled line by line and handed to the logger as one message each, so that every
    // row carries the same prefix as the rest of the library's output
    std::ostringstream header{};
    header << std::left << std::setw(9) << "Backend" << std::setw(6) << "Prec." << std::right << std::setw(14)
           << "Best [ms]" << std::setw(14) << "Mean [ms]" << std::setw(16) << "Best [us/pt]" << std::setw(16)
           << "Mean [us/pt]";
    if (!skipAccuracyCheck) {
        header << std::setw(14) << "MaxRelErr V" << std::setw(14) << "MaxRelErr a";
    }
    header << std::setw(16) << "Checksum";
    POLYHEDRAL_GRAVITY_LOG_INFO("{}", header.str());

    for (const Measurement &measurement: measurements) {
        std::ostringstream row{};
        row << std::left << std::setw(9) << measurement.backend << std::setw(6) << measurement.precision
            << std::right << std::scientific << std::setprecision(4) << std::setw(14) << measurement.best
            << std::setw(14) << measurement.mean << std::setw(16) << perPoint(measurement.best) << std::setw(16)
            << perPoint(measurement.mean);
        if (!skipAccuracyCheck) {
            row << std::setw(14) << measurement.worstPotential << std::setw(14) << measurement.worstAcceleration;
        }
        row << std::setw(16) << measurement.checksum;
        POLYHEDRAL_GRAVITY_LOG_INFO("{}", row.str());
    }
    return 0;
}
