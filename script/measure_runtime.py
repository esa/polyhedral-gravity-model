#!/usr/bin/env python3
from polyhedral_gravity import evaluate, Polyhedron, PolyhedronIntegrity, GravityEvaluable, ComputePrecision, ComputeBackend
import polyhedral_gravity
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
import timeit
import argparse
from loguru import logger
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

# The ways of calling the model, ordered from the most discouraged to the recommended one, which is also
# the order in which they appear within a backend group and the direction in which their colour saturates.
#
# `evaluate` called once with all N points is deliberately absent: GravityModel::evaluate constructs a
# GravityEvaluable, calls it once, and drops it again, so that configuration measures the same thing as
# `GravityEvaluable 1 x N` plus one construction. The remaining three are genuinely different workloads.
CALL_STYLE_MANY_SINGLE = "N x 1"
CALL_STYLE_ONE_BULK = "1 x N"

OPTIONS: List[Tuple[str, str]] = [
    ("evaluate", CALL_STYLE_MANY_SINGLE),
    ("GravityEvaluable", CALL_STYLE_MANY_SINGLE),
    ("GravityEvaluable", CALL_STYLE_ONE_BULK),
]

PRECISIONS: List[Tuple[str, ComputePrecision]] = [
    ("FP32", ComputePrecision.FLOAT32),
    ("FP64", ComputePrecision.FLOAT64),
]

BACKENDS: List[Tuple[str, ComputeBackend]] = [
    ("CPU_SERIAL", ComputeBackend.CPU_SERIAL),
    ("CPU_PARALLEL", ComputeBackend.CPU_PARALLEL),
    ("GPU_PARALLEL", ComputeBackend.GPU_PARALLEL),
]

# Colour carries both remaining levels at once: the hue is the precision, blue for FLOAT32 and red for
# FLOAT64, and the saturation is the option. The pale bars are the way of calling the model one should not
# use, the saturated ones the way one should, so that the recommendation is legible without the legend.
PRECISION_COLORMAPS = {"FP32": plt.get_cmap("Blues"), "FP64": plt.get_cmap("Reds")}
OPTION_SATURATIONS = [0.38, 0.62, 0.92]

OPTION_SHORT_LABELS = {
    ("evaluate", CALL_STYLE_MANY_SINGLE): "evaluate\n$N \\times 1$",
    ("GravityEvaluable", CALL_STYLE_MANY_SINGLE): "Evaluable\n$N \\times 1$",
    ("GravityEvaluable", CALL_STYLE_ONE_BULK): "Evaluable\n$1 \\times N$",
}


def bar_color(option: Tuple[str, str], precision_name: str):
    """The colour of one bar, i.e. the precision's hue at the option's saturation."""
    return PRECISION_COLORMAPS[precision_name](OPTION_SATURATIONS[OPTIONS.index(option)])


@dataclass(frozen=True)
class Configuration:
    """One point of the Cartesian product, i.e. one bar of the plot."""
    interface: str
    call_style: str
    backend_name: str
    precision_name: str


def is_gpu_available() -> bool:
    """Whether this build of the library has a GPU backend.

    ``__parallelization__`` lists the Kokkos execution spaces the module was compiled with, for example
    ``"Serial, OpenMP, Cuda"``. Anything beyond the two host spaces is a device backend.
    """
    spaces = {space.strip().lower() for space in polyhedral_gravity.__parallelization__.split(",")}
    return bool(spaces - {"serial", "openmp", "threads"})


def build_runner(
        configuration: Configuration,
        polyhedron: Polyhedron,
        computation_points: np.ndarray,
        backend: ComputeBackend,
        precision: ComputePrecision,
) -> Callable[[], None]:
    """Returns a nullary callable running exactly the workload of one configuration once.

    The ``GravityEvaluable`` is constructed outside the returned callable, which is the whole point of that
    interface: it keeps the mesh and its caches on the device across calls, whereas ``evaluate`` rebuilds
    them on every invocation.
    """
    if configuration.interface == "evaluate":
        return lambda: [evaluate(polyhedron, point, backend, precision) for point in computation_points]

    evaluable = GravityEvaluable(polyhedron, backend, precision)
    if configuration.call_style == CALL_STYLE_MANY_SINGLE:
        return lambda: [evaluable(point) for point in computation_points]
    return lambda: evaluable(computation_points)


def run_time_measurements(sample_size: int, mesh_files: List[str], repeats: int) -> Dict[Configuration, float]:
    """Measures the Cartesian product of the three options, both precisions, and every available backend.

    Every configuration is timed ``repeats`` times and the fastest run is kept, as the C++ benchmark driver
    does. The mean is the wrong statistic here: everything the operating system does to the process while it
    is being measured can only make a run slower, so the minimum is the closest estimate of what the code
    costs. Without this, `evaluate N x 1` on CPU_PARALLEL -- a thousand OpenMP-backed evaluables constructed
    back to back -- was seen to vary by a factor of three between runs, while every other configuration
    reproduced to within a few percent.

    Returns:
        a mapping from configuration to runtime per computation point in microseconds
    """
    polyhedron = Polyhedron(
        polyhedral_source=mesh_files,
        density=1.0,
        integrity_check=PolyhedronIntegrity.DISABLE,
    )

    # Generate the random cartesian computation points, the same ones for every configuration
    computation_points = np.random.uniform(-2, 2, (sample_size, 3))

    backends = BACKENDS if is_gpu_available() else [entry for entry in BACKENDS if entry[0] != "GPU_PARALLEL"]
    if len(backends) != len(BACKENDS):
        logger.warning("No GPU backend in this build, skipping GPU_PARALLEL")

    results: Dict[Configuration, float] = dict()
    for backend_name, backend in backends:
        logger.info("##########################################################")
        logger.info(f"Backend {backend_name}")
        for interface, call_style in OPTIONS:
            for precision_name, precision in PRECISIONS:
                configuration = Configuration(interface, call_style, backend_name, precision_name)
                runner = build_runner(configuration, polyhedron, computation_points, backend, precision)

                # One untimed run, so that neither the Kokkos runtime's initialization nor the first
                # upload of the mesh to the device is charged to the configuration which happens to run first
                runner()

                timings = []
                for _ in range(repeats):
                    start_time = timeit.default_timer()
                    runner()
                    timings.append(timeit.default_timer() - start_time)

                total_time = min(timings)
                delta = total_time / sample_size * 1e6
                results[configuration] = delta
                logger.info(
                    f"{interface:>16s} {call_style:>5s} {precision_name} | "
                    f"Best {total_time:8.3f} s | {delta:10.3f} us per point | "
                    f"spread {max(timings) / total_time:4.2f}x"
                )
    logger.info("##########################################################")
    return results


def create_plot(runtime_results: Dict[Configuration, float], sample_size: int) -> None:
    """Creates the grouped bar chart of the runtime results.

    The backend is the group on the x-axis and the option is annotated underneath it; colour carries
    the remaining two levels, hue for the precision and saturation for the option. The y-axis is
    logarithmic because a serial FLOAT64 run and a GPU FLOAT32 one are three orders of magnitude apart.
    """
    backend_names = [name for name, _ in BACKENDS if any(key.backend_name == name for key in runtime_results)]
    bars_per_group = len(OPTIONS) * len(PRECISIONS)
    group_width = 0.86
    bar_width = group_width / bars_per_group

    fig, ax = plt.subplots(figsize=(4.2 * len(backend_names), 5.5))

    option_positions, option_labels_below = [], []
    for group_index, backend_name in enumerate(backend_names):
        for option_index, option in enumerate(OPTIONS):
            for precision_index, (precision_name, _) in enumerate(PRECISIONS):
                slot = option_index * len(PRECISIONS) + precision_index
                position = group_index + (slot - (bars_per_group - 1) / 2) * bar_width
                value = runtime_results[Configuration(*option, backend_name, precision_name)]
                ax.bar(
                    position, value, width=bar_width * 0.9,
                    color=bar_color(option, precision_name),
                    edgecolor="black", linewidth=0.5, zorder=3,
                )
                ax.text(
                    position, value * 1.08, f"{value:.1f}",
                    ha="center", va="bottom", rotation=90, fontsize=6, zorder=4,
                )
            # The option's label goes between the FP32 and the FP64 bar of the pair it names
            pair_center = option_index * len(PRECISIONS) + (len(PRECISIONS) - 1) / 2
            option_positions.append(group_index + (pair_center - (bars_per_group - 1) / 2) * bar_width)
            option_labels_below.append(OPTION_SHORT_LABELS[option])

    ax.set_yscale("log")
    ax.grid(True, which="both", axis="y", linestyle=":", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

    # Headroom for the rotated value labels and the two legends above the tallest bar. On a logarithmic
    # axis a constant factor is what a constant amount of space is, hence the multiplication.
    values = list(runtime_results.values())
    ax.set_ylim(bottom=min(values) / 3.0, top=max(values) * 25.0)

    # A thin rule between the backend groups, so that the eye does not have to rely on the gap alone
    for group_index in range(1, len(backend_names)):
        ax.axvline(group_index - 0.5, color="grey", linewidth=0.8, zorder=1)

    # The backend is the only x tick; its label is pushed down to leave room for the option annotations,
    # which are drawn as free text rather than as minor ticks -- matplotlib suppresses a minor tick that
    # coincides with a major one, which is exactly where the middle option of every group sits.
    ax.set_xticks(np.arange(len(backend_names)))
    ax.set_xticklabels(backend_names, fontweight="bold")
    ax.tick_params(axis="x", which="major", length=0, pad=34)
    ax.set_xlim(-0.5, len(backend_names) - 0.5)
    for position, label in zip(option_positions, option_labels_below):
        ax.text(position, -0.015, label, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=6.5)

    ax.set_ylabel(r"Runtime per Point $[\mu s]$")
    ax.set_title(
        f"Runtime, Sample Size $N = {sample_size}$, Interface v{polyhedral_gravity.__version__}\n"
        f"Parallelization {polyhedral_gravity.__parallelization__}, "
        f"Fast Math {'ON' if polyhedral_gravity.__fast_math__ else 'OFF'}",
        fontsize=10,
    )

    # One legend for the whole colour scheme: a row per option, showing that option's FP32 and FP64 swatch
    # side by side, so that both the hue and the saturation are read off the same key
    option_handles = [
        tuple(
            Patch(facecolor=bar_color(option, precision_name), edgecolor="black", linewidth=0.5)
            for precision_name, _ in PRECISIONS
        )
        for option in OPTIONS
    ]
    option_labels = [OPTION_SHORT_LABELS[option].replace("\n", " ") for option in OPTIONS]
    header = "  ".join(name for name, _ in PRECISIONS)
    ax.legend(
        handles=option_handles, labels=option_labels,
        title=f"{header}   (pale = discouraged, saturated = recommended)",
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.4)},
        handlelength=3.0, loc="upper right", fontsize=8, title_fontsize=8, framealpha=0.95,
    )

    fig.tight_layout()
    fig.savefig("runtime_measurements.png", dpi=300)
    logger.info("Wrote runtime_measurements.png")


def main():
    parser = argparse.ArgumentParser(
        description="Command line interface for benchmarking the Python interface of the polyhedral gravity model "
                    "to handle sample size, mesh file inputs, and optional plotting.",
    )
    parser.add_argument(
        '-s', '--sample-size',
        type=int,
        default=1000,
        help="Specify the size of the sample. Defaults to 1000",
    )
    parser.add_argument(
        '-m', '--mesh-files',
        type=str,
        nargs='+',
        default=["mesh/Eros.node", "mesh/Eros.face"],
        help="Input mesh file(s). Provide one or more file paths separated by a space.",
    )
    parser.add_argument(
        '-r', '--repeats',
        type=int,
        default=3,
        help="How often to time each configuration; the fastest run is reported. Defaults to 3",
    )
    parser.add_argument(
        '-p', '--plot',
        action='store_true',
        help="Option to create and display a plot. Use this flag to enable plotting. Defaults to False.",
    )
    args = parser.parse_args()

    logger.info("Benchmarking the Polyhedral Gravity Model")
    logger.info("##########################################################")
    logger.info("Benchmarking Configuration:")
    logger.info(f"Sample Size:       {args.sample_size}")
    logger.info(f"Mesh files:        {args.mesh_files}")
    logger.info(f"Repeats:           {args.repeats}")
    logger.info(f"Plotting Results:  {'Yes' if args.plot else 'No'}")
    logger.info("##########################################################")
    logger.info("Polyhedral Gravity Model Information:")
    logger.info(f"Version:                 {polyhedral_gravity.__version__}")
    logger.info(f"Host Compiler:           {polyhedral_gravity.__host_compiler__}")
    logger.info(f"Device Compiler:         {polyhedral_gravity.__device_compiler__}")
    logger.info(f"Parallelization Backend: {polyhedral_gravity.__parallelization__}")
    logger.info(f"Commit Hash:             {polyhedral_gravity.__commit__}")
    logger.info(f"Fast Math Level:         {'Yes' if polyhedral_gravity.__fast_math__ else 'No'}")
    logger.info(f"Logging Level:           {polyhedral_gravity.__logging__}")
    logger.info("##########################################################")
    results = run_time_measurements(args.sample_size, args.mesh_files, args.repeats)
    if args.plot:
        logger.info("Plotting Results")
        create_plot(results, args.sample_size)


if __name__ == "__main__":
    main()
