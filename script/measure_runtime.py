#!python3
from polyhedral_gravity import evaluate, Polyhedron, PolyhedronIntegrity, GravityEvaluable
import polyhedral_gravity
import numpy as np
import matplotlib.pyplot as plt
import timeit
import argparse
from loguru import logger
from typing import Dict


def run_time_measurements(sample_size: int, mesh_files: list[str]) -> Dict[str, float]:
    """Returns the RuntimeMeasurements for the four configurations as mapping.

    Returns:
        a mapping from configuration name to runtime per point in microseconds
    """
    results = dict()
    polyhedron = Polyhedron(
        polyhedral_source=mesh_files,
        density=1.0,
        integrity_check=PolyhedronIntegrity.DISABLE,
    )

    # Generate 1000 random cartesian points
    computation_points = np.random.uniform(-2, 2, (sample_size, 3))

    start_time = timeit.default_timer()
    for i in range(sample_size):
        evaluate(polyhedron, computation_points[i])
    end_time = timeit.default_timer()

    total_time = end_time - start_time
    delta = total_time / sample_size * 1e6
    logger.info("##########################################################")
    logger.info(f"Benchmarking a {sample_size} times one point with evaluate()")
    logger.info(f"Total Runtime: {total_time:.3f} seconds")
    logger.info(f"Time taken:    {delta:.3f} microseconds per point")
    results[f"evaluate \n ${sample_size} \\times 1$ point"] = delta

    start_time = timeit.default_timer()
    evaluate(polyhedron, computation_points)
    end_time = timeit.default_timer()

    total_time = end_time - start_time
    delta = total_time / sample_size * 1e6
    logger.info("##########################################################")
    logger.info(f"Benchmarking one time {sample_size} points with evaluate()")
    logger.info(f"Total Runtime: {total_time:.3f} seconds")
    logger.info(f"Time taken:    {delta:.3f} microseconds per point")
    results[f"evaluate \n $1 \\times {sample_size}$ points"] = delta

    evaluable = GravityEvaluable(polyhedron)

    start_time = timeit.default_timer()
    for i in range(sample_size):
        evaluable(computation_points[i])
    end_time = timeit.default_timer()

    total_time = end_time - start_time
    delta = total_time / sample_size * 1e6
    logger.info("##########################################################")
    logger.info(f"Benchmarking a {sample_size} times one point with GravityEvaluable")
    logger.info(f"Total Runtime: {total_time:.3f} seconds")
    logger.info(f"Time taken:    {delta:.3f} microseconds per point")
    results[f"GravityEvaluable \n ${sample_size} \\times 1$ point"] = delta

    evaluable = GravityEvaluable(polyhedron)

    start_time = timeit.default_timer()
    evaluable(computation_points)
    end_time = timeit.default_timer()

    total_time = end_time - start_time
    delta = total_time / sample_size * 1e6
    logger.info("##########################################################")
    logger.info(f"Benchmarking one time {sample_size} points with GravityEvaluable")
    logger.info(f"Total Runtime: {total_time:.3f} seconds")
    logger.info(f"Time taken:    {delta:.3f} microseconds per point")
    logger.info("##########################################################")
    results[f"GravityEvaluable \n $1 \\times {sample_size}$ points"] = delta
    return results


def create_plot(runtime_results: Dict[str, float], sample_size: int) -> None:
    """Creates a bar chart plot with given runtime results dictionary"""
    # Legend for X-axis
    configurations = list(runtime_results.keys())
    # Figure of runtime measurements per call (in microseconds)
    runtime_per_call = list(runtime_results.values())

    # Create a numpy range equal to configurations list length
    x_pos = np.arange(len(configurations))

    fig, ax = plt.subplots(figsize=(7, 4))

    # Now, use the cmap to pick colors for each bar
    ax.bar(x_pos, runtime_per_call, align='center', color=["darkorange", "gold", "navy", "blue"])
    ax.grid(True)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(configurations)
    ax.set_ylabel(r'Runtime per Point $[\mu s]$')
    ax.set_title(f'Runtime, Sample Size of ${sample_size}$, Interface v{polyhedral_gravity.__version__}, Parallelization {polyhedral_gravity.__parallelization__}')

    # Save the figure
    plt.savefig("runtime_measurements.png", dpi=300)


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
    logger.info(f"Plotting Results:  {'Yes' if args.plot else 'No'}")
    logger.info("##########################################################")
    logger.info("Polyhedral Gravity Model Information:")
    logger.info(f"Version:                 {polyhedral_gravity.__version__}")
    logger.info(f"Parallelization Backend: {polyhedral_gravity.__parallelization__}")
    logger.info(f"Commit Hash:             {polyhedral_gravity.__commit__}")
    logger.info(f"Logging Level:           {polyhedral_gravity.__logging__}")
    logger.info("##########################################################")
    results = run_time_measurements(args.sample_size, args.mesh_files)
    if args.plot:
        logger.info("Plotting Results")
        create_plot(results, args.sample_size)


if __name__ == "__main__":
    main()
