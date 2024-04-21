#!python3
from polyhedral_gravity import evaluate, Polyhedron, PolyhedronIntegrity, GravityEvaluable
import numpy as np
import matplotlib.pyplot as plt
import mesh_utility
import timeit
import pickle
from typing import Dict


def run_time_measurements(sample_size: int = 1000) -> Dict[str, float]:
    """Returns the RuntimeMeasurements for the four configurations as mapping.

    Returns:
        a mapping from configuration name to runtime per point in microseconds
    """
    results = dict()
    vertices, faces = mesh_utility.read_pk_file("./mesh/Eros.pk")
    vertices, faces = np.array(vertices), np.array(faces)
    polyhedron = Polyhedron(
        polyhedral_source=(vertices, faces),
        density=1.0,
        integrity_check=PolyhedronIntegrity.DISABLE,
    )

    # Generate 1000 random cartesian points
    computation_points = np.random.uniform(-2, 2, (sample_size, 3))

    start_time = timeit.default_timer()
    for i in range(sample_size):
        evaluate(polyhedron, computation_points[i])
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / sample_size * 1e6
    # Print the time in milliseconds
    print("######## Old Single-Point ########")
    print(f"--> {sample_size} times 1 point with old interface")
    print(f"--> Time taken: {delta:.3f} microseconds per point")
    results[f"evaluate \n ${sample_size} \\times 1$ point"] = delta

    start_time = timeit.default_timer()
    evaluate(polyhedron, computation_points)
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / sample_size * 1e6
    # Print the time in milliseconds
    print("######## Old Multi-Point ########")
    print(f"--> 1 time {sample_size} points with old interface")
    print(f"--> Time taken: {delta:.3f} microseconds per point")
    results[f"evaluate \n $1 \\times {sample_size}$ points"] = delta

    evaluable = GravityEvaluable(polyhedron)

    start_time = timeit.default_timer()
    for i in range(sample_size):
        evaluable(computation_points[i])
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / sample_size * 1e6
    # Print the time in milliseconds
    print("######## GravityEvaluable (Single-Point) #########")
    print(f"--> {sample_size} times 1 point with GravityEvaluable")
    print(f"--> Time taken: {delta:.3f} microseconds per point")
    results[f"GravityEvaluable \n ${sample_size} \\times 1$ point"] = delta

    evaluable = GravityEvaluable(polyhedron)

    start_time = timeit.default_timer()
    evaluable(computation_points)
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / sample_size * 1e6
    # Print the time in milliseconds
    print("######## GravityEvaluable (Multi-Point) ########")
    print(f"--> 1 time {sample_size} points with GravityEvaluable")
    print(f"--> Time taken: {delta:.3f} microseconds per point")
    results[f"GravityEvaluable \n $1 \\times {sample_size}$ points"] = delta
    return results


def create_plot(runtime_results: Dict[str, float], sample_size: int = 1000) -> None:
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
    ax.set_ylabel('Runtime per Point $[\mu s]$')
    ax.set_title(f'Runtime Measurements (Sample Size = ${sample_size}$) of the Python Interface v3.0')

    # Save the figure
    plt.savefig("runtime_measurements.png", dpi=300)


if __name__ == '__main__':
    results = run_time_measurements()
    create_plot(results)