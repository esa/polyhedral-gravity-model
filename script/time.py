#!python3
from polyhedral_gravity import evaluate, Polyhedron, PolyhedronIntegrity, GravityEvaluable
import numpy as np
import matplotlib.pyplot as plt
import mesh_utility
import timeit
import pickle


def main():
    vertices, faces = mesh_utility.read_pk_file("/Users/schuhmaj/Programming/polyhedral-gravity-model/script/mesh/Eros.pk")
    vertices, faces = np.array(vertices), np.array(faces)
    polyhedron = Polyhedron(
        polyhedral_source=(vertices, faces),
        density=1.0,
        integrity_check=PolyhedronIntegrity.DISABLE,
    )

    # Generate 1000 random cartesian points
    N = 1000
    computation_points = np.random.uniform(-2, 2, (N, 3))

    ########################################################################################################################
    # OLD
    ########################################################################################################################
    start_time = timeit.default_timer()
    for i in range(N):
        evaluate(polyhedron, computation_points[i])
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / N * 1e6
    # Print the time in milliseconds
    print("######## Old Single-Point ########")
    print(f"--> {N} times 1 point with old interface")
    print(f"--> Time taken: {delta:.3f} microseconds per point")

    start_time = timeit.default_timer()
    evaluate(polyhedron, computation_points)
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / N * 1e6
    # Print the time in milliseconds
    print("######## Old Multi-Point ########")
    print("--> 1 time N points with old interface")
    print(f"--> Time taken: {delta:.3f} microseconds per point")

    ########################################################################################################################
    # NEW
    ########################################################################################################################
    evaluable = GravityEvaluable(polyhedron)

    start_time = timeit.default_timer()
    for i in range(N):
        evaluable(computation_points[i])
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / N * 1e6
    # Print the time in milliseconds
    print("######## GravityEvaluable (Single-Point) #########")
    print(f"--> {N} times 1 point with GravityEvaluable")
    print(f"--> Time taken: {delta:.3f} microseconds per point")

    evaluable = GravityEvaluable(polyhedron)

    start_time = timeit.default_timer()
    evaluable(computation_points)
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / N * 1e6
    # Print the time in milliseconds
    print("######## GravityEvaluable (Multi-Point) ########")
    print(f"--> 1 time {N} points with GravityEvaluable")
    print(f"--> Time taken: {delta:.3f} microseconds per point")


    ########################################################################################################################
    # PICKLE SUPPORT
    ########################################################################################################################


    with open("evaluable.pk", "wb") as f:
        pickle.dump(evaluable, f, pickle.HIGHEST_PROTOCOL)

    with open("evaluable.pk", "rb") as f:
        evaluable2 = pickle.load(f)


    start_time = timeit.default_timer()
    evaluable2(computation_points)
    end_time = timeit.default_timer()

    delta = (end_time - start_time) / N * 1e6
    # Print the time in milliseconds
    print("######## Pickle ########")
    print(f"--> 1 time {N} points with GravityEvaluable")
    print(f"--> Time taken: {delta:.3f} microseconds per point")


if __name__ == '__main__':
    main()