#!python3
import polyhedral_gravity
import numpy as np
import matplotlib.pyplot as plt
import mesh_utility
import timeit

vertices, faces = mesh_utility.read_pk_file("/Users/schuhmaj/Programming/polyhedral-gravity-model/script/mesh/Eros.pk")
vertices, faces = np.array(vertices), np.array(faces)

# Generate 1000 random cartesian points
N = 1000
computation_points = np.random.uniform(-2, 2, (N, 3))

DENSITY = 1.0

########################################################################################################################
# OLD
########################################################################################################################
start_time = timeit.default_timer()
for i in range(N):
    polyhedral_gravity.evaluate((vertices, faces), DENSITY, computation_points[i])
end_time = timeit.default_timer()

delta = (end_time - start_time) / N * 1e6
# Print the time in milliseconds
print("######## Old Single-Point ########")
print(f"--> {N} times 1 point with old interface")
print(f"--> Time taken: {delta:.3f} microseconds per point")

start_time = timeit.default_timer()
polyhedral_gravity.evaluate((vertices, faces), DENSITY, computation_points)
end_time = timeit.default_timer()

delta = (end_time - start_time) / N * 1e6
# Print the time in milliseconds
print("######## Old Multi-Point ########")
print("--> 1 time N points with old interface")
print(f"--> Time taken: {delta:.3f} microseconds per point")

########################################################################################################################
# NEW
########################################################################################################################
evaluable = polyhedral_gravity.GravityEvaluable((vertices, faces), DENSITY)

start_time = timeit.default_timer()
for i in range(N):
    evaluable(computation_points[i])
end_time = timeit.default_timer()

delta = (end_time - start_time) / N * 1e6
# Print the time in milliseconds
print("######## GravityEvaluable (Single-Point) #########")
print(f"--> {N} times 1 point with GravityEvaluable")
print(f"--> Time taken: {delta:.3f} microseconds per point")

evaluable = polyhedral_gravity.GravityEvaluable((vertices, faces), DENSITY)

start_time = timeit.default_timer()
evaluable(computation_points)
end_time = timeit.default_timer()

delta = (end_time - start_time) / N * 1e6
# Print the time in milliseconds
print("######## GravityEvaluable (Multi-Point) ########")
print(f"--> 1 time {N} points with GravityEvaluable")
print(f"--> Time taken: {delta:.3f} microseconds per point")


