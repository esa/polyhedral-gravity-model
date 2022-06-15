import analytical_cube_gravity
import numpy as np
from tqdm import tqdm

X = np.arange(-5, 5.5, 0.5)
Y = np.arange(-5, 5.5, 0.5)
Z = np.arange(-5, 5.5, 0.5)

computation_points = np.array(np.meshgrid(X, Y, Z)).T.reshape(-1, 3)

print(len(computation_points))

with open("analytic_cube_solution.txt", "w") as file:
    for x, y, z in tqdm(computation_points):
        potential = analytical_cube_gravity.evaluate_potential(x, y, z)
        acceleration = analytical_cube_gravity.evaluate_acceleration(x, y, z)
        print("{:.2f} {:.2f} {:.2f} {} {} {} {}".format(x, y, z, potential, acceleration[0],
                                                        acceleration[1], acceleration[2]), file=file)
