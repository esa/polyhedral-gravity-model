import numpy as np
import polyhedral_gravity as gravity_model

import mesh_plotting

cube_vertices = np.array([
    [-1, -1, -1],
    [1, -1, -1],
    [1, 1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [1, 1, 1],
    [-1, 1, 1]
])

cube_faces = np.array([
    [1, 3, 2],
    [0, 3, 1],
    [0, 1, 5],
    [0, 5, 4],
    [0, 7, 3],
    [0, 4, 7],
    [1, 2, 6],
    [1, 6, 5],
    [2, 3, 6],
    [3, 7, 6],
    [4, 5, 6],
    [4, 6, 7]
])

DENSITY = 1.0
VALUES = np.arange(-2, 2.01, 0.01)

########################################################################################################################
# Triangulation
########################################################################################################################

mesh_plotting.plot_triangulation(cube_vertices, cube_faces, "Triangulation of a Cube",
                                 "../figures/cube/cube_triangulation.png")

########################################################################################################################
# Plot of the potential and Acceleration in XY Plane (Z = 0)
########################################################################################################################
computation_points = np.array(np.meshgrid(VALUES, VALUES, [0])).T.reshape(-1, 3)
gravity_results = gravity_model.evaluate(cube_vertices, cube_faces, DENSITY, computation_points)

potentials = -1 * np.array([i[0] for i in gravity_results])
potentials = potentials.reshape((len(VALUES), len(VALUES)))

X = computation_points[:, 0].reshape(len(VALUES), -1)
Y = computation_points[:, 1].reshape(len(VALUES), -1)

mesh_plotting.plot_grid_2d(X, Y, potentials, "Potential of a Cube ($z=0$)", "../figures/cube/cube_potential_2d.png",
                           plot_rectangle=True)
mesh_plotting.plot_grid_3d(X, Y, potentials, "Potential of a Cube ($z=0$)", "../figures/cube/cube_potential_3d.png")

# We just want a slice of x, y values for z = 0
accelerations = np.array([i[1][:] for i in gravity_results])
acc_xy = np.delete(accelerations, 2, 1)
acc_xy = -1 * acc_xy

mesh_plotting.plot_quiver(X, Y, acc_xy, "Acceleration in $x$ and $y$ direction for $z=0$",
                          "../figures/cube/cube_field.png", plot_rectangle=True)
