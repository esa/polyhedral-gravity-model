import numpy as np
from polyhedral_gravity import Polyhedron, evaluate, PolyhedronIntegrity
import mesh_plotting
import mesh_utility


# Meshes from https://github.com/darioizzo/geodesyNets/tree/master/3dmeshes
vertices, faces = mesh_utility.read_pk_file("../mesh/eros.pk")
vertices_lp, faces_lp = mesh_utility.read_pk_file("../mesh/eros_lp.pk")
vertices, faces = np.array(vertices), np.array(faces)
DENSITY = 1.0
VALUES = np.arange(-2, 2.01, 0.01)

eros = Polyhedron(
    polyhedral_source=(vertices, faces),
    density=DENSITY,
    integrity_check=PolyhedronIntegrity.DISABLE,
)

########################################################################################################################
# Triangulation
########################################################################################################################

mesh_plotting.plot_triangulation(vertices_lp, faces_lp, "Triangulation of Eros",
                                 "../figures/eros/eros_triangulation.png")

########################################################################################################################
# Plot of the potential and Acceleration in XY Plane (Z = 0)
########################################################################################################################
computation_points = np.array(np.meshgrid(VALUES, VALUES, [0])).T.reshape(-1, 3)
gravity_results = evaluate(eros, computation_points)

potentials = -1 * np.array([i[0] for i in gravity_results])
potentials = potentials.reshape((len(VALUES), len(VALUES)))

X = computation_points[:, 0].reshape(len(VALUES), -1)
Y = computation_points[:, 1].reshape(len(VALUES), -1)

mesh_plotting.plot_grid_2d(X, Y, potentials, "Potential of Eros XY-Plane ($z=0$)",
                           "../figures/eros/eros_potential_2d_xy.png", vertices=vertices, coordinate=2)
mesh_plotting.plot_grid_3d(X, Y, potentials, "Potential of Eros XY-Plane ($z=0$)",
                           "../figures/eros/eros_potential_3d_xy.png")

# We just want a slice of x, y values for z = 0
accelerations = np.array([i[1][:] for i in gravity_results])
acc_xy = np.delete(accelerations, 2, 1)
acc_xy = -1 * acc_xy

mesh_plotting.plot_quiver(X, Y, acc_xy, "Acceleration in $x$ and $y$ direction for $z=0$",
                          "../figures/eros/eros_field_xy.png", vertices=vertices, coordinate=2)

########################################################################################################################
# Plot of the potential and Acceleration in XZ Plane (Y = 0)
########################################################################################################################
computation_points = np.array(np.meshgrid(VALUES, [0], VALUES)).T.reshape(-1, 3)
gravity_results = evaluate(eros, computation_points)

potentials = -1 * np.array([i[0] for i in gravity_results])
potentials = potentials.reshape((len(VALUES), len(VALUES)))

X = computation_points[:, 0].reshape(len(VALUES), -1)
Z = computation_points[:, 2].reshape(len(VALUES), -1)

mesh_plotting.plot_grid_2d(X, Z, potentials,
                           "Potential of Eros XZ-Plane ($y=0$)", "../figures/eros/eros_potential_2d_xz.png",
                           labels=('$x$', '$z$'), vertices=vertices, coordinate=1)
mesh_plotting.plot_grid_3d(X, Z, potentials,
                           "Potential of Eros XZ-Plane ($y=0$)", "../figures/eros/eros_potential_3d_xz.png",
                           labels=('$x$', '$z$'))

# We just want a slice of x, y values for z = 0
accelerations = np.array([i[1][:] for i in gravity_results])
acc_xy = np.delete(accelerations, 1, 1)
acc_xy = -1 * acc_xy

mesh_plotting.plot_quiver(X, Z, acc_xy, "Acceleration in $x$ and $z$ direction for $y=0$",
                          "../figures/eros/eros_field_xz.png", labels=('$x$', '$z$'),
                          vertices=vertices, coordinate=1)

########################################################################################################################
# Plot of the potential and Acceleration in YZ Plane (X = 0)
########################################################################################################################
computation_points = np.array(np.meshgrid([0], VALUES, VALUES)).T.reshape(-1, 3)
gravity_results = evaluate(eros, computation_points)

potentials = -1 * np.array([i[0] for i in gravity_results])
potentials = potentials.reshape((len(VALUES), len(VALUES)))

X = computation_points[:, 1].reshape(len(VALUES), -1)
Z = computation_points[:, 2].reshape(len(VALUES), -1)

mesh_plotting.plot_grid_2d(Y, Z, potentials,
                           "Potential of Eros YZ-Plane ($x=0$)", "../figures/eros/eros_potential_2d_yz.png",
                           labels=('$y$', '$z$'), vertices=vertices, coordinate=0)
mesh_plotting.plot_grid_3d(Y, Z, potentials,
                           "Potential of Eros YZ-Plane ($x=0$)", "../figures/eros/eros_potential_3d_yz.png",
                           labels=('$y$', '$z$'))

# We just want a slice of x, y values for z = 0
accelerations = np.array([i[1][:] for i in gravity_results])
acc_xy = np.delete(accelerations, 0, 1)
acc_xy = -1 * acc_xy

mesh_plotting.plot_quiver(Y, Z, acc_xy, "Acceleration in $y$ and $z$ direction for $x=0$",
                          "../figures/eros/eros_field_yz.png", labels=('$y$', '$z$'),
                          vertices=vertices, coordinate=0)

