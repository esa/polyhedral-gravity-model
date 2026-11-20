import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def plot_triangulation(vertices, faces, title, filename=None):
    """Plots the triangulation of mesh

    Args:
        vertices: (N, 3)-array-like of the vertices
        faces: (M, 3)-array-like of the triangular faces
        title (str): of the plot
        filename (str): the path/ filename

    Returns:
        None

    """
    print("Plotting Triangulation")
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(5, 5))
    ax.plot_trisurf(vertices[:, 0],
                    vertices[:, 1],
                    vertices[:, 2],
                    triangles=faces,
                    edgecolor=[[0, 0, 0]],
                    linewidth=1.0,
                    alpha=0.5)
    ax.set_title(title)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$z$")


    if filename is not None:
        fig.savefig(filename, dpi=300)


def plot_grid_2d(X, Y, z, title, filename=None, labels=("$x$", "$y$"), limits=(-2, 2), plot_rectangle=False, vertices=None,
                 coordinate=None):
    """Plots the potential in 2D for a given cross-section.

    Args:
        X: first coordinate array
        Y: second coordinate array
        z: potential array
        title (str): the title of the plot
        filename (str): the filename/ path
        labels (str, str): the axes label names
        limits (float, float): the limits of the axes
        plot_rectangle (bool, optional): plot a rectangle centered on the origin
        vertices: the vertices of the polyhedron, if given, plot the polyhedron inside
        coordinate: the coordinate of the cross-section (lazy parameter, could be avoided)

    Returns:
        None

    """
    print("Plotting Grid 2D")
    fig, ax = plt.subplots(figsize=(5, 5))

    Z = z.reshape((len(X), len(Y)))

    surf = ax.contourf(X, Y, Z, cmap=cm.viridis)

    if plot_rectangle:
        rect = patches.Rectangle((-1, -1), 2, 2, linewidth=1, edgecolor='r', facecolor='none')
        ax.add_patch(rect)

    if vertices is not None and coordinate is not None:
        plot_polygon(ax, vertices, coordinate)

    left, right = limits
    ax.set_xlim(left, right)
    ax.set_ylim(left, right)

    ax.axis('equal')

    xl, yl = labels
    ax.set_title(title)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)

    fig.colorbar(surf, aspect=5, orientation='vertical', shrink=0.5)

    if filename is not None:
        fig.savefig(filename, dpi=300)


def plot_grid_3d(X, Y, z, title, filename=None, labels=("$x$", "$y$")):
    """Plots the potential in 3D

    Args:
        X: first coordinate array
        Y: second coordinate array
        z: potential array
        title (str): the title of the plot
        filename (str): the filename/ path
        labels (str, str): the axes label names

    Returns:
        None

    """
    print("Plotting Grid 3D")
    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={"projection": "3d"})

    Z = z.reshape((len(X), len(Y)))

    surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis,
                           linewidth=0, antialiased=False)

    xl, yl = labels
    ax.set_title(title)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)

    # fig.colorbar(surf, aspect=5, orientation='vertical', shrink=0.5)

    fig.set_dpi(300)

    if filename is not None:
        fig.savefig(filename, dpi=300)


def plot_quiver(X, Y, xy, title, filename=None, labels=("$x$", "$y$"), limits=(-2, 2), plot_rectangle=False, vertices=None,
                coordinate=None):
    """Plots a quiver plot, given the accelerations.

    Args:
        X: first coordinate array
        Y: second coordinate array
        xy: acceleration array
        title (str): the title of the plot
        filename (str): the filename/ path
        labels (str, str): the axes label names
        limits (float, float): the limits of the axes
        plot_rectangle (bool, optional): plot a rectangle centered on the origin
        vertices: the vertices of the polyhedron, if given, plot the polyhedron inside
        coordinate: the coordinate of the cross-section (lazy parameter, could be avoided)

    Returns:
        None

    """
    print("Plotting Quiver")
    fig, ax = plt.subplots(figsize=(5, 5))

    U = np.reshape(xy[:, 0], (len(X), -1))
    V = np.reshape(xy[:, 1], (len(Y), -1))

    ax.quiver(X, Y, U, V, angles='xy', linewidth=0.1, color='b', pivot='mid', units='xy')

    if plot_rectangle:
        rect = patches.Rectangle((-1, -1), 2, 2, linewidth=1, edgecolor='r', facecolor='none')
        ax.add_patch(rect)

    if vertices is not None and coordinate is not None:
        plot_polygon(ax, vertices, coordinate)

    left, right = limits
    ax.set_xlim(left, right)
    ax.set_ylim(left, right)

    ax.axis('equal')

    xl, yl = labels
    ax.set_title(title)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)

    if filename is not None:
        fig.savefig(filename, dpi=300)


# Switch to turn plotting the polygon of, if one is to lazy to refactor code
ACTIVATE = False


def plot_polygon(ax, vertices, coordinate):
    """Plots a polygon fitting to the polyhedron inside the cross-section. Not a very clean solution but nice to
    get some feeling about the polyhedron.

    Args:
        ax: matplotlib ax
        vertices: the vertices
        coordinate (int): the coordinate (0=X, 1=Y, 2=Z)

    Returns:
        None

    """
    if ACTIVATE:
        plane_points = vertices[abs(vertices[:, coordinate]) < 1e10]
        xy = np.delete(plane_points, coordinate, 1)
        polygon = patches.Polygon(xy, closed=True, linewidth=1, edgecolor='r', facecolor='none')
        ax.add_patch(polygon)
