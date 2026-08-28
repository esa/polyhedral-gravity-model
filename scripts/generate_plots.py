#!/usr/bin/env python3
"""Generates the potential/ acceleration figures for the example bodies.

This merges the former ``cube_plots.py``, ``eros_plots.py`` and ``torus_plots.py``, which only
differed in the mesh they loaded and in the planes they sliced. The plotting primitives themselves
live next to the notebooks in ``examples/notebooks/mesh_plotting.py``, since they are shared with
them.

Usage:
    ./generate_plots.py                     # all bodies into figures/<body>/
    ./generate_plots.py cube torus          # only a subset
    ./generate_plots.py --output-dir /tmp   # somewhere else
"""
import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Sequence, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "examples" / "data"

# The plotting helpers are shared with the notebooks and therefore kept next to them
sys.path.insert(0, str(REPO_ROOT / "examples" / "notebooks"))

import mesh_plotting  # noqa: E402
import mesh_utility  # noqa: E402
from polyhedral_gravity import Polyhedron, PolyhedronIntegrity, evaluate  # noqa: E402

DENSITY = 1.0
VALUES = np.arange(-2, 2.01, 0.01)

# The three axis-aligned slices through the origin: name, the two axes spanning the plane, and the
# index of the axis which is held at zero (also the one dropped from the acceleration vectors).
PLANES: Dict[str, Tuple[Tuple[int, int], int]] = {
    "xy": ((0, 1), 2),
    "xz": ((0, 2), 1),
    "yz": ((1, 2), 0),
}
AXIS_LABELS = ("$x$", "$y$", "$z$")

CUBE_VERTICES = np.array([
    [-1, -1, -1],
    [1, -1, -1],
    [1, 1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [1, 1, 1],
    [-1, 1, 1],
])

CUBE_FACES = np.array([
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
    [4, 6, 7],
])


@dataclass(frozen=True)
class Body:
    """One example body, i.e. its mesh and which figures to produce for it."""
    name: str
    pretty_name: str
    # Returns (vertices, faces) of the mesh to evaluate and of the (possibly coarser) mesh to draw
    load: Callable[[], Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    planes: Sequence[str] = field(default=("xy", "xz", "yz"))
    # The cube is drawn as its analytic outline, the mesh bodies as the silhouette of their vertices
    plot_rectangle: bool = False


def _load_cube() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return CUBE_VERTICES, CUBE_FACES, CUBE_VERTICES, CUBE_FACES


def _load_pk(filename: str, low_poly_filename: str = None):
    """Loads a mesh from a .pk file, and optionally a coarser one used for the triangulation plot."""

    def load() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        vertices, faces = mesh_utility.read_pk_file(str(DATA_DIR / filename))
        vertices, faces = np.array(vertices), np.array(faces)
        if low_poly_filename is None:
            return vertices, faces, vertices, faces
        vertices_lp, faces_lp = mesh_utility.read_pk_file(str(DATA_DIR / low_poly_filename))
        return vertices, faces, np.array(vertices_lp), np.array(faces_lp)

    return load


# Meshes from https://github.com/darioizzo/geodesyNets/tree/master/3dmeshes
BODIES: Dict[str, Body] = {
    body.name: body for body in (
        Body("cube", "a Cube", _load_cube, planes=("xy",), plot_rectangle=True),
        Body("eros", "Eros", _load_pk("Eros.pk")),
        Body("torus", "a Torus", _load_pk("torus.pk", "torus_lp.pk")),
    )
}


def plot_plane(polyhedron: Polyhedron, vertices: np.ndarray, body: Body, plane: str, output_dir: Path) -> None:
    """Evaluates the model on the given slice through the origin and writes its three figures."""
    (first, second), fixed = PLANES[plane]
    axes = [VALUES, VALUES, VALUES]
    axes[fixed] = [0]
    computation_points = np.array(np.meshgrid(*axes)).T.reshape(-1, 3)
    gravity_results = evaluate(polyhedron, computation_points)

    # The sign convention of the plots is the one of the potential, not the one of the model
    potentials = -1 * np.array([result[0] for result in gravity_results])
    potentials = potentials.reshape((len(VALUES), len(VALUES)))
    accelerations = np.delete(np.array([result[1][:] for result in gravity_results]), fixed, 1)

    A = computation_points[:, first].reshape(len(VALUES), -1)
    B = computation_points[:, second].reshape(len(VALUES), -1)

    labels = (AXIS_LABELS[first], AXIS_LABELS[second])
    plane_title = f"{plane.upper()}-Plane (${AXIS_LABELS[fixed].strip('$')}=0$)"
    outline = {"plot_rectangle": True} if body.plot_rectangle else {"vertices": vertices, "coordinate": fixed}

    mesh_plotting.plot_grid_2d(A, B, potentials, f"Potential of {body.pretty_name} {plane_title}",
                               str(output_dir / f"{body.name}_potential_2d_{plane}.png"),
                               labels=labels, **outline)
    mesh_plotting.plot_grid_3d(A, B, potentials, f"Potential of {body.pretty_name} {plane_title}",
                               str(output_dir / f"{body.name}_potential_3d_{plane}.png"),
                               labels=labels)
    mesh_plotting.plot_quiver(A, B, accelerations,
                              f"Acceleration in {labels[0]} and {labels[1]} direction "
                              f"for ${AXIS_LABELS[fixed].strip('$')}=0$",
                              str(output_dir / f"{body.name}_field_{plane}.png"),
                              labels=labels, **outline)


def plot_body(body: Body, output_root: Path) -> None:
    """Writes the triangulation and every requested plane of one body into output_root/<body>."""
    print(f"Plotting {body.name}")
    vertices, faces, vertices_lp, faces_lp = body.load()
    output_dir = output_root / body.name
    output_dir.mkdir(parents=True, exist_ok=True)

    polyhedron = Polyhedron(
        polyhedral_source=(vertices, faces),
        density=DENSITY,
        integrity_check=PolyhedronIntegrity.DISABLE,
    )

    mesh_plotting.plot_triangulation(vertices_lp, faces_lp, f"Triangulation of {body.pretty_name}",
                                     str(output_dir / f"{body.name}_triangulation.png"))

    for plane in body.planes:
        plot_plane(polyhedron, vertices, body, plane, output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "bodies",
        nargs="*",
        metavar="BODY",
        default=list(BODIES),
        help=f"The bodies to plot, any of {', '.join(BODIES)}. Defaults to all of them.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=REPO_ROOT / "figures",
        help="Directory the figures are written to, one sub-directory per body (default: <repo>/figures).",
    )
    args = parser.parse_args()

    unknown = [name for name in args.bodies if name not in BODIES]
    if unknown:
        parser.error(f"unknown body/ bodies {', '.join(unknown)}; choose from {', '.join(BODIES)}")

    for name in args.bodies:
        plot_body(BODIES[name], args.output_dir)


if __name__ == "__main__":
    main()
