from typing import Tuple, List, Union
from polyhedral_gravity import Polyhedron, GravityEvaluable, evaluate, PolyhedronIntegrity, NormalOrientation
import numpy as np
import pickle
import pytest
from pathlib import Path

CUBE_VERTICES = np.array([
    [-1, -1, -1],
    [1, -1, -1],
    [1, 1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [1, 1, 1],
    [-1, 1, 1]
])

CUBE_FACES_OUTWARDS = np.array([
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

CUBE_FACES_INWARDS = CUBE_FACES_OUTWARDS[:, [1, 0, 2]]
CUBE_FACES_OUTWARDS_MAJOR = CUBE_FACES_OUTWARDS.copy()
CUBE_FACES_OUTWARDS_MAJOR[0] = CUBE_FACES_OUTWARDS_MAJOR[0, [1, 0, 2]]


CUBE_FACES_INWARDS_MAJOR = CUBE_FACES_INWARDS.copy()
CUBE_FACES_INWARDS_MAJOR[2] = CUBE_FACES_INWARDS_MAJOR[2, [1, 0, 2]]
CUBE_FACES_INWARDS_MAJOR[3] = CUBE_FACES_INWARDS_MAJOR[3, [1, 0, 2]]
CUBE_FACES_INWARDS_MAJOR[5] = CUBE_FACES_INWARDS_MAJOR[5, [1, 0, 2]]

DENSITY = 1.0

CUBE_VERTICES_FILE = Path("test/resources/cube.node")
CUBE_FACE_FILE = Path("test/resources/cube.face")

@pytest.mark.parametrize(
    "polyhedral_source,normal_orientation,violating_faces,expected_healed", [
        ((CUBE_VERTICES, CUBE_FACES_OUTWARDS), NormalOrientation.OUTWARDS, [], CUBE_FACES_OUTWARDS),
        ((CUBE_VERTICES, CUBE_FACES_INWARDS), NormalOrientation.INWARDS, [], CUBE_FACES_INWARDS),
        ((CUBE_VERTICES, CUBE_FACES_OUTWARDS_MAJOR), NormalOrientation.OUTWARDS, [0], CUBE_FACES_OUTWARDS),
        ((CUBE_VERTICES, CUBE_FACES_INWARDS_MAJOR), NormalOrientation.INWARDS, [2, 3, 5], CUBE_FACES_INWARDS),
    ],
    ids=["CubeOutwards", "CubeInwards", "CubeOutwardsMajor", "CubeInwardsMajor"]
)
def test_polyhedron(
        polyhedral_source: Tuple[np.ndarray, np.ndarray],
        normal_orientation: NormalOrientation,
        violating_faces: List[int],
        expected_healed: np.ndarray,
) -> None:
    inverse_orientation = NormalOrientation.OUTWARDS if NormalOrientation.INWARDS == normal_orientation else NormalOrientation.INWARDS

    # The wrong orientation always leads to an error
    with pytest.raises(ValueError):
        Polyhedron(polyhedral_source, DENSITY, inverse_orientation, PolyhedronIntegrity.AUTOMATIC)
        Polyhedron(polyhedral_source, DENSITY, inverse_orientation, PolyhedronIntegrity.VERIFY)
    # At least one error will raise, no error won't raises
    if len(violating_faces) == 0:
        Polyhedron(polyhedral_source, DENSITY, normal_orientation, PolyhedronIntegrity.AUTOMATIC)
        Polyhedron(polyhedral_source, DENSITY, normal_orientation, PolyhedronIntegrity.VERIFY)
    else:
        with pytest.raises(ValueError):
            Polyhedron(polyhedral_source, DENSITY, normal_orientation, PolyhedronIntegrity.AUTOMATIC)
            Polyhedron(polyhedral_source, DENSITY, normal_orientation, PolyhedronIntegrity.VERIFY)

    # The healed version obeys to the majority ordering, Hence majority OUTWARDS --> only OUTWARDS
    polyhedron_heal = Polyhedron(polyhedral_source, DENSITY, normal_orientation, PolyhedronIntegrity.HEAL)
    healed_faces = np.array(polyhedron_heal.faces)
    np.testing.assert_array_almost_equal(healed_faces, expected_healed)

    # If disabled, there should not be any change - the only way to create an invalid polyhedron
    polyhedron_non_modified = Polyhedron(polyhedral_source, DENSITY, normal_orientation, PolyhedronIntegrity.DISABLE)
    non_modified_faces = np.array(polyhedron_non_modified.faces)
    np.testing.assert_array_almost_equal(non_modified_faces, polyhedral_source[1])

