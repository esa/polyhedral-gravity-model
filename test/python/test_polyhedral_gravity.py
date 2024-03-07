from typing import Tuple, List, Union
import polyhedral_gravity
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
    [4, 6, 7]
])

DENSITY = 1.0

CUBE_VERTICES_FILE = Path("test/resources/cube.node")
CUBE_FACE_FILE = Path("test/resources/cube.face")


@pytest.fixture(scope="session")
def reference_solution() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reads the cube reference solution from a file and returns it as a tuple"""
    reference_file_path = Path("test/resources/analytic_cube_solution.txt")
    reference = np.loadtxt(reference_file_path)
    points = reference[:, :3]
    expected_potential = reference[:, 3].flatten()
    expected_acceleration = reference[:, 4:]
    return points, expected_potential, expected_acceleration


@pytest.mark.parametrize(
    "polyhedral_source", [(CUBE_VERTICES, CUBE_FACES), [str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)]],
    ids=["with_arrays", "with_file_input"]
)
def test_polyhedral_gravity(
        polyhedral_source: Union[Tuple[np.ndarray, np.ndarray], List[str]],
        reference_solution: Tuple[np.ndarray, np.ndarray, np.ndarray]) -> None:
    """Checks that the evaluate function the correct results and
    is callable with file/ array inputs.
    """
    points, expected_potential, expected_acceleration = reference_solution
    sol = polyhedral_gravity.evaluate(
        polyhedral_source=polyhedral_source,
        density=DENSITY,
        computation_points=points,
        parallel=True
    )
    potential = np.array([result[0] for result in sol])
    acceleration = np.array([result[1] for result in sol])
    np.testing.assert_array_almost_equal(potential, expected_potential)
    np.testing.assert_array_almost_equal(acceleration, expected_acceleration)


@pytest.mark.parametrize(
    "polyhedral_source", [(CUBE_VERTICES, CUBE_FACES), [str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)]],
    ids=["with_arrays", "with_file_input"]
)
def test_polyhedral_gravity_evaluable(
        polyhedral_source: Union[Tuple[np.ndarray, np.ndarray], List[str]],
        reference_solution: Tuple[np.ndarray, np.ndarray, np.ndarray]) -> None:
    """Checks that the evaluable produces the correct results and
    is instantiable with file/ array inputs.
    """
    points, expected_potential, expected_acceleration = reference_solution
    evaluable = polyhedral_gravity.GravityEvaluable(
        polyhedral_source=polyhedral_source,
        density=DENSITY,
    )
    sol = evaluable(
        computation_points=points,
        parallel=True
    )
    potential = np.array([result[0] for result in sol])
    acceleration = np.array([result[1] for result in sol])
    np.testing.assert_array_almost_equal(potential, expected_potential)
    np.testing.assert_array_almost_equal(acceleration, expected_acceleration)


@pytest.mark.parametrize(
    "polyhedral_source", [(CUBE_VERTICES, CUBE_FACES), [str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)]],
    ids=["with_arrays", "with_file_input"]
)
def test_polyhedral_evaluable_pickle(
        polyhedral_source: Union[Tuple[np.ndarray, np.ndarray], List[str]],
        reference_solution: Tuple[np.ndarray, np.ndarray, np.ndarray],
        tmp_path: Path,) -> None:
    """Tests that the evaluable can be pickled and unpicked and that the results
    are still correct afterward (i.e. that the internal cache is correctly
    pickled and unpicked as well).
    """
    points, expected_potential, expected_acceleration = reference_solution
    initial_evaluable = polyhedral_gravity.GravityEvaluable(
        polyhedral_source=polyhedral_source,
        density=DENSITY
    )
    pickle_output = tmp_path.joinpath("evaluable.pk")
    with open(pickle_output, "wb") as f:
        pickle.dump(initial_evaluable, f, pickle.HIGHEST_PROTOCOL)

    with open(pickle_output, "rb") as f:
        read_evaluable = pickle.load(f)

    sol = read_evaluable(
        computation_points=points,
        parallel=True
    )
    potential = np.array([result[0] for result in sol])
    acceleration = np.array([result[1] for result in sol])
    np.testing.assert_array_almost_equal(potential, expected_potential)
    np.testing.assert_array_almost_equal(acceleration, expected_acceleration)
