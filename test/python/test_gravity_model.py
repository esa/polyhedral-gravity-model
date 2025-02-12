from typing import Tuple, List, Union
from polyhedral_gravity import Polyhedron, GravityEvaluable, evaluate, PolyhedronIntegrity, NormalOrientation, MetricUnit
import numpy as np
import pickle
import pytest
from pathlib import Path
from functools import lru_cache

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

CUBE_FACES_INVERTED = CUBE_FACES[:, [1, 0, 2]]

DENSITY = 1.0

CUBE_VERTICES_FILE = Path("test/resources/cube.node")
CUBE_FACE_FILE = Path("test/resources/cube.face")


@lru_cache(maxsize=10)
def reference_solution(density: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reads the cube reference solution from a file and returns it as a tuple"""
    reference_file_path = Path(f"test/resources/analytic_cube_solution_density{int(density)}.txt")
    reference = np.loadtxt(reference_file_path, skiprows=1)
    points = reference[:, :3]
    expected_potential = reference[:, 3].flatten()
    expected_acceleration = reference[:, 4:]
    return points, expected_potential, expected_acceleration


@pytest.mark.parametrize(
    "polyhedral_source,normal_orientation,density", [
        ((CUBE_VERTICES, CUBE_FACES), NormalOrientation.OUTWARDS, 1.0),
        ((CUBE_VERTICES, CUBE_FACES_INVERTED), NormalOrientation.INWARDS, 1.0),
        ([str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)], NormalOrientation.OUTWARDS, 1.0),
        ((CUBE_VERTICES, CUBE_FACES), NormalOrientation.OUTWARDS, 42.0),
        ((CUBE_VERTICES, CUBE_FACES_INVERTED), NormalOrientation.INWARDS, 42.0),
        ([str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)], NormalOrientation.OUTWARDS, 42.0)
    ],
    ids=["with_arrays_outwards01", "with_arrays_inwards01", "with_file_input01",
         "with_arrays_outwards42", "with_arrays_inwards42", "with_file_input42"]
)
def test_polyhedral_gravity(
        polyhedral_source: Union[Tuple[np.ndarray, np.ndarray], List[str]],
        normal_orientation: NormalOrientation,
        density: float
) -> None:
    """Checks that the evaluate function the correct results and
    is callable with file/ array inputs.
    """
    points, expected_potential, expected_acceleration = reference_solution(density)
    polyhedron = Polyhedron(
        polyhedral_source=polyhedral_source,
        density=DENSITY,
        normal_orientation=normal_orientation,
        integrity_check=PolyhedronIntegrity.VERIFY,
    )
    sol = evaluate(
        polyhedron=polyhedron,
        computation_points=points,
        parallel=True,
    )
    potential = np.array([result[0] for result in sol])
    acceleration = np.array([result[1] for result in sol])
    np.testing.assert_array_almost_equal(potential, expected_potential)
    np.testing.assert_array_almost_equal(acceleration, expected_acceleration)


@pytest.mark.parametrize(
    "polyhedral_source,normal_orientation,density", [
        ((CUBE_VERTICES, CUBE_FACES), NormalOrientation.OUTWARDS, 1.0),
        ((CUBE_VERTICES, CUBE_FACES_INVERTED), NormalOrientation.INWARDS, 1.0),
        ([str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)], NormalOrientation.OUTWARDS, 1.0),
        ((CUBE_VERTICES, CUBE_FACES), NormalOrientation.OUTWARDS, 42.0),
        ((CUBE_VERTICES, CUBE_FACES_INVERTED), NormalOrientation.INWARDS, 42.0),
        ([str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)], NormalOrientation.OUTWARDS, 42.0)
    ],
    ids=["with_arrays_outwards01", "with_arrays_inwards01", "with_file_input01",
         "with_arrays_outwards42", "with_arrays_inwards42", "with_file_input42"]
)
def test_polyhedral_gravity_evaluable(
        polyhedral_source: Union[Tuple[np.ndarray, np.ndarray], List[str]],
        normal_orientation: NormalOrientation,
        density: float
) -> None:
    """Checks that the evaluable produces the correct results and
    is instantiable with file/ array inputs.
    """
    points, expected_potential, expected_acceleration = reference_solution(density)
    polyhedron = Polyhedron(
        polyhedral_source=polyhedral_source,
        density=DENSITY,
        normal_orientation=normal_orientation,
        integrity_check=PolyhedronIntegrity.VERIFY,
    )
    evaluable = GravityEvaluable(polyhedron=polyhedron)
    sol = evaluable(
        computation_points=points,
        parallel=True,
    )
    potential = np.array([result[0] for result in sol])
    acceleration = np.array([result[1] for result in sol])
    np.testing.assert_array_almost_equal(potential, expected_potential)
    np.testing.assert_array_almost_equal(acceleration, expected_acceleration)


@pytest.mark.parametrize(
    "polyhedral_source,normal_orientation,density", [
        ((CUBE_VERTICES, CUBE_FACES), NormalOrientation.OUTWARDS, 1.0),
        ((CUBE_VERTICES, CUBE_FACES_INVERTED), NormalOrientation.INWARDS, 1.0),
        ([str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)], NormalOrientation.OUTWARDS, 1.0),
        ((CUBE_VERTICES, CUBE_FACES), NormalOrientation.OUTWARDS, 42.0),
        ((CUBE_VERTICES, CUBE_FACES_INVERTED), NormalOrientation.INWARDS, 42.0),
        ([str(CUBE_VERTICES_FILE), str(CUBE_FACE_FILE)], NormalOrientation.OUTWARDS, 42.0)
    ],
    ids=["with_arrays_outwards01", "with_arrays_inwards01", "with_file_input01",
         "with_arrays_outwards42", "with_arrays_inwards42", "with_file_input42"]
)
def test_polyhedral_evaluable_pickle(
        polyhedral_source: Union[Tuple[np.ndarray, np.ndarray], List[str]],
        normal_orientation: NormalOrientation,
        density: float,
        tmp_path: Path,) -> None:
    """Tests that the evaluable can be pickled and unpicked and that the results
    are still correct afterward (i.e. that the internal cache is correctly
    pickled and unpicked as well).
    """
    points, expected_potential, expected_acceleration = reference_solution(density)
    polyhedron = Polyhedron(
        polyhedral_source=polyhedral_source,
        density=DENSITY,
        normal_orientation=normal_orientation,
        integrity_check=PolyhedronIntegrity.DISABLE,
    )
    initial_evaluable = GravityEvaluable(polyhedron=polyhedron)
    pickle_output = tmp_path.joinpath("evaluable.pk")
    with open(pickle_output, "wb") as f:
        pickle.dump(initial_evaluable, f, pickle.HIGHEST_PROTOCOL)

    with open(pickle_output, "rb") as f:
        read_evaluable = pickle.load(f)

    sol = read_evaluable(
        computation_points=points,
        parallel=True,
    )
    potential = np.array([result[0] for result in sol])
    acceleration = np.array([result[1] for result in sol])
    np.testing.assert_array_almost_equal(potential, expected_potential)
    np.testing.assert_array_almost_equal(acceleration, expected_acceleration)


def test_polyhedron_metric() -> None:
    """Tests the metric conversion/ options
    """
    computation_point = np.array([1.0, 0.0, 1.0])
    gravity_constant = 6.67430e-11
    meter_polyhedron = Polyhedron(
        polyhedral_source=(CUBE_VERTICES, CUBE_FACES),
        density=1.0,
        integrity_check=PolyhedronIntegrity.DISABLE,
        metric_unit=MetricUnit.METER,
    )
    meter_evaluable = GravityEvaluable(polyhedron=meter_polyhedron)
    kilometer_polyhedron = Polyhedron(
        polyhedral_source=(CUBE_VERTICES, CUBE_FACES),
        density=1.0,
        integrity_check=PolyhedronIntegrity.DISABLE,
        metric_unit=MetricUnit.KILOMETER,
    )
    kilometer_evaluable = GravityEvaluable(polyhedron=kilometer_polyhedron)
    unitless_polyhedron = Polyhedron(
        polyhedral_source=(CUBE_VERTICES, CUBE_FACES),
        density=1.0,
        integrity_check=PolyhedronIntegrity.DISABLE,
        metric_unit=MetricUnit.UNITLESS,
    )
    unitless_evaluable = GravityEvaluable(polyhedron=unitless_polyhedron)

    assert evaluate(meter_polyhedron, computation_point)[0] * 1e-9 == evaluate(kilometer_polyhedron, computation_point)[0]
    assert evaluate(meter_polyhedron, computation_point)[0]  == evaluate(unitless_polyhedron, computation_point)[0] * gravity_constant
    assert evaluate(kilometer_polyhedron, computation_point)[0] == evaluate(unitless_polyhedron, computation_point)[0] * 1e-9 * gravity_constant

    assert meter_evaluable.output_units == ["m^2/s^2", "m/s^2", "1/s^2"]
    assert kilometer_evaluable.output_units == ["km^2/s^2", "km/s^2", "1/s^2"]
    assert unitless_evaluable.output_units == ["1/s^2", "1/s^2", "1/s^2"]