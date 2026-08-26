import pickle

import numpy as np
import pytest
from polyhedral_gravity import (
    ComputeBackend,
    ComputePrecision,
    GravityEvaluable,
    MetricUnit,
    NormalOrientation,
    Polyhedron,
    PolyhedronIntegrity,
    evaluate,
)

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
    [1, 3, 2], [0, 3, 1], [0, 1, 5], [0, 5, 4], [0, 7, 3], [0, 4, 7],
    [1, 2, 6], [1, 6, 5], [2, 3, 6], [3, 7, 6], [4, 5, 6], [4, 6, 7],
])

# Covers all singularity cases of Tsoulis' algorithm: inside, on a face, on an edge, on a vertex, outside
COMPUTATION_POINTS = np.array([
    [0.0, 0.0, 0.0],
    [0.5, 0.25, 0.75],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
    [2.0, 0.0, 0.0],
    [-3.0, 2.0, 1.5],
    [10.0, 10.0, 10.0],
])

CPU_BACKENDS = [ComputeBackend.CPU_SERIAL, ComputeBackend.CPU_PARALLEL]


def gpu_available() -> bool:
    """Whether this build of the module has a GPU backend compiled in."""
    try:
        polyhedron = _cube()
        evaluate(
            polyhedron=polyhedron,
            computation_points=[0.0, 0.0, 0.0],
            backend=ComputeBackend.GPU_PARALLEL,
        )
        return True
    except RuntimeError:
        return False


def _cube() -> Polyhedron:
    return Polyhedron(
        polyhedral_source=(CUBE_VERTICES, CUBE_FACES),
        density=1.0,
        normal_orientation=NormalOrientation.OUTWARDS,
        integrity_check=PolyhedronIntegrity.DISABLE,
        metric_unit=MetricUnit.UNITLESS,
    )


def _unpack(solutions):
    """Splits a list of (potential, acceleration, tensor) triplets into three arrays."""
    return (
        np.array([solution[0] for solution in solutions]),
        np.array([solution[1] for solution in solutions]),
        np.array([solution[2] for solution in solutions]),
    )


def test_module_reports_its_execution_spaces():
    """The module tells the user which Kokkos execution spaces it was compiled with."""
    import polyhedral_gravity

    assert isinstance(polyhedral_gravity.__parallelization__, str)
    assert "Serial" in polyhedral_gravity.__parallelization__


@pytest.mark.parametrize("backend", CPU_BACKENDS)
def test_cpu_backends_agree(backend):
    """Every CPU backend runs the identical kernel, so they must agree with each other."""
    evaluable = GravityEvaluable(polyhedron=_cube())
    reference = _unpack(evaluable(computation_points=COMPUTATION_POINTS, backend=ComputeBackend.CPU_SERIAL))
    actual = _unpack(evaluable(computation_points=COMPUTATION_POINTS, backend=backend))
    for actual_values, expected_values in zip(actual, reference):
        np.testing.assert_allclose(actual_values, expected_values, rtol=0.0, atol=1e-13)


@pytest.mark.skipif(not gpu_available(), reason="This build has no GPU backend")
def test_gpu_backend_agrees_with_cpu():
    """The GPU runs the identical kernel and may only differ by the reduction's reassociation."""
    evaluable = GravityEvaluable(polyhedron=_cube())
    reference = _unpack(evaluable(computation_points=COMPUTATION_POINTS, backend=ComputeBackend.CPU_PARALLEL))
    actual = _unpack(evaluable(computation_points=COMPUTATION_POINTS, backend=ComputeBackend.GPU_PARALLEL))
    for actual_values, expected_values in zip(actual, reference):
        np.testing.assert_allclose(actual_values, expected_values, rtol=0.0, atol=1e-10)


@pytest.mark.skipif(gpu_available(), reason="This build has a GPU backend")
def test_gpu_backend_raises_without_gpu():
    """Requesting the GPU without a GPU backend must fail loudly, not silently compute elsewhere."""
    evaluable = GravityEvaluable(polyhedron=_cube())
    with pytest.raises(RuntimeError):
        evaluable(computation_points=[0.0, 0.0, 0.0], backend=ComputeBackend.GPU_PARALLEL)
    with pytest.raises(RuntimeError):
        evaluate(
            polyhedron=_cube(),
            computation_points=COMPUTATION_POINTS,
            backend=ComputeBackend.GPU_PARALLEL,
        )


def test_float32_approximates_float64():
    """Single precision only reproduces a few significant digits, but must stay in the right ballpark."""
    float64 = GravityEvaluable(polyhedron=_cube(), precision=ComputePrecision.FLOAT64)
    float32 = GravityEvaluable(polyhedron=_cube(), precision=ComputePrecision.FLOAT32)
    assert float64.precision == ComputePrecision.FLOAT64
    assert float32.precision == ComputePrecision.FLOAT32

    expected = _unpack(float64(computation_points=COMPUTATION_POINTS))
    actual = _unpack(float32(computation_points=COMPUTATION_POINTS))
    for actual_values, expected_values in zip(actual, expected):
        np.testing.assert_allclose(actual_values, expected_values, rtol=0.0, atol=1e-4)


def test_defaults_are_cpu_parallel_and_float64():
    """The defaults must reproduce what an explicit CPU_PARALLEL/ FLOAT64 evaluation computes."""
    evaluable = GravityEvaluable(polyhedron=_cube())
    assert evaluable.precision == ComputePrecision.FLOAT64
    expected = _unpack(evaluable(computation_points=COMPUTATION_POINTS, backend=ComputeBackend.CPU_PARALLEL))
    actual = _unpack(evaluable(computation_points=COMPUTATION_POINTS))
    for actual_values, expected_values in zip(actual, expected):
        np.testing.assert_array_equal(actual_values, expected_values)


def test_single_point_matches_multi_point():
    """The multi point evaluation uses a team policy and has to agree with the single point one."""
    evaluable = GravityEvaluable(polyhedron=_cube())
    multi = evaluable(computation_points=COMPUTATION_POINTS)
    for index, point in enumerate(COMPUTATION_POINTS):
        single = evaluable(computation_points=point)
        np.testing.assert_allclose(single[0], multi[index][0], rtol=0.0, atol=1e-13)
        np.testing.assert_allclose(single[1], multi[index][1], rtol=0.0, atol=1e-13)
        np.testing.assert_allclose(single[2], multi[index][2], rtol=0.0, atol=1e-13)


def test_free_function_matches_evaluable():
    """The free evaluate function is a thin wrapper and must not change the result."""
    evaluable = GravityEvaluable(polyhedron=_cube())
    expected = _unpack(evaluable(computation_points=COMPUTATION_POINTS, backend=ComputeBackend.CPU_SERIAL))
    actual = _unpack(evaluate(
        polyhedron=_cube(),
        computation_points=COMPUTATION_POINTS,
        backend=ComputeBackend.CPU_SERIAL,
        precision=ComputePrecision.FLOAT64,
    ))
    for actual_values, expected_values in zip(actual, expected):
        np.testing.assert_allclose(actual_values, expected_values, rtol=0.0, atol=1e-15)


def test_pickle_round_trip_keeps_precision(tmp_path):
    """Pickling a GravityEvaluable must preserve the precision and the computed values."""
    original = GravityEvaluable(polyhedron=_cube(), precision=ComputePrecision.FLOAT32)
    pickle_output = tmp_path / "evaluable.pk"
    with open(pickle_output, "wb") as file:
        pickle.dump(original, file, pickle.HIGHEST_PROTOCOL)
    with open(pickle_output, "rb") as file:
        restored = pickle.load(file)

    assert restored.precision == ComputePrecision.FLOAT32
    expected = _unpack(original(computation_points=COMPUTATION_POINTS))
    actual = _unpack(restored(computation_points=COMPUTATION_POINTS))
    for actual_values, expected_values in zip(actual, expected):
        np.testing.assert_array_equal(actual_values, expected_values)


def test_empty_polyhedron_is_rejected():
    """A polyhedron without faces cannot be evaluated and must be rejected, not silently return zeros."""
    empty = Polyhedron(
        polyhedral_source=(np.empty((0, 3)), np.empty((0, 3), dtype=int)),
        density=1.0,
        normal_orientation=NormalOrientation.OUTWARDS,
        integrity_check=PolyhedronIntegrity.DISABLE,
        metric_unit=MetricUnit.UNITLESS,
    )
    with pytest.raises(ValueError):
        GravityEvaluable(polyhedron=empty)
