"""Tests for the compute backend selector of the polyhedral gravity model.

The OpenCL backend is a preference rather than a demand: where it is unavailable -- a build without
OpenCL support, a machine without an OpenCL device, or (for FLOAT64) a device without cl_khr_fp64,
which includes every Apple Silicon GPU -- the evaluation falls back to the CPU. The tests therefore
skip themselves rather than fail when the backend under test did not materialise.
"""
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
    [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
])

CUBE_FACES = np.array([
    [1, 3, 2], [0, 3, 1], [0, 1, 5], [0, 5, 4], [0, 7, 3], [0, 4, 7],
    [1, 2, 6], [1, 6, 5], [2, 3, 6], [3, 7, 6], [4, 5, 6], [4, 6, 7],
])

# The ordinary case plus every singularity case of Tsoulis' algorithm, i.e. P' inside a face,
# on an edge, and on a vertex
COMPUTATION_POINTS = np.array([
    [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [5.0, 5.0, 5.0], [-3.0, 1.5, 0.25],
    [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [-1.0, -1.0, -1.0],
])


@pytest.fixture
def cube() -> Polyhedron:
    """A unitless cube, so that the compared magnitudes are not scaled by the gravitational constant"""
    return Polyhedron(
        polyhedral_source=(CUBE_VERTICES, CUBE_FACES),
        density=1.0,
        normal_orientation=NormalOrientation.OUTWARDS,
        integrity_check=PolyhedronIntegrity.DISABLE,
        metric_unit=MetricUnit.UNITLESS,
    )


def _opencl_evaluable(polyhedron: Polyhedron, precision: ComputePrecision) -> GravityEvaluable:
    """Returns an evaluable running on OpenCL, or skips the test if that is not possible here"""
    evaluable = GravityEvaluable(polyhedron=polyhedron, backend=ComputeBackend.OPENCL, precision=precision)
    if evaluable.compute_backend != ComputeBackend.OPENCL:
        pytest.skip(f"No OpenCL device supporting {precision} is available")
    return evaluable


def _assert_results_close(actual, expected, rtol: float) -> None:
    """Compares a (potential, acceleration, tensor) triplet, with an absolute floor around zero"""
    for actual_part, expected_part in zip(actual, expected):
        np.testing.assert_allclose(actual_part, expected_part, rtol=rtol, atol=rtol)


def test_default_backend_is_opencl_or_falls_back(cube):
    """The default backend must never raise, whatever the machine offers"""
    # pybind11 enum values are not singletons, hence the comparison by value rather than by identity
    evaluable = GravityEvaluable(polyhedron=cube)
    assert evaluable.compute_backend in (ComputeBackend.CPU, ComputeBackend.OPENCL)
    assert evaluable.compute_precision == ComputePrecision.FLOAT64


def test_cpu_backend_is_reported(cube):
    evaluable = GravityEvaluable(polyhedron=cube, backend=ComputeBackend.CPU)
    assert evaluable.compute_backend == ComputeBackend.CPU


@pytest.mark.parametrize("precision, rtol", [
    (ComputePrecision.FLOAT64, 1e-10),
    (ComputePrecision.FLOAT32, 1e-4),
])
def test_opencl_matches_cpu_backend(cube, precision, rtol):
    """Both OpenCL precisions have to reproduce the CPU backend, which is pinned against Tsoulis"""
    opencl_evaluable = _opencl_evaluable(cube, precision)
    cpu_evaluable = GravityEvaluable(polyhedron=cube, backend=ComputeBackend.CPU)

    for point in COMPUTATION_POINTS:
        _assert_results_close(opencl_evaluable(point), cpu_evaluable(point), rtol)


def test_opencl_multi_point_matches_single_point(cube):
    evaluable = _opencl_evaluable(cube, ComputePrecision.FLOAT32)

    batched = evaluable(COMPUTATION_POINTS)
    assert len(batched) == len(COMPUTATION_POINTS)
    for result, point in zip(batched, COMPUTATION_POINTS):
        _assert_results_close(result, evaluable(point), 0.0)


def test_free_evaluate_accepts_the_backend(cube):
    """The free evaluate function has to offer the same selector as the GravityEvaluable"""
    if GravityEvaluable(
        polyhedron=cube, backend=ComputeBackend.OPENCL, precision=ComputePrecision.FLOAT32
    ).compute_backend != ComputeBackend.OPENCL:
        pytest.skip("No OpenCL device is available")

    point = COMPUTATION_POINTS[3]
    actual = evaluate(
        polyhedron=cube,
        computation_points=point,
        backend=ComputeBackend.OPENCL,
        precision=ComputePrecision.FLOAT32,
    )
    expected = evaluate(polyhedron=cube, computation_points=point, backend=ComputeBackend.CPU)
    _assert_results_close(actual, expected, 1e-4)


def test_pickle_round_trip_preserves_the_backend(cube):
    """Pickling has to carry the backend and precision so that the restored evaluable behaves alike"""
    for evaluable in (
        GravityEvaluable(polyhedron=cube, backend=ComputeBackend.CPU),
        GravityEvaluable(polyhedron=cube),
    ):
        restored = pickle.loads(pickle.dumps(evaluable))
        assert restored.compute_backend == evaluable.compute_backend
        assert restored.compute_precision == evaluable.compute_precision
        _assert_results_close(restored(COMPUTATION_POINTS[3]), evaluable(COMPUTATION_POINTS[3]), 1e-14)
