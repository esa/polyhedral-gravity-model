import gc

import numpy as np
import pytest

import polyhedral_gravity
from polyhedral_gravity import MetricUnit, Polyhedron, PolyhedronIntegrity, evaluate

CUBE_VERTICES = np.array([
    [-1, -1, -1],
    [1, -1, -1],
    [1, 1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
    [1, -1, 1],
    [1, 1, 1],
    [-1, 1, 1]
], dtype=np.float64)

CUBE_FACES = np.array([
    [1, 3, 2], [0, 3, 1], [0, 1, 5], [0, 5, 4], [0, 7, 3], [0, 4, 7],
    [1, 2, 6], [1, 6, 5], [2, 3, 6], [3, 7, 6], [4, 5, 6], [4, 6, 7]
], dtype=np.uint64)

COMPUTATION_POINT = [0.5, 0.25, 0.75]


def _make_cube(vertices, faces, integrity_check=PolyhedronIntegrity.DISABLE):
    """Creates the unitless reference cube from the given vertices and faces."""
    return Polyhedron(
        polyhedral_source=(vertices, faces),
        density=1.0,
        integrity_check=integrity_check,
        metric_unit=MetricUnit.UNITLESS,
    )


def _address(array):
    """Returns the address of an array's first element."""
    return array.__array_interface__["data"][0]


def test_matching_dtypes_are_not_copied():
    """A C-contiguous float64/ uint64 mesh becomes the polyhedron's memory as it is."""
    polyhedron = _make_cube(CUBE_VERTICES, CUBE_FACES)
    assert _address(polyhedron.vertices) == _address(CUBE_VERTICES)
    assert _address(polyhedron.faces) == _address(CUBE_FACES)


def test_mesh_properties_are_read_only_views():
    """The mesh is handed back as a view of the polyhedron's own memory, which must not be written."""
    polyhedron = _make_cube(CUBE_VERTICES, CUBE_FACES)
    assert isinstance(polyhedron.vertices, np.ndarray)
    assert polyhedron.vertices.dtype == np.float64
    assert polyhedron.vertices.shape == (8, 3)
    assert not polyhedron.vertices.flags.writeable
    assert polyhedron.faces.dtype == np.uint64
    assert polyhedron.faces.shape == (12, 3)
    assert not polyhedron.faces.flags.writeable
    with pytest.raises(ValueError):
        polyhedron.vertices[0][0] = 42.0


@pytest.mark.parametrize("vertex_dtype", [np.float32, np.float64])
@pytest.mark.parametrize("face_dtype", [np.int32, np.int64, np.uint32, np.uint64])
def test_every_dtype_gives_the_same_result(vertex_dtype, face_dtype):
    """Whatever dtype the arrays have, they are converted into the mesh's own representation."""
    expected = evaluate(_make_cube(CUBE_VERTICES, CUBE_FACES), COMPUTATION_POINT)
    actual = evaluate(
        _make_cube(CUBE_VERTICES.astype(vertex_dtype), CUBE_FACES.astype(face_dtype)),
        COMPUTATION_POINT,
    )
    assert np.allclose(actual[0], expected[0])
    assert np.allclose(actual[1], expected[1])
    assert np.allclose(actual[2], expected[2])


@pytest.mark.parametrize("vertices", [
    CUBE_VERTICES.tolist(),
    np.asfortranarray(CUBE_VERTICES),
    np.repeat(CUBE_VERTICES, 2, axis=0)[::2],
])
def test_non_contiguous_input_is_converted(vertices):
    """Anything which is not a C-contiguous array of the right dtype costs one conversion."""
    expected = evaluate(_make_cube(CUBE_VERTICES, CUBE_FACES), COMPUTATION_POINT)
    actual = evaluate(_make_cube(vertices, CUBE_FACES), COMPUTATION_POINT)
    assert np.allclose(actual[0], expected[0])
    assert np.allclose(actual[1], expected[1])


def test_polyhedron_outlives_its_source_arrays():
    """The polyhedron keeps a reference to the arrays whose memory it uses."""
    vertices = CUBE_VERTICES.copy()
    faces = CUBE_FACES.copy()
    polyhedron = _make_cube(vertices, faces)
    expected = evaluate(polyhedron, COMPUTATION_POINT)

    del vertices, faces
    gc.collect()

    assert np.array_equal(polyhedron.vertices, CUBE_VERTICES)
    assert np.allclose(evaluate(polyhedron, COMPUTATION_POINT)[0], expected[0])


def test_healing_does_not_modify_the_given_faces():
    """Healing the vertex ordering must never write into the caller's array."""
    faces = CUBE_FACES.copy()
    faces[0][0], faces[0][1] = faces[0][1], faces[0][0]
    unhealed = faces.copy()

    healed = _make_cube(CUBE_VERTICES, faces, integrity_check=PolyhedronIntegrity.HEAL)

    assert np.array_equal(faces, unhealed)
    assert np.array_equal(healed.faces, CUBE_FACES)


def test_one_based_indexing_is_shifted():
    """A mesh which never references vertex zero is taken to be indexed from one."""
    expected = evaluate(_make_cube(CUBE_VERTICES, CUBE_FACES), COMPUTATION_POINT)
    actual = evaluate(_make_cube(CUBE_VERTICES, CUBE_FACES + 1), COMPUTATION_POINT)
    assert np.allclose(actual[0], expected[0])


@pytest.mark.parametrize("source", [
    "single-string",
    (np.zeros((4, 2)), CUBE_FACES),
    (CUBE_VERTICES, np.zeros((4, 4), dtype=np.uint64)),
    (CUBE_VERTICES,),
    (CUBE_VERTICES, CUBE_FACES, CUBE_FACES),
])
def test_invalid_sources_are_rejected(source):
    """A polyhedral source which is neither files nor a pair of (N, 3) arrays is an error."""
    with pytest.raises((ValueError, TypeError, RuntimeError)):
        Polyhedron(
            polyhedral_source=source,
            density=1.0,
            integrity_check=PolyhedronIntegrity.DISABLE,
            metric_unit=MetricUnit.UNITLESS,
        )


def test_compiler_attributes_are_exposed():
    """The compilers which translated the host and the device code are embedded at compile time."""
    assert isinstance(polyhedral_gravity.__host_compiler__, str)
    assert polyhedral_gravity.__host_compiler__ != ""
    assert isinstance(polyhedral_gravity.__device_compiler__, str)
    assert polyhedral_gravity.__device_compiler__ != ""
    gpu_available = any(
        space in polyhedral_gravity.__parallelization__ for space in ("Cuda", "HIP", "SYCL")
    )
    assert (polyhedral_gravity.__device_compiler__ != "None") == gpu_available
