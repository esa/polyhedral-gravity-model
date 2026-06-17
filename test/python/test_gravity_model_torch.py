"""
Validate polyhedral_gravity.torch against the C++ polyhedral_gravity reference.

All test classes run against every body defined below:
  - Unit cube            (vertices at ±1, 12 triangular faces)
  - Regular tetrahedron  (inscribed in the ±1 cube, edge length 2√2)
  - Right-angle tetrahedron (legs of length 2 along coordinate axes)
  - Generic tetrahedron  (random vertices, no special symmetries)

The C++ evaluate() is the ground truth; we require relative errors below 1e-10.
Tests are skipped automatically when PyTorch is not installed.
"""

import math

import numpy as np
import pytest

torch = pytest.importorskip("torch", reason="PyTorch not installed")

from polyhedral_gravity import Polyhedron, NormalOrientation, evaluate
from polyhedral_gravity.torch import evaluate as torch_evaluate

REL_TOL = 1e-10

# ---------------------------------------------------------------------------
# Body definitions
# ---------------------------------------------------------------------------

_CUBE = dict(
    id="cube",
    vertices=np.array([
        [-1, -1, -1], [1, -1, -1], [1,  1, -1], [-1,  1, -1],
        [-1, -1,  1], [1, -1,  1], [1,  1,  1], [-1,  1,  1],
    ], dtype=np.float64),
    faces=np.array([
        [1, 3, 2], [0, 3, 1], [0, 1, 5], [0, 5, 4],
        [0, 7, 3], [0, 4, 7], [1, 2, 6], [1, 6, 5],
        [2, 3, 6], [3, 7, 6], [4, 5, 6], [4, 6, 7],
    ], dtype=np.int64),
    exterior_points=np.array([
        [-5, -5, -5], [3, 0, 0], [0, 4, -2], [2, 2, 2],
    ], dtype=np.float64),
    grad_point=np.array([[0., 4., -2.]], dtype=np.float64),
    edge_singular_points=np.array([
        [0.5, 0.5, -5.0], [-0.5, -0.5, 5.0], [0.0, 0.0, -3.0],
    ], dtype=np.float64),
    vertex_singular_points=np.array([
        [1.0, 1.0, -5.0], [-1.0, 1.0, 4.0], [1.0, -1.0, 3.5],
    ], dtype=np.float64),
)

_REGULAR_TET = dict(
    id="regular_tet",
    vertices=np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=np.float64),
    faces=np.array([[0, 1, 2], [1, 3, 2], [0, 2, 3], [0, 3, 1]], dtype=np.int64),
    exterior_points=np.array([
        [5, 0, 0], [-3, -3, 0], [2, 2, -2], [0, -4, 1], [-2, 3, -2],
    ], dtype=np.float64),
    grad_point=np.array([[5., 0., 0.]], dtype=np.float64),
    edge_singular_points=np.array([[2.0, 2.0, -3.0]], dtype=np.float64),
    vertex_singular_points=np.array([[0.0, 2.0, -2.0]], dtype=np.float64),
)

_RIGHT_TET = dict(
    id="right_tet",
    vertices=np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, 0, 2]], dtype=np.float64),
    faces=np.array([[0, 2, 1], [0, 3, 2], [0, 1, 3], [1, 2, 3]], dtype=np.int64),
    exterior_points=np.array([
        [-2, 1, 1], [1, -2, 1], [1, 1, -2], [3, 3, 3], [-1, -1, -1],
    ], dtype=np.float64),
    grad_point=np.array([[3., 3., 3.]], dtype=np.float64),
    edge_singular_points=np.array([[1.0, 1.0, -2.0]], dtype=np.float64),
    vertex_singular_points=np.array([[0.0, 2.0, -2.0]], dtype=np.float64),
)

_GENERIC_TET = dict(
    id="generic_tet",
    vertices=np.array([
        [-1.76691274, -1.08331734, -2.65127211],
        [-2.22567756, -0.14244291, -3.09891698],
        [-1.99141967, -0.38573326, -3.18724320],
        [-2.93988158, -1.39857477, -4.18338299],
    ], dtype=np.float64),
    faces=np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], dtype=np.int64),
    exterior_points=np.array([
        [-3.19803032, -1.76076925, -3.49138362],
        [-3.91751185, -0.54205321, -2.81119925],
        [-3.18486445, -2.05842822, -3.22881105],
        [-2.95495711, -3.35812511, -3.98328743],
    ], dtype=np.float64),
    grad_point=np.array([[-3.76421772, -0.30497729, -1.18500307]], dtype=np.float64),
    edge_singular_points=None,
    vertex_singular_points=None,
)

ALL_BODIES = [_CUBE, _REGULAR_TET, _RIGHT_TET, _GENERIC_TET]
BODIES_WITH_SINGULARITIES = [b for b in ALL_BODIES if b["edge_singular_points"] is not None]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rel_err_pot(mine, ref):
    return np.abs(mine - ref) / (np.abs(ref) + 1e-30)


def _rel_err_acc(mine, ref):
    return np.linalg.norm(mine - ref, axis=-1) / (np.linalg.norm(ref, axis=-1) + 1e-30)


def _esa_batch(poly, pts):
    Vs, gs = [], []
    for p in pts:
        V, g, _ = evaluate(poly, p)
        Vs.append(V)
        gs.append(g)
    return np.array(Vs), np.array(gs)


# ---------------------------------------------------------------------------
# Shared fixture — parametrised over all bodies
# ---------------------------------------------------------------------------

def _make_body(b):
    esa_poly = Polyhedron(
        (b["vertices"], b["faces"]),
        density=1.0,
        normal_orientation=NormalOrientation.OUTWARDS,
    )
    return {
        **b,
        "verts_t": torch.tensor(b["vertices"], dtype=torch.float64),
        "faces_t": torch.tensor(b["faces"], dtype=torch.long),
        "esa_poly": esa_poly,
    }


@pytest.fixture(
    params=[pytest.param(b, id=b["id"]) for b in ALL_BODIES],
    scope="module",
)
def body(request):
    return _make_body(request.param)


@pytest.fixture(
    params=[pytest.param(b, id=b["id"]) for b in BODIES_WITH_SINGULARITIES],
    scope="module",
)
def body_with_singularities(request):
    return _make_body(request.param)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestAccuracy:
    """Point-by-point accuracy for generic exterior query points."""

    def test_exterior_points(self, body):
        pts = body["exterior_points"]
        V, g, _ = torch_evaluate(
            body["verts_t"], body["faces_t"], 1.0,
            torch.tensor(pts, dtype=torch.float64),
        )
        V_np = V.detach().numpy()
        g_np = g.detach().numpy()

        for i, p in enumerate(pts):
            esa_V, esa_g, _ = evaluate(body["esa_poly"], p)
            assert _rel_err_pot(V_np[i], esa_V) < REL_TOL, \
                f"body={body['id']} point={p}: potential rel error too large"
            assert _rel_err_acc(g_np[i], np.array(esa_g)) < REL_TOL, \
                f"body={body['id']} point={p}: acceleration rel error too large"


class TestSingularityCases:
    """Points whose projection onto a face lands on an edge or vertex."""

    def test_edge_singularity(self, body_with_singularities):
        body = body_with_singularities
        V, g, _ = torch_evaluate(
            body["verts_t"], body["faces_t"], 1.0,
            torch.tensor(body["edge_singular_points"], dtype=torch.float64),
        )
        for i, p in enumerate(body["edge_singular_points"]):
            esa_V, esa_g, _ = evaluate(body["esa_poly"], p)
            assert _rel_err_pot(V.detach().numpy()[i], esa_V) < REL_TOL
            assert _rel_err_acc(g.detach().numpy()[i], np.array(esa_g)) < REL_TOL

    def test_vertex_singularity(self, body_with_singularities):
        body = body_with_singularities
        V, g, _ = torch_evaluate(
            body["verts_t"], body["faces_t"], 1.0,
            torch.tensor(body["vertex_singular_points"], dtype=torch.float64),
        )
        for i, p in enumerate(body["vertex_singular_points"]):
            esa_V, esa_g, _ = evaluate(body["esa_poly"], p)
            assert _rel_err_pot(V.detach().numpy()[i], esa_V) < REL_TOL
            assert _rel_err_acc(g.detach().numpy()[i], np.array(esa_g)) < REL_TOL


class TestBatched:
    """Batch evaluation on 50 points on a sphere of radius 10."""

    @pytest.fixture(scope="class")
    def sphere_points(self):
        theta = np.linspace(0, 2 * math.pi, 10, endpoint=False)
        phi = np.linspace(math.pi / 6, 5 * math.pi / 6, 5)
        return np.array([
            [10 * math.sin(p) * math.cos(t),
             10 * math.sin(p) * math.sin(t),
             10 * math.cos(p)]
            for p in phi for t in theta
        ], dtype=np.float64)

    def test_potential_batch(self, body, sphere_points):
        V, _, __ = torch_evaluate(
            body["verts_t"], body["faces_t"], 1.0,
            torch.tensor(sphere_points, dtype=torch.float64),
        )
        ref_V, _ = _esa_batch(body["esa_poly"], sphere_points)
        err = _rel_err_pot(V.detach().numpy(), ref_V)
        assert err.max() < REL_TOL, \
            f"body={body['id']}: max potential rel error = {err.max():.2e}"

    def test_acceleration_batch(self, body, sphere_points):
        _, g, __ = torch_evaluate(
            body["verts_t"], body["faces_t"], 1.0,
            torch.tensor(sphere_points, dtype=torch.float64),
        )
        _, ref_g = _esa_batch(body["esa_poly"], sphere_points)
        err = _rel_err_acc(g.detach().numpy(), ref_g)
        assert err.max() < REL_TOL, \
            f"body={body['id']}: max acceleration rel error = {err.max():.2e}"


class TestGradients:
    """Autodiff correctness."""

    def test_gradient_wrt_vertices_nonzero(self, body):
        verts = body["verts_t"].clone().requires_grad_(True)
        q = torch.tensor(body["exterior_points"][:1], dtype=torch.float64)
        V, _, __ = torch_evaluate(verts, body["faces_t"], 1.0, q)
        V.sum().backward()
        assert verts.grad is not None and verts.grad.abs().max() > 0

    def test_gradient_wrt_scalar_density_nonzero(self, body):
        density = torch.tensor(1.0, dtype=torch.float64, requires_grad=True)
        q = torch.tensor(body["exterior_points"][:1], dtype=torch.float64)
        V, _, __ = torch_evaluate(body["verts_t"], body["faces_t"], density, q)
        V.sum().backward()
        assert density.grad is not None and density.grad.abs() > 0

    def test_gradient_consistent_with_finite_differences(self, body):
        """∂V/∂vertices must agree with central finite differences to 1e-6."""
        verts = body["verts_t"].clone().requires_grad_(True)
        q = torch.tensor(body["grad_point"], dtype=torch.float64)
        V, _, __ = torch_evaluate(verts, body["faces_t"], 1.0, q)
        V.backward()
        grad = verts.grad.clone()

        eps = 1e-5
        fd_grad = torch.zeros_like(grad)
        base = body["verts_t"]
        for i in range(base.shape[0]):
            for j in range(3):
                vp = base.clone(); vp[i, j] += eps
                vm = base.clone(); vm[i, j] -= eps
                Vp, _, __ = torch_evaluate(vp, body["faces_t"], 1.0, q)
                Vm, _, __ = torch_evaluate(vm, body["faces_t"], 1.0, q)
                fd_grad[i, j] = (Vp - Vm) / (2 * eps)

        rel = (grad - fd_grad).abs() / (fd_grad.abs() + 1e-30)
        assert rel.max().item() < 1e-6, \
            f"body={body['id']}: max FD rel error = {rel.max():.2e}"


class TestTensor:
    """Gradient tensor accuracy against the C++ reference."""

    def test_tensor_accuracy(self, body):
        pts = body["exterior_points"]
        _, _, T = torch_evaluate(
            body["verts_t"], body["faces_t"], 1.0,
            torch.tensor(pts, dtype=torch.float64),
        )
        T_np = T.detach().numpy()

        for i, p in enumerate(pts):
            _, _, esa_T = evaluate(body["esa_poly"], p)
            ref = np.array(esa_T)
            # Normalise by the problem scale (max component) so near-zero diagonal
            # components at high-symmetry points don't blow up the relative error.
            scale = np.abs(ref).max() + 1e-30
            err = np.abs(T_np[i] - ref) / scale
            assert err.max() < REL_TOL, \
                f"body={body['id']} point={p}: tensor rel error too large: {err.max():.2e}"

    def test_tensor_scaling(self, body):
        q = torch.tensor(body["exterior_points"][:2], dtype=torch.float64)
        _, _, T1 = torch_evaluate(body["verts_t"], body["faces_t"], 1.0, q)
        _, _, Tr = torch_evaluate(body["verts_t"], body["faces_t"], 2.0, q)
        assert torch.allclose(Tr, 2.0 * T1, rtol=1e-12)


class TestDensityScaling:
    """Potential and acceleration scale linearly with density."""

    @pytest.mark.parametrize("rho", [0.5, 2.0, 1000.0])
    def test_linear_scaling(self, body, rho):
        q = torch.tensor(body["exterior_points"][:2], dtype=torch.float64)
        V1, g1, _ = torch_evaluate(body["verts_t"], body["faces_t"], 1.0, q)
        Vr, gr, _ = torch_evaluate(body["verts_t"], body["faces_t"], rho, q)
        assert torch.allclose(Vr, rho * V1, rtol=1e-12)
        assert torch.allclose(gr, rho * g1, rtol=1e-12)

