"""
Differentiable polyhedral gravity model - PyTorch interface.

Implements the Tsoulis (2012) line-integral formula in pure PyTorch, making
the gravitational potential and acceleration differentiable with respect to
vertex positions and density. Supports CPU and CUDA execution.
"""

try:
    import torch
    from torch import Tensor
except ImportError as _err:
    raise ImportError(
        "polyhedral_gravity.torch requires PyTorch. Install it with:  pip install torch"
    ) from _err

import math

G_SI: float = 6.67430e-11  # m^3 kg^-1 s^-2


def evaluate(
    vertices: Tensor,
    faces: Tensor,
    density: float,
    computation_points: Tensor,
    gravitational_constant: float = G_SI,
) -> tuple[Tensor, Tensor, Tensor]:
    """
    Gravitational potential, acceleration, and gradient tensor for a polyhedron,
    differentiable w.r.t. vertex positions and density.

    Args:
        vertices:                (N, 3) vertex positions [m].
        faces:                   (F, 3) triangle vertex indices, integer dtype.
        density:                 Constant density [kg/m^3].
        computation_points:      (Q, 3) evaluation positions [m].
        gravitational_constant:  Gravitational constant, defaults to 6.67430e-11.

    Returns:
        potential:    (Q,)   gravitational potential [m^2/s^2].
        acceleration: (Q, 3) gravitational acceleration [m/s^2].
        tensor:       (Q, 6) gravity gradient tensor [1/s^2],
                             components ordered as [Vxx, Vyy, Vzz, Vxy, Vxz, Vyz].
    """

    eps = 1e-30

    # face vertex coords: (F, 3, 3)
    face_vertices = vertices[faces]  # [f, vertex_idx, xyz]

    # shift so each query point is at the origin: (Q, F, 3, 3)
    face_vertices_rel = face_vertices[None] - computation_points[:, None, None, :]

    v0 = face_vertices_rel[:, :, 0]  # (Q, F, 3)
    v1 = face_vertices_rel[:, :, 1]
    v2 = face_vertices_rel[:, :, 2]

    # ------------------------------------------------------------------ #
    #  Edge vectors G_pq and outward unit normal (face_normal)
    # ------------------------------------------------------------------ #
    edge0 = v1 - v0  # (Q, F, 3)
    edge1 = v2 - v1
    edge2 = v0 - v2

    edges = torch.stack([edge0, edge1, edge2], dim=2)  # (Q, F, 3, 3)  [q,f,seg,xyz]
    segment_start_points = torch.stack(
        [v0, v1, v2], dim=2
    )  # (Q, F, 3, 3)  segment start vertices

    cross_01 = torch.linalg.cross(edge0, edge1)  # (Q, F, 3)
    face_normal = cross_01 / cross_01.norm(dim=-1, keepdim=True).clamp(min=eps)

    # ------------------------------------------------------------------ #
    #  Segment unit normals n_pq = (G_pq x face_normal) / |...|
    # ------------------------------------------------------------------ #
    face_normal_expanded = face_normal[:, :, None].expand_as(edges)  # (Q, F, 3, 3)
    n_pq = torch.linalg.cross(edges, face_normal_expanded)
    n_pq = n_pq / n_pq.norm(dim=-1, keepdim=True).clamp(min=eps)

    # ------------------------------------------------------------------ #
    #  Signed plane distance, plane normal orientation, projection P'
    # ------------------------------------------------------------------ #
    hp_signed = (face_normal * v0).sum(dim=-1)  # (Q, F)   = N_p . v0
    h_p = hp_signed.abs()
    sigma_p = hp_signed.sign()  # sigma_p

    plane_projection = (
        hp_signed[..., None] * face_normal
    )  # (Q, F, 3)  P' = hp_signed * face_normal

    # ------------------------------------------------------------------ #
    #  Segment normal orientations sigma_pq and edge distances h_pq
    #
    #  h_pq_signed = -(n_pq . (P' - v_start)) is the signed distance from
    #  P' to the edge line within the face plane.  Deriving both sigma_pq
    #  and h_pq from this single expression keeps d(n_pq)/d(v) in the
    #  gradient graph; the equivalent norm-based formula (||P''-P'||) shares
    #  the same value but its graph misses that contribution.
    # ------------------------------------------------------------------ #
    plane_projection_expanded = plane_projection[:, :, None].expand_as(
        segment_start_points
    )  # (Q, F, 3, 3)
    h_pq_raw = -(n_pq * (plane_projection_expanded - segment_start_points)).sum(
        dim=-1
    )  # (Q, F, 3)

    # Threshold relative to edge length: floating-point noise in
    # normalised-vector arithmetic can produce |h_pq_raw| / edge_scale ~ 1e-16
    # for algebraically-zero distances. Zeroing these before .sign() prevents
    # misclassifying the singularity case (all_inside vs on_edge/on_vertex),
    # which would corrupt singularity_scalar by 2-4x.
    edge_length = edges.norm(dim=-1)  # (Q, F, 3)
    edge_scale = edge_length.amax(dim=-1, keepdim=True)  # (Q, F, 1)
    h_pq_signed = torch.where(
        h_pq_raw.abs() < 1e-10 * edge_scale,
        torch.zeros_like(h_pq_raw),
        h_pq_raw,
    )
    sigma_pq = h_pq_signed.sign()  # (Q, F, 3)
    h_pq = h_pq_signed.abs()  # (Q, F, 3)

    # ------------------------------------------------------------------ #
    #  Projection parameter t along each edge (needed for s1, s2)
    # ------------------------------------------------------------------ #
    edge_length_sq = (edges * edges).sum(dim=-1).clamp(min=eps)  # (Q, F, 3)
    t_pq = ((plane_projection_expanded - segment_start_points) * edges).sum(
        dim=-1
    ) / edge_length_sq  # (Q, F, 3)

    # ------------------------------------------------------------------ #
    #  Distances from P (origin) to segment endpoints
    # ------------------------------------------------------------------ #
    segment_end_points = torch.roll(
        segment_start_points, -1, dims=2
    )  # v_end = v_{(q+1)%3}
    l1 = segment_start_points.norm(dim=-1)  # (Q, F, 3)  |v_start|
    l2 = segment_end_points.norm(dim=-1)  # (Q, F, 3)  |v_end|

    # signed 1-D distances along the edge: s1 = -t*|G|,  s2 = (1-t)*|G|
    s1 = -t_pq * edge_length
    s2 = (1.0 - t_pq) * edge_length

    # ------------------------------------------------------------------ #
    #  Transcendental expressions LN_pq and AN_pq  (Tsoulis Eqs 14-15)
    # ------------------------------------------------------------------ #
    ln_num = (s2 + l2).clamp(min=eps)
    ln_den = (s1 + l1).clamp(min=eps)
    log_term = torch.log(ln_num / ln_den)

    h_p_e = h_p[:, :, None].expand_as(h_pq)  # (Q, F, 3)
    degenerate = (h_p_e < 1e-15) | (h_pq < 1e-15)
    arctan_term = torch.atan(h_p_e * s2 / (h_pq * l2).clamp(min=eps)) - torch.atan(
        h_p_e * s1 / (h_pq * l1).clamp(min=eps)
    )
    arctan_term = torch.where(degenerate, torch.zeros_like(arctan_term), arctan_term)

    # ------------------------------------------------------------------ #
    #  Per-face sums S1, S2, singularity correction  (Tsoulis Eqs 11-12)
    # ------------------------------------------------------------------ #
    sum1 = (h_pq_signed * log_term).sum(dim=-1)  # (Q, F)
    sum2 = (sigma_pq * arctan_term).sum(dim=-1)  # (Q, F)
    sum1_tensor = (n_pq * log_term[..., None]).sum(
        dim=-2
    )  # (Q, F, 3)  for gradient tensor

    # Singularity correction - three cases detected from the sigma_pq pattern:
    #   all +1        -> P' strictly inside face -> -2*pi*h_p
    #   two +1, one 0 -> P' on a face edge       -> -pi*h_p
    #   one +1, two 0 -> P' at a face vertex     -> -theta_v*h_p
    n_pos = (sigma_pq > 0).sum(dim=-1)  # (Q, F)
    n_zero = (sigma_pq == 0).sum(dim=-1)  # (Q, F)

    all_inside = n_pos == 3
    on_edge = (n_pos == 2) & (n_zero == 1)
    on_vertex = (n_pos == 1) & (n_zero == 2)

    edge_unit_vector = edges / edge_length[..., None].clamp(min=eps)  # (Q, F, 3, 3)
    cos_t0 = -(edge_unit_vector[..., 0, :] * edge_unit_vector[..., 2, :]).sum(
        dim=-1
    )  # (Q, F)
    cos_t1 = -(edge_unit_vector[..., 1, :] * edge_unit_vector[..., 0, :]).sum(dim=-1)
    cos_t2 = -(edge_unit_vector[..., 2, :] * edge_unit_vector[..., 1, :]).sum(dim=-1)
    theta_0 = torch.acos(cos_t0.clamp(-1.0, 1.0))
    theta_1 = torch.acos(cos_t1.clamp(-1.0, 1.0))
    theta_2 = torch.acos(cos_t2.clamp(-1.0, 1.0))
    theta_v = (
        theta_2 * (sigma_pq[..., 0] > 0).to(h_p.dtype)
        + theta_0 * (sigma_pq[..., 1] > 0).to(h_p.dtype)
        + theta_1 * (sigma_pq[..., 2] > 0).to(h_p.dtype)
    )  # (Q, F)

    singularity_scalar = torch.where(
        all_inside,
        -2.0 * math.pi * h_p,
        torch.where(
            on_edge,
            -math.pi * h_p,
            torch.where(on_vertex, -theta_v * h_p, torch.zeros_like(h_p)),
        ),
    )

    singularity_vector = torch.where(
        all_inside[..., None],
        -2.0 * math.pi * sigma_p[..., None] * face_normal,
        torch.where(
            on_edge[..., None],
            -math.pi * sigma_p[..., None] * face_normal,
            torch.where(
                on_vertex[..., None],
                -theta_v[..., None] * sigma_p[..., None] * face_normal,
                torch.zeros_like(face_normal),
            ),
        ),
    )  # (Q, F, 3)

    face_sum = sum1 + h_p * sum2 + singularity_scalar  # (Q, F)  per-face scalar sum

    # ------------------------------------------------------------------ #
    #  Sum over faces -> potential and acceleration
    # ------------------------------------------------------------------ #
    potential_terms = sigma_p * h_p * face_sum  # (Q, F)
    g_contributions = face_normal * face_sum[..., None]  # (Q, F, 3)

    potential = (gravitational_constant * density / 2.0) * potential_terms.sum(dim=-1)
    acceleration = -gravitational_constant * density * g_contributions.sum(dim=-2)

    # ------------------------------------------------------------------ #
    #  Gradient tensor  (Tsoulis Eq. 13)
    # ------------------------------------------------------------------ #
    tensor_subterm = (
        sum1_tensor + face_normal * (sigma_p * sum2)[..., None] + singularity_vector
    )  # (Q, F, 3)

    # diagonal components: Vxx, Vyy, Vzz  -> face_normal[i] * tensor_subterm[i]
    diag = face_normal * tensor_subterm  # (Q, F, 3)
    # off-diagonal: Vxy=N[0]*sub[1], Vxz=N[0]*sub[2], Vyz=N[1]*sub[2]
    offdiag_normal = torch.stack(
        [face_normal[..., 0], face_normal[..., 0], face_normal[..., 1]], dim=-1
    )
    offdiag_subterm = torch.stack(
        [tensor_subterm[..., 1], tensor_subterm[..., 2], tensor_subterm[..., 2]], dim=-1
    )
    offdiag = offdiag_normal * offdiag_subterm  # (Q, F, 3)

    tensor = (
        gravitational_constant
        * density
        * torch.cat([diag, offdiag], dim=-1).sum(dim=-2)
    )  # (Q, 6)

    return potential, acceleration, tensor
