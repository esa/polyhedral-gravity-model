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
        "polyhedral_gravity.torch requires PyTorch. "
        "Install it with:  pip install torch"
    ) from _err

import math

G_SI: float = 6.67430e-11  # m^3 kg^-1 s^-2


def evaluate(
    vertices: Tensor,
    faces: Tensor,
    density: float,
    computation_points: Tensor,
    G: float = G_SI,
) -> tuple[Tensor, Tensor, Tensor]:
    """
    Gravitational potential, acceleration, and gradient tensor for a polyhedron,
    differentiable w.r.t. vertex positions and density.

    Args:
        vertices:           (N, 3) vertex positions [m].
        faces:              (F, 3) triangle vertex indices, integer dtype.
        density:            Constant density [kg/m^3].
        computation_points: (Q, 3) evaluation positions [m].
        G:                  Gravitational constant, defaults to 6.67430e-11.

    Returns:
        potential:    (Q,)   gravitational potential [m^2/s^2].
        acceleration: (Q, 3) gravitational acceleration [m/s^2].
        tensor:       (Q, 6) gravity gradient tensor [1/s^2],
                             components ordered as [Vxx, Vyy, Vzz, Vxy, Vxz, Vyz].
    """

    EPS = 1e-30

    # face vertex coords: (F, 3, 3)
    v = vertices[faces]  # [f, vertex_idx, xyz]

    # shift so each query point is at the origin: (Q, F, 3, 3)
    v_p = v[None] - computation_points[:, None, None, :]

    v0 = v_p[:, :, 0]  # (Q, F, 3)
    v1 = v_p[:, :, 1]
    v2 = v_p[:, :, 2]

    # ------------------------------------------------------------------ #
    #  Edge vectors G_pq and outward unit normal N_p
    # ------------------------------------------------------------------ #
    G0 = v1 - v0  # (Q, F, 3)
    G1 = v2 - v1
    G2 = v0 - v2

    G_segs = torch.stack([G0, G1, G2], dim=2)   # (Q, F, 3, 3)  [q,f,seg,xyz]
    V_segs = torch.stack([v0, v1, v2], dim=2)   # (Q, F, 3, 3)  segment start vertices

    cross_01 = torch.linalg.cross(G0, G1)       # (Q, F, 3)
    N_p = cross_01 / cross_01.norm(dim=-1, keepdim=True).clamp(min=EPS)

    # ------------------------------------------------------------------ #
    #  Segment unit normals n_pq = (G_pq x N_p) / |...|
    # ------------------------------------------------------------------ #
    N_p_e = N_p[:, :, None].expand_as(G_segs)   # (Q, F, 3, 3)
    n_pq = torch.linalg.cross(G_segs, N_p_e)
    n_pq = n_pq / n_pq.norm(dim=-1, keepdim=True).clamp(min=EPS)

    # ------------------------------------------------------------------ #
    #  Signed plane distance, plane normal orientation, projection P'
    # ------------------------------------------------------------------ #
    hp_signed = (N_p * v0).sum(dim=-1)          # (Q, F)   = N_p . v0
    h_p = hp_signed.abs()
    sigma_p = hp_signed.sign()                   # sigma_p

    Pp = hp_signed[..., None] * N_p             # (Q, F, 3)  P' = hp_signed * N_p

    # ------------------------------------------------------------------ #
    #  Segment normal orientations sigma_pq and edge distances h_pq
    #
    #  h_pq_signed = -(n_pq . (P' - v_start)) is the signed distance from
    #  P' to the edge line within the face plane.  Deriving both sigma_pq
    #  and h_pq from this single expression keeps d(n_pq)/d(v) in the
    #  gradient graph; the equivalent norm-based formula (||P''-P'||) shares
    #  the same value but its graph misses that contribution.
    # ------------------------------------------------------------------ #
    Pp_e = Pp[:, :, None].expand_as(V_segs)                          # (Q, F, 3, 3)
    h_pq_raw = -(n_pq * (Pp_e - V_segs)).sum(dim=-1)                 # (Q, F, 3)

    # Threshold relative to edge length: floating-point noise in
    # normalised-vector arithmetic can produce |h_pq_raw| / edge_scale ~ 1e-16
    # for algebraically-zero distances. Zeroing these before .sign() prevents
    # misclassifying the singularity case (all_inside vs on_edge/on_vertex),
    # which would corrupt singA by 2-4x.
    G_norm = G_segs.norm(dim=-1)                                      # (Q, F, 3)
    edge_scale = G_norm.amax(dim=-1, keepdim=True)                    # (Q, F, 1)
    h_pq_signed = torch.where(
        h_pq_raw.abs() < 1e-10 * edge_scale,
        torch.zeros_like(h_pq_raw),
        h_pq_raw,
    )
    sigma_pq = h_pq_signed.sign()                                     # (Q, F, 3)
    h_pq = h_pq_signed.abs()                                          # (Q, F, 3)

    # ------------------------------------------------------------------ #
    #  Projection parameter t along each edge (needed for s1, s2)
    # ------------------------------------------------------------------ #
    G_sq = (G_segs * G_segs).sum(dim=-1).clamp(min=EPS)       # (Q, F, 3)
    t_pq = ((Pp_e - V_segs) * G_segs).sum(dim=-1) / G_sq      # (Q, F, 3)

    # ------------------------------------------------------------------ #
    #  Distances from P (origin) to segment endpoints
    # ------------------------------------------------------------------ #
    V_segs_end = torch.roll(V_segs, -1, dims=2)   # v_end = v_{(q+1)%3}
    l1 = V_segs.norm(dim=-1)                       # (Q, F, 3)  |v_start|
    l2 = V_segs_end.norm(dim=-1)                   # (Q, F, 3)  |v_end|

    # signed 1-D distances along the edge: s1 = -t*|G|,  s2 = (1-t)*|G|
    s1 = -t_pq * G_norm
    s2 = (1.0 - t_pq) * G_norm

    # ------------------------------------------------------------------ #
    #  Transcendental expressions LN_pq and AN_pq  (Tsoulis Eqs 14-15)
    # ------------------------------------------------------------------ #
    ln_num = (s2 + l2).clamp(min=EPS)
    ln_den = (s1 + l1).clamp(min=EPS)
    LN = torch.log(ln_num / ln_den)

    h_p_e = h_p[:, :, None].expand_as(h_pq)       # (Q, F, 3)
    degenerate = (h_p_e < 1e-15) | (h_pq < 1e-15)
    AN = torch.atan(h_p_e * s2 / (h_pq * l2).clamp(min=EPS)) \
       - torch.atan(h_p_e * s1 / (h_pq * l1).clamp(min=EPS))
    AN = torch.where(degenerate, torch.zeros_like(AN), AN)

    # ------------------------------------------------------------------ #
    #  Per-face sums S1, S2, singularity correction  (Tsoulis Eqs 11-12)
    # ------------------------------------------------------------------ #
    sum1 = (h_pq_signed * LN).sum(dim=-1)          # (Q, F)
    sum2 = (sigma_pq * AN).sum(dim=-1)             # (Q, F)
    sum1_tensor = (n_pq * LN[..., None]).sum(dim=-2)  # (Q, F, 3)  for gradient tensor

    # singA correction - three cases detected from the sigma_pq pattern:
    #   all +1        -> P' strictly inside face -> -2*pi*h_p
    #   two +1, one 0 -> P' on a face edge       -> -pi*h_p
    #   one +1, two 0 -> P' at a face vertex     -> -theta_v*h_p
    n_pos  = (sigma_pq > 0).sum(dim=-1)    # (Q, F)
    n_zero = (sigma_pq == 0).sum(dim=-1)   # (Q, F)

    all_inside = n_pos == 3
    on_edge    = (n_pos == 2) & (n_zero == 1)
    on_vertex  = (n_pos == 1) & (n_zero == 2)

    G_hat = G_segs / G_norm[..., None].clamp(min=EPS)              # (Q, F, 3, 3)
    cos_t0 = -(G_hat[..., 0, :] * G_hat[..., 2, :]).sum(dim=-1)   # (Q, F)
    cos_t1 = -(G_hat[..., 1, :] * G_hat[..., 0, :]).sum(dim=-1)
    cos_t2 = -(G_hat[..., 2, :] * G_hat[..., 1, :]).sum(dim=-1)
    theta_0 = torch.acos(cos_t0.clamp(-1.0, 1.0))
    theta_1 = torch.acos(cos_t1.clamp(-1.0, 1.0))
    theta_2 = torch.acos(cos_t2.clamp(-1.0, 1.0))
    theta_v = (theta_2 * (sigma_pq[..., 0] > 0).to(h_p.dtype) +
               theta_0 * (sigma_pq[..., 1] > 0).to(h_p.dtype) +
               theta_1 * (sigma_pq[..., 2] > 0).to(h_p.dtype))    # (Q, F)

    sing_A = torch.where(all_inside, -2.0 * math.pi * h_p,
             torch.where(on_edge,    -math.pi        * h_p,
             torch.where(on_vertex,  -theta_v        * h_p,
                                      torch.zeros_like(h_p))))

    sing_B = torch.where(all_inside[..., None], -2.0 * math.pi * sigma_p[..., None] * N_p,
             torch.where(on_edge[..., None],    -math.pi        * sigma_p[..., None] * N_p,
             torch.where(on_vertex[..., None],  -theta_v[..., None] * sigma_p[..., None] * N_p,
                                                 torch.zeros_like(N_p))))  # (Q, F, 3)

    A_p = sum1 + h_p * sum2 + sing_A              # (Q, F)  per-face scalar sum

    # ------------------------------------------------------------------ #
    #  Sum over faces -> potential and acceleration
    # ------------------------------------------------------------------ #
    V_contributions = sigma_p * h_p * A_p          # (Q, F)
    g_contributions = N_p * A_p[..., None]         # (Q, F, 3)

    potential    = (G * density / 2.0) * V_contributions.sum(dim=-1)
    acceleration = -G * density * g_contributions.sum(dim=-2)

    # ------------------------------------------------------------------ #
    #  Gradient tensor  (Tsoulis Eq. 13)
    # ------------------------------------------------------------------ #
    subSum = sum1_tensor + N_p * (sigma_p * sum2)[..., None] + sing_B   # (Q, F, 3)

    # diagonal components: Vxx, Vyy, Vzz  -> N_p[i] * subSum[i]
    diag = N_p * subSum                                                  # (Q, F, 3)
    # off-diagonal: Vxy=N[0]*sub[1], Vxz=N[0]*sub[2], Vyz=N[1]*sub[2]
    off_N   = torch.stack([N_p[..., 0],    N_p[..., 0],    N_p[..., 1]],    dim=-1)
    off_sub = torch.stack([subSum[..., 1], subSum[..., 2], subSum[..., 2]], dim=-1)
    offdiag = off_N * off_sub                                            # (Q, F, 3)

    tensor = G * density * torch.cat([diag, offdiag], dim=-1).sum(dim=-2)  # (Q, 6)

    return potential, acceleration, tensor
