"""Step-1 test: variable-D(z) solver, face harmonic means, and the mixed top BC.

B1 (regression): a UNIFORM D passed as a depth array must give exactly the same
    result as the scalar-D path (confirms the face-harmonic-mean machinery
    collapses to the constant case).
B2 (layered physics): a plume straddling a high-D / low-D interface must spread
    faster into the high-D layer, so its vertical centre of mass drifts toward
    the high-D side. Solution must stay bounded in [0, C0].
B3 (mixed top BC): a continuous top source patch into a clean aquifer must grow
    a plume that penetrates downward and advects downstream, with 0 <= C <=
    C_source everywhere and clean (~0) water far upstream.

Run from the repo root:  PYTHONPATH=. python examples/test_solver_layered.py
"""

import numpy as np

from sim.solver import cn_adi_step
from sim.analytic import gaussian_spill


def make_grid():
    xmin, xmax, zmin, zmax = 0.0, 2.0, 0.0, 2.0
    dx = dz = 0.02
    nx = int(round((xmax - xmin) / dx)) + 1
    nz = int(round((zmax - zmin) / dz)) + 1
    xx = np.linspace(xmin, xmax, nx)
    zz = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xx, zz, indexing="ij")
    return xx, zz, X, Z, dx, dz, nx, nz


def test_B1_uniform_array_matches_scalar():
    xx, zz, X, Z, dx, dz, nx, nz = make_grid()
    D, vx, dt, nt = 0.01, 0.1, 0.05, 20
    C0 = gaussian_spill(X, Z, 0.0, 1.0, 0.7, 1.0, 0.1, D, vx)

    C_scalar = C0.copy()
    C_array = C0.copy()
    Darr = np.full(nz, D)
    varr = np.full(nz, vx)
    for _ in range(nt):
        C_scalar = cn_adi_step(C_scalar, D, vx, dx, dz, dt)
        C_array = cn_adi_step(C_array, Darr, varr, dx, dz, dt)
    err = np.max(np.abs(C_scalar - C_array))
    print(f"[B1] uniform array vs scalar : max diff = {err:.2e}")
    assert err < 1e-12, "array-D path differs from scalar-D path"


def test_B2_layered_centroid_drift():
    xx, zz, X, Z, dx, dz, nx, nz = make_grid()
    z_interface = 1.0
    ## . . Upper layer (small z, nearer the surface) is high-D; lower layer low-D
    D = np.where(zz < z_interface, 0.02, 0.002)
    vx = 0.0                                   # . . isolate dispersion
    dt, nt = 0.05, 30

    C0 = gaussian_spill(X, Z, 0.0, 1.0, 1.0, z_interface, 0.1, 0.01, 0.0)
    C = C0.copy()
    for _ in range(nt):
        C = cn_adi_step(C, D, vx, dx, dz, dt)

    ## . . Vertical centre of mass (should drift toward the high-D upper layer)
    z_cm0 = np.sum(Z * C0) / np.sum(C0)
    z_cm = np.sum(Z * C) / np.sum(C)
    print(f"[B2] z centroid  {z_cm0:.4f} -> {z_cm:.4f} (interface {z_interface})")
    print(f"[B2] C range = [{C.min():.3e}, {C.max():.3e}]")
    assert z_cm < z_cm0 - 1e-3, "centroid did not drift toward the high-D layer"
    assert C.min() > -1e-6 and C.max() <= 1.0 + 1e-6, "solution left [0, C0]"


def test_B3_top_source_patch():
    xx, zz, X, Z, dx, dz, nx, nz = make_grid()
    D, vx, dt, nt = 0.02, 0.1, 0.05, 60
    Csrc = 1.0

    ## . . Source patch on the top boundary: 0.6 <= x <= 1.0
    source_mask = (xx >= 0.6) & (xx <= 1.0)

    C = np.zeros((nx, nz))                     # . . clean aquifer
    for _ in range(nt):
        C = cn_adi_step(C, D, vx, dx, dz, dt,
                        source_mask=source_mask, source_val=Csrc)

    j_shallow = 3                              # . . a few cells below the surface
    print(f"[B3] peak just below patch = {C[source_mask, j_shallow].max():.3f}")
    print(f"[B3] C range = [{C.min():.3e}, {C.max():.3e}]")
    ## . . Plume exists below the patch, stays within [0, Csrc], clean far upstream
    assert C[source_mask, j_shallow].max() > 0.1, "no plume below the source"
    assert C.min() > -1e-6 and C.max() <= Csrc + 1e-6, "solution left [0, Csrc]"
    assert C[0, :].max() < 1e-9, "upstream inflow boundary is not clean"
    ## . . Downstream (right of patch) should see more than upstream (left of patch)
    assert C[xx > 1.0, :].sum() > C[xx < 0.6, :].sum(), "plume did not advect downstream"


def main():
    test_B1_uniform_array_matches_scalar()
    test_B2_layered_centroid_drift()
    test_B3_top_source_patch()
    print("Layered / variable-D / mixed-BC tests passed.")


if __name__ == "__main__":
    main()
