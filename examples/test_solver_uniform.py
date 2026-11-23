"""Step-1 test: the constant-coefficient CN-ADI solver vs the analytic Gaussian.

Uniform medium, Gaussian spill away from the boundaries (so the C=0 Dirichlet
box matches the infinite-domain analytic solution). Two cases:
    1. vx = 0  (pure dispersion)
    2. vx != 0 (adds advection; plume centre should translate at vx * t)

Reports L2 and max error and the Pe / Cr diagnostics.
Run from the repo root:  PYTHONPATH=. python examples/test_solver_uniform.py
"""

import numpy as np

from sim.solver import cn_adi_step, cell_diagnostics
from sim.analytic import gaussian_spill


def run_case(vx, label):
    ## . . Grid (depth z and horizontal x), kept fine so Pe, Cr stay small
    xmin, xmax = 0.0, 2.0
    zmin, zmax = 0.0, 2.0
    dx = dz = 0.02
    nx = int(round((xmax - xmin) / dx)) + 1
    nz = int(round((zmax - zmin) / dz)) + 1
    xx = np.linspace(xmin, xmax, nx)
    zz = np.linspace(zmin, zmax, nz)
    X, Z = np.meshgrid(xx, zz, indexing="ij")           # . . X, Z are C[ix, iz]

    ## . . Medium + spill parameters
    D = 0.01
    C0, sigma0 = 1.0, 0.1
    x0, z0 = 0.7, 1.0
    dt = 0.05
    nt = 40
    t_end = nt * dt

    Pe, Cr = cell_diagnostics(D, vx, dx, dt)

    ## . . Initial condition = analytic field at t = 0
    C = gaussian_spill(X, Z, 0.0, C0, x0, z0, sigma0, D, vx)

    ## . . March in time
    for _ in range(nt):
        C = cn_adi_step(C, D, vx, dx, dz, dt)

    ## . . Compare to the analytic solution at t_end
    C_exact = gaussian_spill(X, Z, t_end, C0, x0, z0, sigma0, D, vx)
    err = C - C_exact
    l2 = np.sqrt(np.mean(err ** 2))
    linf = np.max(np.abs(err))

    print(f"--- {label} (vx = {vx}) ---")
    print(f"    Pe = {Pe:.3f}   Cr = {Cr:.3f}")
    print(f"    L2 error  = {l2:.3e}")
    print(f"    max error = {linf:.3e}")
    print(f"    peak numeric/analytic = {C.max():.4f} / {C_exact.max():.4f}")

    ## . . Numerical centre of mass should sit at x0 + vx*t_end
    cx_num = np.sum(X * C) / np.sum(C)
    print(f"    x centroid numeric/expected = {cx_num:.4f} / {x0 + vx * t_end:.4f}")
    return linf


def main():
    e0 = run_case(0.0, "pure dispersion")
    e1 = run_case(0.1, "advection + dispersion")
    assert e0 < 2e-3, "pure-dispersion error too large"
    assert e1 < 2e-3, "advection error too large"
    print("Uniform-medium solver tests passed.")


if __name__ == "__main__":
    main()
