"""Crank-Nicolson ADI (Peaceman-Rachford) solver for 2D advection-dispersion.

    dC/dt + vx dC/dx = d/dx( D dC/dx ) + d/dz( D dC/dz )

This module currently implements the CONSTANT-COEFFICIENT case (uniform D,
constant vx, horizontal flow only) with a fixed-concentration (Dirichlet) box.
It is a direct re-expression of the instructor's
`Solve_2D_Convection_Diffusion_homogeneous_Thomas` (class_notes_05) in our
C[ix, iz] layout, solving every implicit sweep with the batched Thomas solver.

Array layout: C[ix, iz], axis 0 = x (horizontal), axis 1 = z (depth).

One full step is two half-steps (Peaceman-Rachford):
    half 1: implicit in x (advection + x-diffusion), explicit in z
    half 2: implicit in z (z-diffusion), explicit in x
Each half-step is a batch of independent tridiagonal systems, one per grid line.
The 2*dt-scaled "factor of 4" algebra is exactly the course derivation.
"""

import numpy as np

from .thomas import thomas_batch


def cn_adi_step(C, D, vx, dx, dz, dt):
    """Advance the concentration field one time step dt.

    Boundaries are held at fixed (Dirichlet) values taken from the edges of the
    input C (the verification box uses C = 0 on all four sides). Only interior
    nodes are solved for.

    input:
        C  : concentration at time n, shape (nx, nz)  [C[ix, iz]]
        D  : dispersion coefficient (scalar) [m^2/day]
        vx : horizontal seepage velocity (scalar) [m/day]
        dx, dz : grid spacing [m]
        dt : time step [days]
    output:
        Cnew : concentration at time n+1, shape (nx, nz)
    """
    nx, nz = C.shape

    ## . . Dimensionless numbers (vz = 0, so there is no z-advection)
    AX = D * dt / (dx * dx)   # . . x diffusion number  alpha_x
    AZ = D * dt / (dz * dz)   # . . z diffusion number  alpha_z
    CX = vx * dt / dx         # . . x Courant number     C_x

    ## ------------------------------------------------------------------
    ## Half-step 1: implicit in x, explicit in z  ->  C* at n+1/2
    ## ------------------------------------------------------------------
    Cstar = C.copy()   # . . keeps the (fixed) boundary values in place

    ## . . Tridiagonal coefficients in x (constant along every line)
    ni = nx - 2        # . . number of interior x-nodes
    nl = nz - 2        # . . number of interior z-lines (one system each)
    aX = np.full((nl, ni), -(CX + 2 * AX))
    bX = np.full((nl, ni), 4 * (1 + AX))
    cX = np.full((nl, ni), (CX - 2 * AX))

    ## . . Right-hand side: explicit z-diffusion at time n, for interior nodes.
    ## . . Indexed [x, z] then transposed to [line=z, pos=x] for thomas_batch.
    Cx = C[1:nx - 1, :]                                   # (ni, nz)
    dX = (2 * AZ * (Cx[:, 0:nz - 2] + Cx[:, 2:nz])
          + 4 * (1 - AZ) * Cx[:, 1:nz - 1]).T             # (nl, ni)

    ## . . Fold the known Dirichlet x-boundaries into the end equations
    dX[:, 0] += (CX + 2 * AX) * C[0, 1:nz - 1]
    dX[:, -1] += -(CX - 2 * AX) * C[nx - 1, 1:nz - 1]

    Cstar[1:nx - 1, 1:nz - 1] = thomas_batch(aX, bX, cX, dX).T

    ## ------------------------------------------------------------------
    ## Half-step 2: implicit in z, explicit in x  ->  C at n+1
    ## ------------------------------------------------------------------
    Cnew = Cstar.copy()

    ## . . Tridiagonal coefficients in z (constant along every line here)
    mi = nz - 2        # . . number of interior z-nodes
    ml = nx - 2        # . . number of interior x-lines
    aZ = np.full((ml, mi), -2 * AZ)
    bZ = np.full((ml, mi), 4 * (1 + AZ))
    cZ = np.full((ml, mi), -2 * AZ)

    ## . . Right-hand side: explicit advection + x-diffusion at n+1/2.
    ## . . Already indexed [line=x, pos=z], no transpose needed.
    dZ = (-(CX - 2 * AX) * Cstar[2:nx, 1:nz - 1]
          + 4 * (1 - AX) * Cstar[1:nx - 1, 1:nz - 1]
          + (CX + 2 * AX) * Cstar[0:nx - 2, 1:nz - 1])    # (ml, mi)

    ## . . Fold the known Dirichlet z-boundaries (top/bottom) into the ends
    dZ[:, 0] += 2 * AZ * Cstar[1:nx - 1, 0]
    dZ[:, -1] += 2 * AZ * Cstar[1:nx - 1, nz - 1]

    Cnew[1:nx - 1, 1:nz - 1] = thomas_batch(aZ, bZ, cZ, dZ)

    return Cnew


def cell_diagnostics(D, vx, dx, dt):
    """Return the cell Peclet and Courant numbers used to judge a run.

    Pe = |vx| dx / D   (central differences oscillate when Pe > 2)
    Cr = |vx| dt / dx  (CN is stable for any dt but loses accuracy for large Cr)
    """
    Pe = np.inf if D == 0 else abs(vx) * dx / D
    Cr = abs(vx) * dt / dx
    return Pe, Cr
