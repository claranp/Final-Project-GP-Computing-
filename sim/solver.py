"""Crank-Nicolson ADI (Peaceman-Rachford) solver for 2D advection-dispersion.

    dC/dt + vx(z) dC/dx = d/dx( D(z) dC/dx ) + d/dz( D(z) dC/dz )

with horizontal flow only (vz = 0) and a dispersion coefficient that varies with
depth (layered aquifer). This is a generalisation of the instructor's
`Solve_2D_Convection_Diffusion_homogeneous_Thomas` (class_notes_05) for:
  * variable D(z), kept in divergence form with FACE HARMONIC MEANS,
  * a per-layer horizontal velocity vx(z),
  * the project boundary conditions (see below).

Array layout: C[ix, iz], axis 0 = x (horizontal), axis 1 = z (depth).
D and vx are given per depth node (length nz); a scalar is broadcast.

Boundary conditions (project physics, applied every step):
  * Left   (x = 0)   : C = 0            (clean upstream inflow, Dirichlet)
  * Right  (x = L)   : dC/dx = 0        (zero-gradient outflow)
  * Bottom (z = zmax): dC/dz = 0        (impermeable base, no-flux)
  * Top    (z = 0)   : C = C_source on the source patch (Dirichlet),
                       dC/dz = 0 elsewhere (no-flux)

One full step is two half-steps (Peaceman-Rachford), each a batch of tridiagonal
systems solved with the hand-written Thomas algorithm:
    half 1: implicit in x (advection + x-diffusion), explicit in z
    half 2: implicit in z (z-diffusion),             explicit in x
The 2*dt-scaled "factor of 4" algebra is exactly the course derivation.

How the physics enters the tridiagonal coefficients:
  * Because D, vx depend only on depth, along any fixed-z x-line they are
    CONSTANT, so the x-sweep coefficients are constant along each line (but
    differ line to line). The advection term is central-differenced and sits on
    the x sub/super-diagonals as +/- Cx.
  * The z-sweep coefficients VARY along the line: the d/dz(D dC/dz) stencil uses
    the harmonic mean of neighbouring D values at the cell faces j+/-1/2, which
    keeps the diffusive flux continuous across a layer contrast.
"""

import numpy as np

from .thomas import thomas_batch


def _as_depth_array(value, nz):
    """Return value as a length-nz float array (broadcast a scalar)."""
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return np.full(nz, float(arr))
    if arr.shape != (nz,):
        raise ValueError(f"expected scalar or length-{nz} array, got {arr.shape}")
    return arr


def _apply_boundaries(C, source_mask, source_val):
    """Enforce the project boundary conditions on the edges of C in place.

    x-boundaries first, then z-boundaries so the top/bottom no-flux copies pick
    up the already-updated corner columns.
    """
    ## . . Left: clean inflow (Dirichlet 0); Right: zero-gradient outflow
    C[0, :] = 0.0
    C[-1, :] = C[-2, :]

    ## . . Bottom: no-flux (copy the node just above the base)
    C[:, -1] = C[:, -2]

    ## . . Top: no-flux everywhere, then overwrite the source patch with C_source
    C[:, 0] = C[:, 1]
    if source_mask is not None:
        C[source_mask, 0] = source_val


def cn_adi_step(C, D, vx, dx, dz, dt, source_mask=None, source_val=0.0):
    """Advance the concentration field one time step dt.

    input:
        C  : concentration at time n, shape (nx, nz)  [C[ix, iz]]
        D  : dispersion coefficient, scalar or (nz,) [m^2/day]
        vx : horizontal seepage velocity, scalar or (nz,) [m/day]
        dx, dz : grid spacing [m]
        dt : time step [days]
        source_mask : bool array (nx,), True on the top source patch, or None
        source_val  : source concentration on the patch [mg/L]
    output:
        Cnew : concentration at time n+1, shape (nx, nz)
    """
    nx, nz = C.shape
    Dz = _as_depth_array(D, nz)
    Vz = _as_depth_array(vx, nz)

    ## . . Per-depth dimensionless numbers (constant along each x-line)
    AX = Dz * dt / (dx * dx)      # . . x diffusion number  alpha_x(z)
    CX = Vz * dt / dx             # . . x Courant number      C_x(z)

    ## . . Face dispersion via HARMONIC MEAN between neighbouring depth cells.
    ## . . Dface[k] is the value at the face between nodes k and k+1.
    Dface = 2.0 * Dz[:-1] * Dz[1:] / (Dz[:-1] + Dz[1:])
    ## . . z diffusion numbers at the faces just above (+) and below (-) each
    ## . . interior node j = 1..nz-2:
    AZp = Dface[1:nz - 1] * dt / (dz * dz)   # . . face j+1/2  (toward larger j)
    AZm = Dface[0:nz - 2] * dt / (dz * dz)   # . . face j-1/2  (toward smaller j)

    ## ------------------------------------------------------------------
    ## Half-step 1: implicit in x, explicit in z  ->  C* at n+1/2
    ## ------------------------------------------------------------------
    _apply_boundaries(C, source_mask, source_val)   # . . so explicit z read is correct
    Cstar = C.copy()

    ni = nx - 2        # . . interior x-nodes (unknowns per line)
    nl = nz - 2        # . . interior z-lines (one tridiagonal system each)

    ## . . Coefficients: per-line constants (depend on layer j), broadcast over i
    AXl = AX[1:nz - 1][:, None]        # . . shape (nl, 1)
    CXl = CX[1:nz - 1][:, None]
    aX = np.broadcast_to(-(CXl + 2 * AXl), (nl, ni)).copy()
    bX = np.broadcast_to(4 * (1 + AXl), (nl, ni)).copy()
    cX = np.broadcast_to((CXl - 2 * AXl), (nl, ni)).copy()

    ## . . RHS: explicit z-diffusion at n (face means), indexed [x, z] -> [z, x]
    Cxint = C[1:nx - 1, :]                                  # (ni, nz)
    dX = (2 * AZm[None, :] * Cxint[:, 0:nz - 2]
          + (4 - 2 * AZp - 2 * AZm)[None, :] * Cxint[:, 1:nz - 1]
          + 2 * AZp[None, :] * Cxint[:, 2:nz]).T            # (nl, ni)

    ## . . Left Dirichlet folds into the first equation (C[0,:] = 0 -> no effect)
    dX[:, 0] += (CX[1:nz - 1] + 2 * AX[1:nz - 1]) * C[0, 1:nz - 1]
    ## . . Right zero-gradient: C*[nx-1] = C*[nx-2] (the last unknown), fold c -> b
    bX[:, -1] += cX[:, -1]
    cX[:, -1] = 0.0

    Cstar[1:nx - 1, 1:nz - 1] = thomas_batch(aX, bX, cX, dX).T
    _apply_boundaries(Cstar, source_mask, source_val)

    ## ------------------------------------------------------------------
    ## Half-step 2: implicit in z, explicit in x  ->  C at n+1
    ## ------------------------------------------------------------------
    Cnew = Cstar.copy()

    mi = nz - 2        # . . interior z-nodes (unknowns per line)
    ml = nx - 2        # . . interior x-lines

    ## . . Coefficients: vary along the line via face harmonic means, same for
    ## . . every x-line, so broadcast the (mi,) columns over the ml lines.
    aZ = np.broadcast_to(-2 * AZm, (ml, mi)).copy()
    bZ = np.broadcast_to(4 + 2 * AZp + 2 * AZm, (ml, mi)).copy()
    cZ = np.broadcast_to(-2 * AZp, (ml, mi)).copy()

    ## . . RHS: explicit advection + x-diffusion at n+1/2, indexed [x, z]
    AXi = AX[1:nz - 1][None, :]        # . . shape (1, mi)
    CXi = CX[1:nz - 1][None, :]
    dZ = (-(CXi - 2 * AXi) * Cstar[2:nx, 1:nz - 1]
          + 4 * (1 - AXi) * Cstar[1:nx - 1, 1:nz - 1]
          + (CXi + 2 * AXi) * Cstar[0:nx - 2, 1:nz - 1])    # (ml, mi)

    ## . . Top boundary (first unknown, j = 1):
    ## . .   source column  -> C[i,0] = C_source known, add to RHS
    ## . .   no-flux column -> C[i,0] = C[i,1] (the unknown), fold a -> b
    top_src = np.zeros(ml, dtype=bool)
    if source_mask is not None:
        top_src = source_mask[1:nx - 1]
    dZ[top_src, 0] += 2 * AZm[0] * source_val          # . . = -a * C_source
    bZ[~top_src, 0] += aZ[~top_src, 0]                 # . . no-flux fold
    ## . . Bottom boundary (last unknown): no-flux C[i,nz-1] = C[i,nz-2], fold c -> b
    bZ[:, -1] += cZ[:, -1]

    Cnew[1:nx - 1, 1:nz - 1] = thomas_batch(aZ, bZ, cZ, dZ)
    _apply_boundaries(Cnew, source_mask, source_val)

    return Cnew


def cell_diagnostics(D, vx, dx, dt):
    """Return the worst-case cell Peclet and Courant numbers for a run.

    Pe = |vx| dx / D   (central differences oscillate when Pe > 2)
    Cr = |vx| dt / dx  (CN is stable for any dt but loses accuracy for large Cr)
    Accepts scalars or depth arrays; reduces to the worst (largest) cell.
    """
    D = np.asarray(D, dtype=float)
    vx = np.asarray(vx, dtype=float)
    Dmin = np.min(D)
    vmax = np.max(np.abs(vx))
    Pe = np.inf if Dmin == 0 else vmax * dx / Dmin
    Cr = vmax * dt / dx
    return float(Pe), float(Cr)
