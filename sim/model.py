"""Problem set-up and the time-stepping driver.

This module turns user inputs (domain, grid spacing, layer table, source) into
the depth-dependent fields D(z), vx(z) and the arrays the solver needs, and
provides `run_simulation`, the driver that marches the solution in time by
repeatedly calling `sim.solver.cn_adi_step`. The web worker and the examples
use this driver; the verification notebook calls the same `cn_adi_step`
directly with uniform inputs, so there is only ever ONE solver.

Array layout everywhere: C[ix, iz], axis 0 = x (horizontal), axis 1 = z (depth).
"""

from dataclasses import dataclass

import numpy as np

from .solver import cn_adi_step, cell_diagnostics


@dataclass
class Grid:
    """Regular x-z grid. z points downward from the surface (z = 0)."""
    length: float      # . . domain length in x [m]
    depth: float       # . . domain depth in z [m]
    dx: float          # . . x spacing [m]
    dz: float          # . . z spacing [m]
    nx: int
    nz: int
    xx: np.ndarray     # . . x coordinates (nx,)
    zz: np.ndarray     # . . z coordinates (nz,)
    X: np.ndarray      # . . meshgrid X, shape (nx, nz)
    Z: np.ndarray      # . . meshgrid Z, shape (nx, nz)


def build_grid(length, depth, dx, dz):
    """Build a regular grid from domain size and spacing."""
    nx = int(round(length / dx)) + 1
    nz = int(round(depth / dz)) + 1
    xx = np.linspace(0.0, length, nx)
    zz = np.linspace(0.0, depth, nz)
    X, Z = np.meshgrid(xx, zz, indexing="ij")
    return Grid(length, depth, dx, dz, nx, nz, xx, zz, X, Z)


def darcy_velocity(K, gradient, porosity):
    """Seepage velocity from Darcy's law: v = K * i / n.

    K : hydraulic conductivity [m/day], gradient i [-], porosity n [-].
    Lets the interface ask for K, i, n instead of velocity (not both).
    """
    return K * gradient / porosity


def layer_fields(zz, layers):
    """Map a layered-aquifer table onto per-depth-node D(z) and vx(z) arrays.

    input:
        zz     : depth coordinates (nz,) [m]
        layers : list of dicts, top-to-bottom, each with keys
                 'thickness' [m], 'velocity' [m/day], 'dispersion' [m^2/day]
    output:
        D        : dispersion at each depth node (nz,)
        vx       : horizontal velocity at each depth node (nz,)
        interfaces : depths of the internal layer boundaries (for plotting)
    """
    D = np.zeros_like(zz)
    vx = np.zeros_like(zz)

    bottoms = np.cumsum([L["thickness"] for L in layers])   # . . bottom depth of each layer
    for k, layer in enumerate(layers):
        top = 0.0 if k == 0 else bottoms[k - 1]
        if k == len(layers) - 1:
            mask = zz >= top                                # . . last layer runs to the base
        else:
            mask = (zz >= top) & (zz < bottoms[k])
        D[mask] = layer["dispersion"]
        vx[mask] = layer["velocity"]

    interfaces = bottoms[:-1]                               # . . drop the domain base
    return D, vx, interfaces


def source_mask(xx, x_center, width):
    """Boolean mask over x-nodes marking the top source patch (width centred at x_center)."""
    return np.abs(xx - x_center) <= width / 2.0


def receptor_index(grid, x_recv, z_recv):
    """Nearest grid indices (ix, iz) to a receptor point (house/well)."""
    ix = int(np.argmin(np.abs(grid.xx - x_recv)))
    iz = int(np.argmin(np.abs(grid.zz - z_recv)))
    return ix, iz


def run_simulation(grid, D, vx, dt, total_time, src_mask, source_conc,
                   leak_duration=None, C0=None, receptor=None, snapshot_times=None):
    """March the advection-dispersion solution in time with the CN-ADI solver.

    input:
        grid          : Grid from build_grid
        D, vx         : depth fields (nz,) from layer_fields (or scalars)
        dt            : time step [days]
        total_time    : total simulated time [days]
        src_mask      : bool (nx,) top source patch, or None for no source
        source_conc   : source concentration on the patch [mg/L]
        leak_duration : days the source leaks; None = leaks for the whole run.
                        After it stops, the patch reverts to no-flux ground.
        C0            : initial field (nx, nz); None -> clean aquifer (zeros)
        receptor      : (ix, iz) for the breakthrough curve, or None
        snapshot_times: list of times [days] to store concentration maps
    output:
        dict with 'times', 'breakthrough', 'snapshots', 'final', 'Pe', 'Cr'
    """
    nx, nz = grid.nx, grid.nz
    C = np.zeros((nx, nz)) if C0 is None else C0.copy()

    nt = int(round(total_time / dt))
    times = np.arange(nt + 1) * dt

    ## . . Diagnostics the user needs to judge the run
    Pe, Cr = cell_diagnostics(D, vx, grid.dx, dt)

    ## . . Storage: breakthrough curve at the receptor, and requested snapshots
    breakthrough = np.full(nt + 1, np.nan)
    if receptor is not None:
        breakthrough[0] = C[receptor]
    snapshots = {}
    snapshot_times = snapshot_times or []

    def maybe_store(step_time, field):
        for ts in snapshot_times:
            if abs(step_time - ts) < dt / 2 and ts not in snapshots:
                snapshots[ts] = field.copy()

    maybe_store(0.0, C)

    for n in range(1, nt + 1):
        t = times[n]
        ## . . Is the source still leaking at this step?
        leaking = src_mask is not None and (leak_duration is None or t <= leak_duration + 1e-9)
        step_mask = src_mask if leaking else None
        step_val = source_conc if leaking else 0.0

        C = cn_adi_step(C, D, vx, grid.dx, grid.dz, dt,
                        source_mask=step_mask, source_val=step_val)

        if receptor is not None:
            breakthrough[n] = C[receptor]
        maybe_store(t, C)

    return {
        "times": times,
        "breakthrough": breakthrough,
        "snapshots": snapshots,
        "final": C,
        "Pe": Pe,
        "Cr": Cr,
    }
