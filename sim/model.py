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


def run_from_params(p):
    """Build and run a simulation from a plain parameter dict (UI adapter).

    This is a thin, UI-agnostic adapter: it takes a dictionary of inputs, calls
    the same `layer_fields` / `run_simulation` as everything else, and returns a
    dictionary of plain Python numbers and lists (no NumPy arrays) so it
    converts cleanly to a JavaScript object in the browser.

    Expected keys in `p`:
        length, depth, dx, dz, dt, total_time          (floats)
        velocity_mode : 'velocity' or 'darcy'
        layers : list of row dicts, top-to-bottom, with 'thickness',
                 'dispersion', and either 'velocity' (velocity mode) or
                 'K', 'gradient', 'porosity' (darcy mode)
        source   : {'x', 'width', 'conc', 'leak_duration' (None = continuous)}
        receptor : {'x', 'z'}
        nframes  : number of concentration-map frames to return (default 40)

    Returns a dict with grid info, animation frames (flattened row-major
    [ix*nz + iz]), the breakthrough curve, layer interfaces, Pe/Cr, and the peak
    concentration (for a shared colour scale).
    """
    grid = build_grid(p["length"], p["depth"], p["dx"], p["dz"])

    ## . . Build the layer table, converting Darcy inputs to velocity if needed
    mode = p.get("velocity_mode", "velocity")
    layers = []
    for row in p["layers"]:
        if mode == "darcy":
            v = darcy_velocity(row["K"], row["gradient"], row["porosity"])
        else:
            v = row["velocity"]
        layers.append({"thickness": row["thickness"],
                       "velocity": v,
                       "dispersion": row["dispersion"]})
    D, vx, interfaces = layer_fields(grid.zz, layers)

    src = source_mask(grid.xx, p["source"]["x"], p["source"]["width"])
    recv = receptor_index(grid, p["receptor"]["x"], p["receptor"]["z"])

    ## . . Choose evenly spaced snapshot frames ON the time-step grid so they
    ## . . match exactly (no tolerance issues).
    nt = int(round(p["total_time"] / p["dt"]))
    all_times = np.arange(nt + 1) * p["dt"]
    nframes = int(p.get("nframes", 40))
    idx = np.unique(np.round(np.linspace(0, nt, nframes)).astype(int))
    snap_times = [float(all_times[k]) for k in idx]

    res = run_simulation(
        grid, D, vx, p["dt"], p["total_time"],
        src_mask=src, source_conc=p["source"]["conc"],
        leak_duration=p["source"].get("leak_duration"),
        receptor=recv, snapshot_times=snap_times,
    )

    ## . . Order the frames and find the global peak for a shared colour scale
    frame_times = sorted(res["snapshots"].keys())
    frame_arrays = [res["snapshots"][t] for t in frame_times]
    peak = max((float(a.max()) for a in frame_arrays), default=0.0)

    ## . . Breakthrough curve: replace any NaN (defensive) with None for JS
    bt = res["breakthrough"]
    breakthrough = [None if np.isnan(v) else float(v) for v in bt]

    return {
        "nx": int(grid.nx), "nz": int(grid.nz),
        "length": float(grid.length), "depth": float(grid.depth),
        "interfaces": [float(z) for z in interfaces],
        "frame_times": [float(t) for t in frame_times],
        "frames": [a.flatten().tolist() for a in frame_arrays],
        "times": res["times"].tolist(),
        "breakthrough": breakthrough,
        "receptor": {"x": float(grid.xx[recv[0]]), "z": float(grid.zz[recv[1]])},
        "source": {"x": float(p["source"]["x"]), "width": float(p["source"]["width"])},
        "Pe": float(res["Pe"]), "Cr": float(res["Cr"]),
        "peak": peak,
    }
