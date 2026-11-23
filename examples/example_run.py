"""Worked example: a contaminant leak into a three-layer aquifer.

Runs the full model through sim.model.run_simulation (which drives the CN-ADI
solver) and prints the key results. The expected output is listed at the bottom
and in the README, so a fresh clone can confirm the code works.

Run from anywhere:  python examples/example_run.py
"""

import os
import sys

import numpy as np

## . . Make `import sim` work no matter where this is launched from
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sim.model import (build_grid, layer_fields, source_mask,
                       receptor_index, run_simulation)


def main():
    ## ------------------------------------------------------------------
    ## 1. Domain and grid  (100 m long, 40 m deep aquifer cross-section)
    ## ------------------------------------------------------------------
    grid = build_grid(length=100.0, depth=40.0, dx=1.0, dz=1.0)

    ## ------------------------------------------------------------------
    ## 2. Layered aquifer: top-to-bottom (thickness, velocity, dispersion)
    ##    A slow, low-dispersion clay layer sits between two sandy layers.
    ## ------------------------------------------------------------------
    layers = [
        {"thickness": 10.0, "velocity": 0.50, "dispersion": 2.0},   # upper sand
        {"thickness": 10.0, "velocity": 0.05, "dispersion": 0.2},   # clay
        {"thickness": 20.0, "velocity": 0.30, "dispersion": 1.0},   # lower sand
    ]
    D, vx, interfaces = layer_fields(grid.zz, layers)

    ## ------------------------------------------------------------------
    ## 3. Source and receptor
    ##    A 10 m wide leak at x = 20 m on the surface, 100 mg/L, for 50 days.
    ##    A well sits 70 m downstream, 5 m deep.
    ## ------------------------------------------------------------------
    src = source_mask(grid.xx, x_center=20.0, width=10.0)
    recv = receptor_index(grid, x_recv=70.0, z_recv=5.0)

    ## ------------------------------------------------------------------
    ## 4. Run for 200 days with dt = 0.5 day
    ## ------------------------------------------------------------------
    res = run_simulation(
        grid, D, vx, dt=0.5, total_time=200.0,
        src_mask=src, source_conc=100.0, leak_duration=50.0,
        receptor=recv, snapshot_times=[50.0, 200.0],
    )

    ## ------------------------------------------------------------------
    ## 5. Report
    ## ------------------------------------------------------------------
    bt = res["breakthrough"]
    t_peak = res["times"][np.nanargmax(bt)]
    print(f"Layer interfaces at z = {interfaces} m")
    print(f"Cell Peclet Pe = {res['Pe']:.3f}  (warn if > 2)")
    print(f"Cell Courant Cr = {res['Cr']:.3f}")
    print(f"Receptor (ix,iz) = {recv}")
    print(f"Peak concentration at receptor = {np.nanmax(bt):.3f} mg/L at t = {t_peak:.1f} days")
    print(f"Final-field max concentration  = {res['final'].max():.3f} mg/L")


if __name__ == "__main__":
    main()


## Expected output (also in the README):
## Layer interfaces at z = [10. 20.] m
## Cell Peclet Pe = 0.300  (warn if > 2)
## Cell Courant Cr = 0.250
## Receptor (ix,iz) = (70, 5)
## Peak concentration at receptor = 25.959 mg/L at t = 124.5 days
## Final-field max concentration  = 18.419 mg/L
