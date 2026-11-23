"""Step-4 test: the run_from_params UI adapter.

Checks that a plain parameter dict (the shape the web UI sends) produces a
clean, JSON-friendly result, and that velocity mode and Darcy mode agree when
K*i/n reproduces the same velocity.

Run from the repo root:  PYTHONPATH=. python examples/test_run_from_params.py
"""

import numpy as np

from sim.model import run_from_params


BASE = dict(
    length=100.0, depth=40.0, dx=1.0, dz=1.0, dt=0.5, total_time=100.0,
    source={"x": 20.0, "width": 10.0, "conc": 100.0, "leak_duration": 50.0},
    receptor={"x": 70.0, "z": 5.0},
    nframes=25,
)


def test_velocity_mode_shapes():
    p = dict(BASE, velocity_mode="velocity", layers=[
        {"thickness": 10.0, "velocity": 0.5, "dispersion": 2.0},
        {"thickness": 10.0, "velocity": 0.05, "dispersion": 0.2},
        {"thickness": 20.0, "velocity": 0.3, "dispersion": 1.0},
    ])
    r = run_from_params(p)

    ## . . Every frame must be a flat list of length nx*nz
    assert len(r["frames"]) == len(r["frame_times"])
    assert all(len(f) == r["nx"] * r["nz"] for f in r["frames"])
    ## . . Breakthrough has one value per time step and stays in [0, conc]
    assert len(r["breakthrough"]) == len(r["times"])
    vals = [v for v in r["breakthrough"] if v is not None]
    assert min(vals) > -1e-6 and max(vals) <= 100.0 + 1e-6
    assert r["interfaces"] == [10.0, 20.0]
    assert 0.0 < r["peak"] <= 100.0 + 1e-6
    print(f"[velocity] frames={len(r['frames'])}  Pe={r['Pe']:.3f}  "
          f"Cr={r['Cr']:.3f}  peak={r['peak']:.2f}")


def test_darcy_matches_velocity():
    ## . . Pick K, i, n so that v = K*i/n equals the velocities above
    #        0.5  = 5.0 * 0.01 / 0.1 ;  0.05 = 0.5*0.01/0.1 ; 0.3 = 3.0*0.01/0.1
    vel = dict(BASE, velocity_mode="velocity", layers=[
        {"thickness": 10.0, "velocity": 0.5, "dispersion": 2.0},
        {"thickness": 10.0, "velocity": 0.05, "dispersion": 0.2},
        {"thickness": 20.0, "velocity": 0.3, "dispersion": 1.0},
    ])
    dar = dict(BASE, velocity_mode="darcy", layers=[
        {"thickness": 10.0, "K": 5.0, "gradient": 0.01, "porosity": 0.1, "dispersion": 2.0},
        {"thickness": 10.0, "K": 0.5, "gradient": 0.01, "porosity": 0.1, "dispersion": 0.2},
        {"thickness": 20.0, "K": 3.0, "gradient": 0.01, "porosity": 0.1, "dispersion": 1.0},
    ])
    rv = run_from_params(vel)
    rd = run_from_params(dar)
    diff = max(abs(a - b) for a, b in zip(rv["frames"][-1], rd["frames"][-1]))
    print(f"[darcy]  max frame diff vs velocity mode = {diff:.2e}")
    assert diff < 1e-12, "Darcy mode disagrees with equivalent velocity mode"


def main():
    test_velocity_mode_shapes()
    test_darcy_matches_velocity()
    print("run_from_params tests passed.")


if __name__ == "__main__":
    main()
