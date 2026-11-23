// Web Worker: loads Pyodide + numpy, imports the real sim/ package served
// alongside this page, runs a three-layer example, and posts the result back.
// Running in a worker keeps the page responsive while Python works.

// Pinned Pyodide version (see README). numpy is bundled with Pyodide.
importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js");

// The Python files that make up the sim package (served at <site>/sim/).
const SIM_FILES = ["__init__.py", "thomas.py", "solver.py", "analytic.py", "model.py"];

// Python driver: build the model and return a plain dict (numbers + lists) so
// it converts cleanly to a JS object. C is flattened row-major: C[ix*nz + iz].
const PYTHON_RUN = `
import numpy as np
from sim.model import build_grid, layer_fields, source_mask, receptor_index, run_simulation

grid = build_grid(length=100.0, depth=40.0, dx=1.0, dz=1.0)
layers = [
    {"thickness": 10.0, "velocity": 0.50, "dispersion": 2.0},   # upper sand
    {"thickness": 10.0, "velocity": 0.05, "dispersion": 0.2},   # clay
    {"thickness": 20.0, "velocity": 0.30, "dispersion": 1.0},   # lower sand
]
D, vx, interfaces = layer_fields(grid.zz, layers)
src = source_mask(grid.xx, 20.0, 10.0)
recv = receptor_index(grid, 70.0, 5.0)

res = run_simulation(grid, D, vx, dt=0.5, total_time=150.0,
                     src_mask=src, source_conc=100.0, leak_duration=50.0,
                     receptor=recv, snapshot_times=[])

C = res["final"]
{
    "nx": int(grid.nx), "nz": int(grid.nz),
    "length": float(grid.length), "depth": float(grid.depth),
    "C": C.flatten().tolist(),
    "peak": float(C.max()),
    "Pe": float(res["Pe"]), "Cr": float(res["Cr"]),
    "interfaces": [float(z) for z in interfaces],
}
`;

async function init() {
  const pyodide = await loadPyodide();
  await pyodide.loadPackage("numpy");

  // Write the sim package into Pyodide's in-memory filesystem.
  pyodide.FS.mkdir("sim");
  for (const name of SIM_FILES) {
    const resp = await fetch("sim/" + name);
    if (!resp.ok) throw new Error("could not fetch sim/" + name + " (" + resp.status + ")");
    pyodide.FS.writeFile("sim/" + name, await resp.text());
  }
  await pyodide.runPythonAsync("import sys; sys.path.insert(0, '.')");
  return pyodide;
}

const ready = init();

self.onmessage = async () => {
  try {
    const pyodide = await ready;
    const result = await pyodide.runPythonAsync(PYTHON_RUN);
    const data = result.toJs({ dict_converter: Object.fromEntries });
    result.destroy();
    self.postMessage({ ok: true, data });
  } catch (err) {
    self.postMessage({ ok: false, error: String(err) });
  }
};

// Tell the page Pyodide + sim are ready as soon as init resolves.
ready.then(
  () => self.postMessage({ ok: true, ready: true }),
  (err) => self.postMessage({ ok: false, error: String(err) })
);
