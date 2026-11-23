// Web Worker: loads Pyodide + numpy, imports the real sim/ package served
// alongside this page, and runs a simulation from the parameters the UI sends.
// Running in a worker keeps the page responsive while Python works.

// Pinned Pyodide version (see README). numpy is bundled with Pyodide.
importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js");

// The Python files that make up the sim package (served at <site>/sim/).
const SIM_FILES = ["__init__.py", "thomas.py", "solver.py", "analytic.py", "model.py"];

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

self.onmessage = async (e) => {
  try {
    const pyodide = await ready;
    // Pass parameters as a JSON string so null -> None and nested dicts/lists
    // convert cleanly, with no JsProxy surprises.
    pyodide.globals.set("params_json", JSON.stringify(e.data));
    const result = await pyodide.runPythonAsync(`
import json
from sim.model import run_from_params
run_from_params(json.loads(params_json))
`);
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
