# Aquifer Contaminant Transport — 2D Advection–Dispersion (CN-ADI)

Final project for **GPGN 436/536** (Geophysical Computing), Colorado School of Mines.
Student: **Clara Pugh**.

An interactive simulation of a contaminant plume leaking from the ground surface into a
layered aquifer, modelled by the 2D advection–dispersion equation and solved with a
**Crank–Nicolson ADI** scheme written from scratch (no black-box PDE solvers). The same
Python package runs in the browser (via Pyodide) and in the verification notebook.

**Live demo:** https://claranp.github.io/Final-Project-GP-Computing-/

---

## The problem

A contaminant source at the surface (a landfill or chemical plant) leaks into a layered
aquifer. We model how the plume **moves** (advection by groundwater flow) and **spreads**
(dispersion) through a 2D vertical cross-section — `x` horizontal, `z` depth — toward
downstream receptors such as wells or houses.

## Governing equation

The advection–dispersion (convection–diffusion) equation:

```
∂C/∂t + v·∇C = ∇·(D ∇C)
```

- `C(x, z, t)` — concentration [mg/L]
- `v = (vx(z), 0)` — horizontal groundwater velocity, per layer [m/day]
- `D(z)` — dispersion coefficient, per layer [m²/day] (spatially variable → kept in
  divergence form)

## Numerical method

- **Crank–Nicolson ADI (Peaceman–Rachford)** on a regular grid, split into two half-steps:
  - half-step 1: implicit in `x` (advection + x-dispersion), explicit in `z`
  - half-step 2: implicit in `z` (z-dispersion), explicit in `x`
- Each implicit sweep is a tridiagonal system solved by a **hand-written Thomas algorithm**
  (`sim/thomas.py`) — a readable single-line version plus a NumPy-vectorised batched version
  that solves all grid lines at once (the speed path used in the browser).
- **Variable dispersion** `D(z)` is kept in divergence form using **harmonic-mean face
  values**, which keeps the diffusive flux continuous across a layer contrast.
- **Advection** uses central differences in space.
- The scheme is **second-order accurate** in space and time and, being Crank–Nicolson,
  **unconditionally stable** (no time-step limit).

### Boundary and initial conditions

- **Top (z = 0):** fixed-concentration (Dirichlet) source patch of user-set width and
  location, continuous or for a set leak duration; no-flux elsewhere on top.
- **Bottom:** no-flux (impermeable base).
- **Left (upstream):** `C = 0` (clean inflow).
- **Right (downstream):** zero-gradient outflow.
- **Initial condition:** `C = 0` everywhere (or a Gaussian spill, used for verification).

### Diagnostics

- **Cell Péclet number** `Pe = |v|·Δx/D` — central differences oscillate when `Pe > 2`.
  The interface warns (and offers to cancel) before running if `Pe > 2`.
- **Courant number** `Cr = |v|·Δt/Δx` — does not limit stability (CN is unconditionally
  stable) but affects accuracy at sharp fronts.

---

## Repository layout

```
sim/                    standalone simulation package (plain Python + NumPy)
  thomas.py             Thomas algorithm (single-line + batched)
  solver.py             CN-ADI time step + Pe/Cr diagnostics
  model.py              grid, layers → D(z)/vx(z), Darcy, source, run driver + UI adapter
  analytic.py           Gaussian analytic solution (verification only)
web/                    browser interface
  index.html, style.css, app.js   controls + plotting (JS)
  worker.js             loads Pyodide + numpy, imports sim/, runs the model
verification.ipynb      imports the SAME sim/ code; convergence study (saved with outputs)
examples/               test scripts + worked example (expected results below)
requirements.txt        pinned local/verification dependencies
.github/workflows/deploy.yml   GitHub Actions → GitHub Pages
```

---

## Run it locally

Requires Python 3 and the pinned packages:

```bash
pip install -r requirements.txt
```

**Solver + model tests** (run from the repo root):

```bash
PYTHONPATH=. python examples/test_thomas.py           # Thomas vs numpy.linalg.solve
PYTHONPATH=. python examples/test_solver_uniform.py   # vs analytic Gaussian
PYTHONPATH=. python examples/test_solver_layered.py   # variable D, mixed BCs
PYTHONPATH=. python examples/test_run_from_params.py  # UI adapter
```

**Worked example:**

```bash
python examples/example_run.py
```

**Verification notebook:**

```bash
jupyter lab verification.ipynb     # or: jupyter notebook verification.ipynb
```

**The website locally** (assemble `web/` + `sim/` into one directory and serve it, exactly
like the deploy workflow does):

```bash
mkdir -p _site && cp -r web/. _site/ && cp -r sim _site/sim
cd _site && python -m http.server 8000
# then open http://localhost:8000
```

---

## Deployment

A GitHub Actions workflow (`.github/workflows/deploy.yml`) assembles the site (`web/` plus the
`sim/` package) and publishes it to **GitHub Pages** on every push to `main`. In the browser,
**Pyodide `v0.26.4`** (pinned) loads **numpy only** and imports the real `sim/` package — the
exact same code the verification notebook uses. The solver runs in a Web Worker so the page
stays responsive.

---

## Verification

`verification.ipynb` checks the solver against the analytic spreading/translating Gaussian on a
uniform medium (the same `cn_adi_step` the website runs — there is no second solver):

- pure dispersion (`vx = 0`) and advection (`vx ≠ 0`): numeric vs analytic maps, with L2/max
  error; the plume centre of mass advects to exactly `x0 + vx·t`.
- grid-refinement convergence: **observed slope ≈ 2.00 in Δx** and **≈ 2.03 in Δt**,
  confirming second-order accuracy.

The layered / variable-`D` behaviour and the mixed source boundary (no closed form) are covered
by the scripts in `examples/`.

## Expected results (`examples/example_run.py`)

A 100 m × 40 m, three-layer aquifer; a 10 m wide, 100 mg/L surface leak at x = 20 m for 50 days;
a well 70 m downstream at 5 m depth:

```
Layer interfaces at z = [10. 20.] m
Cell Peclet Pe = 0.300  (warn if > 2)
Cell Courant Cr = 0.250
Receptor (ix,iz) = (70, 5)
Peak concentration at receptor = 25.959 mg/L at t = 124.5 days
Final-field max concentration  = 18.419 mg/L
```

---

## References

- Course material (adapted): `class_notes_05_Parabolic_PDEs` and `Lab_04_Parabolic`
  (GPGN 436/536, Dr. Ge Jin & Dr. Jeff Shragge) — the Thomas algorithm and the 2D CN-ADI
  structure. See `AI_USE.md` for exactly what was adapted.
- Thomas algorithm (tridiagonal matrix algorithm):
  https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
- Peaceman, D.W. & Rachford, H.H. (1955), *The numerical solution of parabolic and elliptic
  differential equations*, J. SIAM 3(1), 28–41 (the ADI method).
- Pyodide: https://pyodide.org
