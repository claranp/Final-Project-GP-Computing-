# CLAUDE.md — Aquifer Contaminant Transport (GPGN 436/536 Final Project)

Student: Clara Pugh. Deadline: repo + live GitHub Pages site by **Sun Nov 29, 2026, 11:59 p.m. MT**.
Presentation: Mon Nov 30, 10:36 a.m. (4 min demo, 4 min code + verification walkthrough, 2 min questions answered WITHOUT AI).

The student must be able to explain every line of the simulation code. Favor clear, commented,
course-style code over clever code. When you write or change solver code, explain what you did
and why in plain terms (discretization, stencil, boundary handling, stability/accuracy).
Write code in the same style as the course code in class_notes_05_Parabolic_PDEs copy.ipynb (instructor-created material) and Lab_04_Parabolic copy.ipynb (my lab, could have errors)

Always ask before making assumptions about what the student wants

## Geophysical problem
A contaminant source at the ground surface (landfill, chemical plant) leaks into a layered aquifer.
Model how the plume moves (advection by groundwater flow) and spreads (dispersion) through a
2D vertical cross-section (x horizontal, z depth) toward downstream receptors (houses/wells).

## Governing PDE (advection–dispersion / convection–diffusion)
    ∂C/∂t + v·∇C = ∇·(D ∇C)
- C(x, z, t): concentration [mg/L]
- v = (vx(z), 0): horizontal groundwater velocity per layer [m/day]
- D(z): dispersion coefficient per layer [m²/day] (spatially variable → keep the divergence form)

## Numerical method (must be visible in our own code — no black-box PDE solvers)
- **Crank–Nicolson ADI (Peaceman–Rachford)** on a regular grid:
  - half step 1: implicit in x, explicit in z
  - half step 2: implicit in z, explicit in x
- Every implicit sweep is a tridiagonal system → solve with a **hand-written Thomas algorithm**
  (vectorized over all lines at once with NumPy for speed in Pyodide; also keep a simple
  single-line version that the notebook tests against `numpy.linalg.solve`).
- Advection: central differences in space.
- Variable D: use face values (harmonic mean of neighboring cells) in the
  ∂/∂x(D ∂C/∂x) and ∂/∂z(D ∂C/∂z) stencils.
- Report the diagnostics the user needs to judge the run: cell Péclet number Pe = |v|Δx/D
  (central differences oscillate when Pe > 2) and Courant number Cr = |v|Δt/Δx
  (CN is stable for any Δt but loses accuracy and can oscillate at sharp fronts when Cr is large).

## Boundary and initial conditions
- Top (z = 0): source as a fixed-concentration (Dirichlet) patch of user-set width and x-location,
  either continuous or for a user-set leak duration; no-flux elsewhere on top.
- Bottom: no-flux (impermeable base).
- Left (upstream): C = 0 (clean inflow).
- Right (downstream): zero-gradient outflow.
- Initial condition: C = 0 everywhere (option: initial Gaussian spill for verification).

## Inputs (all with units)
Domain length/depth [m], Δx, Δz [m], Δt [days], total time [days], source x-location [m],
source width [m], source concentration [mg/L], leak duration [days],
number of layers + per-layer table (thickness [m], seepage velocity [m/day], dispersion [m²/day]),
optionally loaded from CSV.
Button to toggle between inputting seepage velocity vs inputting permeability, gradient, and porosity. If used, compute
velocity from Darcy's law v = K·i/n (hydraulic conductivity K [m/day], gradient i [-],
porosity n [-]) instead of also asking for velocity, so the inputs are not redundant.

## Outputs
- Concentration map C(x, z) at the selected time with colorbar [mg/L], axes in m, layer boundaries drawn
- Breakthrough curve C(t) at a receptor point (house/well location) [mg/L vs days]
- Pe and Cr diagnostics shown in the interface, with a warning when Pe > 2

## Repository layout
    sim/                    # standalone simulation package — NO UI/JS/Pyodide code here
      __init__.py
      model.py              # parameters, grid, layers, D/v fields, BCs, ICs
      thomas.py             # Thomas algorithm (single + batched)
      solver.py             # CN-ADI time stepping
      analytic.py           # analytical solutions used for verification
    web/
      index.html, style.css, app.js   # controls + plotting (JS)
      worker.js             # loads Pyodide + numpy, imports sim/, runs the model, posts results
    verification.ipynb      # imports sim/ — never a separate solver; saved with outputs
    examples/example_run.py # working example + expected results (also in README)
    requirements.txt        # numpy, matplotlib, jupyter (pinned versions)
    .github/workflows/deploy.yml   # GitHub Actions → GitHub Pages (copies web/ + sim/ into site)
    README.md               # problem, PDE, method, how to run locally, live link, expected results
    AI_USE.md               # AI tools used and how; references and adapted course code

## Rules
- `sim/` must be importable and runnable with plain Python/NumPy (no Pyodide-specific code).
- Pyodide runs in a Web Worker so the page stays responsive; pin the Pyodide version and list
  loaded packages in README. Keep dependencies to numpy only in the browser.
- The verification notebook imports the SAME `sim/` code the website runs.
- Never write a second solver just for verification.
- Every plot: labeled axes with units, title, colorbar units.
- Commit after each working step with clear messages.

## Verification plan (verification.ipynb)
Uniform medium, single layer, initial Gaussian spill away from boundaries:
    C(x,z,t) = C0 · σ0²/(σ0² + 2Dt) · exp(-[(x - x0 - v t)² + (z - z0)²] / (2(σ0² + 2Dt)))
1. v = 0 (pure diffusion — like Lab 4 part 1, parabolic PDE)
2. v ≠ 0 (tests the advection term)
For each: plots of numerical vs analytical, error map, L2 and max error;
convergence vs Δx (fixed small Δt) and vs Δt (fixed fine grid) on log-log axes with
slope-2 reference lines; short discussion of agreement, Pe/Cr effects, boundary effects, limitations.

## Build order
1. sim/ solver + Thomas, tested in a script
2. verification notebook (convergence plots)
3. minimal Pyodide worker page deployed via GitHub Actions (prove deployment early)
4. full interface (layers table, CSV, plots, diagnostics)
5. README, AI_USE.md, example, fresh-clone test
