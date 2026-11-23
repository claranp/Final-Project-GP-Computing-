# AI use statement

This document records how AI tools were used in this project and which course material was
adapted, in line with the course policy that I must be able to explain every line of the
simulation code myself.

## Tool used

- **Claude Code** (Anthropic), model **Claude Opus 4.8**, used as a pair-programming assistant
  from the terminal in this repository.

## How it was used

- **Deriving the scheme.** I had the assistant walk through the Crank–Nicolson ADI
  discretisation for the variable-`D` advection–dispersion equation — how the divergence-form
  dispersion term with harmonic-mean face values and the central-difference advection term enter
  the tridiagonal coefficients of the `x`- and `z`-sweeps — before any code was written. I
  checked this against the course derivation.
- **Writing the `sim/` package.** The assistant drafted `thomas.py`, `solver.py`, `model.py`,
  and `analytic.py` in the style of the course notebooks (same commenting conventions and
  variable names where sensible). I reviewed each function so I can explain the discretisation,
  the Thomas sweeps, the boundary handling (ghost-node elimination for the Neumann/zero-gradient
  boundaries, Dirichlet source patch), and the stability/accuracy behaviour.
- **Verification.** It helped build `verification.ipynb` and the `examples/` test scripts that
  compare the numerical solution to the analytic Gaussian and measure the Δx/Δt convergence
  slopes. The notebook imports the same `sim/` code the website runs — no separate solver.
- **Web interface and deployment.** It wrote the browser UI (`web/`), the Pyodide Web Worker,
  and the GitHub Actions Pages workflow. This front-end JavaScript/HTML/CSS is supporting
  infrastructure, not the numerical method being graded.
- **Debugging.** It helped fix issues found during review (e.g. computing the cell Péclet number
  per depth layer rather than mixing `min(D)` with `max(|v|)` across layers, colorbar label
  rendering, page layout).

## Course code adapted

The numerical core follows the course material directly:

- **Thomas algorithm.** `sim/thomas.py` is based on `TDMAsolver` from
  `class_notes_05_Parabolic_PDEs` (Dr. Jeff Shragge). I rewrote it to be non-destructive (it
  copies its inputs) and added `thomas_batch`, a NumPy-vectorised version that solves all ADI
  grid lines at once.
- **2D CN-ADI structure.** `sim/solver.py` generalises the instructor's
  `Solve_2D_Convection_Diffusion_homogeneous_Thomas` (same notebook). The course version assumes
  constant `D`, constant velocity, and Dirichlet boundaries on all four sides. I extended it to
  depth-dependent `D(z)` in divergence form with harmonic-mean faces, a per-layer horizontal
  velocity `vx(z)`, and this project's mixed boundary conditions.
- **Lab 4** (`Lab_04_Parabolic`, my own lab) was the reference for the overall solve-and-compare
  workflow and the Gaussian verification approach.

## External references

- Thomas / tridiagonal matrix algorithm — https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
- Peaceman–Rachford ADI method — Peaceman & Rachford (1955), J. SIAM 3(1), 28–41.
- Pyodide documentation — https://pyodide.org
