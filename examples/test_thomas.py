"""Step-1 test: verify the Thomas solvers.

1. thomas() vs numpy.linalg.solve() on random diagonally-dominant systems.
2. thomas_batch() vs a loop over thomas() (must agree to machine precision).

Run from the repo root:  python examples/test_thomas.py
"""

import numpy as np

from sim.thomas import thomas, thomas_batch


def build_dense(a, b, c):
    """Assemble the full tridiagonal matrix from its three diagonals."""
    n = len(b)
    A = np.zeros((n, n))
    for i in range(n):
        A[i, i] = b[i]
        if i > 0:
            A[i, i - 1] = a[i]
        if i < n - 1:
            A[i, i + 1] = c[i]
    return A


def random_system(n, rng):
    """Random tridiagonal system, made diagonally dominant so it is solvable."""
    a = rng.standard_normal(n)
    c = rng.standard_normal(n)
    b = rng.standard_normal(n) + (np.abs(a) + np.abs(c) + 1.0)  # . . dominant
    a[0] = 0.0       # . . unused
    c[-1] = 0.0      # . . unused
    d = rng.standard_normal(n)
    return a, b, c, d


def main():
    rng = np.random.default_rng(0)

    ## . . Test 1: single-line thomas vs dense solve
    max_err = 0.0
    for _ in range(100):
        n = int(rng.integers(3, 50))
        a, b, c, d = random_system(n, rng)
        x_thomas = thomas(a, b, c, d)
        x_dense = np.linalg.solve(build_dense(a, b, c), d)
        max_err = max(max_err, np.max(np.abs(x_thomas - x_dense)))
    print(f"[1] thomas vs numpy.linalg.solve : max error = {max_err:.2e}")
    assert max_err < 1e-9, "thomas disagrees with numpy.linalg.solve"

    ## . . Test 2: batched vs looped single-line
    nlines, n = 37, 25
    a = rng.standard_normal((nlines, n))
    c = rng.standard_normal((nlines, n))
    b = rng.standard_normal((nlines, n)) + (np.abs(a) + np.abs(c) + 1.0)
    a[:, 0] = 0.0
    c[:, -1] = 0.0
    d = rng.standard_normal((nlines, n))

    x_batch = thomas_batch(a, b, c, d)
    x_loop = np.array([thomas(a[k], b[k], c[k], d[k]) for k in range(nlines)])
    batch_err = np.max(np.abs(x_batch - x_loop))
    print(f"[2] thomas_batch vs looped thomas: max error = {batch_err:.2e}")
    assert batch_err < 1e-12, "thomas_batch disagrees with thomas"

    print("All Thomas tests passed.")


if __name__ == "__main__":
    main()
