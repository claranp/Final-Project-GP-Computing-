"""Thomas algorithm (tridiagonal matrix solver).

Solves a tridiagonal linear system

    a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i,     a_0 and c_{n-1} unused,

which is the form that every implicit Crank-Nicolson sweep reduces to (see
class_notes_05 "Solution by Thomas Algorithm"). The algorithm is Gaussian
elimination specialised to a tridiagonal matrix, so it costs O(n) instead of
O(n^3) for a dense solve.

Two versions:
    thomas(a, b, c, d)        - one line at a time, easy to read / explain.
    thomas_batch(a, b, c, d)  - all lines at once, vectorised over NumPy for
                                speed in the browser (Pyodide).
Both return a NEW array and do not modify their inputs.
"""

import numpy as np


def thomas(a, b, c, d):
    """Solve a single tridiagonal system with the Thomas algorithm.

    input:
        a : sub-diagonal   coefficients (length n); a[0]   is not used
        b : main-diagonal  coefficients (length n)
        c : super-diagonal coefficients (length n); c[n-1] is not used
        d : right-hand side             (length n)
    output:
        x : solution vector (length n)
    """
    n = len(d)

    ## . . Work arrays for the modified coefficients (keep inputs intact)
    cp = np.zeros(n)   # . . modified super-diagonal c'
    dp = np.zeros(n)   # . . modified right-hand side d'
    x = np.zeros(n)

    ## . . Forward sweep: eliminate the sub-diagonal row by row
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]      # . . pivot after elimination
        cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m

    ## . . Back substitution: last unknown is known directly, then work up
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]

    return x


def thomas_batch(a, b, c, d):
    """Solve many tridiagonal systems at once (one per row).

    Same algorithm as thomas(), but the inner recurrences are vectorised over
    the first axis so every line is advanced together. This is what we use in
    the ADI sweeps, where each grid line is an independent tridiagonal system.

    input:
        a, b, c, d : arrays of shape (nlines, n). Column 0 of a and column
                     n-1 of c are not used. Coefficients may vary along each
                     line (needed for the variable-D z-sweep).
    output:
        x : solutions, shape (nlines, n)
    """
    nlines, n = d.shape

    cp = np.zeros((nlines, n))
    dp = np.zeros((nlines, n))
    x = np.zeros((nlines, n))

    ## . . Forward sweep (loop over the within-line index, vectorised over lines)
    cp[:, 0] = c[:, 0] / b[:, 0]
    dp[:, 0] = d[:, 0] / b[:, 0]
    for i in range(1, n):
        m = b[:, i] - a[:, i] * cp[:, i - 1]
        cp[:, i] = c[:, i] / m
        dp[:, i] = (d[:, i] - a[:, i] * dp[:, i - 1]) / m

    ## . . Back substitution
    x[:, -1] = dp[:, -1]
    for i in range(n - 2, -1, -1):
        x[:, i] = dp[:, i] - cp[:, i] * x[:, i + 1]

    return x
