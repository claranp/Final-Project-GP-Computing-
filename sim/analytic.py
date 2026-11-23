"""Analytical solutions used only for verification (never a second solver).

Spreading/translating Gaussian for a uniform medium (constant D, constant vx),
on an infinite domain. This is the exact solution of

    dC/dt + vx dC/dx = D ( d2C/dx2 + d2C/dz2 )

for the Gaussian initial condition

    C(x, z, 0) = C0 exp( -[(x - x0)^2 + (z - z0)^2] / (2 sigma0^2) ).
"""

import numpy as np


def gaussian_spill(X, Z, t, C0, x0, z0, sigma0, D, vx):
    """Analytic concentration at coordinates (X, Z) and time t.

    input:
        X, Z  : coordinate arrays (e.g. from np.meshgrid, indexing='ij') [m]
        t     : time [days]
        C0    : initial peak concentration [mg/L]
        x0,z0 : initial plume centre [m]
        sigma0: initial Gaussian width [m]
        D     : dispersion coefficient [m^2/day]
        vx    : horizontal velocity [m/day]
    output:
        C : analytic concentration field, same shape as X and Z
    """
    var = sigma0 ** 2 + 2 * D * t                      # . . variance grows as 2 D t
    amp = C0 * sigma0 ** 2 / var                        # . . peak decays, mass conserved
    xc = x0 + vx * t                                    # . . centre advects with the flow
    return amp * np.exp(-((X - xc) ** 2 + (Z - z0) ** 2) / (2 * var))
