"""Delsarte-Goethals-Seidel LP bound for spherical codes: validation step.

For A(n, theta) with cos(theta) = 1/2 (kissing number), the LP bound is:
  minimize 1 + sum_{k>=1} f_k   over f_k >= 0,
  s.t.  F(u) := 1 + sum_{k>=1} f_k G_k^n(u) <= 0 for all u in [-1, 1/2],
where G_k^n are Gegenbauer polynomials normalized to G_k^n(1) = 1.
Then |C| <= F(1) = 1 + sum f_k.  (Standard form; see e.g. Bachoc-Vallentin
arXiv:math/0608426 Sec. 2, or the Boyvalenkov et al. survey 1507.03631.)

Validation target (dim 5): Odlyzko-Sloane 1979 report 46.345 (survey
value tau_5 <= 46.345); Levenshtein bound 48.

Implementation: discretize [-1, 1/2] finely, solve the LP with scipy
linprog (floats), then do an a-posteriori RIGOROUS-ISH check: evaluate
F on a much finer grid and bound F' to control between-grid excursions.
This validation step is float-based and is labeled as such; it is NOT a
certified bound (Track B rigor comes later, and only matters if we ever
claim a new bound).

Gegenbauer recurrence (normalized, G_k(1)=1), lambda = (n-2)/2:
  G_0 = 1, G_1(u) = u,
  G_k(u) = (2u(k+lambda-1) G_{k-1}(u) - (k-1) G_{k-2}(u)) / (k + 2*lambda - 1)
  -- this is the ultraspherical recurrence rewritten for the
  normalization C_k^lambda(u)/C_k^lambda(1).
"""

import argparse

import numpy as np
from scipy.optimize import linprog


def gegenbauer_normalized(n, kmax, u):
    """G_k^n(u) for k = 0..kmax, normalized G_k(1) = 1; u array."""
    lam = (n - 2) / 2.0
    G = np.zeros((kmax + 1,) + u.shape)
    G[0] = 1.0
    if kmax >= 1:
        G[1] = u
    for k in range(2, kmax + 1):
        G[k] = (2 * u * (k + lam - 1) * G[k - 1]
                - (k - 1) * G[k - 2]) / (k + 2 * lam - 1)
    return G


def delsarte_lp(n, degree, cos_theta=0.5, grid=4001):
    u = np.linspace(-1.0, cos_theta, grid)
    G = gegenbauer_normalized(n, degree, u)
    # variables f_1..f_degree; constraint: sum f_k G_k(u) <= -1
    A_ub = G[1:].T                      # (grid, degree)
    b_ub = -np.ones(grid)
    c = np.ones(degree)                 # minimize sum f_k  (bound = 1 + sum)
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(0, None)] * degree,
                  method="highs")
    assert res.success, res.message
    f = res.x
    bound = 1.0 + f.sum()
    # a-posteriori check on a 100x finer grid
    uf = np.linspace(-1.0, cos_theta, (grid - 1) * 100 + 1)
    Gf = gegenbauer_normalized(n, degree, uf)
    F = 1.0 + f @ Gf[1:]
    worst = F.max()
    return bound, f, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--degrees", default="10,12,16,20,24,30")
    args = ap.parse_args()
    print(f"Delsarte LP bound, dimension {args.n}, angle 60deg "
          f"(float validation, NOT certified)")
    for d in (int(x) for x in args.degrees.split(",")):
        bound, f, worst = delsarte_lp(args.n, d)
        print(f"  degree {d:3d}: bound = {bound:.6f}   "
              f"max F on fine grid = {worst:.2e} (should be ~<= 0)")


if __name__ == "__main__":
    main()
