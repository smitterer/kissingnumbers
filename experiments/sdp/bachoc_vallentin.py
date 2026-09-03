"""Bachoc-Vallentin three-point SDP bound for the kissing number —
SAMPLED float prototype for machinery VALIDATION ONLY.

Formulation: Bachoc & Vallentin, "New upper bounds for kissing numbers
from semidefinite programming" (J. AMS 21 (2008), arXiv:math/0608426),
Theorem 4.2 with the simplification S_k^n(1,1,1) = 0 for k >= 1:

  A(n, theta) <= 1 + min  sum_{k=1..d} a_k + b11 + <F_0, S_0^n(1,1,1)>
  s.t.  [[b11, b12], [b12, b22]] >= 0 (psd),  a_k >= 0,  F_k >= 0,
        sum_k a_k P_k^n(u) + 2 b12 + b22 + 3 sum_k <F_k, S_k^n(u,u,1)> <= -1
                                                for all u in [-1, cos theta],
        b22 + sum_k <F_k, S_k^n(u,v,t)> <= 0   for all (u,v,t) in D',
  D' = {-1 <= u,v,t <= cos theta, 1 + 2uvt - u^2 - v^2 - t^2 >= 0}.

Matrices (their Thm 3.2 with the Remark 3.4 simplification, which
preserves the positivity property (15) and hence validity of the bound):
  Y_k^n[i,j](u,v,t) = u^i v^j Q_k^{n-1}(u,v,t),   0 <= i,j <= d-k,
  Q_k^{n-1}(u,v,t) = sum_{j==k mod 2} p_j (t-uv)^j ((1-u^2)(1-v^2))^{(k-j)/2},
      where p_j are the coefficients of the normalized Gegenbauer
      polynomial P_k^{n-1} (P(1) = 1),
  S_k^n = (1/6) sum over the 6 permutations of (u,v,t) applied to Y_k^n.
  S_0^n(1,1,1) = all-ones matrix in this basis.

RIGOR STATUS: the two polynomial constraints are enforced only on finite
sample grids, so the optimum reported here is in general a LOWER estimate
of the true SDP value and is NOT a valid kissing bound on its own. It is
used solely to validate the formulation against the published value
(n = 5, d = 10: Bachoc-Vallentin report tau_5 <= 45 at d = N = 10).
Certified computations (Track B rigor) would require the SOS/rounding
pipeline; per the project plan they run through ClusteredLowRankSolver.
"""

import argparse
import itertools
import time

import numpy as np


def gegenbauer_coeffs(m, k):
    """Coefficient list (low->high) of normalized Gegenbauer P_k^m,
    P(1) = 1, via the recurrence used in delsarte_lp.py, dimension m."""
    lam = (m - 2) / 2.0
    polys = [np.array([1.0]), np.array([0.0, 1.0])]
    for j in range(2, k + 1):
        a = np.zeros(j + 1)
        a[1:] += 2 * (j + lam - 1) * polys[j - 1]
        a[: j - 1] -= (j - 1) * polys[j - 2]
        polys.append(a / (j + 2 * lam - 1))
    return polys[k] if k > 0 else polys[0]


def q_poly_eval(n, k, U, V, T):
    """Q_k^{n-1}(u,v,t) evaluated on arrays."""
    p = gegenbauer_coeffs(n - 1, k)
    W = (1 - U ** 2) * (1 - V ** 2)
    S = T - U * V
    out = np.zeros(np.broadcast(U, V, T).shape)
    for j in range(k % 2, k + 1, 2):
        if abs(p[j]) > 0:
            out += p[j] * S ** j * W ** ((k - j) // 2)
    return out


def s_matrices(n, d, U, V, T):
    """S_k^n(u,v,t) for k=0..d on m sample triples: list of arrays
    (m, d-k+1, d-k+1)."""
    perms = list(itertools.permutations([U, V, T]))
    out = []
    for k in range(d + 1):
        sz = d - k + 1
        m = np.broadcast(U, V, T).size
        acc = np.zeros((m, sz, sz))
        for (a, b, c) in perms:
            q = q_poly_eval(n, k, a, b, c)
            pa = np.stack([a ** i for i in range(sz)])   # (sz, m)
            pb = np.stack([b ** j for j in range(sz)])
            acc += np.einsum("im,jm,m->mij", pa, pb, q)
        out.append(acc / 6.0)
    return out


def sample_domain(cos_theta, grid):
    """Sample triples in D' with u <= v <= t (S_k are permutation
    symmetric, so this loses nothing)."""
    xs = np.linspace(-1.0, cos_theta, grid)
    U, V, T = np.meshgrid(xs, xs, xs, indexing="ij")
    mask = (U <= V) & (V <= T) & \
        (1 + 2 * U * V * T - U ** 2 - V ** 2 - T ** 2 >= -1e-12)
    return U[mask], V[mask], T[mask]


def solve(n=5, d=10, cos_theta=0.5, ugrid=241, tgrid=28, solver="CLARABEL"):
    import cvxpy as cp

    t0 = time.time()
    # variables
    a = cp.Variable(d, nonneg=True)                 # a_1..a_d
    B = cp.Variable((2, 2), PSD=True)               # [[b11,b12],[b12,b22]]
    Fs = [cp.Variable((d - k + 1, d - k + 1), PSD=True) for k in range(d + 1)]

    cons = []
    # constraint 1 on u-grid
    us = np.linspace(-1.0, cos_theta, ugrid)
    Su = s_matrices(n, d, us, us, np.ones_like(us))
    # P_k^n(u) values
    Pk = np.stack([np.polynomial.polynomial.polyval(
        us, gegenbauer_coeffs(n, k)) for k in range(1, d + 1)])  # (d, m)
    expr1 = Pk.T @ a + 2 * B[0, 1] + B[1, 1]
    three = 0
    for k in range(d + 1):
        three = three + cp.vstack(
            [cp.sum(cp.multiply(Su[k][m], Fs[k])) for m in range(len(us))])
    cons.append(expr1 + 3 * cp.reshape(three, (len(us),), order="C") <= -1)

    # constraint 2 on D' sample
    U, V, T = sample_domain(cos_theta, tgrid)
    St = s_matrices(n, d, U, V, T)
    tot = 0
    for k in range(d + 1):
        tot = tot + cp.vstack(
            [cp.sum(cp.multiply(St[k][m], Fs[k])) for m in range(len(U))])
    cons.append(B[1, 1] + cp.reshape(tot, (len(U),), order="C") <= 0)

    # objective: 1 + sum a + b11 + <F_0, J>
    obj = 1 + cp.sum(a) + B[0, 0] + cp.sum(Fs[0])
    prob = cp.Problem(cp.Minimize(obj), cons)
    prob.solve(solver=solver, verbose=False)
    return prob.value, prob.status, len(us), len(U), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--d", type=int, default=10)
    ap.add_argument("--ugrid", type=int, default=241)
    ap.add_argument("--tgrid", type=int, default=28)
    args = ap.parse_args()
    val, status, nu, nt, dt = solve(args.n, args.d, ugrid=args.ugrid,
                                    tgrid=args.tgrid)
    print(f"BV three-point bound (SAMPLED, float, NOT certified): "
          f"n={args.n} d={args.d}")
    print(f"  status={status}, samples: {nu} (u), {nt} (u,v,t)")
    print(f"  value = {val:.6f}   [published: d=N=10 gives tau_5 <= 45; "
          f"sampled value may dip below the true SDP optimum]")
    print(f"  time = {dt:.1f}s")


if __name__ == "__main__":
    main()
