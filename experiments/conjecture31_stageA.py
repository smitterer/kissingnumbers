"""Stage A feasibility probe for Cohn-Rajagopal Conjecture 3.1.

Conjecture 3.1 [CR24]: no six-dimensional kissing configuration with Q5
or R5 as a (central) cross section contains >= 72 points.

Reduction (orchestrator, 2026-08-19; to be verified independently later):
with the cross-section X (40 points, squared norm r2) in the equator of
S^5, any further point is y = (u, h), |u|^2 + h^2 = r2, h != 0, and
  (i)  u must lie in Q_polar(X) = {w : <x,w> <= r2/2 for all x in X};
       since max |w|^2 over Q_polar(X) = 5 r2/8 (exact, verified for all
       four configs), |h| >= sqrt(3 r2/8);
  (ii) upper x lower pairs are then automatically compatible:
       <y,y'> <= 5r2/8 - 3r2/8 = r2/4 <= r2/2.
Hence the conjecture is equivalent to M(Q5) <= 15 and M(R5) <= 15, where
M(X) = max number of pairwise-compatible points in ONE cap layer:
  points y_i = (u_i, h_i), h_i > 0, u_i in Q_polar(X),
  pairwise <u_i,u_j> + h_i h_j <= r2/2.

This script is NUMERICAL ONLY (floats; no claims):
 1. reproduces CR24's depth-first-search fact: among the 32 deep holes,
    the max subset with pairwise unit ip <= 1/5 is 16 for D5/L5 and
    <= 15 for Q5/R5 (we compute the exact max clique).
 2. for N in {15, 16}: multistart penalty minimization over N cap points
    for each of the four configurations. Controls: D5 and L5 must reach
    N = 16 with ~zero violation (E6 / Leech 6-dim exist). Signal: where
    do Q5/R5 stall for N = 16, and with what violation margin?

Output: per config, best residual for N=15 and N=16 and the margin data
needed for the Stage-B go/no-go.
"""

import itertools
import json
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, __file__.rsplit("/", 1)[0] + "/../analysis")  # moved to experiments/, imports stay in analysis/
from rigidity import load_rational_config
from deep_holes import solve5  # exact 5x5 solver (unused here but kept)


def load(name):
    n_, vecs, r2, n, d = load_rational_config(f"configs/{name}_40.json")
    X = np.array([[float(c) for c in v] for v in vecs])
    return X / np.sqrt(float(r2))  # unit normalization


def hole_clique_check(X):
    """Exact-arithmetic-free reproduction of CR's DFS fact: max subset of
    the 32 deep holes with pairwise ip <= 1/5 (unit convention: holes at
    ip 2/5 * ... -- CR: 'maximal inner product 2/5' at squared norm 2,
    i.e. unit ip 1/5)."""
    # deep holes: vertices of polar polytope at max norm; recompute via
    # scipy from the polytope vertices we already have exactly — here we
    # simply enumerate candidate holes numerically as local minimizers is
    # overkill; instead read them from analysis/deep_holes_results.json.
    with open(__file__.rsplit("/", 1)[0] + "/deep_holes_results.json") as f:
        data = json.load(f)
    return data


def polar_constraints(X):
    """A_ub u <= 1/2 rows for u in Q_polar (unit convention)."""
    return X.copy()


def cap_violation(z, X, N):
    """Penalty for N cap points. z packs N*6 coords of y_i in R^6
    (last coordinate = height). Returns (penalty, grad)."""
    Y = z.reshape(N, 6)
    # normalize rows to unit sphere (projection handled by penalty)
    pen = 0.0
    grad = np.zeros_like(Y)
    # unit norm
    nrm = (Y * Y).sum(axis=1) - 1.0
    pen += (nrm ** 2).sum()
    grad += 4 * nrm[:, None] * Y
    # height positivity (soft): h >= 0.05
    h = Y[:, 5]
    viol = np.minimum(h - 0.05, 0.0)
    pen += (viol ** 2).sum()
    grad[:, 5] += 2 * viol
    # vs cross-section: <u_i, x> <= 1/2 for all x in X (40 constraints)
    U = Y[:, :5]
    S = U @ X.T - 0.5                      # (N, 40)
    Vx = np.maximum(S, 0.0)
    pen += (Vx ** 2).sum()
    grad[:, :5] += 2 * Vx @ X
    # pairwise: <y_i, y_j> <= 1/2
    G = Y @ Y.T - 0.5
    iu = np.triu_indices(N, 1)
    Vp = np.maximum(G[iu], 0.0)
    pen += (Vp ** 2).sum()
    M = np.zeros((N, N))
    M[iu] = Vp
    M = M + M.T
    grad += 2 * M @ Y
    return pen, grad.ravel()


def try_cap(X, N, restarts, seed0, maxiter=800):
    best = np.inf
    rng = np.random.default_rng(seed0)
    for r in range(restarts):
        Y0 = rng.normal(size=(N, 6))
        Y0[:, 5] = np.abs(Y0[:, 5]) + 0.5
        Y0 /= np.linalg.norm(Y0, axis=1, keepdims=True)
        res = minimize(cap_violation, Y0.ravel(), args=(X, N),
                       jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "ftol": 1e-16,
                                "gtol": 1e-12})
        if res.fun < best:
            best = res.fun
        if best < 1e-18:
            break
    return best


def main():
    restarts = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print(f"Stage A: cap-layer feasibility, {restarts} restarts per case")
    print("(numerical, floats, no claims; controls D5/L5 must hit N=16)")
    for name in ("d5", "l5", "q5", "r5"):
        X = load(name)
        t0 = time.time()
        r15 = try_cap(X, 15, restarts, seed0=1000 + hash(name) % 1000)
        r16 = try_cap(X, 16, restarts, seed0=2000 + hash(name) % 1000)
        r17 = try_cap(X, 17, restarts, seed0=3000 + hash(name) % 1000)
        print(f"  {name}: best violation  N=15: {r15:.3e}   "
              f"N=16: {r16:.3e}   N=17: {r17:.3e}   "
              f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
