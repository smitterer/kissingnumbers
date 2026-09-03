"""Validation experiment: reproduce Szollosi's Method-2 run exactly.

Expected (arXiv:2301.08272, Sec. 4): with d = 4, A = {-1, -1/2, 0, 1/2}
(unit normalization), basis B = rows of (1/sqrt(50))*[[5,5,0,0,0],
[5,0,5,0,0],[5,0,0,5,0],[5,0,0,0,5]] and B5 = (1/sqrt(5))*[-1,1,1,1,1]:
the cloud C'_{A,B} contains 78 vectors, the compatibility graph has
exactly four maximum cliques of size omega = 36, two of which extend B to
a 40-point arrangement with the D5 profile and two with the Q5 profile.

Here everything is scaled to squared norm r2 = 50 (angle set scales to
{-50, -25, 0, 25}); the normal is taken as (-1,1,1,1,1) with mu = 5.

Deterministic (no seed). Output: cloud size, clique number, number of
maximum cliques, profile classification of each resulting 40-point
configuration, and checker-format configs written to
experiments/szollosi_repro/ for independent verification.
"""

import itertools
import json
import os
import sys
import time
from collections import Counter
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clique import (Cand, build_graph, check_basis_compatible, clique_to_config,
                    cloud_method2, compatible, dot, max_cliques)

R2 = Fraction(50)
BASIS = [[Fraction(c) for c in row] for row in
         [[5, 5, 0, 0, 0], [5, 0, 5, 0, 0], [5, 0, 0, 5, 0], [5, 0, 0, 0, 5]]]
NORMAL = [Fraction(c) for c in [-1, 1, 1, 1, 1]]
ANGLES = [Fraction(-50), Fraction(-25), Fraction(0), Fraction(25)]

# Unit-normalized profiles of D5 and Q5 (unordered pairs), from
# Cohn-Rajagopal arXiv:2412.00937 Table 2.1.
PROFILES = {
    "D5": {Fraction(-1): 20, Fraction(-1, 2): 240, Fraction(0): 280,
           Fraction(1, 2): 240},
    "Q5": {Fraction(-1): 10, Fraction(-4, 5): 30, Fraction(-1, 2): 180,
           Fraction(-3, 10): 60, Fraction(0): 250, Fraction(1, 5): 10,
           Fraction(1, 2): 240},
}


def full_ip(u, v):
    """Exact inner product of two candidates as (rational, eps, disc-product);
    only usable when the surd part vanishes or disc product is a square."""
    q = dot(u.rat, v.rat)
    e = u.eps * v.eps
    D = u.disc * v.disc
    if e == 0 or D == 0:
        return q
    # exact square root if D is a square of a rational
    num, den = D.numerator, D.denominator
    rn, rd = _isqrt(num), _isqrt(den)
    if rn is not None and rd is not None:
        return q + e * Fraction(rn, rd)
    raise ValueError(f"irrational inner product: {q} + {e}*sqrt({D})")


def _isqrt(n):
    if n < 0:
        return None
    r = int(n ** 0.5)
    for cand in (r - 1, r, r + 1):
        if cand >= 0 and cand * cand == n:
            return cand
    return None


def profile_of(vectors_ip, r2):
    """vectors_ip: function (i, j) -> exact rational ip; n implied 40."""
    c = Counter()
    for i, j in itertools.combinations(range(40), 2):
        c[vectors_ip(i, j) / r2] += 1
    return dict(c)


def main():
    t0 = time.time()
    check_basis_compatible(BASIS, R2)
    cloud = cloud_method2(BASIS, NORMAL, R2, ANGLES)
    print(f"cloud size: {len(cloud)}  (Szollosi: 78)")

    adj = build_graph(cloud, R2)
    omega, cliques = max_cliques(adj, want_all=True)
    print(f"max clique: {omega}  (Szollosi: 36); "
          f"#maximum cliques: {len(cliques)}  (Szollosi: 4)")

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "experiments", "szollosi_repro")
    os.makedirs(outdir, exist_ok=True)
    classification = Counter()
    for ci, clique in enumerate(cliques):
        pts = list(BASIS) + [cloud[k] for k in clique]

        def ip(i, j, pts=pts):
            a, b = pts[i], pts[j]
            if isinstance(a, list) and isinstance(b, list):
                return dot(a, b)
            if isinstance(a, list):
                return dot(a, b.rat)  # basis orthogonal to normal
            if isinstance(b, list):
                return dot(a.rat, b)
            return full_ip(a, b)

        prof = profile_of(ip, R2)
        label = next((k for k, v in PROFILES.items() if v == prof), "OTHER")
        classification[label] += 1
        cfg = clique_to_config(f"szollosi_repro_{ci}_{label}", BASIS, cloud,
                               clique, R2, normal=NORMAL)
        path = os.path.join(outdir, f"clique_{ci}_{label}.json")
        with open(path, "w") as f:
            json.dump(cfg, f, indent=1)
        print(f"  clique {ci}: 4+{len(clique)} = {4 + len(clique)} points, "
              f"profile = {label} -> {path}")
        if label == "OTHER":
            print(f"    !! unexpected profile: "
                  f"{sorted((str(k), v) for k, v in prof.items())}")

    print(f"classification: {dict(classification)}  "
          f"(Szollosi: 2x D5-profile, 2x Q5)")
    print(f"runtime: {time.time() - t0:.2f}s")


if __name__ == "__main__":
    main()
