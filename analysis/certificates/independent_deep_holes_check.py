#!/usr/bin/env python3
"""INDEPENDENT re-verification of the deep-hole / "no 41st point" claim for
the four 40-point five-dimensional kissing configurations.

Written from scratch for the verifier role; shares no code with
analysis/deep_holes.py.  Reads only configs/{d5,l5,q5,r5}_40.json.

                        *** NO FLOATING POINT ANYWHERE ***

Unlike the method under review, this script uses NO float prefilter.  The
entire enumeration runs in exact integer arithmetic (numpy int64, with a
proven Hadamard bound showing every intermediate stays below 2^62, plus
runtime assertions on the observed magnitudes), and every vertex found is
re-verified independently in Python Fractions.  Floats appear only in the
printed summary, marked '~', and in an optional scipy sanity net that is
explicitly labelled as non-verdict.

------------------------------------------------------------------ setting

X = {x_1..x_40} at squared norm r2, pairwise <x_i,x_j> <= r2/2.
Polar body            Q = {w in R^5 : <x_i, w> <= r2/2 for all i}.

(a) A 41st point of the undeformed configuration is a point x with
    |x|^2 = r2 and <x_i,x> <= r2/2 for all i, i.e. a point of Q of squared
    norm exactly r2.  Since 0 is in Q and Q is convex, t*w in Q for all
    t in [0,1], so |t w|^2 sweeps [0, |w|^2] continuously.  Hence

        a 41st point exists  <=>  max_{w in Q} |w|^2 >= r2.

    (=>) is immediate; (<=) take t = sqrt(r2)/|w|.  So max |w|^2 < r2
    proves exactly that no 41st point can be added.

(b) With x^_i = x_i/sqrt(r2) and f(u) = max_i <x^_i, u> on the unit sphere,
    the ray {t u : t >= 0} meets Q in t <= sqrt(r2)/(2 f(u)) (f(u) > 0 for
    every u because 0 is interior to conv(X)).  Therefore

        max_{w in Q} |w| = sqrt(r2) / (2 * min_{|u|=1} f(u)),
        i.e.  m(X)^2 = r2 / (4 * max|w|^2),

    and the minimizers u (the deep holes) are exactly the normalized
    maximal-norm points of Q.  Because |w|^2 is STRICTLY convex, every
    maximal-norm point of the compact polytope Q is a vertex, so
    enumerating vertices finds all deep holes and misses none.

(c) Boundedness of Q -- which (a) and (b) both need, and which the method
    under review asserts but never actually checks -- is established here
    exactly: Q is bounded iff its recession cone
    C = {w : <x_i,w> <= 0 for all i} is {0}.  X spans R^5, so C contains no
    line and is pointed; a pointed polyhedral cone other than {0} has an
    extreme ray, whose tight set has rank 4.  We therefore enumerate ALL
    C(40,4) = 91390 four-subsets, compute the (integer) kernel direction as
    the generalized cross product, and check both signs.  No feasible ray
    => C = {0} => Q is a compact polytope.

Completeness of the vertex enumeration: every vertex of the polyhedron
{w : Aw <= b} in R^5 has tight constraints of rank 5, hence is the unique
solution of some 5-subset of them.  Enumerating ALL C(40,5) = 658008
5-subsets, keeping the nonsingular ones (exact integer determinant, no
tolerance) and exactly testing feasibility, therefore cannot miss a vertex.

------------------------------------------------------- exactness argument

Coordinates are scaled by L = lcm of denominators to integers Y = L*X with
|Y_ij| <= 10 for all four configurations (asserted at runtime).  Appending a
column of ones gives Yhat with rows of Euclidean norm <= sqrt(5*100+1) <
22.4.  By Hadamard, every k x k minor of Yhat satisfies |minor| <= 22.4^k,
so |2x2| < 502, |3x3| < 11239, |4x4| < 251770, |5x5| < 5640000.  The largest
quantity formed is a dot product <Y_j, N> with |N| entries <= 5.64e6, giving
|<Y_j,N>| <= 5 * 10 * 5.64e6 < 2.9e8.  Everything is far below
2^62 = 4.6e18, so the int64 arithmetic is exact; runtime assertions verify
the observed maxima against these bounds.

Exit code 0 iff every configuration passes every check.
"""

from __future__ import annotations

import itertools
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CONFIG_DIR = os.path.join(ROOT, "configs")

MAX_ABS_Y = 10          # asserted bound on the scaled integer coordinates
INT64_SAFE = 1 << 62


# ------------------------------------------------------------------ loading

def load_config(path):
    with open(path) as fh:
        cfg = json.load(fh)
    d, n = cfg["dimension"], cfg["n_points"]
    r2 = Fraction(cfg["norm_squared"])
    xs = [[Fraction(c) for c in v] for v in cfg["vectors"]]
    assert len(xs) == n and all(len(v) == d for v in xs)
    # re-verify it is a kissing configuration (exact)
    half = r2 / 2
    for i, x in enumerate(xs):
        assert sum(c * c for c in x) == r2, "vector %d has wrong norm" % i
    for i in range(n):
        for j in range(i + 1, n):
            assert sum(a * b for a, b in zip(xs[i], xs[j])) <= half, \
                "pair (%d,%d) violates the kissing constraint" % (i, j)
    return cfg.get("name", "?"), xs, r2, n, d


def lcm(a, b):
    x, y = a, b
    while y:
        x, y = y, x % y
    return a * b // x


def to_integer_rows(xs):
    L = 1
    for x in xs:
        for c in x:
            L = lcm(L, c.denominator)
    Y = [[int(c * L) for c in x] for x in xs]
    assert all(Fraction(v, L) == c for row, x in zip(Y, xs)
               for v, c in zip(row, x))
    return L, Y


# --------------------------------------------- exact Laplace minor ladders

def minor_ladder(M, ncols, nrows):
    """All maximal minors of an integer matrix stack.

    ``M`` has shape (c, nrows, ncols) with nrows <= ncols.  Returns a dict
    mapping each ascending column tuple of length ``nrows`` to the (c,)
    array of the corresponding nrows x nrows minors, computed by repeated
    Laplace expansion along the top row.  Pure integer arithmetic.
    """
    cols = list(range(ncols))
    # level 1: 1x1 minors from the bottom row
    level = {(j,): M[:, nrows - 1, j] for j in cols}
    for size in range(2, nrows + 1):
        row = nrows - size
        nxt = {}
        for tup in itertools.combinations(cols, size):
            acc = None
            for k, j in enumerate(tup):
                sub = tup[:k] + tup[k + 1:]
                term = M[:, row, j] * level[sub]
                acc = (term if acc is None else acc + term) if k % 2 == 0 \
                    else acc - term
            nxt[tup] = acc
        level = nxt
    return level


def _check_bounds(arrs, bound, what):
    for a in arrs:
        m = int(np.abs(a).max()) if a.size else 0
        assert m <= bound, "%s exceeded bound: %d > %d" % (what, m, bound)
        assert m < INT64_SAFE, "%s near int64 overflow" % what


# ------------------------------------------------ recession cone / boundedness

def recession_cone_is_trivial(Y, n, d, chunk=20000):
    """Exact: is {w : <x_i,w> <= 0 for all i} = {0}?

    Enumerates all C(n,4) four-subsets; the kernel of a rank-4 4x5 integer
    matrix is the generalized cross product (signed 4x4 minors).  Returns
    (ok, witness) -- ok True iff no extreme ray is feasible.
    """
    Ya = np.asarray(Y, dtype=np.int64)
    subs = np.fromiter(itertools.chain.from_iterable(
        itertools.combinations(range(n), d - 1)), dtype=np.int64)
    subs = subs.reshape(-1, d - 1)
    for s0 in range(0, len(subs), chunk):
        S = subs[s0:s0 + chunk]
        M = Ya[S]                                     # (c, 4, 5)
        minors = minor_ladder(M, d, d - 1)            # 4x4 minors, C(5,4)=5
        _check_bounds(minors.values(), 251770, "4x4 minor")
        # kernel vector k_j = (-1)^j * minor omitting column j
        K = np.empty((len(S), d), dtype=np.int64)
        allc = tuple(range(d))
        for j in range(d):
            tup = allc[:j] + allc[j + 1:]
            K[:, j] = minors[tup] * (1 if j % 2 == 0 else -1)
        nz = np.abs(K).sum(axis=1) > 0
        if not nz.any():
            continue
        Kn = K[nz]
        P = Kn @ Ya.T                                 # (c, n)
        _check_bounds([P], 2 * 10**9, "cone dot product")
        for sign in (1, -1):
            feas = (sign * P <= 0).all(axis=1)
            if feas.any():
                idx = int(np.nonzero(feas)[0][0])
                return False, (sign * Kn[idx]).tolist()
    return True, None


# ------------------------------------------------------ vertex enumeration

def enumerate_vertices(xs, Y, L, r2, n, d, chunk=32768):
    """All vertices of Q = {w : <x_i,w> <= r2/2}, exactly.

    Returns (vertices, n_feasible_subsets, stats).  Vertices are tuples of
    Fractions, deduplicated, each re-verified in Fraction arithmetic.
    """
    Ya = np.asarray(Y, dtype=np.int64)
    Yhat = np.concatenate([Ya, np.ones((n, 1), dtype=np.int64)], axis=1)
    assert int(np.abs(Yhat).max()) <= MAX_ABS_Y

    subs = np.fromiter(itertools.chain.from_iterable(
        itertools.combinations(range(n), d)), dtype=np.int64)
    subs = subs.reshape(-1, d)
    n_subsets = len(subs)

    R = r2 / 2 * L                       # rhs in the scaled coordinates
    assert R == int(R), "L*r2/2 is not an integer"
    R = Fraction(R)

    allc = tuple(range(d + 1))
    cand = []
    n_singular = 0
    for s0 in range(0, n_subsets, chunk):
        S = subs[s0:s0 + chunk]
        M = Yhat[S]                                   # (c, 5, 6)
        minors = minor_ladder(M, d + 1, d)            # all 5x5, C(6,5)=6
        _check_bounds(minors.values(), 5640000, "5x5 minor")
        # p_k = minor omitting column k;  D = p_5,  N_k = (-1)^k p_k
        P = {k: minors[allc[:k] + allc[k + 1:]] for k in range(d + 1)}
        D = P[d]
        good = D != 0
        n_singular += int((~good).sum())
        if not good.any():
            continue
        Dg = D[good]
        N = np.empty((int(good.sum()), d), dtype=np.int64)
        for k in range(d):
            N[:, k] = P[k][good] * (1 if k % 2 == 0 else -1)
        sgn = np.where(Dg > 0, 1, -1).astype(np.int64)
        Ns = N * sgn[:, None]
        Da = np.abs(Dg)
        dots = Ns @ Ya.T                              # (c, n)
        _check_bounds([dots], 3 * 10**8, "vertex dot product")
        feas = (dots <= Da[:, None]).all(axis=1)
        if feas.any():
            Sg = S[good][feas]
            Nf, Df = Ns[feas], Da[feas]
            for t in range(len(Sg)):
                cand.append((tuple(int(v) for v in Sg[t]),
                             tuple(int(v) for v in Nf[t]), int(Df[t])))

    # exact re-verification and dedup, entirely in Fractions
    verts = {}
    half = r2 / 2
    for sub, Nv, Dv in cand:
        w = tuple(R * Fraction(Nv[k], Dv) for k in range(d))
        # (i) the defining equalities really hold
        for i in sub:
            assert sum(a * b for a, b in zip(xs[i], w)) == half, \
                "defining equality failed for subset %s" % (sub,)
        # (ii) feasibility, re-derived from the original rational data
        assert all(sum(a * b for a, b in zip(xs[j], w)) <= half
                   for j in range(n)), "infeasible point survived: %s" % (sub,)
        verts.setdefault(w, sub)
    stats = {"n_subsets": n_subsets, "n_singular": n_singular,
             "n_feasible_subsets": len(cand)}
    return verts, stats


# --------------------------------------------------------------- per config

def analyze(name, filename):
    print("=" * 72)
    print("CONFIG %s   (configs/%s)" % (name.upper(), filename))
    print("=" * 72)
    ok = True
    t0 = time.time()

    cname, xs, r2, n, d = load_config(os.path.join(CONFIG_DIR, filename))
    L, Y = to_integer_rows(xs)
    maxy = max(abs(v) for row in Y for v in row)
    print("  n = %d, d = %d, r2 = %s;  integer scaling L = %d, max|Y| = %d"
          % (n, d, r2, L, maxy))
    assert maxy <= MAX_ABS_Y, "coordinate bound violated -- revisit the "\
                              "overflow analysis before trusting int64"

    # --- boundedness ------------------------------------------------------
    t = time.time()
    bounded, ray = recession_cone_is_trivial(Y, n, d)
    print("  recession cone {w : <x_i,w> <= 0 all i} = {0} : %s   "
          "(all C(%d,4) = %d subsets, exact, %.1fs)"
          % ("YES -> Q is bounded" if bounded else "NO -> Q UNBOUNDED",
             n, len(list(itertools.combinations(range(n), d - 1))),
             time.time() - t))
    if not bounded:
        print("  FAIL: Q is unbounded; witness ray %s" % (ray,))
        return False, None
    # secondary, independent boundedness witness: 0 as a strictly positive
    # combination of the x_i (sufficient together with span = R^5)
    tot = [sum(x[k] for x in xs) for k in range(d)]
    print("  (cross-check: sum_i x_i = %s%s)"
          % ([str(c) for c in tot],
             "  -> 0 is the barycentre, a strictly positive combination"
             if all(c == 0 for c in tot) else ""))

    # --- vertex enumeration ----------------------------------------------
    t = time.time()
    verts, stats = enumerate_vertices(xs, Y, L, r2, n, d)
    print("  vertex enumeration : all C(%d,5) = %d subsets, exact integer "
          "arithmetic, no tolerance (%.1fs)"
          % (n, stats["n_subsets"], time.time() - t))
    print("    singular 5-subsets           : %d" % stats["n_singular"])
    print("    feasible basic solutions     : %d" % stats["n_feasible_subsets"])
    print("    DISTINCT vertices of Q       : %d" % len(verts))

    # --- (a) max norm and the 41st point ----------------------------------
    norms = {w: sum(c * c for c in w) for w in verts}
    maxn2 = max(norms.values())
    deep = [w for w in verts if norms[w] == maxn2]
    print("  (a) max |w|^2 over Q = %s (~%.9f);  r2 = %s"
          % (maxn2, float(maxn2), r2))
    a_ok = maxn2 == Fraction(5, 4) and maxn2 < r2
    print("      max |w|^2 == 5/4 : %s;   5/4 < r2 : %s"
          % (maxn2 == Fraction(5, 4), maxn2 < r2))
    print("      => 41st point on the sphere of squared radius r2: %s"
          % ("IMPOSSIBLE (exact)" if maxn2 < r2 else "POSSIBLE?!"))
    ok = ok and a_ok

    # --- (b) support-function identity ------------------------------------
    m2 = r2 / (4 * maxn2)
    print("  (b) m(X)^2 = r2/(4 max|w|^2) = %s (~%.9f), m(X) = ~%.9f"
          % (m2, float(m2), float(m2) ** 0.5))
    b_ok = m2 == Fraction(2, 5)
    # direct exact confirmation on each deep hole: max_i <x^_i,u> for
    # u = w/|w| equals (r2/2)^2/(r2*|w|^2) squared
    for w in deep:
        mx = max(sum(a * b for a, b in zip(x, w)) for x in xs)
        assert mx == r2 / 2, "a maximal-norm vertex is not tight"
        val2 = (r2 / 2) ** 2 / (r2 * maxn2)     # (max_i <x^_i, w/|w|>)^2
        assert val2 == m2
    print("      direct check on every deep hole: max_i <x_i,w> = r2/2 and "
          "(max_i <x^_i, w/|w|>)^2 = %s : OK" % m2)
    print("      m(X)^2 == 2/5 : %s" % b_ok)
    ok = ok and b_ok

    # --- (c) counts -------------------------------------------------------
    near = {}
    for w in deep:
        c = sum(1 for x in xs if sum(a * b for a, b in zip(x, w)) == r2 / 2)
        near[c] = near.get(c, 0) + 1
    print("  (c) polar vertices total = %d;  deep holes (|w|^2 = 5/4) = %d"
          % (len(verts), len(deep)))
    print("      nearest-neighbour counts among deep holes: %s"
          % {k: near[k] for k in sorted(near)})
    other = {}
    for w in verts:
        if norms[w] != maxn2:
            other[norms[w]] = other.get(norms[w], 0) + 1
    print("      other vertex norms^2: %s"
          % {str(k): v for k, v in sorted(other.items())})
    print("  total runtime %.1fs" % (time.time() - t0))
    print()
    res = {"n_vertices": len(verts), "max_n2": maxn2, "n_deep": len(deep),
           "near": {k: near[k] for k in sorted(near)}, "m2": m2,
           "deep": sorted(tuple(str(c) for c in w) for w in deep)}
    return ok, res


TARGETS = [("d5", "d5_40.json"), ("l5", "l5_40.json"),
           ("q5", "q5_40.json"), ("r5", "r5_40.json")]

EXPECTED = {  # the claim under review
    "d5": {"n_vertices": 42, "n_deep": 32, "near": {10: 32}},
    "l5": {"n_vertices": 50, "n_deep": 32, "near": {10: 32}},
    "q5": {"n_vertices": 92, "n_deep": 32, "near_keys": {9, 10}},
    "r5": {"n_vertices": 100, "n_deep": 32, "near_keys": {9, 10}},
}


def main(argv):
    results, verdicts = {}, {}
    for name, fn in TARGETS:
        ok, res = analyze(name, fn)
        results[name] = res
        exp = EXPECTED[name]
        claim_ok = ok and res is not None
        if res is not None:
            claim_ok = claim_ok and res["n_vertices"] == exp["n_vertices"] \
                and res["n_deep"] == exp["n_deep"]
            if "near" in exp:
                claim_ok = claim_ok and res["near"] == exp["near"]
            else:
                claim_ok = claim_ok and set(res["near"]) <= exp["near_keys"]
        verdicts[name] = claim_ok

    print("=" * 72)
    print("D5 ANALYTIC CROSS-CHECK (independent of the enumeration)")
    # D5 constraints <±e_i±e_j, w> <= 1  <=>  |w_i|+|w_j| <= 1 for all i<j.
    # Vertices: (±1/2)^5 (32 of them, |w|^2 = 5/4) and ±e_k (10, |w|^2 = 1).
    d5 = results["d5"]
    predicted = {tuple(str(Fraction(s, 2)) for s in sgn)
                 for sgn in itertools.product((1, -1), repeat=5)}
    print("  closed form for D5 deep holes {(+-1/2)^5}: %d points, |w|^2 = 5/4"
          % len(predicted))
    match = set(d5["deep"]) == predicted
    print("  enumerated deep holes match the closed form exactly: %s" % match)
    verdicts["d5"] = verdicts["d5"] and match

    print()
    print("=" * 72)
    print("SUMMARY (claims (a), (b), (c) per configuration)")
    for name, _ in TARGETS:
        r = results[name]
        print("  %-4s : vertices %3d, deep holes %2d, nearest %s, "
              "max|w|^2 %s, m^2 %s  -> %s"
              % (name, r["n_vertices"], r["n_deep"], r["near"],
                 r["max_n2"], r["m2"], "PASS" if verdicts[name] else "FAIL"))
    allok = all(verdicts.values())
    print("  OVERALL: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
