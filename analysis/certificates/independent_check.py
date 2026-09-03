#!/usr/bin/env python3
"""INDEPENDENT re-verification of the "all four 40-point 5-dimensional
kissing configurations are infinitesimally jammed" claim.

Written for the verifier role.  It shares NO function body with
analysis/rigidity.py or analysis/stress_lp.py (see IMPLEMENTATION NOTES at
the end of this docstring); it only reads

    configs/{d5,l5,q5,r5}_40.json
    analysis/certificates/{q5,r5}_stress_certificate.json

Everything that influences a verdict is exact (fractions.Fraction / int).
Floats appear only in the printed summary, always marked with '~'.

Setting (CJKT, Geom. Topol. 15 (2011), Sec. 2).  Points x_1..x_n with
<x_i,x_i> = r2 and <x_i,x_j> <= r2/2.  Contacts = pairs with equality.
An infinitesimal flex is v = (v_1..v_n) with

    (T)  <x_i, v_i> = 0                              for every i
    (G)  g_ij(v) := <x_i,v_j> + <x_j,v_i> <= 0       for every contact {i,j}

Trivial flexes: v_i = A x_i with A antisymmetric (dim 10 for d = 5, when the
points span R^5).  Infinitesimally jammed = every flex is trivial.

Certificate checked here, per configuration:

  (S) a self-stress omega with omega_e > 0 on EVERY contact e, meaning that
      for every vertex i,  s_i := sum_{j ~ i} omega_ij x_j  is an exact
      rational multiple of x_i.  Then for any flex v,
          0 >= sum_e omega_e g_e(v) = sum_i <s_i, v_i>
                                    = sum_i lambda_i <x_i, v_i> = 0,
      and since every summand omega_e g_e(v) is <= 0, all of them vanish:
      every flex is an equality flex.  (Note the lambda_i are multiplied by
      <x_i,v_i> = 0, so their signs -- or vanishing -- are irrelevant.)

  (R) the stacked equality system [T rows ; G rows] on the n*d variables has
      rank exactly n*d - 10 over Q.  Since rank over F_p is always <= rank
      over Q (a nonvanishing r x r minor mod p is a nonvanishing minor over
      Q), a rank of 190 mod p gives rank_Q >= 190; the 10 independent
      rotation flexes in the kernel give rank_Q <= 190.  Hence rank_Q = 190
      and the equality kernel is exactly the rotations.

  (S) + (R): every flex is an equality flex, and every equality flex is a
  rotation => infinitesimally jammed => jammed (CJKT Thm 2.2).

IMPLEMENTATION NOTES -- what is done differently from analysis/rigidity.py
and analysis/stress_lp.py (the method under review):

  * the equality system is assembled from the adjacency lists of the contact
    graph as sparse dict-of-rows and scaled by ONE global denominator of the
    configuration (rigidity.py builds dense Fraction rows from the contact
    list and clears denominators row by row);
  * the ten rotation flexes are obtained by multiplying every point with the
    explicit basis matrices E_ab - E_ba of so(d) (rigidity.py writes the two
    nonzero coordinates of A x directly);
  * parallelism s_i || x_i is decided by the vanishing of every 2x2 minor
    s_i[k] x_i[l] - s_i[l] x_i[k], with no division (stress_lp.py computes
    lambda_i = <s_i,x_i>/r2 and tests s_i == lambda_i x_i); s_i is
    accumulated by walking adjacency lists with an edge -> weight map;
  * the dimension of the self-stress space is obtained by duality from the
    rank of the equality system, dim S = (m + n) - rank[T; G] (stress_lp.py
    builds the matrix of tangential components and computes its nullspace);
  * rank over F_p is computed by forward elimination only, with different
    primes, and confirmed by an exact fraction-free Bareiss rank over Q that
    the method under review does not have.

Exit code 0 iff every configuration passes every check.
"""

from __future__ import annotations

import itertools
import json
import math
import operator
import os
import sys
import time
from collections import Counter
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # .../kissingnumbers
CONFIG_DIR = os.path.join(ROOT, "configs")

# Primes chosen independently of analysis/rigidity.py (which uses 1000003).
PRIMES = [2147483647, 998244353]

TRIVIAL_DIM = 10          # dim so(5)


# ----------------------------------------------------------------- loading

def load_config(path):
    with open(path) as fh:
        cfg = json.load(fh)
    d = cfg["dimension"]
    n = cfg["n_points"]
    r2 = Fraction(cfg["norm_squared"])
    xs = [[Fraction(c) for c in v] for v in cfg["vectors"]]
    if len(xs) != n or any(len(v) != d for v in xs):
        raise ValueError("%s: shape mismatch" % path)
    return cfg.get("name", "?"), xs, r2, n, d


def dot(u, v):
    """Exact inner product: sum of the coordinatewise products."""
    return sum(map(operator.mul, u, v))


def common_denominator(xs):
    """Least common multiple of all coordinate denominators (a positive int)."""
    return math.lcm(*[c.denominator for x in xs for c in x])


# --------------------------------------------------- recomputed contact set

def contacts_of(xs, r2):
    """All pairs i<j with <x_i,x_j> exactly r2/2, plus a sanity re-check that
    the configuration really is a kissing configuration (norms and <= r2/2)."""
    n = len(xs)
    half = r2 / 2
    for i, x in enumerate(xs):
        if dot(x, x) != r2:
            raise ValueError("vector %d does not have squared norm r2" % i)
    out = []
    for i, j in itertools.combinations(range(n), 2):
        ip = dot(xs[i], xs[j])
        if ip > half:
            raise ValueError("pair (%d,%d) violates <x_i,x_j> <= r2/2" % (i, j))
        if ip == half:
            out.append((i, j))
    return out


def adjacency(contacts, n):
    """Neighbour sets of the contact graph."""
    nbrs = [set() for _ in range(n)]
    for i, j in contacts:
        nbrs[i].add(j)
        nbrs[j].add(i)
    return nbrs


# --------------------------------------------------------- stress checking

def check_stress(xs, contacts, omega, r2):
    """Exact check of a self-stress, by 2x2 minors.

    Returns (ok, info).  ok is True iff every omega_e > 0 and, for every
    vertex i, s_i = sum_{j~i} omega_ij x_j is parallel to x_i.  Parallelism
    of two vectors s, x in Q^d with x != 0 holds iff all 2x2 minors
    s[k] x[l] - s[l] x[k] vanish; this test involves no division.  The
    scalar lambda_i with s_i = lambda_i x_i is then read off from one
    nonzero coordinate of x_i (and re-checked on every coordinate).
    info carries the lambdas and the multiplicities of omega and lambda.
    """
    n, d = len(xs), len(xs[0])
    if len(omega) != len(contacts):
        return False, {"error": "omega/contacts length mismatch"}
    nonpositive = [e for e, w in enumerate(omega) if w <= 0]
    if nonpositive:
        return False, {"error": "omega not strictly positive",
                       "bad_edges": nonpositive[:10]}
    weight = {}
    for (i, j), w in zip(contacts, omega):
        weight[(i, j)] = w
        weight[(j, i)] = w
    nbrs = adjacency(contacts, n)
    lambdas = []
    for i, x in enumerate(xs):
        s = [Fraction(0)] * d
        for j in nbrs[i]:
            w = weight[(i, j)]
            for k in range(d):
                s[k] += w * xs[j][k]
        for k, l in itertools.combinations(range(d), 2):
            if s[k] * x[l] - s[l] * x[k] != 0:
                return False, {"error": "s_i not parallel to x_i (2x2 minor)",
                               "vertex": i, "minor": (k, l),
                               "s_i": [str(c) for c in s]}
        k0 = next(k for k in range(d) if x[k] != 0)     # x_i != 0 (norm r2 > 0)
        lam = s[k0] / x[k0]
        if any(s[k] != lam * x[k] for k in range(d)):
            return False, {"error": "lambda inconsistent across coordinates",
                           "vertex": i}
        lambdas.append(lam)
    return True, {"lambdas": lambdas,
                  "omega_min": min(omega), "omega_max": max(omega),
                  "omega_mult": Counter(omega),
                  "lambda_mult": Counter(lambdas)}


def fmt_mult(counter):
    """'value xcount, ...' in increasing order of the value."""
    return ", ".join("%s x%d" % (v, c) for v, c in sorted(counter.items()))


# ------------------------------------------------------- equality system

def equality_system(xs, contacts, scale):
    """The equality system [T ; G] as sparse integer rows {column: entry}.

    Assembled from the contact graph: for every vertex i one tangency row
    {(i,k): x_i[k]}, and for every neighbour j > i of i one contact row
    {(j,k): x_i[k], (i,k): x_j[k]}, column index (i,k) -> i*d + k.  Every
    entry is multiplied by the same positive integer `scale` (a common
    denominator of the configuration), which changes neither the row space
    nor the rank.  Rows are dicts, so only the (at most 2d) nonzero entries
    of each row are stored.
    """
    n, d = len(xs), len(xs[0])
    nbrs = adjacency(contacts, n)

    def col(i, k):
        return i * d + k

    rows = []
    for i, x in enumerate(xs):
        rows.append({col(i, k): x[k] * scale for k in range(d) if x[k]})
    for i in range(n):
        for j in sorted(nbrs[i]):
            if j < i:
                continue                  # each contact once, as (i, j), i < j
            row = {}
            for k in range(d):
                if xs[i][k]:
                    row[col(j, k)] = xs[i][k] * scale
                if xs[j][k]:
                    row[col(i, k)] = xs[j][k] * scale
            rows.append(row)
    out = []
    for row in rows:
        r = {}
        for c, v in row.items():
            assert v.denominator == 1, "scale is not a common denominator"
            r[c] = int(v)
        out.append(r)
    return out


def dense(sparse_rows, ncols):
    """Expand sparse rows to dense integer lists (for the rank routines)."""
    return [[row.get(c, 0) for c in range(ncols)] for row in sparse_rows]


def rank_mod_p(int_rows, p):
    """Rank over F_p by Gaussian elimination (forward elimination only)."""
    rows = [[c % p for c in r] for r in int_rows]
    m = len(rows)
    ncols = len(rows[0])
    rank = 0
    for col in range(ncols):
        piv = None
        for r in range(rank, m):
            if rows[r][col]:
                piv = r
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        prow = rows[rank]
        inv = pow(prow[col], p - 2, p)
        if inv:
            prow[:] = [(c * inv) % p for c in prow]
        for r in range(rank + 1, m):
            f = rows[r][col]
            if f:
                rr = rows[r]
                rows[r] = [(a - f * b) % p for a, b in zip(rr, prow)]
        rank += 1
        if rank == m:
            break
    return rank


def rank_exact_Q(rows):
    """Exact rank over Q by fraction-free (Bareiss) forward elimination.

    Integer-preserving, so no Fraction blow-up; the pivot division is exact
    by the Bareiss identity.  Used as a redundant confirmation of the mod-p
    rank certificate.
    """
    mat = [r[:] for r in rows]
    m = len(mat)
    ncols = len(mat[0])
    rank = 0
    prev_piv = 1
    for col in range(ncols):
        piv = None
        for r in range(rank, m):
            if mat[r][col]:
                piv = r
                break
        if piv is None:
            continue
        mat[rank], mat[piv] = mat[piv], mat[rank]
        prow = mat[rank]
        pval = prow[col]
        for r in range(rank + 1, m):
            rr = mat[r]
            f = rr[col]
            if f:
                new = [pval * a - f * b for a, b in zip(rr, prow)]
            else:
                new = [pval * a for a in rr]
            # Bareiss: every entry must be exactly divisible by the previous
            # pivot.  Assert it rather than trusting the identity silently.
            if any(c % prev_piv for c in new):
                raise AssertionError("Bareiss division not exact -- aborting")
            mat[r] = [c // prev_piv for c in new]
        prev_piv = pval
        rank += 1
        if rank == m:
            break
    return rank


# ---------------------------------------------------------- rotation flexes

def so_basis(d):
    """The d(d-1)/2 matrices E_ab - E_ba (a < b), a basis of so(d)."""
    basis = []
    for a in range(d):
        for b in range(a + 1, d):
            A = [[0] * d for _ in range(d)]
            A[a][b] = 1
            A[b][a] = -1
            basis.append(A)
    return basis


def matvec(A, x):
    return [dot(row, x) for row in A]


def rotation_flexes(xs):
    """The trivial flexes v_i = A x_i, one per basis matrix A of so(d)."""
    d = len(xs[0])
    out = []
    for A in so_basis(d):
        flex = []
        for x in xs:
            flex.extend(matvec(A, x))
        out.append(flex)
    return out


def check_rotations(int_rows, flexes):
    """Every rotation flex satisfies every equality row EXACTLY over Q, and
    the flexes are linearly independent (exact rank over Q)."""
    for ri, row in enumerate(int_rows):
        for fi, flex in enumerate(flexes):
            val = sum(Fraction(c) * f for c, f in zip(row, flex) if c)
            if val != 0:
                return False, {"error": "rotation flex violates equality row",
                               "row": ri, "flex": fi, "value": str(val)}
    # exact independence: clear denominators, Bareiss rank over Q
    int_flexes = []
    for flex in flexes:
        den = math.lcm(*[Fraction(c).denominator for c in flex])
        int_flexes.append([int(Fraction(c) * den) for c in flex])
    r = rank_exact_Q(int_flexes)
    return r == len(flexes), {"rank_of_rotations": r, "n_rotations": len(flexes)}


# --------------------------------------------------------------- per-config

def analyze(name, config_file, cert_file=None, exact_rank=True):
    print("=" * 72)
    print("CONFIG %s   (%s)" % (name.upper(), config_file))
    print("=" * 72)
    ok = True

    cname, xs, r2, n, d = load_config(os.path.join(CONFIG_DIR, config_file))
    print("  name in file        : %s" % cname)
    print("  n = %d, d = %d, r2 = %s" % (n, d, r2))

    # --- 1. contacts, recomputed from scratch -------------------------------
    contacts = contacts_of(xs, r2)
    nbrs = adjacency(contacts, n)
    deg = [len(s) for s in nbrs]
    print("  contacts (<x_i,x_j> = r2/2 exactly) : %d" % len(contacts))
    print("  contact-graph degrees               : min %d, max %d%s"
          % (min(deg), max(deg),
             "  (regular)" if min(deg) == max(deg) else "  (IRREGULAR)"))
    if 2 * len(contacts) != sum(deg):
        print("  FAIL: degree sum inconsistent")
        ok = False

    # --- 2. stress certificate ---------------------------------------------
    if cert_file is None:
        omega = [Fraction(1)] * len(contacts)
        stress_source = "uniform omega_e = 1 (no certificate file)"
        cert_contacts = contacts
    else:
        with open(os.path.join(HERE, cert_file)) as fh:
            cert = json.load(fh)
        stress_source = "%s (verdict=%r, stress_dim=%s)" % (
            cert_file, cert.get("verdict"), cert.get("stress_dim"))
        cert_contacts = [tuple(e) for e in cert["contacts"]]
        omega = [Fraction(w) for w in cert["omega"]]
        if cert.get("name") != cname:
            print("  FAIL: certificate name %r != config name %r"
                  % (cert.get("name"), cname))
            ok = False
        # exact match of the contact lists, as ordered lists AND as sets
        if len(cert_contacts) != len(set(cert_contacts)):
            print("  FAIL: certificate contact list has duplicates")
            ok = False
        if any(not (0 <= i < j < n) for i, j in cert_contacts):
            print("  FAIL: certificate contains a malformed pair")
            ok = False
        if set(cert_contacts) != set(contacts):
            missing = sorted(set(contacts) - set(cert_contacts))
            extra = sorted(set(cert_contacts) - set(contacts))
            print("  FAIL: certificate contact set != recomputed contact set "
                  "(missing %d, extra %d)" % (len(missing), len(extra)))
            ok = False
        else:
            print("  certificate contact set == recomputed contact set "
                  "(%d pairs, %s order)"
                  % (len(cert_contacts),
                     "same" if cert_contacts == contacts else "different"))
        if len(omega) != len(cert_contacts):
            print("  FAIL: omega length != contacts length")
            ok = False

    print("  stress source       : %s" % stress_source)
    s_ok, s_info = check_stress(xs, cert_contacts, omega, r2)
    if not s_ok:
        print("  FAIL: stress check: %s" % s_info)
        ok = False
    else:
        lams = s_info["lambdas"]
        wmin, wmax = s_info["omega_min"], s_info["omega_max"]
        print("  stress omega_e > 0 on ALL %d contacts : YES" % len(omega))
        print("    omega min = %s (~%.6f), max = %s (~%.6f), distinct values %d"
              % (wmin, float(wmin), wmax, float(wmax), len(set(omega))))
        print("    s_i = lambda_i x_i exactly for all %d vertices; "
              "lambda values %s"
              % (n, sorted({str(l) for l in lams})))
        print("    omega multiplicities : %s" % fmt_mult(s_info["omega_mult"]))
        print("    lambda multiplicities: %s" % fmt_mult(s_info["lambda_mult"]))
        print("    (lambda signs are irrelevant to the argument: they multiply "
              "<x_i,v_i> = 0)")

    # --- 2b. the points must span R^d, else so(d) would not act faithfully
    scale = common_denominator(xs)
    span = rank_exact_Q([[int(c * scale) for c in x] for x in xs])
    print("  span of the points  : rank %d (need %d)" % (span, d))
    if span != d:
        print("  FAIL: points do not span R^%d" % d)
        ok = False

    # --- 3. rank of the equality system (computed here, printed below) ------
    m = len(contacts)
    nvars = n * d
    target = nvars - TRIVIAL_DIM
    int_rows = dense(equality_system(xs, contacts, scale), nvars)
    rank_lines = ["  equality system     : %d rows x %d columns (%d tangency + %d "
                  "contact)" % (len(int_rows), nvars, n, len(contacts))]
    ranks = {}
    for p in PRIMES:
        t0 = time.time()
        ranks[p] = rank_mod_p(int_rows, p)
        rank_lines.append("    rank mod %d = %d   (%.1fs)"
                          % (p, ranks[p], time.time() - t0))
    rank_ok = all(r == target for r in ranks.values())
    if not rank_ok:
        rank_lines.append("  FAIL: rank mod p != %d" % target)
        ok = False
    rq = None
    if exact_rank:
        t0 = time.time()
        rq = rank_exact_Q(int_rows)
        rank_lines.append("    rank over Q (Bareiss, exact integer) = %d   (%.1fs)"
                          % (rq, time.time() - t0))
        if rq != target:
            rank_lines.append("  FAIL: exact rank over Q = %d != %d" % (rq, target))
            ok = False
            rank_ok = False

    # --- 2c. informational: exact dimension of the whole self-stress space,
    #        by duality.  S = {omega : G^T omega in im(T^T)} is the projection
    #        of ker[G^T | -T^T] onto the omega coordinates; the projection is
    #        injective because T has full row rank n (its rows have disjoint
    #        supports and every x_i != 0), so
    #            dim S = (m + n) - rank[G^T | -T^T] = (m + n) - rank[T ; G].
    rank_T = sum(1 for x in xs if any(x))
    if rank_T != n:
        print("  FAIL: some x_i is the zero vector")
        ok = False
    if rq is not None:
        rank_Q = rq
    elif rank_ok:
        rank_Q = target          # rank_p = target and rotations give rank_Q <= target
    else:
        rank_Q = None
    if rank_Q is not None and rank_T == n:
        print("  self-stress space   : dim %d (exact, over Q; = (m+n) - rank[T;G] "
              "= (%d+%d) - %d)" % (m + n - rank_Q, m, n, rank_Q))
    else:
        print("  self-stress space   : dim undetermined (rank not certified)")

    for line in rank_lines:
        print(line)

    # --- 4. rotation flexes -------------------------------------------------
    flexes = rotation_flexes(xs)
    r_ok, r_info = check_rotations(int_rows, flexes)
    if not r_ok:
        print("  FAIL: rotation flexes: %s" % r_info)
        ok = False
    else:
        print("  rotation flexes     : all %d satisfy every equality row "
              "exactly over Q; exact rank = %d (independent)"
              % (r_info["n_rotations"], r_info["rank_of_rotations"]))

    if ok:
        print("  => equality kernel has dim %d - %d = %d = span(rotations)"
              % (nvars, target, TRIVIAL_DIM))
        print("  => every flex is an equality flex (positive self-stress) and "
              "every equality flex is trivial")
        print("  RESULT: INFINITESIMALLY JAMMED  (hence jammed, CJKT Thm 2.2)")
    else:
        print("  RESULT: CHECKS FAILED")
    print()
    return ok


TARGETS = [
    ("d5", "d5_40.json", None),
    ("l5", "l5_40.json", None),
    ("q5", "q5_40.json", "q5_stress_certificate.json"),
    ("r5", "r5_40.json", "r5_stress_certificate.json"),
]


def main(argv):
    exact_rank = "--no-exact-rank" not in argv
    results = {}
    for name, cfg, cert in TARGETS:
        results[name] = analyze(name, cfg, cert, exact_rank=exact_rank)
    print("=" * 72)
    print("SUMMARY")
    for name, ok in results.items():
        print("  %-4s : %s" % (name, "PASS" if ok else "FAIL"))
    allok = all(results.values())
    print("  OVERALL: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
