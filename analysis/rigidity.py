"""Exact rigidity analysis of kissing configurations (CJKT-style).

Setting (Cohn-Jiao-Kumar-Torquato, "Rigidity of spherical codes", Geom.
Topol. 15 (2011); their Sec. 2): a configuration x_1..x_n on the sphere of
squared radius r2, with the kissing constraints <x_i,x_j> <= r2/2. A
first-order (infinitesimal) flex is v_1..v_n with

    (T)  <x_i, v_i> = 0                      for all i            (stay on sphere)
    (G)  <x_i, v_j> + <x_j, v_i> <= 0        for all contacts i<j (contacts may
                                              only open; contact = pair at
                                              exactly r2/2)

Trivial flexes are v_i = A x_i for antisymmetric A (dimension 10 in R^5,
since the configurations span R^5). The configuration is *infinitesimally
jammed* if every flex is trivial; infinitesimally jammed implies jammed
(Connelly; Roth-Whiteley; see CJKT Thm 2.2).

Certificate used here (all exact rational arithmetic):

1. Uniform self-stress: if s_i := sum of contact neighbors of x_i is
   parallel to x_i for every i, then for every flex v,
       sum_{contacts} (<x_i,v_j> + <x_j,v_i>) = sum_i <s_i, v_i>
                                              = sum_i lambda_i <x_i,v_i> = 0,
   and since each contact term is <= 0, every term is 0: all flexes
   preserve all contacts to first order (are equality flexes).
2. Equality nullspace: the nullspace of the stacked equality system
   [T; G(=)] contains the 10 rotation flexes; if its rank is 190
   (= 200 - 10) the nullspace is exactly the rotations.
   Rank certificate: for an integer matrix, rank over F_p <= rank over Q;
   exhibiting rank 190 mod one prime proves rank >= 190, and the rotation
   flexes prove rank <= 190. Hence rank = 190 exactly.

1 + 2 together: infinitesimally jammed, hence jammed.

If the uniform stress fails or the rank is deficient, this script reports
that and does NOT conclude; the LP / second-order machinery is separate.

Only rational configurations are supported (all committed configs are
rational). Everything runs in fractions.Fraction / exact integers.
"""

import argparse
import itertools
import json
import sys
from fractions import Fraction

PRIMES = [1000003, 1000033, 1000037]  # for the mod-p rank certificate


def load_rational_config(path):
    with open(path) as f:
        cfg = json.load(f)
    vecs = [[Fraction(c) for c in v] for v in cfg["vectors"]]
    r2 = Fraction(cfg["norm_squared"])
    n, d = cfg["n_points"], cfg["dimension"]
    assert len(vecs) == n and all(len(v) == d for v in vecs)
    for v in vecs:
        assert sum(c * c for c in v) == r2, "norm mismatch (run verify/ first)"
    return cfg["name"], vecs, r2, n, d


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def contact_graph(vecs, r2):
    """Pairs at inner product exactly r2/2."""
    n = len(vecs)
    return [(i, j) for i, j in itertools.combinations(range(n), 2)
            if dot(vecs[i], vecs[j]) == r2 / 2]


def uniform_stress_check(vecs, contacts):
    """Is sum of contact neighbors of x_i parallel to x_i, for every i?

    Returns (ok, details) where details lists per-vertex (degree, lambda)
    with lambda s.t. s_i = lambda * x_i, or the first failing vertex.
    """
    n = len(vecs)
    nbrs = [[] for _ in range(n)]
    for i, j in contacts:
        nbrs[i].append(j)
        nbrs[j].append(i)
    lambdas = []
    for i in range(n):
        s = [sum(vecs[j][k] for j in nbrs[i]) for k in range(len(vecs[i]))]
        # solve s = lam * x_i exactly: use any nonzero coordinate of x_i
        pivot = next(k for k, c in enumerate(vecs[i]) if c != 0)
        lam = s[pivot] / vecs[i][pivot]
        if any(s[k] != lam * vecs[i][k] for k in range(len(s))):
            return False, {"failing_vertex": i, "s_i": [str(c) for c in s]}
        lambdas.append((len(nbrs[i]), lam))
    return True, {"degree_lambda_pairs": sorted(set(lambdas))}


def equality_system_int(vecs, contacts, r2):
    """Integer matrix of the equality system [T; G] in variables v (n*d).

    Rows: <x_i, v_i> = 0 for each i; <x_i,v_j> + <x_j,v_i> = 0 per contact.
    Entries are scaled to integers by the common denominator.
    """
    n, d = len(vecs), len(vecs[0])
    rows = []
    for i in range(n):
        row = [Fraction(0)] * (n * d)
        for k in range(d):
            row[i * d + k] = vecs[i][k]
        rows.append(row)
    for i, j in contacts:
        row = [Fraction(0)] * (n * d)
        for k in range(d):
            row[j * d + k] += vecs[i][k]
            row[i * d + k] += vecs[j][k]
        rows.append(row)
    # scale each row to integers (common denominator per row)
    int_rows = []
    for row in rows:
        den = 1
        for c in row:
            if c != 0:
                den = den * c.denominator // _gcd(den, c.denominator)
        int_rows.append([int(c * den) for c in row])
    return int_rows


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def rank_mod_p(int_rows, p):
    """Rank of an integer matrix over F_p (row-reduction, O(m n rank))."""
    rows = [[c % p for c in row] for row in int_rows]
    ncols = len(rows[0])
    rank = 0
    pivot_row = 0
    for col in range(ncols):
        sel = next((r for r in range(pivot_row, len(rows)) if rows[r][col]), None)
        if sel is None:
            continue
        rows[pivot_row], rows[sel] = rows[sel], rows[pivot_row]
        inv = pow(rows[pivot_row][col], p - 2, p)
        prow = rows[pivot_row]
        prow[:] = [(c * inv) % p for c in prow]
        for r in range(len(rows)):
            if r != pivot_row and rows[r][col]:
                f = rows[r][col]
                rr = rows[r]
                rows[r] = [(a - f * b) % p for a, b in zip(rr, prow)]
        pivot_row += 1
        rank += 1
        if pivot_row == len(rows):
            break
    return rank


def rotation_flexes_check(vecs, int_rows):
    """Verify the 10 so(5) flexes v_i = A x_i satisfy the equality system,
    and that they are linearly independent (rank 10 mod a prime)."""
    d = len(vecs[0])
    flexes = []
    for a in range(d):
        for b in range(a + 1, d):
            # A = E_ab - E_ba
            flex = []
            for x in vecs:
                v = [Fraction(0)] * d
                v[a] = x[b]
                v[b] = -x[a]
                flex.extend(v)
            flexes.append(flex)
    # each flex must satisfy every equality row exactly (over Q)
    for row in int_rows:
        for flex in flexes:
            assert sum(Fraction(c) * f for c, f in zip(row, flex) if c) == 0
    # independence: rank of the 10 x (n*d) matrix
    den = 1
    for flex in flexes:
        for c in flex:
            if c != 0:
                den = den * c.denominator // _gcd(den, c.denominator)
    int_flexes = [[int(c * den) for c in flex] for flex in flexes]
    return rank_mod_p(int_flexes, PRIMES[0]) == len(flexes), len(flexes)


def analyze(path):
    name, vecs, r2, n, d = load_rational_config(path)
    nvars = n * d
    contacts = contact_graph(vecs, r2)
    result = {"name": name, "n": n, "d": d, "n_contacts": len(contacts),
              "n_vars": nvars}

    stress_ok, stress_info = uniform_stress_check(vecs, contacts)
    result["uniform_stress"] = stress_ok
    result["stress_info"] = {
        k: [(deg, str(lam)) for deg, lam in v] if k == "degree_lambda_pairs" else v
        for k, v in stress_info.items()
    }

    int_rows = equality_system_int(vecs, contacts, r2)
    rot_ok, n_rot = rotation_flexes_check(vecs, int_rows)
    result["rotation_flexes_in_nullspace_and_independent"] = rot_ok
    result["n_trivial_flexes"] = n_rot

    target_rank = nvars - n_rot
    ranks = {}
    for p in PRIMES:
        ranks[p] = rank_mod_p(int_rows, p)
        if ranks[p] == target_rank:
            break
    result["rank_mod_p"] = ranks
    result["target_rank"] = target_rank
    rank_certified = rot_ok and any(rk == target_rank for rk in ranks.values())
    result["equality_nullspace_is_exactly_rotations"] = rank_certified

    if stress_ok and rank_certified:
        result["verdict"] = "INFINITESIMALLY JAMMED (hence jammed)"
    elif not stress_ok:
        result["verdict"] = ("UNIFORM STRESS FAILS - inconclusive here; "
                             "needs LP-based flex analysis")
    else:
        result["verdict"] = ("EQUALITY NULLSPACE EXCEEDS ROTATIONS - "
                             "nontrivial equality flex exists; needs exact "
                             "nullspace basis + second-order analysis")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("configs", nargs="+")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    for path in args.configs:
        r = analyze(path)
        if args.json:
            print(json.dumps(r))
            continue
        print(f"== {r['name']} ==")
        print(f"  points {r['n']}, contacts (ip = r2/2): {r['n_contacts']}")
        print(f"  uniform stress (sum of contact nbrs parallel to x_i): "
              f"{'YES' if r['uniform_stress'] else 'NO'}  {r['stress_info']}")
        print(f"  rotation flexes valid+independent: "
              f"{r['rotation_flexes_in_nullspace_and_independent']} "
              f"(dim {r['n_trivial_flexes']})")
        print(f"  rank of equality system mod p: {r['rank_mod_p']} "
              f"(target {r['target_rank']})")
        print(f"  VERDICT: {r['verdict']}")
        print()


if __name__ == "__main__":
    main()
