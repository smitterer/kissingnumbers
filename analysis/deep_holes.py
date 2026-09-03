"""Exact deep-hole analysis of kissing configurations via polar-vertex
enumeration.

Question: for a configuration X of 40 vectors at squared norm r2 (unit
vectors after dividing by sqrt(r2)), can a 41st point be added without
deforming X, and how deep are the holes?

Formulation: a 41st point exists iff there is w with <x_i, w> <= r2/2 for
all i and |w|^2 = r2. Let Q = {w in R^5 : <x_i, w> <= r2/2 for all i}.
Boundedness of Q is CHECKED exactly below: sum_i x_i = 0 with X spanning
R^5 puts the origin in the interior of conv(X), which makes Q compact
(independent verification also confirmed recession cone {0} by exhaustive
ray enumeration). Then max_{w in Q} |w|^2 is attained at a vertex of Q
(strict convexity: every maximizer is a vertex, so vertex enumeration
finds ALL deep holes).

  - If max vertex norm^2 < r2: NO 41st point exists (exact impossibility).
  - The unit min-max inner product is m(X) = sqrt(r2) / (2 |w*|), i.e.
    m(X)^2 = r2 / (4 max|w|^2), attained at u* = w*/|w*| — the deep holes
    of X are exactly the maximal-norm vertices of Q, normalized.

Completeness: every vertex of Q is the unique solution of <x_i, w> = r2/2
for some 5-subset of linearly independent constraints (basic feasible
solution). We enumerate ALL C(40,5) = 658,008 subsets. Float linear
algebra (batched numpy) is used only to PREFILTER with loose tolerances
(near-singular systems are skipped only below det threshold; feasibility
tolerance is generous, so any true vertex survives the filter); every
surviving candidate is then recomputed and verified in exact rational
arithmetic, and deduplicated exactly. The reported max is over the exact
verified vertex set.

Outputs per configuration: vertex count of Q, exact max |w|^2, exact
m(X)^2, deep-hole count and their nearest-neighbor counts, verdict on the
41st point. By default nothing is written; with --save the results (exact
rational strings only, no floats, no timings) are written to
analysis/deep_holes_results.json.

Completeness of the vertex list -- that the float prefilter above loses no
vertex -- is certified independently by
analysis/certificates/independent_deep_holes_check.py, which enumerates the
same C(40,5) subsets in exact integer arithmetic and uses no tolerance.

Validation target: for D5, Cohn-Rajagopal (arXiv:2412.00937, Sec. 3)
report 32 deep holes at unit inner product 2*sqrt(2/5)/2 = sqrt(2/5),
i.e. m(D5)^2 = 2/5, with 10 nearest neighbors each.
"""

import itertools
import json
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from rigidity import dot, load_rational_config

DET_TOL = 1e-10
FEAS_TOL = 1e-7


def solve5(rows, rhs):
    """Exact solution of a 5x5 rational system; None if singular."""
    d = len(rows)
    a = [list(map(Fraction, rows[i])) + [Fraction(rhs[i])] for i in range(d)]
    for col in range(d):
        sel = next((r for r in range(col, d) if a[r][col] != 0), None)
        if sel is None:
            return None
        a[col], a[sel] = a[sel], a[col]
        inv = 1 / a[col][col]
        a[col] = [c * inv for c in a[col]]
        for r in range(d):
            if r != col and a[r][col] != 0:
                f = a[r][col]
                a[r] = [x - f * y for x, y in zip(a[r], a[col])]
    return [a[i][d] for i in range(d)]


def analyze(path):
    import numpy as np

    t0 = time.time()
    name, vecs, r2, n, d = load_rational_config(path)
    # boundedness certificate: sum_i x_i = 0 exactly and X spans R^5
    # => 0 in int conv(X) => Q compact
    assert all(sum(v[k] for v in vecs) == 0 for k in range(d)), \
        "sum of points != 0; boundedness of Q not certified by this route"
    Xf = np.array([[float(c) for c in v] for v in vecs])
    half = float(r2) / 2.0

    subsets = np.array(list(itertools.combinations(range(n), d)),
                       dtype=np.int32)
    cand_subsets = []
    chunk = 40000
    for s0 in range(0, len(subsets), chunk):
        S = subsets[s0:s0 + chunk]
        A = Xf[S]                              # (m, 5, 5)
        dets = np.linalg.det(A)
        ok = np.abs(dets) > DET_TOL
        if not ok.any():
            continue
        Sok = S[ok]
        W = np.linalg.solve(A[ok], np.full((int(ok.sum()), d, 1), half))[..., 0]
        feas = (W @ Xf.T <= half + FEAS_TOL).all(axis=1)
        for idx in np.nonzero(feas)[0]:
            cand_subsets.append(tuple(int(i) for i in Sok[idx]))

    # exact verification + dedup
    verts = {}
    for sub in cand_subsets:
        w = solve5([vecs[i] for i in sub], [r2 / 2] * d)
        if w is None:
            continue
        key = tuple(w)
        if key in verts:
            continue
        if all(dot(vecs[j], w) <= r2 / 2 for j in range(n)):
            verts[key] = w
    assert verts, "no vertices found - unexpected"

    norms = {key: dot(list(key), list(key)) for key in verts}
    max_n2 = max(norms.values())
    deep = [verts[k] for k, v in norms.items() if v == max_n2]
    m2 = r2 / (4 * max_n2)
    nearest_counts = sorted({sum(1 for j in range(n)
                                 if dot(vecs[j], w) == r2 / 2) for w in deep})

    result = {
        "name": name,
        "n_polar_vertices": len(verts),
        "max_vertex_norm_sq": str(max_n2),
        "minmax_unit_ip_sq": str(m2),
        "n_deep_holes": len(deep),
        "deep_hole_nearest_counts": nearest_counts,
        "deep_holes_scaled": [[str(c) for c in w] for w in deep[:64]],
        "no_41st_point": bool(max_n2 < r2),
    }
    runtime_s = time.time() - t0
    print(f"== {name} ==")
    print(f"  polar polytope Q vertices (exact, deduped): {len(verts)}")
    print(f"  max |w|^2 over Q = {max_n2}  (r2 = {r2}; 41st point needs |w|^2 >= r2)")
    print(f"  min-max unit ip m(X) = sqrt({m2}) ~ {float(m2) ** 0.5:.6f}")
    print(f"  deep holes: {len(deep)}, nearest-neighbor counts {nearest_counts}")
    print(f"  41st point without deforming: "
          f"{'IMPOSSIBLE (exact: max |w|^2 < r2)' if result['no_41st_point'] else 'POSSIBLE?!'}")
    print(f"  runtime: {runtime_s:.1f}s")
    return result


def main():
    save = "--save" in sys.argv[1:]
    out = []
    for path in sys.argv[1:]:
        if path == "--save":
            continue
        out.append(analyze(path))
        print()
    if save:
        with open(__file__.rsplit("/", 1)[0] + "/deep_holes_results.json", "w") as f:
            json.dump(out, f, indent=1)
        print("saved analysis/deep_holes_results.json")
    else:
        print("(nothing written; pass --save to write analysis/deep_holes_results.json)")


if __name__ == "__main__":
    main()
