"""Decide strict-flex existence for configurations where the uniform stress
fails, via exact self-stress analysis (Stiemke alternative).

Background (continues analysis/rigidity.py): with the equality nullspace
already certified to be exactly the 10 rotations, the configuration is
infinitesimally jammed iff there is NO flex v with Tv = 0, Gv <= 0 and
Gv != 0. By Stiemke's lemma applied in the tangency nullspace, exactly one
of the following holds:

  (a) there exists a strictly positive self-stress: omega in R^{contacts},
      omega_e > 0 for all e, such that for every vertex i,
          sum_{j in contacts(i)} omega_ij x_j  is parallel to  x_i
      -> no strict flex -> (with rank 190) INFINITESIMALLY JAMMED -> jammed.
  (b) there exists a flex v strictly opening at least one contact
      -> first-order unjamming direction (second-order analysis needed).

Method, all verdict-relevant steps in exact rational arithmetic:
 1. Build the stress matrix M (per vertex, the component of
    sum_j omega_ij x_j orthogonal to x_i must vanish; lambda eliminated
    exactly): M omega = 0 defines the stress space.
 2. Exact nullspace basis B of M via Fraction Gaussian elimination.
 3. Search for strictly positive omega in span(B): float LP (Chebyshev-style
    margin maximization) as a *heuristic direction finder only*, then
    rationalize and verify omega > 0 and M omega = 0 exactly.
 4. If the float LP says the max margin is 0 (no positive stress), find the
    strict flex from the LP dual / primal flex LP, rationalize, and verify
    exactly: Tv = 0, Gv <= 0, Gv != 0.
Either way the final certificate is exact; the floats never decide.
"""

import itertools
import json
import sys
from fractions import Fraction

from rigidity import (contact_graph, dot, equality_system_int,
                      load_rational_config)


def stress_matrix(vecs, contacts, r2):
    """Rows: for each vertex i, each coordinate k of
    sum_j omega_ij (x_j - (<x_i,x_j>/r2) x_i) = 0. Columns = contacts."""
    n, d = len(vecs), len(vecs[0])
    m = len(contacts)
    rows = [[Fraction(0)] * m for _ in range(n * d)]
    for e, (i, j) in enumerate(contacts):
        for (a, b) in ((i, j), (j, i)):
            cab = dot(vecs[a], vecs[b]) / r2
            for k in range(d):
                rows[a * d + k][e] += vecs[b][k] - cab * vecs[a][k]
    return [r for r in rows if any(c != 0 for c in r)]


def exact_nullspace(rows, ncols):
    """Nullspace basis of a rational matrix, Fraction RREF."""
    mat = [row[:] for row in rows]
    pivots = {}
    prow = 0
    for col in range(ncols):
        sel = next((r for r in range(prow, len(mat)) if mat[r][col] != 0), None)
        if sel is None:
            continue
        mat[prow], mat[sel] = mat[sel], mat[prow]
        inv = 1 / mat[prow][col]
        mat[prow] = [c * inv for c in mat[prow]]
        for r in range(len(mat)):
            if r != prow and mat[r][col] != 0:
                f = mat[r][col]
                mat[r] = [a - f * b for a, b in zip(mat[r], mat[prow])]
        pivots[col] = prow
        prow += 1
        if prow == len(mat):
            break
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [Fraction(0)] * ncols
        v[fc] = Fraction(1)
        for pc, pr in pivots.items():
            v[pc] = -mat[pr][fc]
        basis.append(v)
    return basis


def find_positive_stress(basis, m):
    """Float LP: max t s.t. omega = B c, omega_e >= t, sum omega = m.
    Returns (t_opt, c_opt) as floats, or (None, None) if LP fails."""
    import numpy as np
    from scipy.optimize import linprog
    k = len(basis)
    B = np.array([[float(basis[b][e]) for b in range(k)] for e in range(m)])
    # variables: c (k), t (1). max t -> min -t
    cobj = np.zeros(k + 1)
    cobj[-1] = -1.0
    A_ub = np.hstack([-B, np.ones((m, 1))])       # t - (Bc)_e <= 0
    b_ub = np.zeros(m)
    A_eq = np.hstack([B.sum(axis=0)[None, :], np.zeros((1, 1))])
    b_eq = np.array([float(m)])
    res = linprog(cobj, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(None, None)] * (k + 1), method="highs")
    if not res.success:
        return None, None
    return res.x[-1], res.x[:-1]


def rationalize(x, max_den=10**6):
    return Fraction(x).limit_denominator(max_den)


def verify_positive_stress(basis, c_rat, vecs, contacts, r2):
    """Exact check: omega = B c > 0 componentwise and stress condition holds
    (recomputed from scratch, not via the basis)."""
    m = len(contacts)
    omega = [sum(basis[b][e] * c_rat[b] for b in range(len(basis)))
             for e in range(m)]
    if any(w <= 0 for w in omega):
        return None
    n, d = len(vecs), len(vecs[0])
    s = [[Fraction(0)] * d for _ in range(n)]
    for e, (i, j) in enumerate(contacts):
        for k in range(d):
            s[i][k] += omega[e] * vecs[j][k]
            s[j][k] += omega[e] * vecs[i][k]
    for i in range(n):
        lam = dot(s[i], vecs[i]) / r2
        if any(s[i][k] != lam * vecs[i][k] for k in range(d)):
            return None
    return omega


def find_strict_flex_float(vecs, contacts, r2):
    """Float LP for a strict flex: Tv=0, (Gv)_e <= 0, sum (Gv)_e = -1.
    Returns float v or None."""
    import numpy as np
    from scipy.optimize import linprog
    n, d = len(vecs), len(vecs[0])
    nv = n * d
    X = [[float(c) for c in v] for v in vecs]
    T = np.zeros((n, nv))
    for i in range(n):
        T[i, i * d:(i + 1) * d] = X[i]
    G = np.zeros((len(contacts), nv))
    for e, (i, j) in enumerate(contacts):
        G[e, i * d:(i + 1) * d] = X[j]
        G[e, j * d:(j + 1) * d] = X[i]
    A_eq = np.vstack([T, G.sum(axis=0)[None, :]])
    b_eq = np.zeros(n + 1)
    b_eq[-1] = -1.0
    res = linprog(np.zeros(nv), A_ub=G, b_ub=np.zeros(len(contacts)),
                  A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)] * nv,
                  method="highs")
    return res.x if res.success else None


def verify_strict_flex(v_rat, vecs, contacts, r2):
    """Exact check: <x_i,v_i>=0 all i; contact derivatives <=0; some < 0.
    Returns (ok, n_opening, openings)."""
    n, d = len(vecs), len(vecs[0])
    for i in range(n):
        if dot(vecs[i], v_rat[i]) != 0:
            return False, 0, None
    openings = []
    for (i, j) in contacts:
        g = dot(vecs[i], v_rat[j]) + dot(vecs[j], v_rat[i])
        if g > 0:
            return False, 0, None
        if g < 0:
            openings.append(((i, j), g))
    return True, len(openings), openings


def analyze(path):
    name, vecs, r2, n, d = load_rational_config(path)
    contacts = contact_graph(vecs, r2)
    m = len(contacts)
    print(f"== {name} ==")
    M = stress_matrix(vecs, contacts, r2)
    basis = exact_nullspace(M, m)
    print(f"  contacts: {m}; stress-space dim (exact): {len(basis)}")
    if not basis:
        print("  stress space trivial -> no positive stress -> strict flex "
              "exists (find it via LP)")
    else:
        t_opt, c_opt = find_positive_stress(basis, m)
        print(f"  float LP max min-stress margin: {t_opt}")
        if t_opt is not None and t_opt > 1e-9:
            for max_den in (10**3, 10**6, 10**12):
                c_rat = [rationalize(c, max_den) for c in c_opt]
                omega = verify_positive_stress(basis, c_rat, vecs, contacts, r2)
                if omega is not None:
                    lo = min(omega)
                    print(f"  EXACT strictly positive self-stress verified "
                          f"(min weight {lo} = ~{float(lo):.4f}, "
                          f"den cap {max_den}).")
                    print("  VERDICT: with rank-190 certificate from "
                          "rigidity.py => INFINITESIMALLY JAMMED (hence "
                          "jammed)")
                    return {"name": name, "verdict": "jammed",
                            "stress_dim": len(basis),
                            "omega": [str(w) for w in omega],
                            "contacts": contacts}
            print("  rationalization failed - refine LP/denominators")
            return {"name": name, "verdict": "undecided"}
    # no positive stress: find and verify a strict flex
    vfl = find_strict_flex_float(vecs, contacts, r2)
    if vfl is None:
        print("  float flex LP infeasible?! inconsistent with stress LP - "
              "investigate")
        return {"name": name, "verdict": "undecided"}
    for max_den in (10**3, 10**6, 10**9, 10**12):
        v_rat = [[rationalize(vfl[i * d + k], max_den) for k in range(d)]
                 for i in range(n)]
        ok, n_open, openings = verify_strict_flex(v_rat, vecs, contacts, r2)
        if ok and n_open > 0:
            print(f"  EXACT strict first-order flex verified: {n_open} "
                  f"contacts strictly opening (den cap {max_den}).")
            print("  VERDICT: NOT infinitesimally jammed - first-order "
                  "unjamming direction exists (second-order analysis next)")
            return {"name": name, "verdict": "strict_flex",
                    "flex": [[str(c) for c in vi] for vi in v_rat],
                    "n_opening": n_open,
                    "opening_edges": [[e, str(g)] for e, g in openings]}
    print("  flex rationalization failed - refine")
    return {"name": name, "verdict": "undecided"}


def main():
    out = {}
    for path in sys.argv[1:]:
        if path == "--save":
            continue
        out[path] = analyze(path)
        print()
    if "--save" in sys.argv:
        import os
        dest = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "certificates")
        os.makedirs(dest, exist_ok=True)
        for k, v in out.items():
            cert = dict(v)
            if "contacts" in cert:
                cert["contacts"] = [list(e) for e in cert["contacts"]]
            path = os.path.join(dest, f"{v['name']}_stress_certificate.json")
            with open(path, "w") as f:
                json.dump(cert, f, indent=1)
            print(f"saved {path}")


if __name__ == "__main__":
    main()
