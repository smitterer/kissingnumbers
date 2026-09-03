"""Szollosi-style multiangular-cloud construction + exact max-clique search.

Reference: F. Szollosi, "A note on five dimensional kissing arrangements,"
arXiv:2301.08272 (2023), Section 3. Three methods; here Methods 1-3:

  Method 1 (cloud in R^d): fix a rational basis B (d x d, rows at common
  squared norm r2, pairwise compatible) and a finite angle set A of allowed
  inner products (all <= r2/2). Candidates are the v in span(B) = R^d with
  <v, B_i> = t_i for a tuple t in A^d and <v,v> = r2:  v = t G^{-1} B where
  G = B B^T. All exact rational.

  Method 2 (extend a basis of R^d to R^{d+1}): B is d x (d+1) of rank d,
  n a vector orthogonal to all rows of B with <n,n> = r2. For t in A^d,
  v0 = t G^{-1} B has squared norm q = t G^{-1} t^T; if q <= r2 the two
  candidates are v0 +- sqrt((r2-q)/r2) * n. Coordinates live in Q(sqrt(m));
  inner products between candidates are q_uv + eps_u eps_v sqrt(du dv)/r2
  with rational q_uv, du = r2 - q_u; compatibility (<= r2/2) is decided
  EXACTLY by isolating the radical and squaring with sign analysis.

  Method 3 (precomputed vectors): any finite candidate list; here used with
  unions of clouds.

Compatibility graph: vertices = candidates, edges between u != v with
<u,v> <= r2/2. A clique of size k plus the d (or d+1) basis vectors of a
compatible basis yields a kissing arrangement of size d + k (Methods 1-2,
when basis vectors are pairwise compatible and compatible-by-construction
with all cloud vectors) or k (Method 3).

Max clique: exact branch-and-bound with greedy-coloring upper bounds
(Ostergard-flavored). Guaranteed exact maximum (no heuristics), suitable
for clouds up to a few thousand vertices.

Everything verdict-relevant is exact (Fraction / integer surd arithmetic).
"""

import itertools
import json
from fractions import Fraction


# ---------------------------------------------------------------------------
# exact linear algebra over Q

def mat_inv(m):
    """Exact inverse of a rational matrix (Gauss-Jordan)."""
    n = len(m)
    a = [row[:] + [Fraction(int(i == j)) for j in range(n)]
         for i, row in enumerate(m)]
    for col in range(n):
        piv = next(r for r in range(col, n) if a[r][col] != 0)
        a[col], a[piv] = a[piv], a[col]
        inv = 1 / a[col][col]
        a[col] = [c * inv for c in a[col]]
        for r in range(n):
            if r != col and a[r][col] != 0:
                f = a[r][col]
                a[r] = [x - f * y for x, y in zip(a[r], a[col])]
    return [row[n:] for row in a]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def mat_vec(m, v):
    return [dot(row, v) for row in m]


# ---------------------------------------------------------------------------
# candidates: rational part + coefficient along a fixed surd direction
# vector value = rat + eps*sqrt(disc)*n_unitized, disc rational >= 0.
# We store (coords_rational, eps, disc) with the convention that the actual
# inner product of two candidates u, v is:
#   <u.rat, v.rat> + u.eps*v.eps*sqrt(u.disc*v.disc)/r2   (n has norm^2 r2)

class Cand:
    __slots__ = ("rat", "eps", "disc", "tag")

    def __init__(self, rat, eps, disc, tag):
        self.rat = tuple(rat)
        self.eps = eps          # -1, 0, +1
        self.disc = disc        # Fraction >= 0; 0 iff eps == 0
        self.tag = tag

    def key(self):
        return (self.rat, self.eps, self.disc)


def ip_leq(q, e, s, bound):
    """Decide exactly whether q + e*sqrt(s) <= bound, for rational q, s>=0,
    e in {-1,0,1}."""
    if e == 0 or s == 0:
        return q <= bound
    rhs = bound - q
    if e > 0:
        # sqrt(s) <= rhs  <=>  rhs >= 0 and s <= rhs^2
        return rhs >= 0 and s <= rhs * rhs
    # -sqrt(s) <= rhs  <=>  rhs >= 0 or (rhs < 0 and s >= rhs^2)
    return rhs >= 0 or s >= rhs * rhs


def compatible(u, v, r2):
    # ip(u, v) = <u.rat, v.rat> + eps_u*eps_v*sqrt(disc_u*disc_v)
    # (disc = r2 - |v0|^2; the normal's normalization cancels)
    q = dot(u.rat, v.rat)
    return ip_leq(q, u.eps * v.eps, u.disc * v.disc, r2 / 2)


def exact_ip_str(u, v, r2):
    """Human-readable exact inner product of two candidates (divide by r2
    for the unit-normalized value)."""
    q = dot(u.rat, v.rat)
    e = u.eps * v.eps
    D = u.disc * v.disc
    if e == 0 or D == 0:
        return str(q)
    return f"{q} {'+' if e > 0 else '-'} sqrt({D})"


# ---------------------------------------------------------------------------
# cloud construction

def cloud_method1(basis, r2, angles):
    """All v in span(basis)=R^d with <v,b_i> in angles and <v,v> = r2."""
    G = [[dot(a, b) for b in basis] for a in basis]
    Ginv = mat_inv(G)
    d = len(basis)
    out, seen = [], set()
    for t in itertools.product(angles, repeat=d):
        c = mat_vec(Ginv, list(t))            # coefficients
        if dot(c, list(t)) != r2:             # <v,v> = t G^-1 t^T
            continue
        v = [sum(c[i] * basis[i][k] for i in range(d))
             for k in range(len(basis[0]))]
        cand = Cand(v, 0, Fraction(0), tuple(t))
        if cand.key() not in seen:
            seen.add(cand.key())
            out.append(cand)
    return out


def cloud_method2(basis, normal, r2, angles):
    """Rank-d basis in R^{d+1}, normal orthogonal to all its rows (any
    rational normalization mu = <normal,normal>). Candidates
    v0 +- sqrt((r2-q)/mu)*normal with q = |v0|^2 <= r2; stored with
    disc = r2 - q so that cross inner products are
    <u.rat,v.rat> + eps_u*eps_v*sqrt(disc_u*disc_v)."""
    for b in basis:
        assert dot(b, normal) == 0, "normal not orthogonal to basis"
    G = [[dot(a, b) for b in basis] for a in basis]
    Ginv = mat_inv(G)
    d = len(basis)
    out, seen = [], set()
    for t in itertools.product(angles, repeat=d):
        c = mat_vec(Ginv, list(t))
        q = dot(c, list(t))                   # squared norm of v0
        if q > r2:
            continue
        v0 = [sum(c[i] * basis[i][k] for i in range(d))
              for k in range(len(basis[0]))]
        disc = r2 - q
        if disc == 0:
            cands = [Cand(v0, 0, Fraction(0), tuple(t))]
        else:
            cands = [Cand(v0, +1, disc, tuple(t)),
                     Cand(v0, -1, disc, tuple(t))]
        for cand in cands:
            if cand.key() not in seen:
                seen.add(cand.key())
                out.append(cand)
    return out


def check_basis_compatible(basis, r2):
    for a, b in itertools.combinations(basis, 2):
        assert dot(a, b) <= r2 / 2, "basis vectors not pairwise compatible"
    for a in basis:
        assert dot(a, a) == r2, "basis vector with wrong norm"


# ---------------------------------------------------------------------------
# exact max clique (branch and bound with greedy coloring bound)

def build_graph(cands, r2):
    n = len(cands)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if compatible(cands[i], cands[j], r2):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return adj


def max_cliques(adj, want_all=False):
    """Exact maximum clique via BnB with greedy coloring. Returns
    (size, list_of_cliques) — all maximum cliques if want_all."""
    n = len(adj)
    best = [0, []]

    def color_bound(cand_mask):
        """Greedy coloring upper bound + ordered vertex list."""
        order, bounds = [], []
        classes = []
        m = cand_mask
        while m:
            v = m.bit_length() - 1  # highest set bit
            m &= ~(1 << v)
            for ci, cmask in enumerate(classes):
                if not (adj[v] & cmask):
                    classes[ci] |= 1 << v
                    break
            else:
                classes.append(1 << v)
        for ci, cmask in enumerate(classes):
            mm = cmask
            while mm:
                v = mm.bit_length() - 1
                mm &= ~(1 << v)
                order.append(v)
                bounds.append(ci + 1)
        return order, bounds

    def expand(cur, cur_size, cand_mask):
        order, bounds = color_bound(cand_mask)
        for idx in range(len(order) - 1, -1, -1):
            v = order[idx]
            if cur_size + bounds[idx] < best[0] or (
                    not want_all and cur_size + bounds[idx] == best[0]):
                return
            new_mask = cand_mask & adj[v]
            cur.append(v)
            if cur_size + 1 > best[0]:
                best[0] = cur_size + 1
                best[1] = [cur[:]]
            elif cur_size + 1 == best[0] and want_all:
                if new_mask == 0:
                    best[1].append(cur[:])
            if new_mask:
                expand(cur, cur_size + 1, new_mask)
            cur.pop()
            cand_mask &= ~(1 << v)

    expand([], 0, (1 << n) - 1)
    # deduplicate (want_all can record subsets of later-found maxima)
    seen, out = set(), []
    for cl in best[1]:
        if len(cl) == best[0] and tuple(sorted(cl)) not in seen:
            seen.add(tuple(sorted(cl)))
            out.append(sorted(cl))
    return best[0], out


# ---------------------------------------------------------------------------

def frac_str(f):
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def cand_coord_strings(cand, r2):
    """Exact coordinate strings for a candidate (checker format)."""
    out = []
    for k, c in enumerate(cand.rat):
        out.append(frac_str(c))
    return out  # NB: only valid for eps == 0 candidates


def clique_to_config(name, basis, cands, clique, r2, normal=None):
    """Assemble basis + clique candidates into a checker-format config dict.
    Handles eps != 0 by emitting sqrt() expressions along `normal`:
    v = rat + eps*sqrt(disc/mu)*normal, mu = <normal,normal>."""
    vectors = [[frac_str(c) for c in b] for b in basis]
    for ci in clique:
        cand = cands[ci]
        if cand.eps == 0:
            vectors.append([frac_str(c) for c in cand.rat])
        else:
            assert normal is not None
            mu = dot(normal, normal)
            row = []
            for k in range(len(cand.rat)):
                base = frac_str(cand.rat[k])
                if normal[k] == 0:
                    row.append(base)
                else:
                    mult = frac_str(Fraction(cand.eps) * normal[k])
                    row.append(f"({base}) + ({mult})*sqrt(({frac_str(cand.disc)})/({frac_str(mu)}))")
            vectors.append(row)
    return {
        "name": name,
        "dimension": len(basis[0]),
        "n_points": len(vectors),
        "source": "search/clique.py output; verify with verify/exact_check.py",
        "norm_squared": frac_str(r2),
        "vectors": vectors,
    }
