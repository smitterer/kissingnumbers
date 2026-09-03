"""Rigorous branch-and-bound for the cap-layer capacity M(X) — exact
integer-lattice engine.

Goal: prove M(X) <= T-1 (no T pairwise-compatible points in the cap layer
over X), where the cap layer over X (40 points, squared norm r2 = 2) is
  C(X) = {(u, h) : u in Q_polar(X), h = sqrt(2 - |u|^2)},
compatibility: g(u,u') = <u,u'> + h h' <= 1.

Verified inputs used: max |u|^2 over Q_polar(X) = 5/4 (exact, independently
verified) — hence h >= sqrt(3/4) and the reduction lemma (blind-verified)
makes M(X) <= 15 for X in {Q5, R5} equivalent to Cohn-Rajagopal
Conjecture 3.1 (in full), and M(X) <= 12 the sharp target.

EXACTNESS: every verdict-relevant quantity is an integer computation.
Box bounds lie on the dyadic grid (9/S units with S = 2^j): coordinate
value = c * 9/S for integer c. All tests below are integer inequalities;
numpy int64 is exact as long as magnitudes stay < 2^63 (asserted: the
largest intermediate is 81*(dot+hh') <= ~8.4*S^2, requiring S <= 2^30;
refinement depth is capped accordingly). Square roots enter only via
integer isqrt with explicit floor direction (under-approximation of h,
which is the sound direction for a LOWER bound of g). No floats anywhere.

Soundness of the relaxation (necessary conditions only):
- a box is kept whenever it COULD intersect Q_polar (min over box of each
  constraint <= bound) — coverage is never lost;
- boxes A, B get an edge unless the exact lower bound
    dotmin(A,B) + hlo(A) hlo(B) > 1
  holds — so any T compatible points yield a T-clique of (distinct,
  because every box is refined until self-incompatible) boxes;
- k-core pruning removes only boxes with < T-1 potential neighbors.

Refinement strategy: after k-core pruning, ALL surviving boxes are split
(region refinement), so the active region shrinks geometrically; the
exact clique test runs each round and PROVED is declared when either the
core is empty or no T-clique exists.

Regression required: X = D5 with T = 17 must prove M(D5) <= 16 while
16-cliques keep surviving (E6 layer exists); same for L5.
"""

import argparse
import itertools
import json
import math
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0] + "/../analysis")  # moved to experiments/, imports stay in analysis/
from rigidity import load_rational_config

MAX_J = 30  # S = 2^j <= 2^30 keeps 81*(dot + h*h') < 2^63 (see module doc)


class Cover:
    """Box cover at scale S = 2^j. Bounds are int arrays lo, hi with
    coordinate value = c * 9 / S. All boxes share one scale."""

    def __init__(self, j, lo, hi):
        self.j = j
        self.S = 1 << j
        self.lo = lo.astype(np.int64)
        self.hi = hi.astype(np.int64)

    def n(self):
        return len(self.lo)


def initial_cover(init_split):
    """Uniform split of [-9/8, 9/8]^5: grid step (9/4)/init_split.
    Choose j so that step = 9/S * k with integer k: S = 4*init_split
    (then step corresponds to integer width 1... use width w = S//(4? )).
    Simplest: S = 2^j with 2^j = 4 * init_split (init_split power of 2);
    box integer bounds run over [-S/8, S/8] in steps of S/(4*init_split)."""
    assert init_split & (init_split - 1) == 0, "init_split must be 2^k"
    j = int(math.log2(4 * init_split))
    S = 1 << j
    w = S // (4 * init_split) * 2  # integer width per box; S/8*2/init...
    # bounds: from -S/8 to S/8 in init_split steps => width = (S/4)/init_split
    w = (S // 4) // init_split
    edges = [-(S // 8) + k * w for k in range(init_split + 1)]
    los, his = [], []
    for cell in itertools.product(range(init_split), repeat=5):
        los.append([edges[c] for c in cell])
        his.append([edges[c + 1] for c in cell])
    return Cover(j, np.array(los), np.array(his))


def feasible_mask(cov, Xint, xden):
    """Boxes that could intersect Q_polar. Xint = X * xden (int).
    Constraint <x, u> <= 1: min over box of 9*<Xint, c> <= xden*S."""
    lo, hi, S = cov.lo, cov.hi, cov.S
    keep = np.ones(cov.n(), dtype=bool)
    for x in Xint:
        m = np.zeros(cov.n(), dtype=np.int64)
        for k in range(5):
            m += np.where(x[k] >= 0, x[k] * lo[:, k], x[k] * hi[:, k])
        keep &= (9 * m) <= (xden * S)
    # norm: min |u|^2 over box <= 5/4  <=>  4*81*n2min <= 5*S^2
    n2min = np.zeros(cov.n(), dtype=np.int64)
    for k in range(5):
        a, b = lo[:, k], hi[:, k]
        c = np.where(a > 0, a, np.where(b < 0, b, 0))
        n2min += c * c
    keep &= (4 * 81 * n2min) <= (5 * S * S)
    return keep


def h_lo(cov):
    """Integer under-approx of h*S/9 per box: h^2 = 2 - min(n2max, 5/4).
    n2max units: (9/S)^2 * n2max_int; h^2*(S/9)^2 = (2*S^2 - 81*n2)/81
    capped: n2_eff = min(81*n2max_int, ceil? use max with 3/4 floor):
    2 - cap(...) >= 3/4 always. hlo = isqrt((2*S^2 - n2cap)//81) is a
    floor -> under-approximation. Sound (we need h lower bounds)."""
    lo, hi, S = cov.lo, cov.hi, cov.S
    n2max = np.zeros(cov.n(), dtype=np.int64)
    for k in range(5):
        n2max += np.maximum(lo[:, k] ** 2, hi[:, k] ** 2)
    n2cap = np.minimum(81 * n2max, 5 * S * S // 4)  # floor cap is sound:
    # capping n2 DOWN can only increase 2 - n2, i.e. raise hlo — check:
    # points in P satisfy 81*|c|^2 <= 5S^2/4 exactly? |u|^2 <= 5/4 means
    # 81*n2 <= 5S^2/4; integer floor of the cap keeps hlo an
    # under-approx of sqrt(2 - min(|u|^2, 5/4)) since we never lower n2
    # below the true bound: floor(5S^2/4) >= any valid 81*n2? No —
    # SOUNDNESS: h = sqrt(2 - |u|^2) with |u|^2 <= min(n2max_val, 5/4).
    # We need hlo <= min over box of h = sqrt(2 - max feasible |u|^2)
    # where max feasible |u|^2 <= min(n2max_val, 5/4). Using the cap
    # min(81*n2max, floor(5S^2/4)) >= 81 * (max feasible |u|^2 scaled)
    # ... floor only shrinks the cap, shrinking 2-n2's subtrahend and
    # RAISING hlo, which would be UNSOUND if the true |u|^2 could exceed
    # floor(5S^2/4)/81. Since 81*|u_int-units|^2 is an integer and
    # |u|^2 <= 5/4 exactly, 81*n2_true <= 5S^2/4, and as an integer
    # 81*n2_true <= floor(5S^2/4). Sound.
    hsq = (2 * S * S - n2cap) // 81
    return np.array([math.isqrt(int(v)) for v in hsq], dtype=np.int64)


def pairwise_edges(cov, hlo, block=2048):
    """Boolean adjacency: edge unless provably incompatible:
    81*(dotmin + hlo_i*hlo_j) > S^2  => no edge.
    dotmin per coord: min of the four products of interval endpoints."""
    n, S = cov.n(), cov.S
    lo, hi = cov.lo, cov.hi
    adj = np.zeros((n, n), dtype=bool)
    S2 = np.int64(S) * np.int64(S)
    for i0 in range(0, n, block):
        i1 = min(n, i0 + block)
        dm = np.zeros((i1 - i0, n), dtype=np.int64)
        for k in range(5):
            a0 = lo[i0:i1, k, None]
            a1 = hi[i0:i1, k, None]
            b0 = lo[None, :, k]
            b1 = hi[None, :, k]
            p = np.minimum(np.minimum(a0 * b0, a0 * b1),
                           np.minimum(a1 * b0, a1 * b1))
            dm += p
        g = 81 * (dm + hlo[i0:i1, None] * hlo[None, :])
        adj[i0:i1] = g <= S2
    return adj


def self_ok_mask(cov, hlo):
    """Boxes provably unable to host two points: 81*(dotmin(B,B)+hlo^2) > S^2."""
    n, S = cov.n(), cov.S
    lo, hi = cov.lo, cov.hi
    dm = np.zeros(n, dtype=np.int64)
    for k in range(5):
        a, b = lo[:, k], hi[:, k]
        dm += np.minimum(np.minimum(a * a, a * b), b * b)
    return 81 * (dm + hlo * hlo) > np.int64(S) * np.int64(S)


def split_all(cov, mask=None):
    """Split every (masked) box along its widest coordinate; scale doubles
    if needed (widths become odd)."""
    lo, hi, j = cov.lo, cov.hi, cov.j
    if mask is None:
        mask = np.ones(cov.n(), dtype=bool)
    widths = hi - lo
    wk = np.argmax(widths, axis=1)
    need_double = ((hi[np.arange(cov.n()), wk] +
                    lo[np.arange(cov.n()), wk]) % 2 != 0).any()
    if need_double or True:
        # always double the scale to keep midpoints integral
        assert j + 1 <= MAX_J, "refinement depth exceeds int64-safe scale"
        lo, hi, j = lo * 2, hi * 2, j + 1
    los, his = [], []
    for i in range(cov.n()):
        if not mask[i]:
            los.append(lo[i]); his.append(hi[i])
            continue
        k = wk[i]
        mid = (lo[i, k] + hi[i, k]) // 2
        l1, h1 = lo[i].copy(), hi[i].copy()
        l2, h2 = lo[i].copy(), hi[i].copy()
        h1[k] = mid
        l2[k] = mid
        los.extend([l1, l2]); his.extend([h1, h2])
    return Cover(j, np.array(los), np.array(his))


def kcore(adj, T):
    """Iterated (T-1)-degree pruning. Returns kept index array."""
    alive = np.ones(len(adj), dtype=bool)
    while True:
        deg = (adj & alive[None, :]).sum(axis=1) - 1  # exclude self if set
        drop = alive & (deg < (T - 1))
        if not drop.any():
            return np.nonzero(alive)[0]
        alive &= ~drop


def clique_geq(adj_bool, T, budget=5_000_000):
    """Exact: contains a clique of size >= T? bitset BnB with coloring."""
    n = len(adj_bool)
    rows = []
    for i in range(n):
        r = 0
        for j in np.nonzero(adj_bool[i])[0]:
            if j != i:
                r |= 1 << int(j)
        rows.append(r)
    nodes = [0]
    res = {"found": False, "witness": None, "over": False}

    def color_bound(cand):
        order, bounds, classes = [], [], []
        m = cand
        while m:
            v = m.bit_length() - 1
            m &= ~(1 << v)
            for ci, cm in enumerate(classes):
                if not (rows[v] & cm):
                    classes[ci] |= 1 << v
                    break
            else:
                classes.append(1 << v)
        for ci, cm in enumerate(classes):
            mm = cm
            while mm:
                v = mm.bit_length() - 1
                mm &= ~(1 << v)
                order.append(v)
                bounds.append(ci + 1)
        return order, bounds

    def expand(cur, cand):
        nodes[0] += 1
        if nodes[0] > budget:
            res["over"] = True
            return True
        if len(cur) >= T:
            res["found"], res["witness"] = True, cur[:]
            return True
        order, bounds = color_bound(cand)
        for idx in range(len(order) - 1, -1, -1):
            if len(cur) + bounds[idx] < T:
                return False
            v = order[idx]
            if expand(cur + [v], cand & rows[v]):
                return True
            cand &= ~(1 << v)
        return False

    expand([], (1 << n) - 1)
    if res["over"]:
        return None, None
    return res["found"], res["witness"]


def prove(config_path, T, init_split=8, max_rounds=40, max_boxes=400_000,
          clique_budget=800_000):
    name, vecs, r2, npts, d = load_rational_config(config_path)
    assert r2 == Fraction(2)
    xden = 1
    for v in vecs:
        for c in v:
            xden = xden * c.denominator // math.gcd(xden, c.denominator)
    Xint = np.array([[int(c * xden) for c in v] for v in vecs],
                    dtype=np.int64)

    cov = initial_cover(init_split)
    keep = feasible_mask(cov, Xint, xden)
    cov = Cover(cov.j, cov.lo[keep], cov.hi[keep])
    print(f"[{name}] T={T}: initial {cov.n()} feasible boxes at scale 2^{cov.j}", flush=True)
    t0 = time.time()

    for rnd in range(1, max_rounds + 1):
        # ensure all boxes self-incompatible
        while True:
            hlo = h_lo(cov)
            ok = self_ok_mask(cov, hlo)
            if ok.all():
                break
            cov = split_all(cov, mask=~ok)
            keep = feasible_mask(cov, Xint, xden)
            cov = Cover(cov.j, cov.lo[keep], cov.hi[keep])
            if cov.n() > max_boxes:
                return {"name": name, "T": T, "verdict": "ABORT-boxcount",
                        "boxes": int(cov.n())}
        hlo = h_lo(cov)
        adj = pairwise_edges(cov, hlo)
        np.fill_diagonal(adj, False)
        kept = kcore(adj, T)
        print(f"[{name}] T={T} round {rnd}: kcore {len(kept)}/{cov.n()} at 2^{cov.j}", flush=True)
        if len(kept) == 0:
            dt = time.time() - t0
            print(f"[{name}] T={T}: PROVED (empty core) round {rnd} [{dt:.0f}s]")
            return {"name": name, "T": T, "verdict": "PROVED",
                    "rounds": rnd, "seconds": round(dt, 1), "j": cov.j}
        cov = Cover(cov.j, cov.lo[kept], cov.hi[kept])
        adj = adj[np.ix_(kept, kept)]
        found, wit = clique_geq(adj, T, budget=clique_budget)
        dt = time.time() - t0
        if found is None:
            print(f"[{name}] T={T} round {rnd}: clique budget exceeded on "
                  f"{cov.n()} boxes [{dt:.0f}s] — refining anyway", flush=True)
        elif not found:
            print(f"[{name}] T={T}: PROVED (no {T}-clique, {cov.n()} boxes) "
                  f"round {rnd} [{dt:.0f}s]")
            return {"name": name, "T": T, "verdict": "PROVED",
                    "rounds": rnd, "boxes": int(cov.n()),
                    "seconds": round(dt, 1), "j": cov.j}
        else:
            print(f"[{name}] T={T} round {rnd}: {cov.n()} boxes at 2^{cov.j}, "
                  f"{T}-clique survives [{dt:.0f}s]", flush=True)
        # region refinement: split every surviving box
        cov = split_all(cov)
        keep = feasible_mask(cov, Xint, xden)
        cov = Cover(cov.j, cov.lo[keep], cov.hi[keep])
        if cov.n() > max_boxes:
            return {"name": name, "T": T, "verdict": "ABORT-boxcount",
                    "boxes": int(cov.n()), "rounds": rnd}
    return {"name": name, "T": T, "verdict": "INCONCLUSIVE-maxrounds",
            "boxes": int(cov.n())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--T", type=int, required=True)
    ap.add_argument("--init-split", type=int, default=8)
    ap.add_argument("--max-rounds", type=int, default=40)
    ap.add_argument("--max-boxes", type=int, default=400_000)
    args = ap.parse_args()
    res = prove(args.config, args.T, init_split=args.init_split,
                max_rounds=args.max_rounds, max_boxes=args.max_boxes)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
