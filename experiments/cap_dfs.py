"""Tuple branch-and-prune prover for the cap-layer capacity M(X).

Proves M(X) <= T-1: no T points y_i = (u_i, sqrt(2-|u_i|^2)), u_i in
P = Q_polar(X), with pairwise <u_i,u_j> + h_i h_j <= 1.

Architecture (third design; see cap_bnb.py header for exactness
conventions — unchanged: dyadic integer boxes, integer interval tests,
integer isqrt height bounds, no floats on the verdict path):

STATE = a multiset of T boxes ("slots"), kept as a lexicographically
sorted tuple. Invariant: any T-point configuration whose points are
covered by the slots' boxes remains covered in exactly one child branch.

Expansion: pick the physically widest slot, split it along its widest
coordinate into c1, c2 (children filtered by the polytope feasibility
test); branch into (slots \\ {b}) ∪ {c1} and ∪ {c2} and — because several
POINTS may lie in the SAME slot box — also ∪ {c1, c2} \\ {b, some other
occurrence?}  -- NO: to keep the multiset semantics simple and sound we
give each POINT its own slot from the start (T slots all equal to a root
box covering P). Splitting a slot sends ITS point into c1 or c2: exactly
two children. Points in other slots are unaffected. Sorting + memo on
the sorted tuple collapses permutation-equivalent states (sound because
slots are interchangeable).

Pruning (all sound):
- pair prune: if some slot pair (including equal boxes, via the
  independent-points dotmin) is provably incompatible, the branch dies;
- feasibility prune: a child box that provably misses P kills that child
  (its point has nowhere to go IN THAT BRANCH; the sibling branch covers
  the rest);
- memoization: sorted state seen before => skip (bounded LRU).

Termination of a branch at the scale floor (all slots tiny, all pairs
still possibly compatible) => near-feasible witness: verdict
INCONCLUSIVE with slot centers reported. Full search with no witness =>
PROVED.

Regression: D5/T=17 must be PROVED; D5/T=16 must return a witness near
the E6 layer.
"""

import argparse
import itertools
import json
import math
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0] + "/../analysis")  # moved to experiments/, imports stay in analysis/
from rigidity import load_rational_config

MAX_J = 30


def rescale(box, j):
    bj, lo, hi = box
    f = 1 << (j - bj)
    return (j, tuple(c * f for c in lo), tuple(c * f for c in hi))


def dotmin_pair(A, B):
    _, alo, ahi = A
    _, blo, bhi = B
    s = 0
    for a0, a1, b0, b1 in zip(alo, ahi, blo, bhi):
        s += min(a0 * b0, a0 * b1, a1 * b0, a1 * b1)
    return s


def h_lo(box):
    j, lo, hi = box
    S = 1 << j
    n2 = sum(max(a * a, b * b) for a, b in zip(lo, hi))
    n2cap = min(81 * n2, 5 * S * S // 4)
    return math.isqrt((2 * S * S - n2cap) // 81)


def compatible(boxA, boxB):
    j = max(boxA[0], boxB[0])
    A = rescale(boxA, j) if boxA[0] < j else boxA
    B = rescale(boxB, j) if boxB[0] < j else boxB
    S = 1 << j
    return 81 * (dotmin_pair(A, B) + h_lo(A) * h_lo(B)) <= S * S


def feasible(box, Xint, xden):
    j, lo, hi = box
    S = 1 << j
    n2min = 0
    for a, b in zip(lo, hi):
        c = a if a > 0 else (b if b < 0 else 0)
        n2min += c * c
    if 4 * 81 * n2min > 5 * S * S:
        return False
    for x in Xint:
        m = 0
        for xk, a, b in zip(x, lo, hi):
            m += xk * a if xk >= 0 else xk * b
        if 9 * m > xden * S:
            return False
    return True


def width_phys(box):
    j, lo, hi = box
    return max(h - l for l, h in zip(lo, hi)) / (1 << j)


def split(box):
    j, lo, hi = box
    if any((h - l) < 2 for l, h in zip(lo, hi)):
        box = rescale(box, j + 1)
        j, lo, hi = box
        assert j <= MAX_J, "scale floor exceeded"
    widths = [h - l for l, h in zip(lo, hi)]
    k = widths.index(max(widths))
    mid = (lo[k] + hi[k]) // 2
    lo1, hi1 = list(lo), list(hi)
    lo2, hi2 = list(lo), list(hi)
    hi1[k] = mid
    lo2[k] = mid
    return (j, tuple(lo1), tuple(hi1)), (j, tuple(lo2), tuple(hi2))


def prove(config_path, T, jfloor_width=2e-4, node_budget=500_000_000,
          memo_cap=4_000_000, progress_every=60.0):
    name, vecs, r2, npts, d = load_rational_config(config_path)
    assert r2 == Fraction(2)
    xden = 1
    for v in vecs:
        for c in v:
            xden = xden * c.denominator // math.gcd(xden, c.denominator)
    Xint = [[int(c * xden) for c in v] for v in vecs]

    root = (3, tuple([-1] * 5), tuple([1] * 5))  # [-9/8, 9/8]^5 at j=3
    state0 = tuple([root] * T)
    stats = {"nodes": 0, "pruned_pair": 0, "pruned_memo": 0,
             "t0": time.time(), "last": time.time(), "maxdepth": 0}
    memo = set()
    witness = []

    def pair_ok(state):
        for a, b in itertools.combinations(state, 2):
            if not compatible(a, b):
                return False
        return True

    def dfs(state, depth):
        stats["nodes"] += 1
        stats["maxdepth"] = max(stats["maxdepth"], depth)
        if stats["nodes"] > node_budget:
            raise RuntimeError("node-budget")
        now = time.time()
        if now - stats["last"] > progress_every:
            stats["last"] = now
            widest = max(width_phys(b) for b in state)
            print(f"[{name}] T={T} nodes={stats['nodes']} depth={depth} "
                  f"widest={widest*9:.4f} pairprunes={stats['pruned_pair']} "
                  f"memo={stats['pruned_memo']} [{now-stats['t0']:.0f}s]",
                  flush=True)
        # widest slot
        wi = max(range(T), key=lambda i: width_phys(state[i]))
        if width_phys(state[wi]) * 9 <= jfloor_width:
            centers = [[9.0 * (l + h) / 2 / (1 << j)
                        for l, h in zip(lo_, hi_)]
                       for (j, lo_, hi_) in state]
            witness.append(centers)
            return True
        c1, c2 = split(state[wi])
        rest = state[:wi] + state[wi + 1:]
        for c in (c1, c2):
            if not feasible(c, Xint, xden):
                continue
            ok = True
            for b in rest:
                if not compatible(c, b):
                    stats["pruned_pair"] += 1
                    ok = False
                    break
            if ok and not compatible(c, c):
                # a slot must host one point; self-compat not required
                pass
            if not ok:
                continue
            child = tuple(sorted(rest + (c,)))
            if child in memo:
                stats["pruned_memo"] += 1
                continue
            if len(memo) < memo_cap:
                memo.add(child)
            if dfs(child, depth + 1):
                return True
        return False

    try:
        found = dfs(state0, 0)
    except RuntimeError as e:
        return {"name": name, "T": T, "verdict": f"ABORT-{e}",
                "nodes": stats["nodes"]}
    dt = time.time() - stats["t0"]
    out = {"name": name, "T": T, "nodes": stats["nodes"],
           "pruned_pair": stats["pruned_pair"],
           "pruned_memo": stats["pruned_memo"],
           "maxdepth": stats["maxdepth"], "seconds": round(dt, 1)}
    if not found:
        print(f"[{name}] T={T}: PROVED M <= {T-1}  ({stats['nodes']} nodes, "
              f"{dt:.0f}s)", flush=True)
        out["verdict"] = "PROVED"
    else:
        print(f"[{name}] T={T}: INCONCLUSIVE — near-feasible tuple at scale "
              f"floor ({dt:.0f}s)", flush=True)
        out["verdict"] = "INCONCLUSIVE-witness"
        out["witness_centers"] = witness[0]
        with open(f"analysis/dfs_witness_{name}_T{T}.json", "w") as f:
            json.dump(out, f, indent=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--T", type=int, required=True)
    ap.add_argument("--jfloor-width", type=float, default=2e-4)
    args = ap.parse_args()
    res = prove(args.config, args.T, jfloor_width=args.jfloor_width)
    print(json.dumps({k: v for k, v in res.items()
                      if k != "witness_centers"}))


if __name__ == "__main__":
    main()
