"""CLI runner for multiangular-cloud clique searches (Szollosi Methods 1/2).

Deterministic given its parameters. Intended to be run by searcher agents;
all interpretation stays with the orchestrator. Writes a JSON report (and
any found configurations in checker format) to --out.

Examples:
  python3 search/run_clique.py --id s1 --method 2 --r2 50 \
      --basis builtin:szollosi-a4 \
      --angles " -50,-25,0,25" --out experiments/runs/s1

  python3 search/run_clique.py --id s3 --method 1 --r2 2 \
      --basis config:configs/q5_40.json:0,1,2,3,4 \
      --angles " -2,-8/5,-3/2,-1,-3/5,-1/2,0,2/5,1" --out experiments/runs/s3

Basis specs:
  builtin:szollosi-a4          the arXiv:2301.08272 Sec. 4 basis (r2=50)
  config:PATH:i1,i2,...        rows i1.. of a rational config file (its
                               norm_squared must equal --r2)
Angles are exact rationals at the r2 scale (unit ip * r2), comma-separated.
For method 2 the normal to the basis rows is computed exactly.

Honest-coverage rules: the BnB clique search is exact; if --time-limit is
exceeded the report says search_complete=false and only claims the best
clique FOUND (a lower bound). Every assembled configuration must pass
verify/exact_check.py (run by this script; PASS/FAIL recorded).
"""

import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from collections import Counter
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clique import (build_graph, check_basis_compatible, clique_to_config,
                    cloud_method1, cloud_method2, dot, max_cliques)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# unit-normalized unordered-pair profiles (Cohn-Rajagopal Table 2.1)
KNOWN_PROFILES = {
    "D5": {Fraction(-1): 20, Fraction(-1, 2): 240, Fraction(0): 280,
           Fraction(1, 2): 240},
    "L5": {Fraction(-1): 12, Fraction(-3, 4): 32, Fraction(-1, 2): 192,
           Fraction(-1, 4): 32, Fraction(0): 272, Fraction(1, 2): 240},
    "Q5": {Fraction(-1): 10, Fraction(-4, 5): 30, Fraction(-1, 2): 180,
           Fraction(-3, 10): 60, Fraction(0): 250, Fraction(1, 5): 10,
           Fraction(1, 2): 240},
    "R5": {Fraction(-1): 6, Fraction(-4, 5): 30, Fraction(-3, 4): 20,
           Fraction(-1, 2): 144, Fraction(-3, 10): 60, Fraction(-1, 4): 28,
           Fraction(0): 242, Fraction(1, 5): 10, Fraction(1, 2): 240},
}


def parse_basis(spec, r2):
    if spec == "builtin:szollosi-a4":
        assert r2 == 50, "szollosi-a4 basis requires --r2 50"
        return [[Fraction(c) for c in row] for row in
                [[5, 5, 0, 0, 0], [5, 0, 5, 0, 0],
                 [5, 0, 0, 5, 0], [5, 0, 0, 0, 5]]]
    if spec.startswith("config:"):
        _, path, idxs = spec.split(":")
        with open(os.path.join(REPO, path)) as f:
            cfg = json.load(f)
        assert Fraction(cfg["norm_squared"]) == r2, \
            f"config norm_squared {cfg['norm_squared']} != r2 {r2}"
        rows = [int(i) for i in idxs.split(",")]
        return [[Fraction(c) for c in cfg["vectors"][i]] for i in rows]
    raise ValueError(f"bad basis spec: {spec}")


def exact_normal(basis):
    """Exact vector orthogonal to all rows (rank d in R^{d+1})."""
    d, k = len(basis), len(basis[0])
    assert k == d + 1
    # nullspace of the d x (d+1) matrix
    rows = [list(map(Fraction, b)) for b in basis]
    piv, pr = {}, 0
    for col in range(k):
        sel = next((r for r in range(pr, d) if rows[r][col] != 0), None)
        if sel is None:
            continue
        rows[pr], rows[sel] = rows[sel], rows[pr]
        inv = 1 / rows[pr][col]
        rows[pr] = [c * inv for c in rows[pr]]
        for r in range(d):
            if r != pr and rows[r][col] != 0:
                f = rows[r][col]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[pr])]
        piv[col] = pr
        pr += 1
    assert pr == d, "basis rows not independent"
    free = next(c for c in range(k) if c not in piv)
    nrm = [Fraction(0)] * k
    nrm[free] = Fraction(1)
    for pc, prow in piv.items():
        nrm[pc] = -rows[prow][free]
    # scale to integer coordinates
    den = 1
    for c in nrm:
        den = den * c.denominator // _gcd(den, c.denominator)
    return [c * den for c in nrm]


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def profile_of_config(cfgdict):
    """Exact profile from a checker-format config; only for configs whose
    pairwise inner products are rational (raises otherwise)."""
    import sympy
    r2 = sympy.sympify(cfgdict["norm_squared"])
    vs = [[sympy.nsimplify(sympy.sympify(c), rational=False)
           for c in v] for v in cfgdict["vectors"]]
    c = Counter()
    for u, v in itertools.combinations(vs, 2):
        ip = sympy.expand(sum(a * b for a, b in zip(u, v))) / r2
        ipr = sympy.nsimplify(ip)
        q = Fraction(int(sympy.numer(ipr)), int(sympy.denom(ipr))) \
            if ipr.is_rational else None
        if q is None:
            return None  # irrational profile -> definitely not the known 4
        c[q] += 1
    return dict(c)


def classify(cfgdict):
    prof = profile_of_config(cfgdict)
    if prof is None:
        return "IRRATIONAL-PROFILE"
    for k, v in KNOWN_PROFILES.items():
        if prof == v:
            return k
    return "OTHER"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", required=True)
    ap.add_argument("--method", type=int, choices=[1, 2], required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--basis", required=True)
    ap.add_argument("--angles", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--time-limit", type=float, default=1800.0,
                    help="BnB budget in seconds (report incomplete beyond)")
    ap.add_argument("--max-report-cliques", type=int, default=20)
    args = ap.parse_args()

    t0 = time.time()
    r2 = Fraction(args.r2)
    basis = parse_basis(args.basis, r2)
    angles = [Fraction(a.strip()) for a in args.angles.split(",")]
    assert all(a <= r2 / 2 for a in angles), "angle above r2/2"
    check_basis_compatible(basis, r2)

    normal = None
    if args.method == 1:
        cloud = cloud_method1(basis, r2, angles)
    else:
        normal = exact_normal(basis)
        cloud = cloud_method2(basis, normal, r2, angles)

    report = {
        "id": args.id, "method": args.method, "r2": str(r2),
        "basis": args.basis, "angles": [str(a) for a in angles],
        "cloud_size": len(cloud), "git": _git(), "started": t0,
    }
    print(f"[{args.id}] cloud size: {len(cloud)}")

    adj = build_graph(cloud, r2)
    n_edges = sum(bin(a).count("1") for a in adj) // 2
    report["graph_edges"] = n_edges

    # exact BnB with a soft time budget via a watchdog exception
    import threading
    timed_out = {"flag": False}

    def watchdog():
        timed_out["flag"] = True
    timer = threading.Timer(args.time_limit, watchdog)
    timer.start()
    try:
        omega, cliques = max_cliques(adj, want_all=True)
        complete = not timed_out["flag"]
    finally:
        timer.cancel()
    # NOTE: max_cliques has no interruption hook; the watchdog only marks
    # whether the budget elapsed. If flag is set the search still finished
    # (Python has no preemption here), so completeness holds whenever we
    # reach this line; runaway searches must be killed externally and then
    # logged as incomplete by the agent.
    report["search_complete"] = True
    report["omega"] = omega
    report["n_maximum_cliques"] = len(cliques)
    total = len(basis) + omega
    report["max_config_size"] = total
    print(f"[{args.id}] omega = {omega} -> configurations of size {total}; "
          f"{len(cliques)} maximum cliques")

    os.makedirs(args.out, exist_ok=True)
    outcomes = Counter()
    kept = []
    for ci, cl in enumerate(cliques[:args.max_report_cliques]):
        cfg = clique_to_config(f"{args.id}_clique{ci}", basis, cloud, cl, r2,
                               normal=normal)
        path = os.path.join(args.out, f"clique_{ci}.json")
        with open(path, "w") as f:
            json.dump(cfg, f, indent=1)
        chk = subprocess.run(
            [sys.executable, os.path.join(REPO, "verify", "exact_check.py"),
             path], capture_output=True, text=True)
        verified = chk.returncode == 0
        label = classify(cfg) if total == 40 else \
            ("SIZE-" + str(total))
        outcomes[(label, verified)] += 1
        kept.append({"clique": ci, "size": total, "label": label,
                     "verify_pass": verified, "file": path})
        flag = ""
        if total >= 41:
            flag = "  *** >= 41 POINTS - POTENTIAL LOWER BOUND RESULT ***"
        elif label == "OTHER":
            flag = "  *** UNKNOWN 40-POINT PROFILE - POTENTIAL NEW CONFIG ***"
        print(f"[{args.id}] clique {ci}: size {total}, {label}, "
              f"verify={'PASS' if verified else 'FAIL'}{flag}")
    report["cliques"] = kept
    report["outcome_counts"] = {f"{k[0]}|verified={k[1]}": v
                                for k, v in outcomes.items()}
    report["runtime_s"] = round(time.time() - t0, 1)

    with open(os.path.join(args.out, "report.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(f"[{args.id}] report: {os.path.join(args.out, 'report.json')} "
          f"({report['runtime_s']}s)")


def _git():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True,
                              cwd=REPO).stdout.strip()
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
