#!/usr/bin/env python3
"""Single entry point for the verification behind

    Rigidity of the known five-dimensional kissing configurations
    and a proof of a conjecture of Cohn and Rajagopal

Runs, in order, every script that backs a computational claim of the paper
(Theorems 1, 2, 6 and Lemmas 4, 8), parses each script's stdout for the
values in EXPECTED below, and asserts equality.  Any mismatch, or a nonzero
exit code of any script, is printed and makes this script exit 1; otherwise
the last line printed is

    ALL CHECKS PASS

Outputs (both regenerated on every run, both git-ignored):
    results/summary.json   every check with expected/computed/pass, and the
                           exit code and wall time of every step
    results/log.txt        full stdout and stderr of every script

Nothing under version control is modified: analysis/deep_holes.py and
analysis/stress_lp.py run without --save, and search/run_szollosi_repro.py
rewrites experiments/szollosi_repro/*.json byte-identically.

Deliberately NOT run here: search/run_clique.py (the discovery search,
10-70 minutes, not needed for any claim) and everything under experiments/
(floating point, exploratory, no expected output).

Runtime: well under two minutes on a laptop.
"""

from __future__ import annotations

import ast
import datetime
import json
import os
import platform
import re
import subprocess
import sys
import time
from fractions import Fraction

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
RESULTS = os.path.join(ROOT, "results")

FOUR = ["configs/%s_40.json" % n for n in ("d5", "l5", "q5", "r5")]
SEVEN = FOUR + ["configs/q5_szollosi_40.json",
                "configs/q5_ext64.json", "configs/r5_ext64.json"]

# ---------------------------------------------------------------------------
# The claims.  Hard-coded; a change here is a change of the paper.
# ---------------------------------------------------------------------------
EXPECTED = {
    # Theorem 1 (jamming): contact graph, stress certificates, equality rank
    "contacts":       {"d5": 240, "l5": 240, "q5": 240, "r5": 240},
    "contact_degree": {"d5": 12, "l5": 12, "q5": 12, "r5": 12},
    "antipodal":      {"D5": 40, "L5": 24, "Q5": 20, "R5": 12},
    "rank":           {"d5": 190, "l5": 190, "q5": 190, "r5": 190},
    "kernel_dim":     {"d5": 10, "l5": 10, "q5": 10, "r5": 10},
    "uniform_stress_fails": {"d5": False, "l5": False, "q5": True, "r5": True},
    "stress_omega": {
        "d5": {"1": 240},
        "l5": {"1": 240},
        "q5": {"20/21": 210, "4/3": 30},
        "r5": {"200/207": 210, "280/207": 18, "220/207": 12},
    },
    "stress_lambda": {
        "d5": {"6": 40},
        "l5": {"6": 40},
        "q5": {"40/7": 20, "44/7": 20},
        "r5": {"400/69": 20, "140/23": 12, "440/69": 8},
    },
    "r5_middle_weight": "220/207",
    # Theorem 2 (deep holes, no 41st point)
    "vertices":    {"d5": 42, "l5": 50, "q5": 92, "r5": 100},
    "max_w2":      {"d5": "5/4", "l5": "5/4", "q5": "5/4", "r5": "5/4"},
    "deep_holes":  {"d5": 32, "l5": 32, "q5": 32, "r5": 32},
    "max_hole_ip": {"D5": "3/4", "L5": "3/4", "Q5": "23/20", "R5": "23/20"},
    # Lemma 4 (hole-graph clique numbers)
    "clique": {"D5": 16, "L5": 16, "Q5": 12, "R5": 12},
    # Lemma 8 (Gegenbauer identity)
    "lemma8_identity": True,
    "lemma8_F1": 16,
    # Theorem 6 (64-point extensions)
    "ext64": {"Q5": "PASS, 64 points", "R5": "PASS, 64 points"},
    # Szollosi reproduction (validation of the search machinery)
    "szollosi": {"cloud": 78, "max_clique": 36, "n_max_cliques": 4,
                 "classification": {"Q5": 2, "D5": 2}},
}

CHECKS = []      # dicts: step, check, expected, computed, pass
STEPS = []       # dicts: step, command, exit_code, wall_time_s


def add(step, check, expected, computed):
    CHECKS.append({"step": step, "check": check, "expected": expected,
                   "computed": computed, "pass": expected == computed})


def grab(pattern, text, conv=str, default=None, flags=0):
    m = re.search(pattern, text, flags)
    if not m:
        return default
    try:
        return conv(m.group(1))
    except (ValueError, SyntaxError):
        return m.group(1)


def blocks(text, marker):
    """Split stdout into {name: text} at lines matching `marker` (1 group)."""
    out, cur = {}, None
    for line in text.splitlines():
        m = re.match(marker, line)
        if m:
            cur = m.group(1)
            out[cur] = []
        elif cur is not None:
            out[cur].append(line)
    return {k: "\n".join(v) for k, v in out.items()}


def mult_dict(line):
    """'20/21 x210, 4/3 x30' -> {'20/21': 210, '4/3': 30}."""
    return {v: int(c) for v, c in re.findall(r"(\S+) x(\d+)", line or "")}


def run(label, argv, log):
    cmd = [PY] + argv
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    wall = round(time.time() - t0, 2)
    STEPS.append({"step": label, "command": " ".join(argv),
                  "exit_code": proc.returncode, "wall_time_s": wall})
    log.write("=" * 78 + "\n")
    log.write("STEP %s\n$ %s\nexit code %d, wall time %.2fs\n"
              % (label, " ".join(cmd), proc.returncode, wall))
    log.write("-" * 78 + "\n" + proc.stdout)
    if proc.stderr:
        log.write("\n--- stderr ---\n" + proc.stderr)
    log.write("\n\n")
    print("  [%5.1fs] %-40s exit %d" % (wall, label, proc.returncode))
    add(label, "exit code", 0, proc.returncode)
    return proc.stdout


# ---------------------------------------------------------------------------
# steps
# ---------------------------------------------------------------------------

def step_pytest(log):
    label = "pytest verify/"
    out = run(label, ["-m", "pytest", "verify/", "-q"], log)
    failed = (grab(r"(\d+) failed", out, int, 0) or 0) + \
             (grab(r"(\d+) error", out, int, 0) or 0)
    passed = grab(r"(\d+) passed", out, int, 0)
    add(label, "tests failed (%s passed)" % passed, 0, failed)


def step_exact_check(log):
    label = "verify/exact_check.py x7"
    out = run(label, ["verify/exact_check.py"] + SEVEN, log)
    # each report: "Configuration : <name>" / "Source : <path>" / ... / "VERDICT: X"
    verdicts = dict(re.findall(r"Source\s+: (\S+)[\s\S]*?VERDICT: (\w+)", out))
    for path in SEVEN:
        add(label, "%s verdict" % os.path.basename(path), "PASS",
            verdicts.get(path, "missing"))


def step_rigidity(log):
    label = "analysis/rigidity.py"
    out = run(label, ["analysis/rigidity.py"] + FOUR, log)
    bl = blocks(out, r"^== (\w+) ==$")
    for n in ("d5", "l5", "q5", "r5"):
        b = bl.get(n, "")
        add(label, "%s contacts" % n, EXPECTED["contacts"][n],
            grab(r"contacts \(ip = r2/2\): (\d+)", b, int))
        add(label, "%s uniform stress fails" % n,
            EXPECTED["uniform_stress_fails"][n],
            {"YES": False, "NO": True}.get(grab(r"uniform stress .*?: (YES|NO)", b)))
        add(label, "%s rotation flexes valid, dim" % n, 10,
            grab(r"rotation flexes valid\+independent: True \(dim (\d+)\)", b, int))
        ranks = grab(r"rank of equality system mod p: (\{[^}]*\})", b,
                     ast.literal_eval, {})
        vals = sorted(set(ranks.values())) if isinstance(ranks, dict) else []
        add(label, "%s rank of equality system" % n, EXPECTED["rank"][n],
            vals[0] if len(vals) == 1 else str(ranks))


def step_stress_lp(log):
    label = "analysis/stress_lp.py"
    out = run(label, ["analysis/stress_lp.py"] + FOUR, log)
    bl = blocks(out, r"^== (\w+) ==$")
    for n in ("d5", "l5", "q5", "r5"):
        b = bl.get(n, "")
        exp_min = str(min(Fraction(k) for k in EXPECTED["stress_omega"][n]))
        add(label, "%s exact positive self-stress, min weight" % n, exp_min,
            grab(r"self-stress verified \(min weight (\S+) =", b))
        add(label, "%s verdict" % n, "jammed",
            "jammed" if re.search(r"VERDICT: .*INFINITESIMALLY JAMMED", b)
            else grab(r"VERDICT: (.*)", b, default="missing"))


def step_independent_check(log):
    label = "analysis/certificates/independent_check.py"
    out = run(label, ["analysis/certificates/independent_check.py"], log)
    bl = blocks(out, r"^CONFIG (\w+) ")
    for n in ("d5", "l5", "q5", "r5"):
        b = bl.get(n.upper(), "")
        add(label, "%s contacts" % n, EXPECTED["contacts"][n],
            grab(r"contacts \(<x_i,x_j> = r2/2 exactly\) : (\d+)", b, int))
        m = re.search(r"contact-graph degrees\s+: min (\d+), max (\d+)", b)
        deg = (int(m.group(1)) if m and m.group(1) == m.group(2)
               else (m.group(0) if m else None))
        add(label, "%s contact degree" % n, EXPECTED["contact_degree"][n], deg)
        om = mult_dict(grab(r"omega multiplicities\s*: (.+)", b))
        lm = mult_dict(grab(r"lambda multiplicities\s*: (.+)", b))
        add(label, "%s stress weight multiplicities" % n,
            EXPECTED["stress_omega"][n], om)
        add(label, "%s lambda multiplicities" % n,
            EXPECTED["stress_lambda"][n], lm)
        if n == "r5":
            keys = sorted(om, key=Fraction)
            add(label, "r5 middle stress weight", EXPECTED["r5_middle_weight"],
                keys[1] if len(keys) == 3 else str(keys))
        add(label, "%s rank over Q (Bareiss)" % n, EXPECTED["rank"][n],
            grab(r"rank over Q \(Bareiss, exact integer\) = (\d+)", b, int))
        add(label, "%s equality kernel dim" % n, EXPECTED["kernel_dim"][n],
            grab(r"equality kernel has dim \d+ - \d+ = (\d+)", b, int))
        add(label, "%s result" % n, "jammed",
            "jammed" if "RESULT: INFINITESIMALLY JAMMED" in b
            else grab(r"RESULT: (.*)", b, default="missing"))
    add(label, "OVERALL", "PASS", grab(r"OVERALL: (\w+)", out))


def step_deep_holes(log):
    label = "analysis/deep_holes.py"
    out = run(label, ["analysis/deep_holes.py"] + FOUR, log)
    bl = blocks(out, r"^== (\w+) ==$")
    for n in ("d5", "l5", "q5", "r5"):
        b = bl.get(n, "")
        add(label, "%s polar vertices" % n, EXPECTED["vertices"][n],
            grab(r"vertices \(exact, deduped\): (\d+)", b, int))
        add(label, "%s max |w|^2" % n, EXPECTED["max_w2"][n],
            grab(r"max \|w\|\^2 over Q = (\S+)", b))
        add(label, "%s deep holes" % n, EXPECTED["deep_holes"][n],
            grab(r"deep holes: (\d+),", b, int))
        add(label, "%s 41st point" % n, "IMPOSSIBLE",
            grab(r"41st point without deforming: (\w+)", b))
    add(label, "no file written", True, "saved" not in out)


def step_independent_deep_holes(log):
    label = "analysis/certificates/independent_deep_holes_check.py"
    out = run(label, ["analysis/certificates/independent_deep_holes_check.py"], log)
    pat = (r"^\s+(d5|l5|q5|r5)\s+: vertices\s+(\d+), deep holes (\d+), nearest .*?, "
           r"max\|w\|\^2 (\S+), m\^2 (\S+)\s+-> (PASS|FAIL)")
    rows = {m[0]: m for m in re.findall(pat, out, re.M)}
    for n in ("d5", "l5", "q5", "r5"):
        r = rows.get(n)
        add(label, "%s polar vertices" % n, EXPECTED["vertices"][n],
            int(r[1]) if r else None)
        add(label, "%s deep holes" % n, EXPECTED["deep_holes"][n],
            int(r[2]) if r else None)
        add(label, "%s max |w|^2" % n, EXPECTED["max_w2"][n], r[3] if r else None)
        add(label, "%s m(X)^2" % n, "2/5", r[4] if r else None)
        add(label, "%s verdict" % n, "PASS", r[5] if r else None)
    add(label, "OVERALL", "PASS", grab(r"OVERALL: (\w+)", out))


def step_lemma4_lemma8(log):
    label = "analysis/verify_lemma4_lemma8.py"
    out = run(label, ["analysis/verify_lemma4_lemma8.py"], log)
    bl = blocks(out, r"^== (.+?) ==$")
    line = re.compile(r"^\s+(.+?)\s+computed=\s*(.+?)\s+expected=\s*(.+?)\s+(PASS|FAIL)\s*$")

    def table(b):
        return {m.group(1): (m.group(2), m.group(4))
                for m in (line.match(l) for l in b.splitlines()) if m}

    def val(t, key, conv=str):
        if key not in t:
            return None
        try:
            return conv(t[key][0])
        except ValueError:
            return t[key][0]

    for n in ("D5", "L5", "Q5", "R5"):
        t = table(bl.get(n, ""))
        add(label, "%s antipodal points" % n, EXPECTED["antipodal"][n],
            val(t, "antipodal points", int))
        add(label, "%s polytope vertices" % n, EXPECTED["vertices"][n.lower()],
            val(t, "polytope vertices", int))
        add(label, "%s max |w|^2" % n, "5/4", val(t, "max |w|^2"))
        add(label, "%s deep holes" % n, 32, val(t, "deep holes", int))
        add(label, "%s max hole inner product" % n, EXPECTED["max_hole_ip"][n],
            val(t, "max hole inner product"))
        add(label, "%s clique number (networkx)" % n, EXPECTED["clique"][n],
            val(t, "clique number (networkx)", int))
        add(label, "%s clique number (own B&B)" % n, EXPECTED["clique"][n],
            val(t, "clique number (own B&B)", int))
    t = table(bl.get("Lemma 8", ""))
    add(label, "Lemma 8 identity F = 125/16 (t+3/5)^2 (t-1/5)",
        EXPECTED["lemma8_identity"],
        val(t, "identity F = 125/16 (t+3/5)^2 (t-1/5)") == "True")
    add(label, "Lemma 8 F(1)", EXPECTED["lemma8_F1"], val(t, "F(1)", int))
    t = table(bl.get("Theorem 6", ""))
    for n in ("Q5", "R5"):
        key = "%s ext64 valid (64 pts, norms, ips)" % n
        got = t.get(key)
        add(label, "%s ext64" % n, EXPECTED["ext64"][n],
            "PASS, 64 points" if got == ("True", "PASS") else str(got))
    # the config-file comparison (Step 1.2 of the submission prep)
    t = table(bl.get("configs/*_40.json vs. rebuilt configurations", ""))
    for n in ("d5", "l5", "q5", "r5"):
        key = "%s_40.json == rebuilt %s (point sets)" % (n, n.upper())
        add("config files vs. rebuilt configurations", key, True,
            t.get(key) == ("True", "PASS"))
    add(label, "OVERALL", "ALL CHECKS PASS", grab(r"OVERALL: (.*)", out))


def step_szollosi(log):
    label = "search/run_szollosi_repro.py"
    out = run(label, ["search/run_szollosi_repro.py"], log)
    e = EXPECTED["szollosi"]
    add(label, "cloud size", e["cloud"], grab(r"cloud size: (\d+)", out, int))
    add(label, "max clique", e["max_clique"], grab(r"max clique: (\d+)", out, int))
    add(label, "number of maximum cliques", e["n_max_cliques"],
        grab(r"#maximum cliques: (\d+)", out, int))
    add(label, "classification", e["classification"],
        grab(r"classification: (\{[^}]*\})", out, ast.literal_eval))


# ---------------------------------------------------------------------------

def git_head():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


def main():
    os.makedirs(RESULTS, exist_ok=True)
    t_start = time.time()
    print("run_all.py: verification of the kissing-configuration paper "
          "(git %s, %s, %s)" % (git_head(), platform.python_implementation()
                                + " " + platform.python_version(), platform.system()))
    with open(os.path.join(RESULTS, "log.txt"), "w") as log:
        log.write("run_all.py log, %s, git %s, %s\n\n"
                  % (datetime.datetime.now().isoformat(timespec="seconds"),
                     git_head(), sys.version))
        for step in (step_pytest, step_exact_check, step_rigidity, step_stress_lp,
                     step_independent_check, step_deep_holes,
                     step_independent_deep_holes, step_lemma4_lemma8,
                     step_szollosi):
            step(log)
    total = round(time.time() - t_start, 1)

    failures = [c for c in CHECKS if not c["pass"]]
    print()
    print("%-40s %-46s %-24s %-24s %s" % ("step", "check", "expected", "computed", "ok"))
    for c in CHECKS:
        print("%-40s %-46s %-24s %-24s %s" % (c["step"][:40], c["check"][:46],
                                              str(c["expected"])[:24],
                                              str(c["computed"])[:24],
                                              "PASS" if c["pass"] else "FAIL"))
    summary = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "git": git_head(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "n_checks": len(CHECKS),
        "n_failed": len(failures),
        "all_pass": not failures,
        "total_wall_time_s": total,
        "steps": STEPS,
        "checks": CHECKS,
    }
    with open(os.path.join(RESULTS, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    print()
    print("%d checks, %d failed, total wall time %.1fs; results/summary.json, "
          "results/log.txt written" % (len(CHECKS), len(failures), total))
    if failures:
        print("FAILED CHECKS:")
        for c in failures:
            print("  %s :: %s : expected %r, computed %r"
                  % (c["step"], c["check"], c["expected"], c["computed"]))
        print("VERIFICATION FAILED")
        return 1
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
