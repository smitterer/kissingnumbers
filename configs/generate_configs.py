"""Regenerate (or, by default, re-check) the configuration JSONs from the
committed source excerpts.

Provenance: coordinates for D5/L5/Q5/R5 come from Cohn & Rajagopal,
arXiv:2412.00937, Table 2.2 (all at squared norm 2); the independent Q5
transcription comes from Szollosi, arXiv:2301.08272, appendix (integers at
squared norm 50). The excerpts under configs/sources/ were extracted with
`pdftotext -layout` from the arXiv PDFs on 2026-08-18.

As a generation-side sanity check (NOT the authoritative verification --
that is verify/exact_check.py), each configuration's unordered-pair
inner-product profile is compared against Table 2.1 of arXiv:2412.00937.

Usage (from the repo root):

    python3 configs/generate_configs.py            # compare mode (default):
        regenerate in memory, compare with the committed JSON files as
        point sets, print PASS/FAIL per file, exit 1 on any FAIL;
        nothing is written.
    python3 configs/generate_configs.py --write    # (re)write the five files.
"""

import itertools
import json
import re
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).parent

# Table 2.1 of arXiv:2412.00937: unordered pairs at each unit inner product.
EXPECTED_PROFILES = {
    "d5": {Fraction(-1): 20, Fraction(-1, 2): 240, Fraction(0): 280,
           Fraction(1, 2): 240},
    "l5": {Fraction(-1): 12, Fraction(-3, 4): 32, Fraction(-1, 2): 192,
           Fraction(-1, 4): 32, Fraction(0): 272, Fraction(1, 2): 240},
    "q5": {Fraction(-1): 10, Fraction(-4, 5): 30, Fraction(-1, 2): 180,
           Fraction(-3, 10): 60, Fraction(0): 250, Fraction(1, 5): 10,
           Fraction(1, 2): 240},
    "r5": {Fraction(-1): 6, Fraction(-4, 5): 30, Fraction(-3, 4): 20,
           Fraction(-1, 2): 144, Fraction(-3, 10): 60, Fraction(-1, 4): 28,
           Fraction(0): 242, Fraction(1, 5): 10, Fraction(1, 2): 240},
}

SOURCES = {
    "d5": "Minimal vectors of the D5 root lattice (Korkine & Zolotareff 1873); coordinates as in Cohn & Rajagopal, arXiv:2412.00937, Table 2.2 (squared norm 2).",
    "l5": "Leech (1967) non-lattice 40-point kissing configuration; coordinates from Cohn & Rajagopal, arXiv:2412.00937, Table 2.2 (squared norm 2).",
    "q5": "Szollosi (2023) Q5 configuration; coordinates from Cohn & Rajagopal, arXiv:2412.00937, Table 2.2 (squared norm 2).",
    "r5": "Cohn & Rajagopal (2024) R5 configuration; coordinates from arXiv:2412.00937, Table 2.2 (squared norm 2).",
}

SZOLLOSI_SOURCE = (
    "Szollosi arXiv:2301.08272 appendix, machine-readable integer form "
    "(rescale 1/(5*sqrt(2)) => squared norm 50). Independent transcription "
    "of Q5 for cross-validation.")


def parse_table_2_2():
    text = (HERE / "sources" / "cohn_rajagopal_2412.00937_table2.2.txt").read_text()
    text = text.replace("−", "-")
    pts = []
    for t in re.findall(r"\(([-0-9.,]+)\)", text):
        parts = t.split(",")
        if len(parts) == 5:
            pts.append(tuple(Fraction(p) for p in parts))
    assert len(pts) == 160, f"expected 160 points in Table 2.2, got {len(pts)}"
    return {"d5": pts[0:40], "l5": pts[40:80], "q5": pts[80:120], "r5": pts[120:160]}


def parse_szollosi_appendix():
    text = (HERE / "sources" / "szollosi_2301.08272_appendix.txt").read_text()
    vecs = [tuple(Fraction(x) for x in m.split(","))
            for m in re.findall(r"\{([-0-9,]+)\}", text)]
    vecs = [v for v in vecs if len(v) == 5]
    assert len(vecs) == 40, f"expected 40 Szollosi vectors, got {len(vecs)}"
    return vecs


def profile(vectors, norm_squared):
    c = Counter()
    for u, v in itertools.combinations(vectors, 2):
        c[sum(a * b for a, b in zip(u, v)) / norm_squared] += 1
    return dict(c)


def frac_str(f):
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def config_dict(name, vectors, norm_squared, source):
    assert len(set(vectors)) == len(vectors) == 40
    for v in vectors:
        assert sum(a * a for a in v) == norm_squared, f"{name}: bad norm for {v}"
    return {
        "name": name,
        "dimension": 5,
        "n_points": 40,
        "source": source,
        "norm_squared": frac_str(Fraction(norm_squared)),
        "vectors": [[frac_str(c) for c in v] for v in vectors],
    }


def serialize(cfg):
    return json.dumps(cfg, indent=1) + "\n"


def build_all():
    """{file name: config dict} for the five files, profiles checked."""
    out = {}
    blocks = parse_table_2_2()
    for name, vecs in blocks.items():
        got = profile(vecs, 2)
        assert got == EXPECTED_PROFILES[name], f"{name}: profile mismatch vs Table 2.1"
        out[f"{name}_40.json"] = config_dict(name, vecs, 2, SOURCES[name])
    sz = parse_szollosi_appendix()
    assert profile(sz, 50) == EXPECTED_PROFILES["q5"], "szollosi q5: profile mismatch"
    out["q5_szollosi_40.json"] = config_dict("q5_szollosi", sz, 50, SZOLLOSI_SOURCE)
    return out


def compare(path, cfg):
    """Committed file vs. regenerated configuration. True iff equal as point
    sets with equal metadata (name, dimension, n_points, norm_squared)."""
    if not path.exists():
        print(f"{path.name}: FAIL (file missing)")
        return False
    committed = json.loads(path.read_text())
    same_points = ({tuple(Fraction(c) for c in v) for v in committed["vectors"]}
                   == {tuple(Fraction(c) for c in v) for v in cfg["vectors"]})
    same_meta = all(committed.get(k) == cfg[k]
                    for k in ("name", "dimension", "n_points", "norm_squared"))
    identical = path.read_text() == serialize(cfg)
    ok = same_points and same_meta and len(committed["vectors"]) == 40
    print(f"{path.name}: {'PASS' if ok else 'FAIL'}  (equal as point sets: "
          f"{same_points}; metadata equal: {same_meta}; file byte-identical to "
          f"regeneration: {identical})")
    return ok


def main(argv):
    write = "--write" in argv
    cfgs = build_all()
    allok = True
    for fname, cfg in cfgs.items():
        path = HERE / fname
        if write:
            path.write_text(serialize(cfg))
            print(f"{fname}: written (40 points, profile matches Table 2.1)")
        else:
            allok = compare(path, cfg) and allok
    if write:
        return 0
    print("OVERALL:", "PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
