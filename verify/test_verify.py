"""Tests for the exact-arithmetic configuration checker.

Run with:  python -m pytest verify/ -q
"""

from __future__ import annotations

import copy
import json
import os
import subprocess
import sys

import pytest
import sympy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import exact_check  # noqa: E402
from exact_check import (  # noqa: E402
    ConfigError,
    check_config,
    check_config_data,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D5_PATH = os.path.join(REPO_ROOT, "configs", "d5_40.json")
CHECKER = os.path.join(REPO_ROOT, "verify", "exact_check.py")


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------


def load_d5() -> dict:
    with open(D5_PATH, "r", encoding="utf-8") as fh:
        return json.load(fh)


def write_config(tmp_path, data: dict, name: str = "cfg.json") -> str:
    path = os.path.join(str(tmp_path), name)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh)
    return path


def run_cli(*args):
    return subprocess.run(
        [sys.executable, CHECKER, *args],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )


def exact_eq(text: str, expected) -> bool:
    """Exact comparison of a reported expression against an expected value."""
    return sympy.simplify(sympy.sympify(text) - sympy.sympify(expected)) == 0


# ----------------------------------------------------------------------
# 1. the reference D5 configuration passes
# ----------------------------------------------------------------------


def test_d5_passes():
    res = check_config(D5_PATH)
    assert res["pass"] is True
    assert res["name"] == "d5"
    assert res["dimension"] == 5
    assert res["n_points"] == 40
    assert res["n_pairs"] == 780
    assert res["norm_failures"] == []
    assert res["violations"] == []
    # worst pair achieves the bound exactly: ip = r^2/2 = 1, margin = 0
    assert exact_eq(res["worst_ip"], 1)
    assert exact_eq(res["margin"], 0)
    assert exact_eq(res["norm_squared"], 2)
    assert exact_eq(res["threshold"], 1)
    i, j = res["worst_pair"]
    assert 0 <= i < j < 40


def test_d5_cli_exit_code_zero_and_reports_worst_pair():
    proc = run_cli("configs/d5_40.json")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "VERDICT: PASS" in proc.stdout
    assert "Worst pair" in proc.stdout
    assert "margin r^2/2-ip = 0" in proc.stdout


def test_json_flag_emits_result_dict():
    proc = run_cli("--json", "configs/d5_40.json")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    payload = json.loads(proc.stdout)
    assert payload["pass"] is True
    assert payload["name"] == "d5"
    assert payload["worst_pair"] is not None
    assert exact_eq(payload["worst_ip"], 1)


# ----------------------------------------------------------------------
# 2. perturbed D5 fails, with the offending pair identified
# ----------------------------------------------------------------------


def test_perturbed_d5_fails(tmp_path):
    original = load_d5()
    vecs = [[sympy.Rational(x) for x in v] for v in original["vectors"]]

    # perturb one coordinate that equals 1:  "1" -> "101/100"
    col = next(k for k, x in enumerate(vecs[0]) if x == 1)
    data = load_d5()
    data["vectors"][0][col] = "101/100"
    data["name"] = "d5-perturbed"
    path = write_config(tmp_path, data)

    # exactly which pairs must now break: those at the bound whose partner
    # has a +1 in the perturbed coordinate (ip becomes 1 + 1/100 > 1)
    expected_pairs = {
        (0, j)
        for j in range(1, len(vecs))
        if sum(a * b for a, b in zip(vecs[0], vecs[j])) == 1 and vecs[j][col] == 1
    }
    assert expected_pairs, "perturbation must break at least one pair"

    res = check_config(path)
    assert res["pass"] is False

    # the perturbed vector no longer has squared norm 2
    assert {f["index"] for f in res["norm_failures"]} == {0}

    offending = {tuple(v["pair"]) for v in res["violations"]}
    assert offending == expected_pairs
    v0 = res["violations"][0]
    assert exact_eq(v0["inner_product"], sympy.Rational(101, 100))
    assert exact_eq(v0["margin"], sympy.Rational(-1, 100))

    proc = run_cli(path)
    assert proc.returncode == 1
    assert "VERDICT: FAIL" in proc.stdout


# ----------------------------------------------------------------------
# 3. D5 plus a duplicated 41st point fails
# ----------------------------------------------------------------------


def test_duplicate_41st_point_fails(tmp_path):
    data = load_d5()
    data["vectors"].append(list(data["vectors"][0]))  # duplicate of vector 0
    data["n_points"] = len(data["vectors"])
    data["name"] = "d5-plus-duplicate"
    path = write_config(tmp_path, data)

    res = check_config(path)
    assert res["pass"] is False
    assert res["n_points"] == 41
    assert res["n_pairs"] == 820
    assert res["norm_failures"] == []  # the duplicate has the right norm
    assert res["violations"], "duplicate point must violate the pair condition"
    assert [0, 40] in [v["pair"] for v in res["violations"]]
    v = next(v for v in res["violations"] if v["pair"] == [0, 40])
    assert exact_eq(v["inner_product"], 2)  # ip = r^2 = 2 > r^2/2 = 1
    assert exact_eq(v["margin"], -1)
    # the worst pair reported is the duplicate pair
    assert exact_eq(res["worst_ip"], 2)
    assert exact_eq(res["margin"], -1)

    proc = run_cli(path)
    assert proc.returncode == 1
    assert "VERDICT: FAIL" in proc.stdout


# ----------------------------------------------------------------------
# 4. norm violation (check b) fails even when all pairs are fine
# ----------------------------------------------------------------------


def test_norm_violation_fails(tmp_path):
    data = {
        "name": "bad-norm",
        "dimension": 5,
        "n_points": 2,
        "norm_squared": "2",
        "vectors": [
            ["1", "1", "0", "0", "0"],
            ["0", "0", "1", "0", "0"],  # squared norm 1, not 2
        ],
    }
    res = check_config_data(data, source="bad-norm")
    assert res["pass"] is False
    assert res["violations"] == []  # ip = 0 <= 1, pairs are fine
    assert [f["index"] for f in res["norm_failures"]] == [1]
    f = res["norm_failures"][0]
    assert exact_eq(f["norm_squared"], 1)
    assert exact_eq(f["expected"], 2)
    assert exact_eq(f["difference"], -1)

    path = write_config(tmp_path, data)
    proc = run_cli(path)
    assert proc.returncode == 1
    assert "NORM VIOLATIONS" in proc.stdout
    assert "VERDICT: FAIL" in proc.stdout


# ----------------------------------------------------------------------
# 5. decimals and floats are rejected outright (exit code 2)
# ----------------------------------------------------------------------


def test_decimal_coordinate_rejected(tmp_path):
    data = {
        "name": "decimal",
        "dimension": 5,
        "n_points": 1,
        "norm_squared": "2",
        "vectors": [["0.5", "1", "0", "0", "0"]],
    }
    with pytest.raises(ConfigError) as exc:
        check_config_data(data, source="decimal")
    assert "decimal" in str(exc.value).lower()

    path = write_config(tmp_path, data)
    proc = run_cli(path)
    assert proc.returncode == 2
    assert "ERROR" in proc.stdout
    assert "VERDICT: PASS" not in proc.stdout


def test_exponent_notation_rejected():
    data = {
        "name": "expnot",
        "dimension": 2,
        "n_points": 1,
        "norm_squared": "2",
        "vectors": [["1e-3", "1"]],
    }
    with pytest.raises(ConfigError):
        check_config_data(data, source="expnot")


def test_json_number_rejected():
    """Coordinates must be strings; a raw JSON number is not accepted."""
    data = {
        "name": "rawnumber",
        "dimension": 2,
        "n_points": 1,
        "norm_squared": "2",
        "vectors": [[1, 1]],
    }
    with pytest.raises(ConfigError):
        check_config_data(data, source="rawnumber")


def test_symbolic_junk_rejected():
    data = {
        "name": "junk",
        "dimension": 2,
        "n_points": 1,
        "norm_squared": "2",
        "vectors": [["x", "1"]],
    }
    with pytest.raises(ConfigError):
        check_config_data(data, source="junk")


# ----------------------------------------------------------------------
# 6. structural errors (exit code 2)
# ----------------------------------------------------------------------


def test_n_points_mismatch_is_structural_error(tmp_path):
    data = load_d5()
    data["n_points"] = 41
    with pytest.raises(ConfigError):
        check_config_data(data, source="mismatch")
    path = write_config(tmp_path, data)
    assert run_cli(path).returncode == 2


def test_dimension_mismatch_is_structural_error(tmp_path):
    data = load_d5()
    data["vectors"][3] = ["1", "1", "0", "0"]  # only 4 coordinates
    with pytest.raises(ConfigError):
        check_config_data(data, source="dim-mismatch")
    path = write_config(tmp_path, data)
    assert run_cli(path).returncode == 2


def test_missing_file_is_error():
    proc = run_cli(os.path.join(REPO_ROOT, "configs", "does_not_exist.json"))
    assert proc.returncode == 2


# ----------------------------------------------------------------------
# 7. irrational coordinates: Q(sqrt 2) configurations
# ----------------------------------------------------------------------


def test_irrational_config_passes():
    """Hand-verified in Q(sqrt 2), norm_squared = 2, threshold = 1.

    x0 = (sqrt2, 0, 0, 0, 0)                 <x0,x0> = 2
    x1 = (sqrt2/2, sqrt2/2, 1, 0, 0)         <x1,x1> = 1/2 + 1/2 + 1 = 2
    x2 = (-sqrt2/2, sqrt2/2, 0, 1, 0)        <x2,x2> = 1/2 + 1/2 + 1 = 2
    <x0,x1> = sqrt2*sqrt2/2 = 1   = r^2/2   (tight, margin 0)
    <x0,x2> = -sqrt2*sqrt2/2 = -1           (margin 2)
    <x1,x2> = -1/2 + 1/2 = 0                (margin 1)
    """
    data = {
        "name": "q-sqrt2-pass",
        "dimension": 5,
        "n_points": 3,
        "norm_squared": "2",
        "vectors": [
            ["sqrt(2)", "0", "0", "0", "0"],
            ["sqrt(2)/2", "sqrt(2)/2", "1", "0", "0"],
            ["-sqrt(2)/2", "sqrt(2)/2", "0", "1", "0"],
        ],
    }
    res = check_config_data(data, source="q-sqrt2-pass")
    assert res["pass"] is True
    assert res["norm_failures"] == []
    assert res["violations"] == []
    assert res["worst_pair"] == [0, 1]
    assert exact_eq(res["worst_ip"], 1)
    assert exact_eq(res["margin"], 0)


def test_irrational_config_with_strictly_positive_irrational_margin_passes():
    """Exercises the exact sign test on an irrational margin.

    norm_squared = 2, threshold = 1.
    x0 = (sqrt2, 0, 0, 0, 0)          <x0,x0> = 2
    x1 = (1/2, sqrt(7)/2, 0, 0, 0)    <x1,x1> = 1/4 + 7/4 = 2
    <x0,x1> = sqrt2/2 (~0.7071) <= 1, margin = 1 - sqrt2/2 (~0.2929) > 0.
    """
    data = {
        "name": "q-sqrt2-strict",
        "dimension": 5,
        "n_points": 2,
        "norm_squared": "2",
        "vectors": [
            ["sqrt(2)", "0", "0", "0", "0"],
            ["1/2", "sqrt(7)/2", "0", "0", "0"],
        ],
    }
    res = check_config_data(data, source="q-sqrt2-strict")
    assert res["pass"] is True
    assert exact_eq(res["worst_ip"], sympy.sqrt(2) / 2)
    assert exact_eq(res["margin"], 1 - sympy.sqrt(2) / 2)
    assert res["margin_float"] == pytest.approx(1 - 2 ** 0.5 / 2)


def test_irrational_config_with_negative_irrational_margin_fails():
    """margin = 1 - sqrt2 (~ -0.4142) < 0 must be detected exactly."""
    data = {
        "name": "q-sqrt2-fail",
        "dimension": 5,
        "n_points": 2,
        "norm_squared": "2",
        "vectors": [
            ["sqrt(2)", "0", "0", "0", "0"],
            ["1", "1", "0", "0", "0"],
        ],
    }
    res = check_config_data(data, source="q-sqrt2-fail")
    assert res["pass"] is False
    assert res["norm_failures"] == []
    assert [v["pair"] for v in res["violations"]] == [[0, 1]]
    assert exact_eq(res["worst_ip"], sympy.sqrt(2))
    assert exact_eq(res["margin"], 1 - sympy.sqrt(2))


def test_rationalised_denominator_coordinates():
    """'1/sqrt(10)'-style entries are handled exactly."""
    data = {
        "name": "inv-sqrt",
        "dimension": 2,
        "n_points": 2,
        "norm_squared": "1",
        "vectors": [
            ["3/sqrt(10)", "1/sqrt(10)"],
            ["-1/sqrt(10)", "3/sqrt(10)"],
        ],
    }
    res = check_config_data(data, source="inv-sqrt")
    assert res["pass"] is True
    assert exact_eq(res["worst_ip"], 0)  # orthogonal
    assert exact_eq(res["margin"], sympy.Rational(1, 2))


# ----------------------------------------------------------------------
# 8. the exact sign machinery itself
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "expr, expected_sign",
    [
        ("0", 0),
        ("1/3", 1),
        ("-7/2", -1),
        ("sqrt(2) - 1", 1),
        ("1 - sqrt(2)", -1),
        ("sqrt(2) + sqrt(3) - sqrt(6)", 1),  # ~0.697
        ("sqrt(2) + sqrt(3) - sqrt(5)", 1),  # ~0.910
        ("sqrt(2)*sqrt(2) - 2", 0),  # exact zero through radicals
        ("sqrt(6)*sqrt(10) - 2*sqrt(15)", 0),  # non-squarefree product
        ("(1+sqrt(5))/2 - 1618/1000", 1),  # golden ratio vs a close rational
        ("(1+sqrt(5))/2 - 1619/1000", -1),
        ("1/sqrt(2) - sqrt(2)/2", 0),
    ],
)
def test_exact_sign(expr, expected_sign):
    fe = exact_check.to_field(exact_check.parse_exact(expr, "test"), "test")
    assert exact_check._fe_sign(fe) == expected_sign
    # cross-check against sympy's own numeric evaluation (sanity only, not
    # part of the verdict path)
    approx = float(sympy.sympify(expr).evalf(50))
    assert expected_sign == (0 if abs(approx) < 1e-30 else (1 if approx > 0 else -1))


def test_zero_test_is_exact_not_numeric():
    """A tiny but non-zero algebraic quantity must not be rounded to zero."""
    expr = "(1+sqrt(5))/2 - 1618033988749894848204586834365638117720/1000000000000000000000000000000000000000"
    fe = exact_check.to_field(exact_check.parse_exact(expr, "test"), "test")
    assert exact_check._fe_sign(fe) == 1  # phi is slightly larger
    assert fe  # non-empty dict => provably non-zero


def test_nested_radical_is_denested_or_reported_undecided():
    """sqrt(3+2*sqrt(2)) = 1 + sqrt(2); the checker must not guess."""
    expr = "sqrt(3+2*sqrt(2)) - 1 - sqrt(2)"
    try:
        fe = exact_check.to_field(exact_check.parse_exact(expr, "test"), "test")
    except exact_check.UndecidedError as exc:
        assert "UNDECIDED" in str(exc)
    else:
        assert exact_check._fe_sign(fe) == 0


# ----------------------------------------------------------------------
# 9. multiple configs on one command line
# ----------------------------------------------------------------------


def test_cli_multiple_configs(tmp_path):
    bad = load_d5()
    bad["vectors"].append(list(bad["vectors"][0]))
    bad["n_points"] = len(bad["vectors"])
    bad_path = write_config(tmp_path, bad, "bad.json")
    proc = run_cli("configs/d5_40.json", bad_path)
    assert proc.returncode == 1
    assert proc.stdout.count("VERDICT:") == 2
    assert "VERDICT: PASS" in proc.stdout
    assert "VERDICT: FAIL" in proc.stdout


def test_check_config_does_not_mutate_input():
    data = load_d5()
    before = copy.deepcopy(data)
    check_config_data(data, source="immutability")
    assert data == before
