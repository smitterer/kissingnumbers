#!/usr/bin/env python3
"""Exact-arithmetic checker for spherical-code configurations.

Reads a configuration in the JSON format documented in ``configs/README.md``
and verifies, in exact arithmetic only:

  (a) structural consistency: ``n_points == len(vectors)`` and
      ``dimension == len(v)`` for every vector ``v``;
  (b) for every i:      <x_i, x_i> == norm_squared   (exactly);
  (c) for every i < j:  <x_i, x_j> <= norm_squared/2 (exactly).

Together (b) and (c) say that ``x_i / r`` are unit vectors with pairwise
inner products at most 1/2, i.e. a kissing configuration.

No floating point value ever influences a verdict.  Floats appear only in
the human-readable report, always prefixed with '~' to mark them as
approximations.

------------------------------------------------------------------------
How decidability is guaranteed
------------------------------------------------------------------------
Every coordinate is parsed with sympy (decimal points and floats are
rejected outright) and then *canonicalised* into an element of the
multiquadratic ring

    R = span_Q { sqrt(d) : d a squarefree positive integer }

represented as a finite dict ``{d: Rational}``.  Because the numbers
``sqrt(d)`` for distinct squarefree ``d`` are linearly independent over Q,
this representation is *faithful*: an element is zero iff its dict is
empty.  Zero-testing is therefore exact and immediate, and it needs no
simplification heuristics at all.

Sign testing is exact too, by descent on the prime support.  Writing
``a = u + sqrt(p) * v`` with ``u, v`` free of the prime ``p``:

  * if u and v have the same (weak) sign, that is the sign of a;
  * otherwise sign(a) = +-sign(u^2 - p*v^2), and ``u^2 - p*v^2`` lives in a
    strictly smaller field, so the recursion terminates.

All arithmetic in this procedure is over Q.  There is no numerical
tolerance anywhere.  If an input expression cannot be brought into this
canonical form (e.g. an unsupported nested radical that sympy cannot
denest), the checker raises :class:`UndecidedError` and reports UNDECIDED
rather than guessing -- it never silently passes.

Exit codes: 0 = PASS, 1 = FAIL, 2 = structural/parse error or UNDECIDED.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import sympy
from sympy import Rational
from sympy.parsing.sympy_parser import parse_expr

__all__ = [
    "check_config",
    "check_config_data",
    "ConfigError",
    "UndecidedError",
    "format_report",
    "main",
]


# ----------------------------------------------------------------------
# errors
# ----------------------------------------------------------------------


class ConfigError(Exception):
    """Structural problem or non-exact / unparseable entry in a config."""


class UndecidedError(Exception):
    """An exact comparison could not be decided.  Never treated as a pass."""


class _Unsupported(Exception):
    """Internal: expression outside the canonical multiquadratic grammar."""


# ----------------------------------------------------------------------
# field elements: dict {squarefree positive int d: nonzero Rational c}
# meaning sum_d c * sqrt(d).   The empty dict is 0.
# ----------------------------------------------------------------------

FieldElem = Dict[int, Rational]

_ZERO: FieldElem = {}
_ONE: FieldElem = {1: Rational(1)}

_SQFREE_CACHE: Dict[int, Tuple[int, int]] = {}


def _sqfree_split(n: int) -> Tuple[int, int]:
    """Return (s, d) with n == s*s*d and d squarefree.  Requires n >= 1."""
    if n < 1:
        raise ValueError("_sqfree_split needs a positive integer")
    cached = _SQFREE_CACHE.get(n)
    if cached is not None:
        return cached
    s = 1
    d = 1
    for prime, exp in sympy.factorint(n).items():
        s *= prime ** (exp // 2)
        if exp % 2:
            d *= prime
    _SQFREE_CACHE[n] = (s, d)
    return s, d


def _fe_from_rational(r: Rational) -> FieldElem:
    return {} if r == 0 else {1: Rational(r)}


def _fe_sqrt_rational(q: Rational) -> FieldElem:
    """Canonical form of sqrt(q) for a non-negative rational q."""
    if q == 0:
        return {}
    if q < 0:
        raise _Unsupported("square root of a negative rational (not real)")
    p, r = q.p, q.q  # q = p/r in lowest terms, both positive
    s, d = _sqfree_split(p * r)  # sqrt(p/r) = sqrt(p*r)/r = (s/r)*sqrt(d)
    return {d: Rational(s, r)}


def _fe_add(a: FieldElem, b: FieldElem) -> FieldElem:
    if not a:
        return dict(b)
    if not b:
        return dict(a)
    out = dict(a)
    for d, c in b.items():
        nc = out.get(d)
        nc = c if nc is None else nc + c
        if nc == 0:
            out.pop(d, None)
        else:
            out[d] = nc
    return out


def _fe_neg(a: FieldElem) -> FieldElem:
    return {d: -c for d, c in a.items()}


def _fe_sub(a: FieldElem, b: FieldElem) -> FieldElem:
    return _fe_add(a, _fe_neg(b))


def _fe_scale(a: FieldElem, k: Rational) -> FieldElem:
    if k == 0:
        return {}
    return {d: c * k for d, c in a.items()}


def _fe_mul(a: FieldElem, b: FieldElem) -> FieldElem:
    if not a or not b:
        return {}
    out: FieldElem = {}
    for d1, c1 in a.items():
        for d2, c2 in b.items():
            if d1 == 1:
                d, coeff = d2, c1 * c2
            elif d2 == 1:
                d, coeff = d1, c1 * c2
            elif d1 == d2:
                d, coeff = 1, c1 * c2 * d1
            else:
                s, d = _sqfree_split(d1 * d2)
                coeff = c1 * c2 * s
            prev = out.get(d)
            coeff = coeff if prev is None else prev + coeff
            if coeff == 0:
                out.pop(d, None)
            else:
                out[d] = coeff
    return out


def _pick_prime(a: FieldElem) -> int:
    """Smallest prime dividing some radicand of ``a`` (a must be non-rational)."""
    primes = set()
    for d in a:
        if d != 1:
            primes.update(sympy.factorint(d).keys())
    if not primes:
        raise _Unsupported("no radical prime found")
    return min(int(p) for p in primes)


def _fe_split_prime(a: FieldElem, p: int) -> Tuple[FieldElem, FieldElem]:
    """Write a = u + sqrt(p)*v with u, v free of the prime p."""
    u: FieldElem = {}
    v: FieldElem = {}
    for d, c in a.items():
        if d % p == 0:
            v[d // p] = c
        else:
            u[d] = c
    return u, v


def _fe_sign(a: FieldElem) -> int:
    """Exact sign of a field element: -1, 0 or +1.  Pure rational arithmetic."""
    if not a:
        return 0
    if len(a) == 1:
        (_d, c), = a.items()
        return -1 if c < 0 else 1  # sqrt(d) > 0 for every d >= 1
    p = _pick_prime(a)
    u, v = _fe_split_prime(a, p)
    su = _fe_sign(u)
    sv = _fe_sign(v)
    if su >= 0 and sv >= 0:
        return 1 if (su or sv) else 0
    if su <= 0 and sv <= 0:
        return -1 if (su or sv) else 0
    # Opposite strict signs: compare u^2 against p*v^2 in the smaller field.
    w = _fe_sub(_fe_mul(u, u), _fe_scale(_fe_mul(v, v), Rational(p)))
    t = _fe_sign(w)
    if t == 0:
        # Impossible for a faithful representation (would make sqrt(p) lie in
        # the smaller field).  Refuse to guess.
        raise UndecidedError(
            "sign descent produced u^2 - %d*v^2 == 0; representation not faithful" % p
        )
    return t if su > 0 else -t


def _fe_is_zero(a: FieldElem) -> bool:
    return not a


def _fe_to_expr(a: FieldElem) -> sympy.Expr:
    """Rebuild a sympy expression (for display) from a canonical element."""
    if not a:
        return sympy.Integer(0)
    terms = []
    for d in sorted(a):
        c = a[d]
        terms.append(c if d == 1 else c * sympy.sqrt(d))
    return sympy.Add(*terms, evaluate=False) if len(terms) > 1 else terms[0]


def _fe_to_float(a: FieldElem) -> float:
    """Approximation for the report only.  Never used for a verdict."""
    total = 0.0
    for d, c in a.items():
        total += float(c) * (1.0 if d == 1 else math.sqrt(d))
    return total


# ----------------------------------------------------------------------
# parsing: string -> sympy expr -> canonical field element
# ----------------------------------------------------------------------

_DECIMAL_RE = re.compile(r"\d\s*[eE]\s*[+-]?\s*\d")

_PARSE_GLOBALS = {
    "sqrt": sympy.sqrt,
    "Integer": sympy.Integer,
    "Rational": sympy.Rational,
    "Float": sympy.Float,
    "Symbol": sympy.Symbol,
}


def _validate_exact(expr: sympy.Basic, where: str) -> None:
    """Reject anything outside {Rational, Add, Mul, Pow(_, Rational)}."""
    for node in sympy.preorder_traversal(expr):
        if isinstance(node, sympy.Float):
            raise ConfigError(
                "%s: floating point value %r is not exact; use an exact "
                "expression (integer, rational, or radical)" % (where, str(node))
            )
        if node.is_Rational:
            continue
        if node.is_Add or node.is_Mul:
            continue
        if node.is_Pow:
            if not node.exp.is_Rational:
                raise ConfigError(
                    "%s: exponent %r is not rational; only rational powers "
                    "(radicals) are allowed" % (where, str(node.exp))
                )
            continue
        raise ConfigError(
            "%s: unsupported term %r (only integers, rationals, radicals and "
            "+ - * / ( ) are allowed)" % (where, str(node))
        )


def parse_exact(text: object, where: str) -> sympy.Expr:
    """Parse one config entry into an exact sympy expression.

    Decimal points, exponent notation and floats are rejected.
    """
    if not isinstance(text, str):
        raise ConfigError(
            "%s: entry must be a JSON *string* holding an exact expression, "
            "got %r (%s)" % (where, text, type(text).__name__)
        )
    s = text.strip()
    if not s:
        raise ConfigError("%s: empty expression" % where)
    if "." in s:
        raise ConfigError(
            "%s: decimal point in %r is not allowed; configurations must be "
            "exact (use e.g. '1/2' instead of '0.5')" % (where, text)
        )
    if _DECIMAL_RE.search(s):
        raise ConfigError(
            "%s: exponent/float notation in %r is not allowed; "
            "configurations must be exact" % (where, text)
        )
    try:
        expr = parse_expr(s, global_dict=dict(_PARSE_GLOBALS), evaluate=True)
    except Exception as exc:  # noqa: BLE001 - any parse failure is a config error
        raise ConfigError("%s: cannot parse %r (%s)" % (where, text, exc)) from exc
    expr = sympy.sympify(expr)
    if expr.free_symbols:
        raise ConfigError(
            "%s: %r contains free symbols %s; only numeric expressions are "
            "allowed" % (where, text, sorted(str(x) for x in expr.free_symbols))
        )
    _validate_exact(expr, where)
    return expr


def _canon(expr: sympy.Expr) -> FieldElem:
    """Canonicalise an exact expression into the multiquadratic ring."""
    if expr.is_Rational:
        return _fe_from_rational(Rational(expr))
    if expr.is_Add:
        out: FieldElem = {}
        for arg in expr.args:
            out = _fe_add(out, _canon(arg))
        return out
    if expr.is_Mul:
        out = dict(_ONE)
        for arg in expr.args:
            out = _fe_mul(out, _canon(arg))
            if not out:
                return {}
        return out
    if expr.is_Pow:
        base, exp = expr.as_base_exp()
        if not exp.is_Rational:
            raise _Unsupported("non-rational exponent %s" % exp)
        if exp.is_Integer:
            n = int(exp)
            b = _canon(base)
            if n < 0:
                b = _fe_inv(b)
                n = -n
            out = dict(_ONE)
            for _ in range(n):
                out = _fe_mul(out, b)
            return out
        num, den = int(exp.p), int(exp.q)
        if den != 2:
            raise _Unsupported("exponent %s: only square roots are supported" % exp)
        if not base.is_Rational:
            raise _Unsupported("nested radical %s" % expr)
        val = Rational(base) ** num
        return _fe_sqrt_rational(val)
    raise _Unsupported("unsupported expression %s" % expr)


def _fe_inv(a: FieldElem) -> FieldElem:
    if not a:
        raise ConfigError("division by zero in a configuration entry")
    if len(a) == 1:
        (d, c), = a.items()
        # 1 / (c * sqrt(d)) = sqrt(d) / (c * d)
        return {d: Rational(1) / (c * d)}
    p = _pick_prime(a)
    conj = {d: (-c if d % p == 0 else c) for d, c in a.items()}
    denom = _fe_mul(a, conj)  # lies in the field without sqrt(p)
    return _fe_mul(conj, _fe_inv(denom))


def to_field(expr: sympy.Expr, where: str) -> FieldElem:
    """Canonicalise, with one sympy-assisted retry for awkward radicals."""
    try:
        return _canon(expr)
    except _Unsupported:
        pass
    try:
        retry = sympy.expand(sympy.radsimp(sympy.sqrtdenest(sympy.expand(expr))))
        return _canon(retry)
    except _Unsupported as exc:
        raise UndecidedError(
            "UNDECIDED: %s: cannot bring %r into an exact canonical form (%s). "
            "Refusing to guess -- no verdict is issued for this configuration."
            % (where, str(expr), exc)
        ) from exc


# ----------------------------------------------------------------------
# the checks
# ----------------------------------------------------------------------


def _load_json(path: str) -> dict:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError as exc:
        raise ConfigError("no such configuration file: %s" % path) from exc
    except json.JSONDecodeError as exc:
        raise ConfigError("%s: invalid JSON (%s)" % (path, exc)) from exc
    if not isinstance(data, dict):
        raise ConfigError("%s: top level must be a JSON object" % path)
    return data


def check_config_data(data: dict, source: str = "<dict>") -> dict:
    """Run all exact checks on an already-loaded configuration dict.

    Returns the result dict described in :func:`check_config`.
    Raises :class:`ConfigError` on structural/parse problems and
    :class:`UndecidedError` if an exact comparison cannot be decided.
    """
    name = data.get("name", "<unnamed>")
    if not isinstance(name, str):
        raise ConfigError("%s: 'name' must be a string" % source)

    for key in ("dimension", "n_points", "norm_squared", "vectors"):
        if key not in data:
            raise ConfigError("%s: missing required field %r" % (source, key))

    dimension = data["dimension"]
    n_points = data["n_points"]
    vectors = data["vectors"]
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension < 1:
        raise ConfigError("%s: 'dimension' must be a positive integer" % source)
    if not isinstance(n_points, int) or isinstance(n_points, bool) or n_points < 0:
        raise ConfigError("%s: 'n_points' must be a non-negative integer" % source)
    if not isinstance(vectors, list):
        raise ConfigError("%s: 'vectors' must be a list" % source)

    # (a) structural
    if n_points != len(vectors):
        raise ConfigError(
            "%s: n_points = %d but 'vectors' has %d entries"
            % (source, n_points, len(vectors))
        )
    for i, vec in enumerate(vectors):
        if not isinstance(vec, list):
            raise ConfigError("%s: vector %d is not a list" % (source, i))
        if len(vec) != dimension:
            raise ConfigError(
                "%s: vector %d has %d coordinates, expected dimension = %d"
                % (source, i, len(vec), dimension)
            )

    ns_expr = parse_exact(data["norm_squared"], "%s: norm_squared" % source)
    ns = to_field(ns_expr, "%s: norm_squared" % source)
    if _fe_sign(ns) <= 0:
        raise ConfigError(
            "%s: norm_squared must be strictly positive (got %s)"
            % (source, _fe_to_expr(ns))
        )
    half = _fe_scale(ns, Rational(1, 2))

    pts: List[List[FieldElem]] = []
    for i, vec in enumerate(vectors):
        row = []
        for k, coord in enumerate(vec):
            where = "%s: vector %d, coordinate %d" % (source, i, k)
            row.append(to_field(parse_exact(coord, where), where))
        pts.append(row)

    # (b) norms
    norm_failures: List[dict] = []
    for i, row in enumerate(pts):
        sq: FieldElem = {}
        for c in row:
            sq = _fe_add(sq, _fe_mul(c, c))
        diff = _fe_sub(sq, ns)
        if not _fe_is_zero(diff):
            norm_failures.append(
                {
                    "index": i,
                    "norm_squared": str(_fe_to_expr(sq)),
                    "norm_squared_float": _fe_to_float(sq),
                    "expected": str(_fe_to_expr(ns)),
                    "expected_float": _fe_to_float(ns),
                    "difference": str(_fe_to_expr(diff)),
                    "difference_float": _fe_to_float(diff),
                }
            )

    # (c) pairwise inner products
    n = len(pts)
    worst_pair: Optional[Tuple[int, int]] = None
    worst_ip: FieldElem = {}
    worst_margin: Optional[FieldElem] = None
    violations: List[dict] = []

    for i in range(n):
        ri = pts[i]
        for j in range(i + 1, n):
            rj = pts[j]
            ip: FieldElem = {}
            for a, b in zip(ri, rj):
                ip = _fe_add(ip, _fe_mul(a, b))
            margin = _fe_sub(half, ip)  # want margin >= 0
            sgn = _fe_sign(margin)  # exact; raises UndecidedError if undecidable
            if worst_margin is None or _fe_sign(_fe_sub(margin, worst_margin)) < 0:
                worst_margin = margin
                worst_ip = ip
                worst_pair = (i, j)
            if sgn < 0:
                violations.append(
                    {
                        "pair": [i, j],
                        "inner_product": str(_fe_to_expr(ip)),
                        "inner_product_float": _fe_to_float(ip),
                        "margin": str(_fe_to_expr(margin)),
                        "margin_float": _fe_to_float(margin),
                    }
                )

    passed = not norm_failures and not violations

    result = {
        "name": name,
        "source": source,
        "pass": passed,
        "dimension": dimension,
        "n_points": n_points,
        "n_pairs": n * (n - 1) // 2,
        "norm_squared": str(_fe_to_expr(ns)),
        "norm_squared_float": _fe_to_float(ns),
        "threshold": str(_fe_to_expr(half)),
        "threshold_float": _fe_to_float(half),
        "worst_pair": list(worst_pair) if worst_pair is not None else None,
        "worst_ip": str(_fe_to_expr(worst_ip)) if worst_pair is not None else None,
        "worst_ip_float": _fe_to_float(worst_ip) if worst_pair is not None else None,
        "margin": str(_fe_to_expr(worst_margin)) if worst_margin is not None else None,
        "margin_float": _fe_to_float(worst_margin) if worst_margin is not None else None,
        "norm_failures": norm_failures,
        "violations": violations,
    }
    return result


def check_config(path: str) -> dict:
    """Check the configuration stored at ``path``.

    Returns a dict with (at least) the fields:
      ``name``, ``pass`` (bool), ``worst_pair`` ([i, j] or None),
      ``worst_ip`` (str of the exact sympy expression), ``margin`` (str of
      the exact r^2/2 - ip for the worst pair), plus ``norm_failures``,
      ``violations`` and float approximations suffixed ``_float``.

    Raises :class:`ConfigError` for structural/parse errors and
    :class:`UndecidedError` if an exact comparison cannot be decided.
    """
    data = _load_json(path)
    return check_config_data(data, source=path)


# ----------------------------------------------------------------------
# reporting / CLI
# ----------------------------------------------------------------------


def _f(x: float) -> str:
    return "~%.12g" % x


def format_report(result: dict) -> str:
    lines = []
    lines.append("Configuration : %s" % result["name"])
    lines.append("Source        : %s" % result["source"])
    lines.append("Dimension     : %d" % result["dimension"])
    lines.append("Points        : %d   (%d pairs)" % (result["n_points"], result["n_pairs"]))
    lines.append(
        "norm_squared  : %s   (%s)"
        % (result["norm_squared"], _f(result["norm_squared_float"]))
    )
    lines.append(
        "threshold r^2/2: %s   (%s)"
        % (result["threshold"], _f(result["threshold_float"]))
    )
    lines.append("")
    if result["worst_pair"] is None:
        lines.append("Worst pair    : (none -- fewer than two points)")
    else:
        i, j = result["worst_pair"]
        lines.append("Worst pair    : (%d, %d)" % (i, j))
        lines.append(
            "  <x_%d, x_%d>      = %s   (%s)"
            % (i, j, result["worst_ip"], _f(result["worst_ip_float"]))
        )
        lines.append(
            "  margin r^2/2-ip = %s   (%s)"
            % (result["margin"], _f(result["margin_float"]))
        )
    lines.append("")

    nf = result["norm_failures"]
    if nf:
        lines.append("NORM VIOLATIONS (check b): %d" % len(nf))
        for f in nf[:20]:
            lines.append(
                "  vector %d: <x,x> = %s (%s), expected %s (%s), difference %s (%s)"
                % (
                    f["index"],
                    f["norm_squared"],
                    _f(f["norm_squared_float"]),
                    f["expected"],
                    _f(f["expected_float"]),
                    f["difference"],
                    _f(f["difference_float"]),
                )
            )
        if len(nf) > 20:
            lines.append("  ... and %d more" % (len(nf) - 20))
    else:
        lines.append("Norm checks (b): all %d vectors have <x,x> = r^2 exactly."
                     % result["n_points"])

    vi = result["violations"]
    if vi:
        lines.append("INNER-PRODUCT VIOLATIONS (check c): %d" % len(vi))
        for v in vi[:20]:
            i, j = v["pair"]
            lines.append(
                "  pair (%d, %d): ip = %s (%s), margin = %s (%s)"
                % (
                    i,
                    j,
                    v["inner_product"],
                    _f(v["inner_product_float"]),
                    v["margin"],
                    _f(v["margin_float"]),
                )
            )
        if len(vi) > 20:
            lines.append("  ... and %d more" % (len(vi) - 20))
    else:
        lines.append(
            "Pair checks (c): all %d pairs satisfy <x_i,x_j> <= r^2/2 exactly."
            % result["n_pairs"]
        )

    lines.append("")
    lines.append("VERDICT: %s" % ("PASS" if result["pass"] else "FAIL"))
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="exact_check.py",
        description="Exact-arithmetic checker for kissing-number configurations.",
    )
    parser.add_argument("configs", nargs="+", metavar="config.json")
    parser.add_argument(
        "--json",
        dest="as_json",
        action="store_true",
        help="emit the result dict(s) as JSON instead of a human-readable report",
    )
    args = parser.parse_args(argv)

    results = []
    exit_code = 0
    for k, path in enumerate(args.configs):
        try:
            result = check_config(path)
        except ConfigError as exc:
            if args.as_json:
                results.append({"source": path, "pass": False, "error": str(exc),
                                "error_type": "ConfigError"})
            else:
                if k:
                    print()
                print("Configuration : %s" % path)
                print("ERROR (structural/parse): %s" % exc)
                print("VERDICT: ERROR")
            exit_code = 2
            continue
        except UndecidedError as exc:
            if args.as_json:
                results.append({"source": path, "pass": False, "error": str(exc),
                                "error_type": "UndecidedError"})
            else:
                if k:
                    print()
                print("Configuration : %s" % path)
                print("UNDECIDED: %s" % exc)
                print("VERDICT: UNDECIDED (no pass is claimed)")
            exit_code = 2
            continue

        results.append(result)
        if not args.as_json:
            if k:
                print()
            print(format_report(result))
        if not result["pass"] and exit_code == 0:
            exit_code = 1

    if args.as_json:
        payload = results[0] if len(results) == 1 else results
        print(json.dumps(payload, indent=2))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
