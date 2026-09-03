# Configuration file format

A configuration is a JSON file with exact coordinates. The format avoids
floating point entirely.

```json
{
  "name": "d5",
  "dimension": 5,
  "n_points": 40,
  "source": "citation for where this configuration comes from",
  "norm_squared": "2",
  "vectors": [
    ["1", "1", "0", "0", "0"],
    ["1", "-1", "0", "0", "0"]
  ]
}
```

- Every coordinate and `norm_squared` is a **string holding an exact
  symbolic expression**: integers, rationals (`"3/4"`), and square roots
  (`"sqrt(2)"`, `"1/sqrt(10)"`, `"(1+sqrt(5))/4"`), composed with
  `+ - * / ( )`. Parseable by `sympy.sympify`. No decimal points allowed.
- Vectors are given at common squared norm `norm_squared` = r² (not
  necessarily 1, to keep coordinates clean). The kissing conditions are:
  - for every i: ⟨x_i, x_i⟩ = r² exactly,
  - for every i < j: ⟨x_i, x_j⟩ ≤ r²/2 exactly
  (equivalent to unit vectors with pairwise inner product ≤ 1/2).
- `n_points` must equal the length of `vectors`; `dimension` the length of
  each vector.

The independent checker for this format lives in `verify/exact_check.py`.
