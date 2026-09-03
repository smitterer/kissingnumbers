# Rigidity of the known five-dimensional kissing configurations — verification repository

Code, exact data and certificates for the paper *Rigidity of the known
five-dimensional kissing configurations and a proof of a conjecture of Cohn
and Rajagopal* (`docs/preprint/kissing5_note.tex`).

Repository: <https://github.com/smitterer/kissingnumbers>, release v1.0
(how to cite: section 9).

## 1. What this repository certifies, and what it does not

The repository certifies, in exact rational arithmetic, the computational
claims of the paper: **Theorem 1** (each of the four known 40-point kissing
configurations D5, L5, Q5, R5 in R⁵ is infinitesimally jammed: a strictly
positive self-stress exists on all 240 contacts and the equality system of
40 tangency and 240 contact equations has rank exactly 190 in 200
unknowns), **Theorem 2** (the polytope P_X = {w : ⟨x,w⟩ ≤ 1 for all x ∈ X}
has 42/50/92/100 vertices with maximal squared norm 5/4, attained at exactly
32 vertices, so no 41st point can be adjoined), **Theorem 6** (explicit
64-point six-dimensional kissing configurations with Q5 and R5 as central
cross sections), **Lemma 4** (the hole graphs have clique number 16, 16,
12, 12) and **Lemma 8** (the Gegenbauer identity
F = 1 + (30/7)G₁ + (25/4)G₂ + (125/28)G₃ = (125/16)(t+3/5)²(t−1/5) with
F(1) = 16). It certifies **nothing about κ(5) or κ(6) themselves**: the
known bounds 40 ≤ κ(5) ≤ 44 and 72 ≤ κ(6) ≤ 77 are untouched. The results
concern only the four named configurations and six-dimensional
configurations that contain Q5 or R5.

## 2. Install and run

```
pip install -r requirements.txt     # sympy, pytest, numpy, scipy, networkx
make verify                         # = python3 run_all.py
```

`make verify` runs every verification script in order, parses their output,
compares every value with the hard-coded `EXPECTED` table in `run_all.py`,
and writes `results/summary.json` and `results/log.txt` (both git-ignored).
Expected runtime: about 25 seconds on a laptop (under two minutes). The
expected final line is

```
ALL CHECKS PASS
```

Any mismatch or any script failure is printed and the exit code is 1.
To run with a virtual environment's interpreter instead of `python3`, pass
the Makefile variable: `make verify PYTHON=<venv>/bin/python`.
`make test` runs only the unit tests of the exact checker; `make clean`
removes `results/` and caches. Nothing under version control is modified by
a run (`git status` stays clean).

## 3. Expected values

Hard-coded in `run_all.py` (`EXPECTED`); these are the numbers stated in the paper.

| quantity | D5 | L5 | Q5 | R5 |
|---|---|---|---|---|
| contacts (pairs at inner product 1) | 240 | 240 | 240 | 240 |
| contact-graph degree | 12 | 12 | 12 | 12 |
| antipodal points (x with −x in X) | 40 | 24 | 20 | 12 |
| rank of the equality system | 190 | 190 | 190 | 190 |
| equality kernel dimension (= rotations) | 10 | 10 | 10 | 10 |
| self-stress weights ω (multiplicity) | 1 ×240 | 1 ×240 | 20/21 ×210, 4/3 ×30 | 200/207 ×210, 280/207 ×18, 220/207 ×12 |
| λ with s_i = λ_i x_i (multiplicity) | 6 ×40 | 6 ×40 | 40/7 ×20, 44/7 ×20 | 400/69 ×20, 140/23 ×12, 440/69 ×8 |
| uniform stress ω ≡ 1 works | yes | yes | **fails** | **fails** |
| vertices of P_X | 42 | 50 | 92 | 100 |
| max \|w\|² over P_X | 5/4 | 5/4 | 5/4 | 5/4 |
| deep holes (vertices with \|w\|² = 5/4) | 32 | 32 | 32 | 32 |
| max inner product between distinct holes | 3/4 | 3/4 | 23/20 | 23/20 |
| hole-graph clique number | 16 | 16 | 12 | 12 |

| quantity | expected |
|---|---|
| Lemma 8: F = (125/16)(t+3/5)²(t−1/5) identity | True |
| Lemma 8: F(1) | 16 |
| Theorem 6: 64-point extensions of Q5 and R5 | PASS, 64 points each |
| Szöllősi reproduction (search machinery) | cloud 78, max clique 36, 4 maximum cliques, classification {Q5: 2, D5: 2} |

### Computed values

Generated from `results/summary.json` by the command shown below, not typed
by hand; every row is one check of `run_all.py`.

<!-- COMPUTED-VALUES-BEGIN -->
Run of 2026-09-03T16:17:19, release v1.0, Python 3.14.4, 150 checks, 0 failed, wall time 24.9s.

| step | check | expected | computed | pass |
|---|---|---|---|---|
| pytest verify/ | exit code | 0 | 0 | PASS |
| pytest verify/ | tests failed (33 passed) | 0 | 0 | PASS |
| verify/exact_check.py x7 | exit code | 0 | 0 | PASS |
| verify/exact_check.py x7 | d5_40.json verdict | PASS | PASS | PASS |
| verify/exact_check.py x7 | l5_40.json verdict | PASS | PASS | PASS |
| verify/exact_check.py x7 | q5_40.json verdict | PASS | PASS | PASS |
| verify/exact_check.py x7 | r5_40.json verdict | PASS | PASS | PASS |
| verify/exact_check.py x7 | q5_szollosi_40.json verdict | PASS | PASS | PASS |
| verify/exact_check.py x7 | q5_ext64.json verdict | PASS | PASS | PASS |
| verify/exact_check.py x7 | r5_ext64.json verdict | PASS | PASS | PASS |
| analysis/rigidity.py | exit code | 0 | 0 | PASS |
| analysis/rigidity.py | d5 contacts | 240 | 240 | PASS |
| analysis/rigidity.py | d5 uniform stress fails | False | False | PASS |
| analysis/rigidity.py | d5 rotation flexes valid, dim | 10 | 10 | PASS |
| analysis/rigidity.py | d5 rank of equality system | 190 | 190 | PASS |
| analysis/rigidity.py | l5 contacts | 240 | 240 | PASS |
| analysis/rigidity.py | l5 uniform stress fails | False | False | PASS |
| analysis/rigidity.py | l5 rotation flexes valid, dim | 10 | 10 | PASS |
| analysis/rigidity.py | l5 rank of equality system | 190 | 190 | PASS |
| analysis/rigidity.py | q5 contacts | 240 | 240 | PASS |
| analysis/rigidity.py | q5 uniform stress fails | True | True | PASS |
| analysis/rigidity.py | q5 rotation flexes valid, dim | 10 | 10 | PASS |
| analysis/rigidity.py | q5 rank of equality system | 190 | 190 | PASS |
| analysis/rigidity.py | r5 contacts | 240 | 240 | PASS |
| analysis/rigidity.py | r5 uniform stress fails | True | True | PASS |
| analysis/rigidity.py | r5 rotation flexes valid, dim | 10 | 10 | PASS |
| analysis/rigidity.py | r5 rank of equality system | 190 | 190 | PASS |
| analysis/stress_lp.py | exit code | 0 | 0 | PASS |
| analysis/stress_lp.py | d5 exact positive self-stress, min weight | 1 | 1 | PASS |
| analysis/stress_lp.py | d5 verdict | jammed | jammed | PASS |
| analysis/stress_lp.py | l5 exact positive self-stress, min weight | 1 | 1 | PASS |
| analysis/stress_lp.py | l5 verdict | jammed | jammed | PASS |
| analysis/stress_lp.py | q5 exact positive self-stress, min weight | 20/21 | 20/21 | PASS |
| analysis/stress_lp.py | q5 verdict | jammed | jammed | PASS |
| analysis/stress_lp.py | r5 exact positive self-stress, min weight | 200/207 | 200/207 | PASS |
| analysis/stress_lp.py | r5 verdict | jammed | jammed | PASS |
| analysis/certificates/independent_check.py | exit code | 0 | 0 | PASS |
| analysis/certificates/independent_check.py | d5 contacts | 240 | 240 | PASS |
| analysis/certificates/independent_check.py | d5 contact degree | 12 | 12 | PASS |
| analysis/certificates/independent_check.py | d5 stress weight multiplicities | {'1': 240} | {'1': 240} | PASS |
| analysis/certificates/independent_check.py | d5 lambda multiplicities | {'6': 40} | {'6': 40} | PASS |
| analysis/certificates/independent_check.py | d5 rank over Q (Bareiss) | 190 | 190 | PASS |
| analysis/certificates/independent_check.py | d5 equality kernel dim | 10 | 10 | PASS |
| analysis/certificates/independent_check.py | d5 result | jammed | jammed | PASS |
| analysis/certificates/independent_check.py | l5 contacts | 240 | 240 | PASS |
| analysis/certificates/independent_check.py | l5 contact degree | 12 | 12 | PASS |
| analysis/certificates/independent_check.py | l5 stress weight multiplicities | {'1': 240} | {'1': 240} | PASS |
| analysis/certificates/independent_check.py | l5 lambda multiplicities | {'6': 40} | {'6': 40} | PASS |
| analysis/certificates/independent_check.py | l5 rank over Q (Bareiss) | 190 | 190 | PASS |
| analysis/certificates/independent_check.py | l5 equality kernel dim | 10 | 10 | PASS |
| analysis/certificates/independent_check.py | l5 result | jammed | jammed | PASS |
| analysis/certificates/independent_check.py | q5 contacts | 240 | 240 | PASS |
| analysis/certificates/independent_check.py | q5 contact degree | 12 | 12 | PASS |
| analysis/certificates/independent_check.py | q5 stress weight multiplicities | {'20/21': 210, '4/3': 30} | {'20/21': 210, '4/3': 30} | PASS |
| analysis/certificates/independent_check.py | q5 lambda multiplicities | {'40/7': 20, '44/7': 20} | {'40/7': 20, '44/7': 20} | PASS |
| analysis/certificates/independent_check.py | q5 rank over Q (Bareiss) | 190 | 190 | PASS |
| analysis/certificates/independent_check.py | q5 equality kernel dim | 10 | 10 | PASS |
| analysis/certificates/independent_check.py | q5 result | jammed | jammed | PASS |
| analysis/certificates/independent_check.py | r5 contacts | 240 | 240 | PASS |
| analysis/certificates/independent_check.py | r5 contact degree | 12 | 12 | PASS |
| analysis/certificates/independent_check.py | r5 stress weight multiplicities | {'200/207': 210, '280/207': 18, '220/207': 12} | {'200/207': 210, '220/207': 12, '280/207': 18} | PASS |
| analysis/certificates/independent_check.py | r5 lambda multiplicities | {'400/69': 20, '140/23': 12, '440/69': 8} | {'400/69': 20, '140/23': 12, '440/69': 8} | PASS |
| analysis/certificates/independent_check.py | r5 middle stress weight | 220/207 | 220/207 | PASS |
| analysis/certificates/independent_check.py | r5 rank over Q (Bareiss) | 190 | 190 | PASS |
| analysis/certificates/independent_check.py | r5 equality kernel dim | 10 | 10 | PASS |
| analysis/certificates/independent_check.py | r5 result | jammed | jammed | PASS |
| analysis/certificates/independent_check.py | OVERALL | PASS | PASS | PASS |
| analysis/deep_holes.py | exit code | 0 | 0 | PASS |
| analysis/deep_holes.py | d5 polar vertices | 42 | 42 | PASS |
| analysis/deep_holes.py | d5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/deep_holes.py | d5 deep holes | 32 | 32 | PASS |
| analysis/deep_holes.py | d5 41st point | IMPOSSIBLE | IMPOSSIBLE | PASS |
| analysis/deep_holes.py | l5 polar vertices | 50 | 50 | PASS |
| analysis/deep_holes.py | l5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/deep_holes.py | l5 deep holes | 32 | 32 | PASS |
| analysis/deep_holes.py | l5 41st point | IMPOSSIBLE | IMPOSSIBLE | PASS |
| analysis/deep_holes.py | q5 polar vertices | 92 | 92 | PASS |
| analysis/deep_holes.py | q5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/deep_holes.py | q5 deep holes | 32 | 32 | PASS |
| analysis/deep_holes.py | q5 41st point | IMPOSSIBLE | IMPOSSIBLE | PASS |
| analysis/deep_holes.py | r5 polar vertices | 100 | 100 | PASS |
| analysis/deep_holes.py | r5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/deep_holes.py | r5 deep holes | 32 | 32 | PASS |
| analysis/deep_holes.py | r5 41st point | IMPOSSIBLE | IMPOSSIBLE | PASS |
| analysis/deep_holes.py | no file written | True | True | PASS |
| analysis/certificates/independent_deep_holes_check.py | exit code | 0 | 0 | PASS |
| analysis/certificates/independent_deep_holes_check.py | d5 polar vertices | 42 | 42 | PASS |
| analysis/certificates/independent_deep_holes_check.py | d5 deep holes | 32 | 32 | PASS |
| analysis/certificates/independent_deep_holes_check.py | d5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/certificates/independent_deep_holes_check.py | d5 m(X)^2 | 2/5 | 2/5 | PASS |
| analysis/certificates/independent_deep_holes_check.py | d5 verdict | PASS | PASS | PASS |
| analysis/certificates/independent_deep_holes_check.py | l5 polar vertices | 50 | 50 | PASS |
| analysis/certificates/independent_deep_holes_check.py | l5 deep holes | 32 | 32 | PASS |
| analysis/certificates/independent_deep_holes_check.py | l5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/certificates/independent_deep_holes_check.py | l5 m(X)^2 | 2/5 | 2/5 | PASS |
| analysis/certificates/independent_deep_holes_check.py | l5 verdict | PASS | PASS | PASS |
| analysis/certificates/independent_deep_holes_check.py | q5 polar vertices | 92 | 92 | PASS |
| analysis/certificates/independent_deep_holes_check.py | q5 deep holes | 32 | 32 | PASS |
| analysis/certificates/independent_deep_holes_check.py | q5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/certificates/independent_deep_holes_check.py | q5 m(X)^2 | 2/5 | 2/5 | PASS |
| analysis/certificates/independent_deep_holes_check.py | q5 verdict | PASS | PASS | PASS |
| analysis/certificates/independent_deep_holes_check.py | r5 polar vertices | 100 | 100 | PASS |
| analysis/certificates/independent_deep_holes_check.py | r5 deep holes | 32 | 32 | PASS |
| analysis/certificates/independent_deep_holes_check.py | r5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/certificates/independent_deep_holes_check.py | r5 m(X)^2 | 2/5 | 2/5 | PASS |
| analysis/certificates/independent_deep_holes_check.py | r5 verdict | PASS | PASS | PASS |
| analysis/certificates/independent_deep_holes_check.py | OVERALL | PASS | PASS | PASS |
| analysis/verify_lemma4_lemma8.py | exit code | 0 | 0 | PASS |
| analysis/verify_lemma4_lemma8.py | D5 antipodal points | 40 | 40 | PASS |
| analysis/verify_lemma4_lemma8.py | D5 polytope vertices | 42 | 42 | PASS |
| analysis/verify_lemma4_lemma8.py | D5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/verify_lemma4_lemma8.py | D5 deep holes | 32 | 32 | PASS |
| analysis/verify_lemma4_lemma8.py | D5 max hole inner product | 3/4 | 3/4 | PASS |
| analysis/verify_lemma4_lemma8.py | D5 clique number (networkx) | 16 | 16 | PASS |
| analysis/verify_lemma4_lemma8.py | D5 clique number (own B&B) | 16 | 16 | PASS |
| analysis/verify_lemma4_lemma8.py | L5 antipodal points | 24 | 24 | PASS |
| analysis/verify_lemma4_lemma8.py | L5 polytope vertices | 50 | 50 | PASS |
| analysis/verify_lemma4_lemma8.py | L5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/verify_lemma4_lemma8.py | L5 deep holes | 32 | 32 | PASS |
| analysis/verify_lemma4_lemma8.py | L5 max hole inner product | 3/4 | 3/4 | PASS |
| analysis/verify_lemma4_lemma8.py | L5 clique number (networkx) | 16 | 16 | PASS |
| analysis/verify_lemma4_lemma8.py | L5 clique number (own B&B) | 16 | 16 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 antipodal points | 20 | 20 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 polytope vertices | 92 | 92 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 deep holes | 32 | 32 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 max hole inner product | 23/20 | 23/20 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 clique number (networkx) | 12 | 12 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 clique number (own B&B) | 12 | 12 | PASS |
| analysis/verify_lemma4_lemma8.py | R5 antipodal points | 12 | 12 | PASS |
| analysis/verify_lemma4_lemma8.py | R5 polytope vertices | 100 | 100 | PASS |
| analysis/verify_lemma4_lemma8.py | R5 max |w|^2 | 5/4 | 5/4 | PASS |
| analysis/verify_lemma4_lemma8.py | R5 deep holes | 32 | 32 | PASS |
| analysis/verify_lemma4_lemma8.py | R5 max hole inner product | 23/20 | 23/20 | PASS |
| analysis/verify_lemma4_lemma8.py | R5 clique number (networkx) | 12 | 12 | PASS |
| analysis/verify_lemma4_lemma8.py | R5 clique number (own B&B) | 12 | 12 | PASS |
| analysis/verify_lemma4_lemma8.py | Lemma 8 identity F = 125/16 (t+3/5)^2 (t-1/5) | True | True | PASS |
| analysis/verify_lemma4_lemma8.py | Lemma 8 F(1) | 16 | 16 | PASS |
| analysis/verify_lemma4_lemma8.py | Q5 ext64 | PASS, 64 points | PASS, 64 points | PASS |
| analysis/verify_lemma4_lemma8.py | R5 ext64 | PASS, 64 points | PASS, 64 points | PASS |
| config files vs. rebuilt configurations | d5_40.json == rebuilt D5 (point sets) | True | True | PASS |
| config files vs. rebuilt configurations | l5_40.json == rebuilt L5 (point sets) | True | True | PASS |
| config files vs. rebuilt configurations | q5_40.json == rebuilt Q5 (point sets) | True | True | PASS |
| config files vs. rebuilt configurations | r5_40.json == rebuilt R5 (point sets) | True | True | PASS |
| analysis/verify_lemma4_lemma8.py | OVERALL | ALL CHECKS PASS | ALL CHECKS PASS | PASS |
| search/run_szollosi_repro.py | exit code | 0 | 0 | PASS |
| search/run_szollosi_repro.py | cloud size | 78 | 78 | PASS |
| search/run_szollosi_repro.py | max clique | 36 | 36 | PASS |
| search/run_szollosi_repro.py | number of maximum cliques | 4 | 4 | PASS |
| search/run_szollosi_repro.py | classification | {'Q5': 2, 'D5': 2} | {'Q5': 2, 'D5': 2} | PASS |
<!-- COMPUTED-VALUES-END -->

Regeneration command (from the repository root, after `make verify`):

```
python3 -c "import json; d=json.load(open('results/summary.json')); print('Run of %s, release v1.0, Python %s, %d checks, %d failed, wall time %ss.\n' % (d['generated'], d['python'], d['n_checks'], d['n_failed'], d['total_wall_time_s'])); print('| step | check | expected | computed | pass |'); print('|---|---|---|---|---|'); [print('| %s | %s | %s | %s | %s |' % (c['step'], c['check'], c['expected'], c['computed'], 'PASS' if c['pass'] else 'FAIL')) for c in d['checks']]"
```

## 4. Map from the paper to the scripts and data

| paper | script(s) that verify it | certificate / data files read |
|---|---|---|
| Theorem 1 (jamming) | `analysis/rigidity.py` (contacts, uniform stress, rank mod p, rotation flexes); `analysis/stress_lp.py` (exact strictly positive self-stress; float LP only proposes it); `analysis/certificates/independent_check.py` (second implementation, see §6; prints the ω and λ multiplicities) | `configs/{d5,l5,q5,r5}_40.json`; `analysis/certificates/{q5,r5}_stress_certificate.json` |
| Theorem 2 (no 41st point, deep holes) | `analysis/deep_holes.py`; `analysis/certificates/independent_deep_holes_check.py` (second implementation, no tolerance) | `configs/*_40.json`; `analysis/deep_holes_results.json` is the *output* (the 32 holes per configuration, exact strings) |
| Theorem 3 (Conjecture 3.1 of Cohn–Rajagopal) | no computation of its own: Lemma 5 + Lemma 10 in the paper | — |
| Lemma 4 (hole-graph clique numbers) | `analysis/verify_lemma4_lemma8.py` (two clique algorithms must agree) | configurations rebuilt from the Cohn–Rajagopal rules and compared, as point sets, with `configs/*_40.json` |
| Lemma 5 (reduction) | proof in the paper; its inputs — 240 contacts, Σx = 0, span R⁵, max \|w\|² = 5/4 — are the outputs of the Theorem 1 and Theorem 2 scripts (see `analysis/certificates/reduction_lemma_check.md`) | — |
| Theorem 6 (64-point extensions) | `verify/exact_check.py` on the shipped coordinates; rebuilt from a 12-clique of the hole graph in `analysis/verify_lemma4_lemma8.py` | `configs/{q5,r5}_ext64.json` |
| Lemma 7 (norm coupling) | proof in the paper | — |
| Lemma 8 (LP bound) and Remark 9 (sharpness) | `analysis/verify_lemma4_lemma8.py` (sympy identity, F(1) = 16, demihypercube inner products {1/5, −3/5}) | — |
| Lemma 10 (equality-case rigidity) | proof in the paper, using Theorem 2 and Lemma 4 | — |
| the configurations themselves | `verify/exact_check.py` (+ `verify/test_verify.py`, 33 tests); `configs/generate_configs.py` compares the files with the source excerpts | `configs/*.json`, `configs/sources/` |
| Szöllősi's published run (validation of the search code) | `search/run_szollosi_repro.py` | writes `experiments/szollosi_repro/` (deterministic) |

Layout: `verify/` exact checker, `configs/` exact data, `analysis/` the
verification scripts and certificates, `search/` the clique search,
`experiments/` exploratory code and logs, `docs/preprint/` the paper source
(LaTeX and PDF).

## 5. Exactness

Every accepting step is exact (Python `fractions.Fraction`, integers, or
sympy rationals and radicals). Floating point occurs only in candidate
prefilters — the vertex prefilter of `analysis/deep_holes.py` and
`analysis/verify_lemma4_lemma8.py`, and the LP proposal of
`analysis/stress_lp.py` — never in a step that accepts a value:

- `analysis/deep_holes.py` and `analysis/verify_lemma4_lemma8.py` solve the
  C(40,5) = 658 008 five-by-five systems in numpy to prefilter vertex
  candidates; every surviving candidate is rationalised, deduplicated as an
  exact tuple, and checked (feasibility and the rank-5 vertex test) in
  `Fraction`, and the count is reported from the exact set. Completeness of
  the vertex list is certified by `independent_deep_holes_check.py`, which
  enumerates all subsets in int64 integer arithmetic with a proven overflow
  bound and uses no tolerance anywhere.
- `analysis/stress_lp.py` uses a scipy linear program to propose a stress
  direction; the proposal is rationalised and the certificate (positivity
  on all 240 contacts, s_i ∥ x_i) is verified exactly before it is accepted,
  and the committed certificates are re-checked by `independent_check.py`
  without any float.

Everything else — the checker (`verify/exact_check.py`, canonical
multiquadratic representation with a decidable sign test), the ranks
(mod p and fraction-free Bareiss over Q), the clique numbers, the
Gegenbauer identity, the Szöllősi cloud (exact surd arithmetic) — is exact.
Floats printed in reports are marked with `~`. No data or certificate file
contains a floating-point number.

## 6. Independence

Three claims are checked by two independently written implementations:

- **Theorem 1.** `analysis/rigidity.py` + `analysis/stress_lp.py` versus
  `analysis/certificates/independent_check.py`. The second builds the
  equality system from the contact graph as sparse dict-of-rows scaled by
  one global denominator (not dense rows scaled per row), generates the ten
  rotation flexes from the explicit basis matrices E_ab − E_ba of so(5),
  tests s_i ∥ x_i by the vanishing of all 2×2 minors (not by a
  λ-multiple), obtains the self-stress dimension 90 by duality from the
  equality rank, and uses different primes plus an exact Bareiss rank over
  Q. An AST comparison (all identifiers, attributes and string constants
  normalised, docstrings stripped) finds **no function body shared** between
  `independent_check.py` and `rigidity.py`/`stress_lp.py`; before
  2026-09-03 two bodies were identical and two more were renamed copies.
- **Theorem 2.** `analysis/deep_holes.py` (float prefilter, exact
  re-verification) versus `independent_deep_holes_check.py` (tolerance-free
  integer enumeration, plus an exhaustive recession-cone check of
  boundedness and a closed-form cross-check for D5).
- **Lemma 4.** `networkx.find_cliques` (Bron–Kerbosch) and an own
  branch-and-bound with greedy-colouring bound, in
  `analysis/verify_lemma4_lemma8.py`; both must return 16/16/12/12.

Checked by one implementation: the Gegenbauer identity of Lemma 8 (sympy
expansion). Theorem 6 is checked by the exact checker on the shipped
coordinates and, separately, rebuilt from a hole clique. The configuration
data are tied together three ways: `configs/generate_configs.py` regenerates
them from the source excerpts, `verify_lemma4_lemma8.py` rebuilds them from
the Cohn–Rajagopal construction rules, and `exact_check.py` verifies the
kissing property directly.

## 7. What is not part of the verification

- `experiments/` — floating-point, exploratory code (the abandoned cap-layer
  branch-and-bound provers, the numerical Stage-A probe, the SDP prototypes
  in `experiments/sdp/`, which alone need `requirements-experiments.txt`,
  i.e. cvxpy) and the run logs. Not part of any claim, may run for hours,
  no documented expected output. See `experiments/README.md`.
- `search/run_clique.py` — the discovery search over multiangular clouds
  (runs s1–s4, 10–70 minutes). It is how the configurations were re-found,
  not how anything is verified; `run_all.py` does not run it. Only the
  deterministic Szöllősi reproduction `search/run_szollosi_repro.py`
  (0.1 s) is part of `make verify`.

## 8. Disclosure

Computer assistance. All code, certificates and the first draft of this
note were produced with Claude Fable 5 (Anthropic) under the author's
direction (with Claude Opus 5 and Claude Sonnet 5 as subagents for the
verifier implementations and literature retrieval, respectively). The
results were then independently re-verified from the committed code in
fresh sessions with Claude Fable 5.1, Gemini Pro and MiniMax M3, and by
the exact-arithmetic checks described in Section 5. The author is
responsible for the mathematics, the repository and this text.

Every result is reproduced by `make verify` from the committed code and
data.

## 9. Citing

Repository: <https://github.com/smitterer/kissingnumbers>. Cite as
`github.com/smitterer/kissingnumbers, release v1.0`.
