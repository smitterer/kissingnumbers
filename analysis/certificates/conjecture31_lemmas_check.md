# Lemmas 4, 8, 10 and Theorem 6 of the paper: where the checks live

This file used to be a prose verification report (a blind re-derivation of
"Lemmas B, C, D" of the project's working notes, i.e. Lemmas 7, 8, 10 of the
paper, together with the inputs (F1) = Theorem 2 and (F2) = Lemma 4). It
cited working scripts (`vc.py`, `vc2.py`, `vc3.py`) that were never
committed, so nothing in it could be re-run; neither the prose nor those
scripts are part of this repository.

Every machine-checkable claim it made is now verified by committed scripts,
all of which `run_all.py` executes and compares against hard-coded expected
values:

| claim (paper numbering) | script | data it reads |
|---|---|---|
| (F1) max \|w\|² over P_X = 5/4, attained at exactly 32 vertices; 42/50/92/100 vertices (Theorem 2) | `analysis/deep_holes.py`; independently `analysis/certificates/independent_deep_holes_check.py` (no tolerance anywhere) | `configs/{d5,l5,q5,r5}_40.json` |
| (F2) hole-graph clique numbers 16/16/12/12 (Lemma 4); max inner product between distinct holes 3/4 (D5, L5) and 23/20 (Q5, R5); the two 16-cliques of D5 are the parity classes | `analysis/verify_lemma4_lemma8.py` — two clique algorithms (`networkx.find_cliques` and an own branch-and-bound) must agree | configurations rebuilt from the Cohn–Rajagopal rules and compared, as point sets, with `configs/*_40.json` |
| Gegenbauer identity F = (125/16)(t+3/5)²(t−1/5) with F = 1 + (30/7)G₁ + (25/4)G₂ + (125/28)G₃, F(1) = 16; demihypercube inner products {1/5, −3/5} (Lemma 8, Remark 9) | `analysis/verify_lemma4_lemma8.py` (sympy, exact) | — |
| 64-point extensions of Q5 and R5 (Theorem 6) | `verify/exact_check.py` on the shipped coordinates; rebuilt independently from a 12-clique of the hole graph in `analysis/verify_lemma4_lemma8.py` | `configs/{q5,r5}_ext64.json` |
| Lemma 7 (norm coupling) and Lemma 10 (equality-case rigidity) | proved in the paper; they use only the verified inputs above | — |
