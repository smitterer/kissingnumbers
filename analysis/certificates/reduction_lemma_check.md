# Reduction lemma (Lemma 5 of the paper): where the checks live

This file used to be a prose verification report (blind protocol) of the
reduction lemma; that prose is not part of this repository.

Lemma 5 is proved in the paper (Section 4). It is not a computation; its
three computational inputs are verified by committed scripts, all run by
`run_all.py`:

| input | script |
|---|---|
| each configuration has exactly 240 contacts and a 12-regular contact graph (used for the forced centrality, part 1) | `analysis/rigidity.py`; independently `analysis/certificates/independent_check.py` |
| Σᵢ xᵢ = 0 and the points span R⁵ (compactness of P_X) | `analysis/deep_holes.py` (assertion), `analysis/certificates/independent_deep_holes_check.py` (recession cone {0} by exhaustive ray enumeration, barycentre 0), `analysis/certificates/independent_check.py` (span rank 5) |
| max \|w\|² over P_X = 5/4 < 2, hence h² ≥ 3/4 for every off-section point (parts 1–3) | `analysis/deep_holes.py`; independently `analysis/certificates/independent_deep_holes_check.py` |
