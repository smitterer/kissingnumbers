# experiments/ — exploratory code and run logs

**Nothing in this directory is part of any claim in the paper.** The
scripts here use floating point, are exploratory, may run for hours, and
have no documented expected output. `run_all.py` does not execute them,
and they are not needed to reproduce Theorems 1, 2, 6 or Lemmas 4, 8.

## Exploratory scripts (moved here from `analysis/` and `sdp/` on 2026-09-03)

| file | what it is | status |
|---|---|---|
| `conjecture31_stageA.py` | numerical (float) multistart optimisation of cap-layer sizes over Q5/R5; feasibility probe from 2026-08-19 | superseded by the analytic proof in the paper |
| `cap_bnb.py` | integer-lattice branch-and-bound for the cap-layer capacity M(X) | abandoned (too slow); route superseded |
| `cap_dfs.py` | tuple branch-and-prune prover for M(X), third design | abandoned; route superseded |
| `sdp/delsarte_lp.py` | Delsarte LP bound for dimension 5, discretised, scipy `linprog` | machinery validation only (reproduces ~46.34) |
| `sdp/bachoc_vallentin.py` | sampled three-point SDP bound prototype, needs `cvxpy` (`requirements-experiments.txt`) | machinery validation only, not a valid bound |

The three moved analysis scripts still import `rigidity.py` / `deep_holes.py`
from `../analysis/`; their `sys.path` line was adjusted when they were moved.

## Run logs and search outputs (unchanged)

- `log.md` — the experiment log (date, commit, parameters, runtime, result)
  kept during the project, including the Lasserre level-2 runs that are
  mentioned nowhere in the paper.
- `runs/` — the outputs of the `search/run_clique.py` cloud searches s3/s4
  (77 maximum cliques each, all size 40, all classified as D5/L5/Q5/R5) are
  not shipped with this repository; `run_clique.py` writes them here when
  re-run, and no claim of the paper depends on them.
- `szollosi_repro/` — ships empty. `search/run_szollosi_repro.py` writes the
  four 40-point configurations here on every `make verify` (deterministic,
  byte-identical between runs).
