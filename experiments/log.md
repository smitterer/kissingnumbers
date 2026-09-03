# Experiment log

Every experiment gets an entry: date, git commit, command/seed/parameters,
runtime, raw result, and where the artifacts live. Interpretation goes into the paper
and the project's working notes, not here.

| date | what | seed/params | runtime | result |
|------|------|-------------|---------|--------|
| 2026-08-18 | Config extraction: D5/L5/Q5/R5 from arXiv:2412.00937 Table 2.2 + Szollosi Q5 (arXiv:2301.08272 appendix) | deterministic (generate_configs.py) | <1s | all 5 configs: 40 distinct points, exact norms, profiles match Table 2.1 (generation-side check; independent verify/ pending) |
| 2026-08-18 | verify/exact_check.py on all 5 configs | deterministic | ~1s | all PASS exactly, worst-pair margin exactly 0 |
| 2026-08-18 | analysis/rigidity.py on d5,l5,q5,r5 | deterministic, primes 1000003.. | 0.8s total | D5: jammed (reproduces CJKT). L5: jammed (uniform stress, rank 190). Q5/R5: uniform stress fails, rank 190 |
| 2026-08-18 | analysis/stress_lp.py on q5,r5 | deterministic (float LP heuristic + exact verify) | 14s | exact strictly positive self-stresses found and verified: min weights 20/21 (Q5), 200/207 (R5) => both infinitesimally jammed, PENDING independent verifier |
| 2026-08-18 | search/run_szollosi_repro.py | deterministic | 0.05s | reproduces arXiv:2301.08272 Sec.4 exactly: cloud 78, omega=36, 4 max cliques, 2x D5-profile + 2x Q5; all 4 configs PASS verify/ |
| 2026-08-18 | analysis/deep_holes.py all 4 configs | deterministic | ~4s | all: max poly vertex norm^2=5/4 < 2, m^2=2/5, 32 deep holes; VERIFIED by independent exact re-enumeration |
| 2026-08-18 | run_clique.py s3_q5cloud(_full), Q5-basis Method-1 cloud, 9 angles | deterministic | 247s + reclass | cloud 285, omega=35 (no 41), 77 max cliques = 1 D5 + 10 L5 + 11 Q5 + 55 R5, all PASS; no fifth profile in this coverage |
| 2026-08-18 | sdp/delsarte_lp.py n=5 | float, grid 4001 | ~2s | 46.33757 (plateau from degree 10); literature 46.34 — machinery validated |
| 2026-08-18 | sdp/bachoc_vallentin.py n=5 d=10 sampled | float, 241+2837 samples | 40s | 45.14 (sampled relaxation; published exact d=10 ~45.4) |
| 2026-08-18 | CLRS three_point_spherical_codes(5, 1/2, d, d), 256 bits | deterministic | 4.3/14.8/97.9s | d2=d3=3: 48.0; 5: 46.3376 (=Delsarte); 7: 46.1535 — trivariate part activates from d3~7; MV 44.9989 needs d3~14 |
| 2026-08-18 | run_clique.py s4_r5cloud(_full), R5-basis Method-1 cloud, 9 angles | deterministic | 526s + reclass | cloud 285, edges 33644, omega=35, 77 max cliques = 1 D5 + 10 L5 + 11 Q5 + 55 R5 (identical to s3); no fifth profile, no 41-point |
| 2026-08-18 | las2 smoke n=4 d1=delta=d2=6, 256 bits, 4 threads | deterministic | 40s, 86 CPU-s, 1.1GB | obj = 26.0 (valid weak bound; sharp value 24 needs full degree) — authors' level-2 pipeline operational locally |
| 2026-08-18 | las2 ladder rung 1: n=5 d=6, n=4 d=8, n=5 d=8 (d1=delta=d2=d) | deterministic | 40/193/161s, 956 CPU-s total, 2.4GB peak | n=5: 48.0 at both d=6,8 (plateau); n=4: 26.0 at d=8 (no change from d=6). Level-2 activates only at higher degree; cost ~x4 per +2 degree |
| 2026-08-18 | CLRS 3pt n=5 d2=d3=9,10, 256 bits, 4 threads | deterministic | 302/584s, 5.4GB peak | 45.312 / 45.163 — monotone convergence toward MV 44.9989 (their truncation 14). Three-point VALIDATION COMPLETE per plan; no further degree pushed (cap reserved for level-2) |
| 2026-08-18 | run_clique.py s1 (7 angles), s2 (9 angles), A4 basis | deterministic | ~70min + 10min, killed | INCOMPLETE: clouds 1086 / 3414 — beyond pure-Python exact pipeline; omega unknown, no claims; next step would be cliquer (Ostergard) integration |
| 2026-08-18 | las2 n=5 d1=delta=d2=10, 256 bits (authors' level-2 code) | deterministic | 1831s solve (~30min), 6.4GB peak | **obj = 46.33757256** — the previously UNREPORTED dim-5 level-2 value at this degree; equals the Delsarte LP plateau. Ladder: 48.0 (d=6,8) -> 46.3376 (d=10). n=4 d=12 calibration cancelled (90min zonal build, low value); n=5 d=12 launched as final rung |
| 2026-08-19 | las2 n=5 d1=delta=d2=12, 256 bits, 6 threads | deterministic | 13977s wall, 12.8 CPU-h, 8.0GB peak | **obj = 46.30520292**. Final ladder: 48.0 / 48.0 / 46.3376 / 46.3052 (d=6/8/10/12). d=14 extrapolates to 140+ CPU-h > cap; ladder closed. ~16 of 48 CPU-h used |
| 2026-08-19 | conjecture31_stageA.py (originally under analysis/, later moved to experiments/): multistart penalty minimisation of cap-layer sizes N over D5, L5, Q5, R5 (Stage A probe for Conjecture 3.1) | float, L-BFGS-B, 200 restarts per case, seeds 1000/2000/3000 + hash(name) | not recorded | controls D5/L5 reach N=16 to machine precision and stall at 17; Q5/R5 reach N=12 to machine precision and stall from N=13 on (residual 7e-3 at 13 rising to 6e-2 at 16); numerical only, no claim. Entry written 2026-09-03 from the Stage-A feasibility notes (not part of this repository), the run's stdout was not kept |
