# Completed arithmetic census and separate historical 11952 controls

All ranks below are **certified lower bounds**, not exact ranks. This report runs no point searches or class-group probes.

Prospective census: 2080 fibres; BASE {'PASS': 2080}; LOCAL {'PASS': 898, 'UNKNOWN_TIMEOUT': 1182}.

## Cohort isolation and replay

The 2,080 rows are labeled `cohort = prospective_broad_2080`, `selection_mode = frozen_prospective` in `prospective_profiles.json`. The historical rows are labeled `cohort = historical_external`, `selection_mode = retrospective_known_high_rank` in `historical_profiles.json`. These are final-analysis views: the original worker schemas/checkpoints remain unedited.

The full 445-row repository inventory roster was frozen before selecting every 11952 fibre with certified lower bound at least 27. All five qualifying certificates were replayed exactly: point membership, independent finite-reduction columns, rational 2-torsion exclusion, specialization isomorphism and point transport. The rank-28 entry is a public-point reproduction; its original local-search bound is 27.

Historical arithmetic calls the **identical unchanged** BASE/LOCAL worker and controller with the census's Sage version, 60/180-second per-worker budgets and memory policy. At most two historical workers run together. No cached factorization or known points are supplied to the workers. Every timeout/error is retained as UNKNOWN; no refills.

[Arithmetic replay](ARITHMETIC_REPLAY.json) checks all five BASE results and all three completed LOCAL results, using integral-monic 2-division fields, certified PARI maximal orders, prime-factor proofs and exact local reductions. This shares Sage/PARI with the worker; it is not a separate full descent. Neither timed-out LOCAL result is filled in.

The hashes and file set of all 9232 existing census files are checked unchanged. This includes the original plan, source snapshots, worker inputs, results, summaries, frozen CLASS control selection and launch state. The original prospective stratification and Spearman diagnostics are reproduced exactly. Historical rows never enter those computations or any population histogram.

## Historical controls: per fibre

| 11952 parameter | Rank LB | BASE | LOCAL | u+n | #Phi_m | #Phi_a | log2 \|D_K\| | Ramified primes | Signature | Root | log2 conductor | Forced g lower* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -2448/11 | 27 | PASS | PASS | 8 | 6 | 0 | 384.677 | 11 | [3, 0] | -1 | 407.222 | 19 |
| 110314/102227 | 28 | PASS | PASS | 10 | 9 | 0 | 422.117 | 9 | [1, 1] | 1 | 460.187 | 18 |
| 2828/2015 | 27 | PASS | PASS | 8 | 5 | 1 | 401.953 | 12 | [3, 0] | -1 | 425.769 | 19 |
| 4286/1881 | 27 | PASS | UNKNOWN_TIMEOUT | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |
| 2012/211 | 27 | PASS | UNKNOWN_TIMEOUT | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN | UNKNOWN |

## Separate retrospective/prospective comparison

Medians use available values only, never zero for UNKNOWN. Groups overlap (the LB-25 fibre also belongs to the tail and full population); they are not independent samples. The historical panel is selected after knowing high lower bounds and is not a prospective cohort.

| Group | Fibres | LOCAL PASS | Median u+n | Median #Phi_m | Median #Phi_a | Median log2 \|D_K\| | Median ramified primes | Median log2 conductor | Median forced g lower* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Historical 11952, LB >=27 | 5 | 3 | 8.000 | 6.000 | 0.000 | 401.953 | 11.000 | 425.769 | 19.000 |
| Broad 11952 @ 921/653, LB 25 | 1 | 1 | 9.000 | 5.000 | 1.000 | 353.753 | 5.000 | 375.469 | 16.000 |
| Prospective LB >=23 tail | 30 | 17 | 7.000 | 4.000 | 1.000 | 376.061 | 5.000 | 390.448 | 16.000 |
| Prospective LB 23 | 25 | 15 | 7.000 | 4.000 | 1.000 | 376.061 | 6.000 | 390.448 | 16.000 |
| Prospective LB 24 | 4 | 1 | 9.000 | 7.000 | 0.000 | 383.011 | 3.000 | 411.834 | 15.000 |
| Prospective LB 25 | 1 | 1 | 9.000 | 5.000 | 1.000 | 353.753 | 5.000 | 375.469 | 16.000 |
| Full prospective population | 2080 | 898 | 5.000 | 3.000 | 1.000 | 405.906 | 6.000 | 420.211 | 12.000 |
| Prospective 11952 only (same-family diagnostic) | 320 | 145 | 6.000 | 3.000 | 1.000 | 408.864 | 6.000 | 423.455 | 12.000 |

The JSON/CSV comparison retains the nonmissing count and range for each metric, BASE/LOCAL status counts, root-number counts and field signatures.

*`forced g lower = max(0, rank_LB-(u+n))` uses the known rank lower bound. It is **not independent equation-derived evidence** and is not a class-group upper bound.

## Unchanged prospective strata

| Rank LB stratum | Fibres | Known u+n | Median u+n | Median log2 \|D_K\| | Median log2 conductor |
| --- | --- | --- | --- | --- | --- |
| 17-18 | 1574 | 649 | 5.000 | 418.938 | 432.646 |
| 19-20 | 351 | 169 | 6.000 | 378.478 | 392.804 |
| 21-22 | 125 | 63 | 6.000 | 371.026 | 381.795 |
| 23 | 25 | 15 | 7.000 | 376.061 | 390.448 |
| 24 | 4 | 1 | 9.000 | 383.011 | 411.834 |
| 25+ | 1 | 1 | 9.000 | 353.753 | 375.469 |

## Unchanged prospective descriptive diagnostics

| Metric | Pairs | Spearman vs final rank LB | Spearman vs follow-up gain |
| --- | --- | --- | --- |
| bk_local_term | 898 | 0.089 | 0.043 |
| log2_abs_minimal_discriminant | 2080 | -0.505 | -0.151 |
| log2_abs_field_discriminant | 898 | -0.535 | -0.146 |
| log2_conductor | 898 | -0.549 | -0.148 |

The completed-case local-term association is weak. Discriminant/conductor sizes show a stronger negative association with discovered lower bound. Neither establishes a rank mechanism: the original search was adaptive, ranks are censored lower bounds, parent composition can confound pooled diagnostics, and LOCAL completion is heavily censored by the fixed timeout. No p-values, exact ranks, full Selmer dimensions, or class-group upper bounds are inferred. A weak local association does not locate a missing signal in g.

`check` byte-rebuilds every analysis output and verifies the untouched census and bound sources.
