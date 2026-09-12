# Censoring and prospective-population ramification audit

Retrospective analysis of the unchanged prospective 2,080-fibre population. Historical controls are excluded. No new arithmetic, searches or class groups.

## LOCAL completion by BASE discriminant-size quintile

Quintile boundaries use BASE information available on every fibre, without rank or LOCAL outcomes. Boundaries: [457.05237123752363, 476.0531992042401, 492.32547351738253, 513.0141645464843].

| Quintile | Fibres | LOCAL PASS | Completion | Median log2 \|Delta\| |
| --- | --- | --- | --- | --- |
| 1 | 416 | 235 | 0.5649 | 440.1903 |
| 2 | 416 | 186 | 0.4471 | 466.1009 |
| 3 | 416 | 190 | 0.4567 | 484.4978 |
| 4 | 416 | 160 | 0.3846 | 500.5551 |
| 5 | 416 | 127 | 0.3053 | 539.5511 |

## Completion by certified lower-bound stratum

| Rank LB | Fibres | LOCAL PASS | Completion | Known ramification >=9 | Worst-case population fraction |
| --- | --- | --- | --- | --- | --- |
| 17-18 | 1574 | 649 | 0.4123 | 121 | [0.07687420584498093, 0.6645489199491741] |
| 19-20 | 351 | 169 | 0.4815 | 28 | [0.07977207977207977, 0.5982905982905983] |
| 21-22 | 125 | 63 | 0.5040 | 9 | [0.072, 0.568] |
| 23 | 25 | 15 | 0.6000 | 2 | [0.08, 0.48] |
| 24 | 4 | 1 | 0.2500 | 0 | [0.0, 0.75] |
| 25+ | 1 | 1 | 1.0000 | 0 | [0.0, 0.0] |

## Parent/family completion

| Family | Fibres | LOCAL PASS | Completion | Median ramified primes (known) |
| --- | --- | --- | --- | --- |
| 074d9 | 320 | 151 | 0.4719 | 6.0000 |
| 07ca9 | 320 | 166 | 0.5188 | 7.0000 |
| 08234 | 320 | 111 | 0.3469 | 6.0000 |
| 08f72 | 320 | 142 | 0.4437 | 6.0000 |
| 103b2 | 320 | 148 | 0.4625 | 6.0000 |
| 11952 | 320 | 145 | 0.4531 | 6.0000 |
| x1092-class1 | 80 | 3 | 0.0375 | 5.0000 |
| x1092-original | 80 | 32 | 0.4000 | 6.0000 |

## All prospective: ramification strata

| Ramified primes | Fibres | LB >=23 | Tail fraction | Median rank LB | Median BASE log2 \|Delta\| |
| --- | --- | --- | --- | --- | --- |
| 0-4 | 168 | 5 | 0.0298 | 17.0000 | 478.6748 |
| 5-6 | 312 | 7 | 0.0224 | 17.0000 | 482.6167 |
| 7-8 | 258 | 3 | 0.0116 | 17.0000 | 476.9977 |
| 9+ | 160 | 2 | 0.0125 | 17.0000 | 472.0977 |
| UNKNOWN | 1182 | 13 | 0.0110 | 17.0000 | 488.8310 |

## Prospective 11952 only: ramification strata

| Ramified primes | Fibres | LB >=23 | Tail fraction | Median rank LB | Median BASE log2 \|Delta\| |
| --- | --- | --- | --- | --- | --- |
| 0-4 | 35 | 2 | 0.0571 | 17.0000 | 476.0594 |
| 5-6 | 51 | 2 | 0.0392 | 18.0000 | 478.5613 |
| 7-8 | 42 | 0 | 0.0000 | 17.0000 | 481.4725 |
| 9+ | 17 | 0 | 0.0000 | 17.0000 | 464.9066 |
| UNKNOWN | 175 | 2 | 0.0114 | 17.0000 | 482.7286 |

## Descriptive diagnostics

{
  "completed_11952_ramification_association": {
    "n": 145,
    "spearman_ramification_vs_rank_LB": -0.13435254534365634
  },
  "completed_case_ramification_association": {
    "n": 898,
    "spearman_ramification_vs_rank_LB": 0.00802020131534974
  },
  "completion_vs_base_size_spearman": -0.1798294980299543,
  "within_family_and_size_rank_association": {
    "centered_rank_association": -0.05442432735614881,
    "informative_strata": 30,
    "pairs": 883
  },
  "within_family_rank_association": {
    "centered_rank_association": -0.03001037633632978,
    "informative_strata": 7,
    "pairs": 895
  }
}


Worst-case fraction intervals assign every unknown to low or high ramification. No missing-at-random assumption, inverse weighting or imputation.

All ranks are certified lower bounds; no exact ranks, p-values, predictive validation, Selmer dimensions or class ranks inferred.

A completion association is evidence of observed censoring structure, not proof of why a particular worker timed out. Conditioning on completed rows, family and coarse size bands does not remove unobserved selection bias. This audit cannot establish or refute a class-group mechanism.
