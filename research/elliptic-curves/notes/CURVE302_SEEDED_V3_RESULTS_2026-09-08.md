# Curve302: completed V3 amplification and transfer results

## Verified outcomes

The [frozen V3 calibration](ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md) completed
generic **17→18→…→31 in 1169 charts**, following its successful fixed-M30
gate (101 charts). This supersedes the pending-run status in the old handoff,
not the immutable V1/V2 experiment records.

The separate two-seed experiment is **COMPLETE_TWO_SEED_AMPLIFIER**.
Both cases have status `PASS_INDEPENDENT_SEEDED_V3_REPLAY` and stop reason
`TARGET_LOWER_BOUND_REACHED`: a certified rank-18 seed subgroup acquires
thirteen further independent directions. The retained
[compact replay summaries](../../artifacts/generated-results/elliptic-curves/curve302_seeded_v3_amplifier_v1.json)
are a verbatim JSON-value export of the local controller summary, not a new
arithmetic replay.

| Certified transition | recovered-strict-02 charts | recovered-strict-03 charts |
| --- | ---: | ---: |
| 18→19 | 7 | 1 |
| 19→20 | 17 | 3 |
| 20→21 | 27 | 17 |
| 21→22 | 1 | 50 |
| 22→23 | 27 | 1 |
| 23→24 | 37 | 120 |
| 24→25 | 50 | 51 |
| 25→26 | 18 | 43 |
| 26→27 | 12 | 127 |
| 27→28 | 290 | 5 |
| 28→29 | 76 | 32 |
| 29→30 | 141 | 95 |
| 30→31 | 152 | 92 |
| Total | 855 | 637 |

Every stage records its rebuilt landscape and separate mod-3/mod-5 rank
certificates, in addition to the incremental mod-2 audits. Each run scores
all extensions of its retained anchors, from 64 cosets at M18 to 262144
at M30. These are not all cosets of the full M_r/2M_r quotient.

## What the seed comparison measures

The [retained-atlas unlock diagnostic](../../artifacts/generated-results/elliptic-curves/curve302_unlock_seed_closure_v2.json)
ties these two seeds at minimax numerator 102542623, scaled height
25.635655750. That is exact dynamic programming on a finite, retrospective
atlas, not a global CVP optimum or a runtime prediction. Its v1 artifact is
preserved; v2 carries the upgraded schema.

The measured chart totals differ despite that tie. This is a descriptive
two-case result, not proof that one seed is generally superior. Known residual
labels are used to supply the initial seed. After its proof is frozen, the
search uses the certified subgroup and frozen V3 policy, not later residual
points as an execution oracle. **This is retrospective known-seed amplification,
not prospective seed selection, a new rank record, or an exact-rank proof.**

## Completed transfer controls

| Experiment | Cases | Executed charts | Certified outcome |
| --- | ---: | ---: | --- |
| Native11952 V3 control | 1 | 82 | 17→17 |
| Warm11952: 41, 72, 186 | 3 | 1508, 1608, 1512 | each 27→27 |
| Determinant1092 V3 pilot | 8 | 82 each; 656 total | each 17→17 |
| Same eight, V4 fresh bootstrap | 8 | 512 each; 4096 total | each 17→17 |

The warm roster and the two eight-fibre panels completed independent replay;
these are completed finite no-gain exposures, not ongoing startup failures.
See the [warm runbook](V3_WARM_RUNBOOK.md),
[V3 pilot](DET1092_V3_EIGHT_PILOT.md), and
[V4 bootstrap](DET1092_V4_WIDE_BOOTSTRAP.md) for frozen populations and budgets.
The V4 experiment scores the full 2^17 quotient but searches only 512 fresh
exact-CVP-certified parities per fibre. Its survivor TSV is not the full quotient.
No automatic reserve expansion or V4-to-V3 cascade was released.

These controls do not prove rank exactly 17/27, arithmetic absence of a jump,
or poor visibility as the unique cause. Successful amplification on302 and
failed bootstrap elsewhere are different conditional observations.

## Evidence and reproduction

All paths below are relative to `research/`; raw checkpoints and failed
startup histories remain unchanged.

- Seeded controller: `artifacts/local/elliptic-curves/curve302-seeded-v3-amplifier-v1/summary.json`.
  Each seed subdirectory retains `protocol.json`, `seed-proof.json`,
  `seeded-verified.json` and `replay-M17/` transcripts.
- [Runner and independent case replay](../cas/run_curve302_seeded_v3_amplifier.py);
  [Sage preflight](../cas/check_curve302_seeded_v3_amplifier_preflight.py).
  Both preflight seeds passed at rank18 before resumption; preflight itself
  executes zero charts.
- V3 controls: `artifacts/local/elliptic-curves/det1092-v3-eight-pilot/summary.json`.
- V4 controls: `artifacts/local/elliptic-curves/det1092-v4-wide-bootstrap-v2/summary.json`.
- Warm completion: `artifacts/local/elliptic-curves/v3-warm-start-overnight-v1/sessions/a6c275cb-e1ba-4c5a-93ba-d0d2a86cf02d/session.json`.
- Native control: `artifacts/local/elliptic-curves/v3-transfer-11952-v3/control-native11952/`.

The compact seeded export records full protocol and terminal SHA256 values and
all thirteen stage summaries for each case. Documentation refresh does not
rerun searches, alter policies, or replace the underlying exact certificates.
