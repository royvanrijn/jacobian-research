# Lower-height six-fibre cohort

The completed 9 September experiment retains subgroup lower bounds
**22, 24, 17, 22, 25, 23**, reviewed on 12 September. These are dated cohort
results; use [MATH_STATUS.json](../../MATH_STATUS.json) for claim authority
and the [inventory](../INVENTORY.md) for subsequent results on these curves.
No exact rank, public novelty, conductor record or general success rate follows.

| Family | Native parameter | Initial seed bound | Seed calls | Final dated bound |
|---|---|---:|---:|---:|
| 074d9 | 25/3 | 18 | 1 | 22 |
| 07ca9 | 142/561 | 18 | 1 | 24 |
| 08234 | 383/388 | 17 | 86 | 17 |
| 08f72 | 254/17 | 18 | 1 | 22 |
| 103b2 | -95/728 | 18 | 1 | 25 |
| 11952 | 238/27 | 18 | 1 | 23 |

The [selector](../cas/select_fresh6_lowheight_r17.py) restricts the same 6,144
saved score rows to parameter height at most 1,024. Unchanged per-family score
order and frozen rational-isomorphism/address exclusions leave 366 eligible
addresses. Validation-band values and measured ranks do not enter selection.
The six curves are pairwise rationally nonisomorphic. Five first-chart M18
successes in this selected cohort do not establish a general height advantage.

The exact-maximum-class seed constructor made 91 point calls without point
timeouts. Limits were 1,800 seconds/2 GiB per search, 5 seconds/1 GiB per map,
and 10 seconds per point call at height 125,000. Historical geometry/map/point
replays and the six standalone native-point certificates passed.
[Seed packets](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_first_m18_v1/result.json),
[original queue](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_first_m18_v1/amplification-queue.json)
and [reservation snapshot](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_reserved_results_v1.json)
remain retained. The 08234 bounded miss gives no rank upper bound.

## Reconcile saved clouds before searching again

All five winning-parent passes completed. Full-cloud audits stopped further
work until the saved points were reconciled and finite independence certified.

| Family | Terminal → reconciled basis | Point calls | Reconciled packet |
|---|---|---:|---|
| 074d9 | 19 → 22 | 1 | [M22](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/074d9-lowheight-001/reconciled.json) |
| 07ca9 | 19 → 20 | 1 | [M20](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/07ca9-lowheight-001/reconciled.json) |
| 08f72 | 19 → 22 | 1 | [M22](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/08f72-lowheight-001/reconciled.json) |
| 103b2 | 21 → 22 | 3 | [M22](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/103b2-lowheight-001/reconciled.json) |
| 11952 | 19 → 23 | 1 | [M23](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/11952-lowheight-001/reconciled.json) |

Seven amplification point calls were followed by saved-witness reconciliation
that added **twelve independent directions across five fibres with zero
additional point searches**. A second exact reconciliation agreed. Retain
all returned points when admission stops at a first gain; a terminal basis
can omit directions already present in its audited cloud.

[Continuation preparation](../cas/prepare_reconciled_v3_continuation.py)
uses the reconciled basis and rechecks its finite certificate. Complementary
banks contain at most sixteen exact maximum classes outside the productive
parent span, checked by both exact CVP solvers. A new basis needs its own
verified landscape; earlier compatible chart receipts remain reusable.

## Completed complementary-parent passes

All five first passes and their historical full replays completed without
point or map timeouts. The [continuation ledger](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/continuation-queue.json)
preserves their inputs. It is not a current scheduling instruction.

| Family | Subgroup before → after | Calls | Result and retained cursor |
|---|---|---:|---|
| 11952 | 23 → 23 | 100 | [Result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/11952-M23-complement/result.json); centre 50 of 636 |
| 074d9 | 22 → 22 | 100 | [Result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/074d9-M22-complement/result.json); centre 50 |
| 08f72 | 22 → 22 | 100 | [Result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/08f72-M22-complement/result.json); centre 50 of 486 |
| 103b2 | 22 → 25 | 100 | [M25 result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/103b2-M22-complement/result.json); centre 38, factor_free, of 682 |
| 07ca9 | 20 → 21 → 24 | 1 | [Run](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/07ca9-M20-complement/result.json) · [M24 reconciliation](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/07ca9-M20-complement/reconciled.json) |

103b2 gained at cumulative calls 3, 14 and 23, then completed 77 calls at M25.
07ca9 stopped at M21 when its cloud audit found more directions; reconciliation
added three points to reach M24 without another point search. Bounded misses
prove no upper bounds. These subgroup certificates do not prove exact rank.

The 11952 and 103b2 conductor audits each hit their 30-second factorization
caps. Exact conductors remain UNKNOWN in those audits. The retained
[11952 bounds](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/11952-conductor-bounds.json)
include local lower bound 4573532994030; [103b2's bounds](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/103b2-conductor-bounds.json)
and both timeout receipts also remain available.

## Replay and history

The [shared retained-rank replay](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md#retained-rank-replay)
checks all subgroup endpoints above without search. It does not rerun the
original full cloud/map/selection replays or the conductor producer. Raw
amplification witnesses remain under `artifacts/local/elliptic-curves/fresh6-lowheight-amplification-v1/`.
The [review receipt](../../archive/repository-cleanup-2026-09-12/fresh6-source-review/REVIEW.json)
preserves the old chronology, timings and next-candidate instructions as history.
