# Second retained six-fibre cohort

The completed 9 September experiment retains subgroup lower bounds
**20, 20, 22, 17, 26, 22**, reviewed on 12 September. These are dated cohort
results; use [MATH_STATUS.json](../../MATH_STATUS.json) for claim authority
and the [inventory](../INVENTORY.md) for subsequent results on these curves.
No exact rank, public novelty, conductor record or general success rate follows.

| Family | Native parameter | Initial seed bound | Seed calls | Final dated bound |
|---|---|---:|---:|---:|
| 074d9 | 88/2551 | 18 | 9 | 20 |
| 07ca9 | -2475/2848 | 18 | 3 | 20 |
| 08234 | 2570/2143 | 18 | 5 | 22 |
| 08f72 | -2433/3479 | 17 | 98 | 17 |
| 103b2 | 2815/1088 | 18 | 7 | 26 |
| 11952 | 327/1403 | 18 | 11 | 22 |

## Selection and first seeds

The [selector](../cas/select_fresh6_retained_r17_second_cohort.py) reuses the
6,144 saved score rows and unchanged per-family score order. Its frozen
exclusions include the first cohort reservation, rationally isomorphic curves
and scheduled addresses. The 5,822 eligible addresses are snapshot-relative,
not independent curves. Validation-band values and measured ranks do not
enter selection; all six selected curves are pairwise rationally nonisomorphic.

Exact maximum-class seed construction made 133 point calls with no point
timeouts. The limits were 1,800 seconds/2 GiB per search, 5 seconds/1 GiB per
map, and 10 seconds per point call at height 125,000. Full geometry/map/point
replays passed historically. Separate finite certificates establish five M18
seeds and one M17 lower bound. A bounded seed miss gives no rank upper bound.

- [Generic M17 packets](../../artifacts/generated-results/elliptic-curves/fresh6_second_generic17_v1/result.json) · [six exact seed packets](../../artifacts/generated-results/elliptic-curves/fresh6_second_first_m18_v1/result.json).
- [Original seed queue](../../artifacts/generated-results/elliptic-curves/fresh6_second_first_m18_v1/amplification-queue.json) · [reservation snapshot](../../artifacts/generated-results/elliptic-curves/fresh6_second_reserved_results_v1.json).

## Certified amplification endpoints

All five winning-parent passes completed. Four branches stayed at M18 after
six calls each; 103b2 reached M20 in eighteen calls. All five complementary
maximum-parent branches then completed. Their retained historical replays
cover rational CVP, maps, points, gain provenance and complete-cloud audits.

| Endpoint | Calls in this pass | Primary evidence |
|---|---:|---|
| 103b2 M18 → M20 | 18 | [Winning-parent result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-fresh-002/result.json) |
| 103b2 M20 → M22 → M23 | 10; reconciliation uses 0 new searches | [Full run](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-complement/result.json) · [M23 packet](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-complement/reconciled.json) |
| 103b2 M23 → M26 | 100; gains at 31, 54, 97 | [M26 result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-M23/result.json) |
| 11952 M18 → M22 | 100 | [M22 result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/11952-complement/result.json) |
| 08234 M18 → M22 | 100 | [M22 result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/08234-complement/result.json) |
| 074d9 M18 → M20 | 100 | [M20 result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/074d9-complement/result.json) |
| 07ca9 M18 → M20 | 100 | [M20 result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/07ca9-complement/result.json) |

103b2's cloud audit stopped the first complementary pass at a retained M22
basis. Replaying its saved points certified M23 without further searching.
The subsequent rebuilt pass reached M26 after the earlier M25 progress
snapshot; the sealed terminal determines its endpoint. These are subgroup
lower bounds, with no exact rank or upper bound. No point or map timeouts
occurred in these amplification passes.

## Completed continuations and conductor limit

Four cached 100-call M26 batches added no direction. Including the gaining
M23 restart, that branch completed 500 calls. Its
[latest retained result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-M26-cached-v4/result.json)
and [suffix](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-M26-cached-v4/suffix-queue.json)
record centre 201, `factor_free`, in the 689-centre M26 landscape; that centre's
preconditioned box had already been searched. Earlier queues link to descendants.

A separate fixed 128-parity sample supplied sixteen norm10 parents, thirteen
outside 103b2's maximum-class hyperplane. Its
[100-call lower-shell pass](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-sampled-shell/result.json)
also stayed at M26, with a [separate suffix](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-sampled-shell/suffix-queue.json).
Generic parity coverage does not predict rational points or specialized gains.
These bounded misses and completed prefixes must not be treated as unsearched
work or rank upper bounds. The [five-branch ledger](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/complement-queue.json)
retains the other branch checkpoints as historical inputs, not a scheduled queue.

The 103b2 conductor audit hit its 30-second/1 GiB cap with a 123-digit residual
cofactor. Exact conductor remains UNKNOWN in that audit. Its local bounds and
timeout are retained under `artifacts/local/elliptic-curves/second-fresh103b2-conductor-v1/`.

## Replay and history

The [shared retained-rank replay](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md#retained-rank-replay)
checks all subgroup endpoints above without search. Its compact witness is
separate from the original full protocol/map replays and conductor producer.
Version checks and cached PASS receipts are not fresh full replays.
The [review receipt](../../archive/repository-cleanup-2026-09-12/fresh6-source-review/REVIEW.json)
preserves the original chronology, timings and obsolete “not started” statements.
The [lower-height cohort](LOWHEIGHT_FRESH6_SEED_COHORT_2026-09-09.md) is a separate
completed experiment.
