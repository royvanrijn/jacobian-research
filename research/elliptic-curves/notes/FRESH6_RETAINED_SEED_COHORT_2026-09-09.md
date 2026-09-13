# First retained six-fibre cohort

These are the completed 9 September experiments, reviewed on 12 September.
The final retained lower bounds are **24, 23, 17, 24, 25, 17** in the family
order below. [MATH_STATUS.json](../../MATH_STATUS.json) supplies claim authority;
the [inventory](../INVENTORY.md) supplies later results on the same curves.
No exact rank, public novelty, conductor record or general success rate follows.

| Family | Native parameter | Initial seed bound | Seed calls | Final dated bound |
|---|---|---:|---:|---:|
| 074d9 | -2707/3437 | 18 | 47 | 24 |
| 07ca9 | 587/1074 | 18 | 1 | 23 |
| 08234 | 2071/4060 | 17 | 86 | 17 |
| 08f72 | 4032/1663 | 18 | 35 | 24 |
| 103b2 | 3889/2172 | 18 | 2 | 25 |
| 11952 | -436/3149 | 17 | 98 | 17 |

## Selection and bounded exposure

The [selector](../cas/select_fresh6_retained_r17_v2.py) uses 6,144 saved H4096
score rows, ordered by selection-band S1, good-prime count, denominator and
signed numerator. Validation-band values, measured ranks and points do not
enter ordering. Frozen exclusions comprise 1,764 equation tuples and 1,146
scheduled addresses; exact rational isomorphism checks leave 5,828 eligible
addresses. These are snapshot-relative address counts, not independent curves
or public novelty. The six selected equations are pairwise nonisomorphic.
The failed v1 selector stopped on a singular exclusion entry; v2 skips it.

[Generic preparation](../cas/prepare_fresh6_generic_seeds.py) specializes the
17 atlas sections and certifies their independence. The first-extra-point
constructor tests exact generic maximum classes, specialized CVP representatives
and bounded pointed boxes, retaining every winning witness. Its 269 point calls
had no timeouts. Limits were 1,800 seconds/2 GiB per search, 5 seconds/1 GiB
per map, and 10 seconds per point call at height 125,000. Bounded misses give
no rank upper bound. The first replay's tuple/list comparison failure was
fixed by JSON normalization; its search and failed receipt remain preserved.

- [Generic M17 packets](../../artifacts/generated-results/elliptic-curves/fresh6_generic17_packets_v1/result.json) · [M17/M18 seed packets](../../artifacts/generated-results/elliptic-curves/fresh6_first_m18_v1/result.json).
- [Original four-seed queue](../../artifacts/generated-results/elliptic-curves/fresh6_first_m18_v1/amplification-queue.json) · [reservation snapshot](../../artifacts/generated-results/elliptic-curves/fresh6_reserved_seed_results_v1.json).

## Certified amplification endpoints

Winning-parent V3 completed all four branches: 074d9, 08f72 and 103b2 stayed
at M18 after six calls each; 07ca9 reached M19 in nine calls. Complementary
maximum-parent banks then completed all four branches. Each certified gain
triggered rebuilding. Their historical full replay receipts cover rational
CVP, maps, points, gain provenance and complete-cloud audits.

| Endpoint | Calls in this pass | Primary evidence |
|---|---:|---|
| 07ca9 M18 → M19 | 9 | [Winning-parent result](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-fresh-001/result.json) |
| 07ca9 M19 → M22 → M23 | 7; reconciliation uses 0 new searches | [Full run](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-exact-complement/result.json) · [M23 packet](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-exact-complement/reconciled.json) |
| 074d9 M18 → M24 | 100 | [M24 result](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/074d9-exact-complement/result.json) |
| 08f72 M18 → M24 | 100 | [M24 result](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/08f72-exact-complement/result.json) |
| 103b2 M18 → M25 | 100 | [M25 result](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement/result.json) |

07ca9's full-cloud audit found a direction beyond the frozen M22 terminal.
[Saved-witness reconciliation](../cas/reconcile_verified_v3_cloud.py) certified
M23 without another point search. Retain all returned points when stopping
on a first gain; reconcile the cloud before preparing another basis.

## Completed misses and retained cursors

These are dated branch checkpoints, not instructions to start or resume work.
Earlier queues point to their descendants; check those links before any new
budget. None of the following misses establishes an upper bound.

| Branch | Completed exposure | Last retained checkpoint |
|---|---|---|
| 07ca9 M23 complementary bank | 100 further calls, no gain | [Result](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-M23/result.json) · [suffix](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-M23/suffix-queue.json) |
| 074d9 M24 cached branch | 200 cumulative calls including the gaining pass | [Cached result](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/074d9-exact-complement-cached-v1/result.json) · [centre 58, preconditioned_full; 621 of 679 centres remain](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/074d9-exact-complement-cached-v1/suffix-queue.json) |
| 103b2 M25 cached branch | 300 cumulative calls including the gaining pass | [Second cached result](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement-cached-v2/result.json) · [centre 126, preconditioned_full; 576 of 702 remain](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement-cached-v2/suffix-queue.json) |
| 103b2 disjoint maximum parents | 100 calls, no gain | [Separate branch](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-disjoint-max/result.json) |
| 103b2 sampled norm10 parents | 100 calls, no gain | [Separate branch](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-sampled-shell/result.json) |
| 08234 sampled parents outside the maximum span | 32 boxes, still M17 | [Separate seed attempt](../../artifacts/generated-results/elliptic-curves/08234_outside_maximum_seed_v1/result.json) |

The [exact parity audit](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2_sampled_parent_span.json)
shows that all 43 maximum classes of 103b2 span dimension 16, with annihilator
mask 45903. Thirteen of sixteen sampled norm10 parents lie outside that
hyperplane; together with the productive parent the bank spans dimension 17.
Different parent masks need not enlarge a previous bank's span. Full generic
parity coverage predicts neither specialized independent gains nor point
existence. The [sampled-parent preparer](../cas/prepare_sampled_continuation_parents.py)
uses a fixed 128-parity sample and both exact CVP solvers; it is not a full
lower-shell census.

## Conductor evidence

The [four bounded audits](../../artifacts/generated-results/elliptic-curves/fresh6_conductor_audit_v1/result.json)
retain exact local data and a second completed execution for 08f72 at 4032/1663.
Its conductor is

```
5010188826156411226384182822992530591002057751473896628642100179406951935900788324249280619427023331755172627345571322858912720052190
```

It exceeds the pinned rank-at-least-24 benchmark for curve548, whose conductor
is 361660950250291954113011326751617907627414070443777168347285410734060239220929310535800.
This comparison concerns that saved catalogue only. The 074d9, 07ca9 and 103b2
audits hit 30-second/1 GiB caps with residual cofactors of 124, 104 and 113
digits. Their exact conductors remain UNKNOWN in these audits.
[The audit script](../cas/audit_fresh6_conductors.sage) performs local reduction
and proof-enabled factorization; it has no read-only replay mode.

## Retained rank replay

From the repository root, `make verify-fresh6-retained-ranks` checks this cohort,
the [second cohort](SECOND_FRESH6_SEED_COHORT_2026-09-09.md) and the
[lower-height cohort](LOWHEIGHT_FRESH6_SEED_COHORT_2026-09-09.md). The
[pinned 660 kB witness](../../artifacts/generated-results/elliptic-curves/fresh6_retained_rank_witnesses_v1.json)
retains 37 subgroup endpoints for 17 claims. The
[checker](../cas/verify_fresh6_retained_ranks.py) reconstructs native sections,
checks exact rational points, recomputes 626 sufficient small-prime blocks
selected from 5,912 saved blocks, and verifies torsion exclusion and pairwise
nonisomorphism within each seed cohort. All 37 passed in 2.87 seconds on
12 September. This shares the established finite-group implementation;
it is not a second independent algorithm or an exact-rank computation.

`python3 research/elliptic-curves/cas/verify_fresh6_retained_ranks.py --check-sources`
compares the projection with 30 hashed original exports, without arithmetic.
The default arithmetic replay needs only the retained witness, pinned atlas
and Python code, with no original worker, persistent cache, CVP or point search.
It does not replay selection, maps, complete clouds, timing or conductor data.

The original [lean](../cas/run_lean_preconditioned_seed_v3.py) and
[complementary](../cas/run_complement_seed_v3.py) drivers have explicit `replay`
modes. After version guards, an existing matching `verified.json` makes them
return without a fresh full replay. Without that receipt, replay may reconstruct
reference CVP landscapes. Preserve these checkpoints; removing them to force
replay is not part of cleanup. Saved-cloud reconciliation separately needs raw
chart/epoch files. All 1,295 inspected bindings matched, and those referenced
files were present locally in the [provenance audit](../../archive/repository-cleanup-2026-09-12/fresh6-source-review/PROVENANCE.json).
They are not all distributed with a fresh checkout. Historical assurance flags
were retained; no full replay or conductor factorization was run by cleanup.

The [review receipt](../../archive/repository-cleanup-2026-09-12/fresh6-source-review/REVIEW.json)
preserves the complete earlier narrative, timings and superseded instructions.
