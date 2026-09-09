# Second fresh six-fibre seed cohort

Five of six newly selected native R17 fibres now have independently certified
M18 seeds. The remaining fibre retains lower bound17 after a bounded miss.
No exact rank, public novelty or record is claimed.

Selection uses the same6144 saved score rows and unchanged per-family order,
with a new frozen exclusion snapshot containing the first cohort reservation.
It selects one additional fibre per family, excluding rationally isomorphic
curves and scheduled addresses in the named snapshots.5822 address rows
remain eligible at selection; this is not a predicted rank density.
Selection and replay take2.096 and2.043seconds. Generic17 preparation and
replay take19.666 and14.261seconds, without point search.

| Family | Parameter | Certified lower bound | Seed point calls |
| --- | --- | ---: | ---: |
| 074d9 | 88/2551 | 18 | 9 |
| 07ca9 | -2475/2848 | 18 | 3 |
| 08234 | 2570/2143 | 18 | 5 |
| 08f72 | -2433/3479 | 17 | 98 |
| 103b2 | 2815/1088 | 18 | 7 |
| 11952 | 327/1403 | 18 | 11 |

The unchanged first-M18 constructor uses exact generic maximum classes,
specialized exact CVP representatives, lean bounded maps and cached finite
admission. It stops at the first certified extra direction. All six search
and geometry/map/point replays pass. Total search time is117.234seconds
for133 point calls, with no point timeouts. Each search is bounded by
1800seconds/2GiB; each map5seconds/1GiB; each point call10seconds at height125000.

A separate pure-Python packet checker reconstructs the generic sections,
rechecks finite independence and pairwise rational nonisomorphism. It passes
in7.874seconds. The five M18 seeds are queued; amplification has not started.
Full winning witnesses are retained, including unprocessed points after the
first gain. No general seed success rate is inferred from this small cohort.

- [Selector](../cas/select_fresh6_retained_r17_second_cohort.py).
- [Generic17 packets](../../artifacts/generated-results/elliptic-curves/fresh6_second_generic17_v1/result.json).
- [Exact seed packets](../../artifacts/generated-results/elliptic-curves/fresh6_second_first_m18_v1/result.json).
- [Independent packet checker](../cas/verify_fresh6_seed_packets.py).
- [Five-seed queue](../../artifacts/generated-results/elliptic-curves/fresh6_second_first_m18_v1/amplification-queue.json).
- [Reservation ledger](../../artifacts/generated-results/elliptic-curves/fresh6_second_reserved_results_v1.json).

Raw evidence: `artifacts/local/elliptic-curves/fresh6-second-selection-v1/`,
`fresh6-second-generic-seeds-v1/` and `fresh6-second-seed-confirmation-v1/`.
The [first cohort](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md) and its unfinished
suffixes remain preserved.

## First productive-parent amplification

All five M18 seeds complete their winning-parent V3 policies and independent
replays. The103b2 fibre at2815/1088 gains two certified directions:18→19 on
its first chart, then19→20 after rebuilding and three more charts. Fourteen
further charts complete its finite policy without a gain. Search takes35.845
seconds; independent replay takes9.638seconds. Its subgroup lower bound is20.
The other four each complete six charts and remain M18. No point or map
timeouts occur. These finite-policy misses establish no rank upper bounds.

| Family | Final certified lower bound | Point calls |
| --- | ---: | ---: |
| 074d9 | 18 | 6 |
| 07ca9 | 18 | 6 |
| 08234 | 18 | 6 |
| 103b2 | 20 | 18 |
| 11952 | 18 | 6 |

Each resulting basis now has a prepared16-parent complementary maximum-class
bank, with all remaining generic minima checked by both exact solvers. The
M20 basis is carried forward for103b2. These complementary searches have not
started. No exact rank, public novelty or record is claimed.

- [M20 run evidence](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-fresh-002/result.json).
- [Five prepared complementary-parent branches](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/complement-queue.json).

Next start the broader-parent pass from the new M20 state, then cover the
remaining four branches with the same bounded protocol.

## 103b2 complementary-parent amplification and M23 reconciliation

At2815/1088, complementary-parent V3 raises M20→M21 after9 calls and
M21→M22 on the next call. The complete cloud audit certifies23 directions,
triggering `ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION`. Search takes
31.348seconds; full independent geometry/map/point/cloud replay passes in
25.974seconds. The frozen M22 terminal and its stop remain preserved.

Arithmetic-only reconciliation of the saved charts extends the retained basis
to **M23**. A standalone finite certificate proves independence. Construction
and a second full reconciliation replay each take4.754seconds, with no
additional point search. A fresh sealed M23 preparation retains the same
complementary parents and is ready for rebuilding. No exact rank, public
novelty or record is claimed. The four other complementary branches remain
prepared but unsearched.

- [Full bounded-run evidence](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-complement/result.json).
- [Reconciled M23 packet](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-complement/reconciled.json).
- [Current continuation queue](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/complement-queue.json).

## Rebuilt103b2 continuation: M23 to M26

At2815/1088, the fresh continuation from the reconciled M23 basis finds
three further certified directions. Epochs23→24,24→25 and25→26 use31,
23 and43 calls respectively. The third gain is on cumulative call97;
after rebuilding, three more calls complete the100-call budget at M26.
No point or map timeouts occur. Search takes235.402seconds; full independent
rational-CVP, exact map/point, gain-provenance and complete-cloud replay
passes in401.103seconds. The longer replay includes four successive bases.
The earlier progress snapshot at M25 preceded the final gain; the sealed
terminal and certificate establish the lower bound26.

The terminal is `CHART_BUDGET_EXHAUSTED`. Its exact M26 continuation cursor
is preserved, with no exact rank, upper bound, public novelty or record claim.
The four other complementary-parent branches remain prepared and unsearched.

- [Certified M26 continuation](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-M23/result.json).
- [Exact M26 suffix cursor](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-M23/suffix-queue.json).

### Cached M26 continuation

The first cached continuation adds100 point calls,200 cumulative from the
reconciled M23 restart, without a new certified direction or point/map timeout.
Search takes94.428seconds; independent cached replay takes11.014seconds.
The lower bound remains26. The exact cursor is centre51, `factor_free`,
so its preconditioned box is already tested and its factor-free box is next.
The prior queue points to this verified descendant. No rank upper bound follows.

- [Cached M26 result](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-M26-cached-v1/result.json).
- [Current M26 cursor](../../artifacts/generated-results/elliptic-curves/fresh6_second_amplification_v1/103b2-M26-cached-v1/suffix-queue.json).
