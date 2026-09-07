# Certified discovery yield across retained MW16 score levels

**Complete. All 2,580 boxes, sixty exact baseline and terminal point proofs,
independent accounting checks and 240 standalone replay stages pass. The
extreme-top arm yields 10 certified added directions, the moderate arm 1,
and the lower arm 0. Top-only yield is 2.7203 times the equal three-arm
portfolio's yield per discovery worker-second in this matched sample.
The queued R17 sweep remains cancelled; no new parameter scan follows.**

## Completed outcome

Each arm contains twenty curves and completes all 860 allocated boxes.
There are no chart timeouts, worker failures, missing CPU markers or
unresolved point certificates. Discovery time includes map construction and
point-worker time; it is summed worker time, not elapsed wall time with two
workers. The initial population scan and shared orchestration costs are not
charged selectively to an arm.

| Arm | Certified added directions | Discovery worker-seconds | Directions / 1,000 discovery seconds | Seconds including verification |
|---|---:|---:|---:|---:|
| Extreme top | 10 | 875.808 | 11.418 | 1,229.967 |
| Moderately strong | 1 | 867.214 | 1.153 | 1,215.689 |
| Further down | 0 | 877.662 | 0 | 1,227.883 |

The equal three-arm portfolio yields 11 directions in 2,620.684 discovery
worker-seconds: 4.197 per 1,000 seconds. Top-only yield is **2.7203 times**
that rate, and **2.7152 times** the portfolio rate when baseline and terminal
verification costs are included. Top remains ahead in each of the five
leave-one-family-out comparisons. The frozen diversification criterion is
not met; the recorded decision is `TOP_FAVOURED_IN_THIS_FINITE_SAMPLE`.
Keep priority on the strongest retained initial-score candidates rather than
spreading this point budget equally across these score levels. This finding
does not authorize another larger parameter sweep.

A subsequent [read-only CPU-cost sensitivity check](../../artifacts/generated-results/elliptic-curves/strata60_gp_cpu_yield_v1.json)
uses the already independently audited GP timing markers:449,847 milliseconds
for top,443,004 for moderate and453,708 for lower. Its exact top-to-portfolio
yield ratio is4488530/1649439, approximately2.72125. This agrees with the
discovery-worker-time conclusion. GP CPU measures enumeration only and
excludes Python admission, map preparation and verification; it is a secondary
cost measure, not a replacement criterion or total-computation estimate.
`audit_strata60_gp_cpu_yield.py --check` replays the calculation.

Only three curves gained certified directions. One top-arm curve supplies
nine of its ten directions, so this remains a small, uneven finite sample.
The overlapping policy estimates are not independent samples. The result
neither proves optimality nor predicts record-tail discovery rates. A zero
certified gain does not prove that a curve has no additional rational points.

The [experiment report](../../artifacts/generated-results/elliptic-curves/retained_mw16_score_strata_experiment_v1.json)
and [independent accounting replay](../../artifacts/generated-results/elliptic-curves/retained_mw16_score_strata_accounting_replay_v1.json)
retain every curve, family, matched triplet, cost and exposure endpoint.
The [standalone replay](../../artifacts/generated-results/elliptic-curves/strata60_mw16_point_portable_replay_v1.json)
passes all 240 stages from a 93,533,321-byte archive.

### New rank-at-least-25 curve

The top arm finds family02 at **t=-32999/14074**, now stable inventory ID
`new-20260906-200`, with 25 independent rational points. Its proved global
minimal Weierstrass equation has coefficients

```
[1, 0, 0,
 -1600850384897451561881895982234065750603549638504054515,
 779215362299612638813009572863885921153122866957851908328254774486756288990448225]
```

The invariant gcd is 361; exact local scaling exclusions prove minimality.
The [minimal-model certificate](../../artifacts/generated-results/elliptic-curves/strata60_high_rank_minimal_proof_v1.json)
and [executable Sage equation and point export](../../artifacts/generated-results/elliptic-curves/new_strata60_high_rank_curves.sage)
pass. All sixty comparison equations are mutually nonisomorphic over Q and
unmatched among 593 pinned catalogue and 1,345 prior equations, as checked
in the [novelty certificate](../../artifacts/generated-results/elliptic-curves/strata60_mw16_novelty_v1.json).
The other gains are one direction each at family01 t=32719/6374 (top) and
family02 t=-134337/196864 (moderate). Final lower bounds are fifty-seven16s,
two17s and one25; they are not exact ranks.

The [V20 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v20.json)
contains 200 distinct curves: eight with bound27, eighteen with26,
thirty-eight with25, fifty with24, forty-six with23 and forty with22.
All point certificates, source bindings, catalogue exclusions, distinct
j-invariants and the CSV pass the
[200-curve replay](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v20_memory_replay_v1.json).
No new rank28/32, exact conductor or universal novelty is asserted.

### Where the successful curve was already retained

A [post-outcome lookup and replay](../../artifacts/generated-results/elliptic-curves/strata60_retained_score_positions_v1.json)
accounts for all sixty matched curves in the existing earlier score tables.
Only the twenty top-arm curves already had extended scores; no new trace
or validation-prime result was calculated for this lookup.

The new25 curve was rank21 in its signed initial-score block, rank70 in the
combined-sign family/band initial ranking, and **rank9 of1024 by the earlier
extended score**, below the original six-finalist cutoff. Its extended score
is154.868682625006 versus155.172108737511 at that cutoff. Thus this discovery
was already available inside the earlier scalar-scored population. This
locates a shortlist limit, not a need for another billion-address scan; it
does not establish an optimal new cutoff or retroactively change either
experiment.

## Frozen question and design

The production question is whether concentrating point-search work at the
extreme score tail yields more certified independent directions per unit
of computation than a portfolio spread across score levels. Better selection
scores alone cannot answer it. The original corrected campaign is a separate
experiment and will not be pooled into the primary comparison: its existing
rank-triggered stop would make realized exposure different.

## Population and prospective allocation

Use the 1,310,720 candidates retained by the completed corrected MW16 scan.
Every one has the same corrected score over the 3,510 selection primes
through 32749, including restored good reductions at 5 and 13. Define strata
using that common score. The later 65521 scores available only for the earlier
top 10,240 candidates do not define or alter this comparison. Consequently
the policy comparison concerns allocation across the common initial-score
strata of this retained pool. It does not directly compare alternatives to
the later 65521-prime finalist-ranking pipeline, or retained candidates with
unretained parameter addresses.

There are twenty matching blocks: five families, two existing height bands
(16384 < H <= 65536 and 65536 < H <= 262144), and two parameter signs.
Each block has 65,536 retained candidates. Rank by decreasing corrected score,
decreasing good-prime count, denominator and signed numerator. In each block,
freeze one matched triplet:

| Arm | One-based rank within the retained block | Fixed candidate pool |
|---|---:|---:|
| Extreme top | 1–128 | All 128 |
| Moderately strong | 1025–8192 | SHA256 sample of 1024 |
| Further down | 16385–65536 | SHA256 sample of 2048 |

This gives twenty curves per arm and sixty overall. A fixed seed orders
candidate pools without point outcomes or validation primes. Exclude rational
isomorphs of the original sixty corrected finalists using their frozen
prospective models, and deduplicate the new sixty. No catalogue or public
rank labels enter selection.

Match within the exact block using the multiplicative Weil height
H(j)=max(|numerator(j)|, denominator(j)) of the reduced rational j-invariant.
The largest and smallest H(j) in a triplet must differ by at most a factor
four; parameter heights must differ by at most a factor two. This is an exact,
model-independent arithmetic-height criterion, not an exact conductor or a
minimal-coefficient-height claim. Also record the actual search-model
coefficient sizes as computational covariates. Select the first feasible
top anchor in seeded order, then a feasible moderate/lower pair in a fixed
deterministic order. Failure to fill a block is reported as a matching failure;
do not widen strata, relax the calipers or search an extra population.

The selection and its independent replay must finish before point work.
The selector reads no point outcomes. It may run only after the corrected
campaign's terminal point, geometry and standalone proof stages finish.

## Identical point exposure and computation accounting

All sixty curves start from exactly the sixteen generic sections, whose
specialized independence is checked. Every curve receives the same 43 generic
parity labels, point height 125000 and ten-second cap per chart: 2,580 possible
boxes. All map files precede all point searches. Use two workers, the same
600-second per-curve cap, and a seeded, counterbalanced triplet order.
There is **no rank-based early stop**, adaptive wave, refill, extra chart or
arm-dependent retry. The original corrected campaign's stop rule is unchanged.

For every curve retain attempted and completed boxes, chart timeouts,
worker failures, elapsed worker seconds, map seconds, exact-verification
seconds and the available GP search CPU milliseconds. Missing CPU time on
censored calls is unavailable, never zero. Completed boxes are counted from
their successful bounded-search transcripts, not inferred from chart count.
Report all allocated curves, including failures. Verify any retained partial
point evidence separately; an incomplete computation is not a zero-rank-gain
or point-absence theorem.

## Endpoints and production decision

The primary rank endpoint is the certified lower-bound gain over the same
specialized sixteen-point subgroup. Verify exact points, transports and
independence; online estimates and Nagao scores are not results. Report:

- summed certified additional directions and their distribution per arm;
- certified directions per discovery worker-second, including map time;
- the same rate including exact verification cost;
- completed boxes, completion fraction, timeouts and censored attempts;
- yield per completed box as a secondary exposure-normalized description;
- matched-triplet and family breakdowns, including height balance.

Compare the top-only policy with an equal three-arm portfolio using the same
observed curves and costs; do not treat the overlapping policy estimates as
independent samples. Also compare top directly with moderate and lower arms.
The initial population and trace-table costs are sunk shared costs for this
retained-pool decision and are reported separately, not charged selectively.

A production recommendation to diversify requires higher pooled certified
gain per discovery time and per completed box, no lower completion fraction,
and persistence of the gain-per-time advantage when each family is omitted
in turn. These are declared finite budget-decision criteria, not a significance
test or a theorem of rank density. Otherwise report an inconclusive or
top-favouring result as the measurements warrant; do not automatically launch
a larger sweep. Low power, censoring or inability to match are valid outcomes.

Primes 65537 through 131071 remain separate validation data. If computed,
they are applied equally to the frozen arms and reported separately; they
never change strata, matching, exposure, the primary endpoint or a failed
decision criterion. No new validation calculation is needed to select this
experiment. Catalogue comparison and any inventory promotion follow terminal
exact proofs and remain separate from the production-policy endpoint.

The stopped R17 controller protocols and their waiting ledgers remain
unaltered. The cancellation journal is
`artifacts/local/elliptic-curves/retained-score-stratification-v1/cancelled-queued-r17-sweep.json`.
The already-started R17 trace-table and byte-verification work has finished;
its downstream parameter, scalar and point controllers remain stopped.

## Execution and replay status

The original corrected campaign and all 182 standalone checks finished
unchanged. The matching selector and exact replay now pass for all twenty
triplets, with no missing block, relaxed caliper or enlarged pool.
No corrected rank outcome entered this comparison's allocation.

The separate [point controller](../cas/finish_strata60_mw16_experiment.py)
has passed its freeze, all sixty maps, all sixty exact generic16 baseline
certificates, all point attempts, seven bounded proof stages per result and the
[accounting report and replay](../cas/report_strata60_mw16_experiment.py).
A matching failure prevents point work. Map or baseline gate failures also
prevent point work and trigger an all-allocation unresolved-outcome report.
Search or proof failures do not remove individual curves from accounting.
No retry or automatic follow-on campaign is scheduled.

The [point worker](../cas/strata60_mw16_pari_batch.py) checkpoints each chart
start and returned transcript separately from point admission. The
[verifier](../cas/verify_strata60_mw16_points.py) can certify a terminal partial
prefix, including exact points returned before interrupted admission.
It replays admission history, generic transports, raw maps and exposure,
then independently checks the retained point cloud modulo 2, 3 and 5.
Missing certificates or unverified completion counts remain null in the
report; they prevent a resolved policy recommendation.

Eight focused regression checks pass for exact matching, missing-result
accounting, partial exposure, computation-normalized yield, completion and
leave-one-family-out robustness. These checks validate machinery, not a
scientific outcome. The sixty curves, maps, exact baseline certificates, point attempts and
terminal proofs all pass; the completed outcome is recorded above. Protocols and ledgers are under
`artifacts/local/elliptic-curves/retained-score-stratification-v1/`;
point artifacts are under `strata60-mw16-pari-v1`.

A separate [accounting audit](../cas/audit_strata60_mw16_accounting.py) has passed
both its terminal build and read-only replay. It reconstructs all sixty curve
costs, chart completion, missing CPU timing, gain bindings, arm/family rates
and the policy criterion directly from saved evidence without importing the
report generator. Three additional regression checks reject zero-filling
unresolved outcomes, loss of timeout cost and corrupted counts or nonfinite
rates. All eleven matching/accounting checks pass. Its local controller is
`retained-score-stratification-v1/accounting-controller`.

The [frozen selected roster](../../artifacts/generated-results/elliptic-curves/retained_mw16_score_strata_selection_v1.json)
contains twenty curves per arm. Realized retained-block rank ranges are
1–104 (top), 1,085–7,763 (moderate), and 18,939–64,639 (lower).
The largest within-triplet j-height ratio is about 2.74659, with median
1.02291; the largest parameter-height ratio is about 1.79313. The exact
rational ratios are retained in the selection certificate. These are balance
diagnostics, not discovery outcomes or reasons to alter the allocation.

The [standalone point bundle builder](../cas/package_strata60_mw16_points.py)
and [isolated verifier](../cas/verify_strata60_mw16_points_portable.py) have
finished after the terminal comparison and independent accounting audit.
All sixty allocations remain in the bundle. Each certified row receives
four isolated checks: admission history, rational geometry/provenance,
mod-2 point independence and mod-3/5 point independence, at most 240 stages.
Unresolved rows are retained without a zero-gain or rank claim. The local
accounting audit and computation logs are embedded as immutable evidence;
the isolated point replay does not rerun search timings or population scores.

The [relocation preflight](../../artifacts/generated-results/elliptic-curves/strata60_mw16_portability_preflight_v1.json)
passes for 116 prospective files under isolated Python in a separate root.
It checks the sixty frozen input bindings and proof imports without copying
or reading any point outcomes. This closes an input-portability gap before
the terminal proof bundle is built; it is not an already-completed point
proof replay. The portable controller is under
`retained-score-stratification-v1/portable-controller`.

A separate [terminal novelty check](../cas/certify_strata60_mw16_novelty.py)
has passed build and replay after the independent accounting result. It compares all sixty
selected equations up to rational isomorphism with the 593 pinned catalogue
equations and 1,345 previously measured equations, including the completed
corrected cohort. Certified rank bounds are attached only through the
separate verified point clouds; unresolved bounds stay null. Catalogue
absence cannot change the score-stratum allocation or production criterion.
The local controller is `retained-score-stratification-v1/novelty-controller`.

The declared decision criterion measures certified-direction yield. It does
not establish which policy most often reaches rank 28 or the record target;
that broader objective remains open even if this finite comparison supports
a change in retained-pool allocation.
