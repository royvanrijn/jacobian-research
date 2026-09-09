# Test parents outside the previously productive span

The [retained-label audit](../cas/audit_curve302_parent_novelty.py) checks two
completed M18→M31 cascades on curve302. In `recovered-strict-02`, all13 winning
generic masks are distinct and independent over F2. In `recovered-strict-03`,
there are11 distinct winning masks, also independent; two later gains reuse
earlier masks. In the first run, all four transitions from M27 onward use
classes outside the span of its earlier winning masks. In the second, two of
those four transitions do.

These are retrospective chart labels from two seeds on one curve. This is
not a new rank replay, a cross-fibre predictor, or proof that an omitted
winning mask is the only way to find its point. It does show that a policy
restricted to earlier productive masks, even their entire F2 span, excludes
many classes used by the known successful mechanism. Pairwise XORs cannot
escape that span.

Evidence: `artifacts/local/elliptic-curves/curve302-parent-novelty-audit-v2/`.
The v1 exploratory report is retained; v2 adds the reproducible source and
explicit terminal-hash verification. The independent twelve-seed302 panel
was inspected read-only and remains under its existing controller.

## Prospective construction on curve48

Curve48's eight retained productive masks span dimension8 inside the native
17-dimensional generic parity space. Their exact generic minimum is12
(scaled norm24 in the doubled Gram matrix). No curve302 points or winning
mask labels are transferred to curve48.

The [bounded preparation](../cas/prepare_sampled_complement_parent_bank_v2.py)
hashes `family/outside-productive-span/index` with SHA256 to generate128
distinct parity masks outside that original span. Each generic minimum is
checked by both the integer and independent rational CVP solvers. The observed
scaled minima are8 (one mask),12 (37),16 (51), and20 (39). All128 checks pass.

Among the39 masks with largest sampled minimum20, ascending mask order first
extends the old span and then fills the bank to16 classes outside the original
span. The old and new masks together span dimension17. This does not mean the
16 selected classes cover all generic parities or prove new rational points.
The resulting norm10 bank is a finite prospective choice, not a complete shell
bank or a demonstrated optimal choice.

Preparation completes in2.986 seconds under120 seconds/1GiB, with zero point
searches. Its data and all exact128-candidate/16-selected checks are retained
under `artifacts/local/elliptic-curves/curve48-complement-preparation-v3/`.

Two unsuccessful preparations remain as evidence. Complete enumeration
through scaled norm24 hits its20-million-node cap in47.470 seconds without a
complete bank. The first128-sample variant requires a norm24 candidate, finds
none and stops in2.465 seconds before producing a bank. The adopted version
explicitly selects the largest sampled norm instead; no norm24 completeness
or impossibility claim is made.

## Exposure and bounded discovery protocol

The [complement selector](../cas/visibility_complement_subset.py) retains V3
extension scoring and exact CVPs on this explicitly frozen16-parent bank.
The [runner](../cas/run_complement_seed_v3.py) uses bounded lean maps, exact
finite-group admission and immediate rebuild after a certified gain. Its
100-invocation budget is shared across adaptive epochs; target lower bound32.
Map limits are5 seconds/1GiB, point limits10 seconds at height125000, with
1800 seconds/3GiB per search and independent replay phase.

It also uses exact signed-permutation box deduplication, including infinity,
as tested in [the coordinate preflight](CURVE92_PRODUCTIVE_V3_SEARCH_2026-09-09.md).
Existing protocols and their deduplication rules remain unchanged.

The [new preflight](../cas/prepare_complement_exposure.py) builds779 centres
in72.800 seconds under300 seconds/1GiB. None overlaps up to sign with the248
distinct centres in the three named histories: the latest300-box productive
cloud and the two earlier49-centre attempts. This is exposure only against
those histories, not a point-existence or rank claim.

## Completed first pass

The first100 actual point invocations complete at50 centres, with zero map
or point timeouts: M27→M27. Search takes87.931 seconds, peak process-tree
RSS159051776 bytes. Full independent replay passes in225.810 seconds,
peak304615424 bytes, including the original rational-CVP landscape, exact
maps, point witnesses, receipt chain and whole-cloud mod2/3/5 certificates.
All three lower bounds remain27. Preparation and preflight are separate costs.

- [Sealed result](../../artifacts/generated-results/elliptic-curves/curve48_complement_v3_v1/result.json).
- Raw evidence: `artifacts/local/elliptic-curves/curve48-complement-v3-discovery-v1/`.
- The [cursor recorder](../cas/queue_verified_complement_pass.py) retains
  centre50, `preconditioned_full`, with729 centres unvisited.
- A [cached continuation adapter](../cas/run_complement_cached_seed_v3.py)
  supports this protocol; its first completed continuation is below.

The budget is exhausted, not the selected policy. The outside-span branch has
now received a prospective bounded test, but no rank28 direction, conductor
improvement, rank upper bound or general amplification prediction follows.

### Curve48 cached continuation

The next100 actual point invocations complete at centres50–99,200 cumulative
in the outside-span branch. M27→M27, zero map or point timeouts. Search takes
90.711 seconds, peak RSS159178752 bytes; replay passes in2.513 seconds,
peak76042240 bytes. It inherits the bound independently verified full CVP
landscape and rechecks the maps, all200 point receipts and the extended cloud.
The mod2, mod3 and mod5 lower bounds remain27.

[Sealed continuation](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v1/result.json).
Raw evidence is in `artifacts/local/elliptic-curves/curve48-complement-cached-v3-discovery-v1/`.
The cursor retains centre100, `preconditioned_full`, with679 centres unvisited;
the first pass's `continuation-head.json` points to this verified descendant.
The separate productive-parent branch has300 completed boxes. No rank28
direction, conductor improvement or upper bound is obtained by this continuation.

## Curve92 transfer

The same bounded preparation and runner now apply to the native MW16 M26
seed on curve92. Its15 selected norm9.5 parents extend the old productive
span from dimension3 to16. All708 centres are new relative to the302 distinct
centres in the two compared completed policies. The first100 boxes replay
without gain or timeout, M26→M26. Search92.466 seconds and full replay247.369
seconds pass. Its first cached continuation completes100 additional boxes,
200 cumulative, with no gain or timeout. Search100.828 seconds and replay
3.052 seconds pass; the earlier independent full CVP proof is inherited.
The second continuation adds100 boxes,300 cumulative, again without gain or
timeout. Search102.124 seconds and replay4.076 seconds pass, leaving558 centres.
See [the curve92 record](CURVE92_PRODUCTIVE_V3_SEARCH_2026-09-09.md) for exact
bindings, exposure limits and its still-conditional conductor payoff.

## Curve48 second complementary-parent cached continuation

After the fresh-seed portfolio passes, the existing native M27 curve48
complementary-parent branch receives100 additional point calls,300 cumulative.
All complete without point/map timeout and without a new certified direction.
Search takes90.002seconds; independent cached replay passes in3.341seconds.
The lower bound remains27. Its exact next cursor is centre150,
`preconditioned_full`, leaving629 of779 centres unvisited. These are bounded
misses, not a rank upper bound. No old search or source is modified.

- [Second cached result](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v2/result.json).
- [Current complementary-parent cursor](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v2/suffix-queue.json).

## Curve71 sampled complementary-parent pass

The native103b2 fibre at3726/881 starts from its verified M27 basis.
A fixed128-parity sample outside its productive-parent span is checked with
both exact CVP solvers. Sixteen parents of doubled norm20 are selected;
preparation completes in2.983seconds without a point search.

The bounded V3 pass completes100 point calls without a new certified
direction or point/map timeout. Search takes152.403seconds; full independent
rational-CVP, map, point and complete-cloud replay passes in206.492seconds.
The lower bound remains27. The terminal is `CHART_BUDGET_EXHAUSTED`, with
its exact unvisited suffix queued. No rank upper bound or record is established.

- [Replayed curve71 complementary-parent pass](../../artifacts/generated-results/elliptic-curves/curve71_sampled_shell_v1/result.json).
- [Complementary-parent continuation cursor](../../artifacts/generated-results/elliptic-curves/curve71_sampled_shell_v1/suffix-queue.json).

### Curve71 cached complementary-parent continuation

The first cached pass completes100 further calls,200 cumulative, without
a new certified direction or point/map timeout. Search takes81.327seconds;
independent cached replay takes2.409seconds. The lower bound remains27.
The next cursor is centre100, `preconditioned_full`, leaving651 of751
centres unvisited. The first-pass queue now points to this verified descendant.
These bounded misses establish no rank upper bound.

- [Cached result](../../artifacts/generated-results/elliptic-curves/curve71_sampled_shell_cached_v1/result.json).
- [Current cursor](../../artifacts/generated-results/elliptic-curves/curve71_sampled_shell_cached_v1/suffix-queue.json).

## Curve48 deeper three-batch complementary-parent continuation

Three sequential cached batches each add100 point calls, with independent
replay between batches and a stop-for-reassessment condition on any gain.
All three complete without gains or point/map timeouts. The lower bound
remains27; total complementary-parent exposure is now600 calls.

| Cached version | Cumulative calls | Search seconds | Replay seconds |
| --- | ---: | ---: | ---: |
| 3 | 400 | 92.110 | 4.230 |
| 4 | 500 | 95.549 | 5.051 |
| 5 | 600 | 90.079 | 5.931 |

The final cursor is centre300, `preconditioned_full`, leaving479 of779
centres unvisited. Every earlier queue points to its verified descendant.
These bounded misses establish no rank upper bound or record.

- [Current result](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v5/result.json).
- [Current cursor](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v5/suffix-queue.json).

## Exact maximum parents now prepared for curve48

The complete R17 maximum-class catalogue now allows direct preparation of
norm12 parents outside curve48's original eight-dimensional productive
span. This avoids repeating the earlier capped enumeration. The new bank
contains sixteen classes of doubled norm24; every generic CVP is checked
with both integer and rational exact solvers. It starts from the independently
replayed M27 terminal of `curve48-complement-cached-v3-discovery-v5`.

Preparation completes in4.719 seconds under120 seconds/1GiB with zero point
searches. Its [sealed preparation](../../artifacts/generated-results/elliptic-curves/curve48_exact_maximum_preparation_v1/result.json)
has completed its first bounded100-call pass (below). These maximum classes are distinct
from the previously sampled norm10 parents, but a larger generic minimum
is not a guarantee of a new point. The old sampled-bank searches and their
suffix receipts remain intact. No new rank is claimed by preparation.

Curve71's corresponding exact maximum preparation also passes in4.278 seconds,
with sixteen norm12 parents outside its original productive span and both
exact CVP solvers agreeing. It uses the verified M27 terminal of
`curve71-sampled-shell-cached-v1`. Its [sealed preparation](../../artifacts/generated-results/elliptic-curves/curve71_exact_maximum_preparation_v1/result.json)
has completed its first bounded100-call pass after curve48. Preparation makes zero
point searches and does not claim a new rank.

## Curve48 exact maximum first pass

The new norm12 bank completes100 point calls in155.549 seconds and full
independent replay in227.707 seconds: M27 remains M27, with zero map or
point timeouts. This is a bounded miss and supplies no rank upper bound.
The final landscape has779 centres. The sealed suffix resumes at centre
50, policy `preconditioned_full`, leaving729
centres in this search order. Exact geometry, map receipts, rational points,
standalone independence and complete-cloud audits all pass replay.

- [Replayed result](../../artifacts/generated-results/elliptic-curves/curve48_exact_maximum_v1/result.json).
- [Exact suffix cursor](../../artifacts/generated-results/elliptic-curves/curve48_exact_maximum_v1/suffix-queue.json).

Curve71's prepared norm12 bank has now completed its first pass (below).
Curve48's completed prefix and older sampled-bank evidence remain preserved
for later continuation. No record is claimed.

## Curve71 exact maximum first pass and overlap audit

Curve71's norm12 bank completes100 point calls in150.451 seconds and full
independent replay in215.165 seconds: M27 remains M27, with zero map or
point timeouts. Its776-centre landscape resumes at centre50, policy
`preconditioned_full`, leaving726 centres. This is a bounded miss, not a
rank upper bound.

- [Replayed result](../../artifacts/generated-results/elliptic-curves/curve71_exact_maximum_v1/result.json).
- [Exact suffix cursor](../../artifacts/generated-results/elliptic-curves/curve71_exact_maximum_v1/suffix-queue.json).

The new [selected-centre overlap checker](../cas/audit_verified_centre_overlap.py)
requires sealed replayed terminals, their bound initial selections, identical
curve models and initial bases, and verifies every selected centre on the
curve. It compares exact rational points up to elliptic negation. Repeated
runs reproduce both saved receipts; a self-comparison control returns all779
curve48 centres as shared.

| Curve | Prior sampled selection | New norm12 selection | Shared up to sign |
| --- | ---: | ---: | ---: |
| 48 | 779 | 779 | 0 |
| 71 | 751 | 776 | 0 |

The [curve48 overlap receipt](../../artifacts/generated-results/elliptic-curves/curve48_exact_maximum_v1/sampled-centre-overlap.json)
and [curve71 overlap receipt](../../artifacts/generated-results/elliptic-curves/curve71_exact_maximum_v1/sampled-centre-overlap.json)
compare only the two named selected sets per curve. They do not certify
that every corresponding quartic box is inequivalent, compare all historical
searches, or predict independent new point directions. Box deduplication
remains a separate execution gate.

Both first passes now have verified cursors. The next bounded continuation
can reuse the prepared exact geometry and completed coverage, beginning with
curve48's norm12 suffix at centre50, then curve71's. No new rank or conductor
record has been obtained from these two passes.

## Cached norm12 passes and a third M27 family

Both norm12 branches complete100 further point calls, reaching200 cumulative
calls each, with M27 unchanged and no map or point timeouts. Independent
cached replay passes in both cases. Each new suffix starts at centre100,
policy `preconditioned_full`:679 centres remain for curve48 and676 for curve71.
The initial suffix records now point to their verified descendants.

| Curve | Additional search seconds | Cached replay seconds | Cumulative calls |
| --- | ---: | ---: | ---: |
| 48 | 93.032 | 2.531 | 200 |
| 71 | 85.169 | 2.361 | 200 |

- [Curve48 cached result](../../artifacts/generated-results/elliptic-curves/curve48_exact_maximum_cached_v1/result.json).
- [Curve71 cached result](../../artifacts/generated-results/elliptic-curves/curve71_exact_maximum_cached_v1/result.json).

The cheap cached replay reuses sealed geometry from the initial full
independent replay and checks the new receipts. These timings measure
different replay scopes; they are not a claim that full geometry verification
has become90 times faster. Both bounded misses leave rank upper bounds open.

For broader M27 exposure, curve40 at family074d9 parameter2818/1535 now has
sixteen exact norm12 parents outside its original productive span. Both exact
CVP solvers agree. Preparation uses the verified M27 terminal from
`curve40-preconditioned-censored-followup-v1`, completes in3.665 seconds and
performs zero point searches. The [sealed preparation](../../artifacts/generated-results/elliptic-curves/curve40_exact_maximum_preparation_v1/result.json)
has now completed its first bounded100-call pass. It does not repair or restart any old
censored search and does not claim a new rank.

## Exact maximum preparation on the known M28 curve188

Curve188 / ICARM619 already has a certified M28 subgroup. Its historical
productive bank has nine classes of doubled norms20/24, and its pairwise
bank stays inside that original span. The complete maximum-class catalogue
now supplies sixteen norm12 parents outside the original productive span.
Both exact generic CVP solvers agree; preparation takes4.993 seconds under
120 seconds/1GiB and performs no point searches.

The [sealed preparation](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_preparation_v1/result.json)
uses the verified M28 terminal of `curve188-pairwise-v3-discovery-v1`.
This is the known M28 input, not a newly discovered rank28 curve. The first bounded100-call pass below tests whether parents outside the old
span can provide an additional certified direction. Public higher-rank points or curve302
winning labels are not imported into its parent selection.

## Curve40 norm12 first pass completed

Curve40 completes100 point calls in172.833 seconds and full independent
replay in227.532 seconds: M27 remains M27, with zero map or point timeouts.
The new landscape has759 centres. Its verified suffix resumes at centre50,
policy `preconditioned_full`, leaving709 centres in this ordering.
This is a bounded miss, not a rank upper bound.

- [Replayed result](../../artifacts/generated-results/elliptic-curves/curve40_exact_maximum_v1/result.json).
- [Exact suffix cursor](../../artifacts/generated-results/elliptic-curves/curve40_exact_maximum_v1/suffix-queue.json).

The [exact centre-overlap audit](../../artifacts/generated-results/elliptic-curves/curve40_exact_maximum_v1/productive-centre-overlap.json)
reproduces on a second run: none of the759 new selected centres occurs,
even up to elliptic negation, among the424
centres of `curve40-bounded-maps-v3-discovery-v1`. This compares the two
named selections only, not all historical or equivalent quartic boxes.

The next prioritized branch is the newly prepared norm12 bank on curve188's
known M28 basis. Its objective is a certified29th direction; preparation
alone does not establish one. All three M27 norm12 branches retain their
verified continuation cursors.

## Curve188 M28 exact maximum first pass

The norm12 bank completes100 point calls in204.753 seconds. Full independent
replay passes in334.225 seconds, including exact geometry, maps, rational
points, standalone independence and full-cloud audits. M28 remains M28;
there are zero map or point timeouts. No29th direction or rank upper bound
is established.

The784-centre landscape has a verified suffix starting at centre50,
policy `preconditioned_full`, leaving734 centres in its ordering.

- [Replayed result](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_v1/result.json).
- [Exact suffix cursor](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_v1/suffix-queue.json).

The exact centre-overlap checker reproduces two audits. None of the784
new centres coincides up to elliptic negation with either the earlier
450 productive centres or1727 pairwise centres, on the
same curve and initial basis. This verifies new selected-centre exposure;
it does not compare all historical quartic boxes or predict a rank gain.

- [Productive overlap receipt](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_v1/productive-centre-overlap.json).
- [Pairwise overlap receipt](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_v1/pairwise-centre-overlap.json).

This known M28 fibre remains the closest active norm12 branch to a certified
29th direction. The next bounded pass resumes its verified suffix using
cached exact geometry, with full saved-cloud reconciliation if required.

## Curve188 cached norm12 continuation to400 calls

Three further100-call passes complete and independently replay. Each
retains M28, with zero map or point timeouts. The old cursors are marked
consumed by their verified descendants, so the first400 calls are not
rescheduled. The latest784-centre landscape resumes at centre200, policy
`preconditioned_full`, leaving584 centres in this finite ordering.

| Cached pass | Cumulative calls | Search seconds | Cached replay seconds |
| --- | ---: | ---: | ---: |
| 1 | 200 | 93.603 | 2.581 |
| 2 | 300 | 89.666 | 3.468 |
| 3 | 400 | 95.615 | 4.407 |

- [200-call result](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v1/result.json).
- [300-call result](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v2/result.json).
- [400-call result](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v3/result.json).
- [Current suffix](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v3/suffix-queue.json).

Cached replay reuses the sealed original exact geometry and verifies the
new retained receipts. The finite-policy misses do not prove exact rank28
or rule out29. Further exposure, if scheduled, starts from the last verified
suffix; the goal of a new certified29th direction remains open.

A fourth cached norm12 pass on curve188 completes100 additional calls in
100.102 seconds and passes independent cached replay in5.191 seconds.
M28 remains M28 after500 cumulative calls, with zero point or map timeouts.
The [500-call result](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v4/result.json)
and [current suffix](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v4/suffix-queue.json)
resume at centre250 of784, policy `preconditioned_full`, leaving534 centres.
The old400-call cursor is marked consumed. No29th direction or rank upper
bound is established by this bounded miss.

The fifth cached norm12 pass on curve188 completes100 additional calls in
96.595 seconds and passes independent cached replay in6.115 seconds.
M28 remains M28 after600 cumulative calls, with zero point or map timeouts.
The [600-call result](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v5/result.json)
and [current suffix](../../artifacts/generated-results/elliptic-curves/curve188_exact_maximum_cached_v5/suffix-queue.json)
resume at centre300 of784, policy `preconditioned_full`, leaving484 centres.
The old500-call cursor is marked consumed. No29th direction or rank upper
bound is established by this bounded miss.
