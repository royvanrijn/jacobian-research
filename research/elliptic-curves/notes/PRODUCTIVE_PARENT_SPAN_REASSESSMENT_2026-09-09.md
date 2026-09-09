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
