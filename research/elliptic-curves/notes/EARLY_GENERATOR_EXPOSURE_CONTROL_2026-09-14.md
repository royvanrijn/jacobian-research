# Early generator exposure: actual accessibility, failed performance gate

**The two-block performance gate fails.** Both arms recover 49 directions,
including 28 later gains, and reach 12 of 24 three-direction targets. Refreshing
uses 904.600762 complete CPU seconds versus 779.756541 for the fixed bank.
All 28 later candidate gains occur on newly enabled anchors and have exact
finite visibility witnesses that survive the generator-deletion and translated-
representative checks below. This demonstrates the representation mechanism on
these controls, without a repeatable later-gain performance advantage.

The [frozen protocol](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/protocol.json)
has SHA256
`4d42c491807a0eef8d32d3f73d4d110e87b86e10138cc1e68d488d52b73e36e5`.
The [supervisor receipts](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/supervision.json)
record all 48 completed attempts and their charged process CPU. All independent
worker replays pass, with no preparation or infrastructure failure. The
[completion record](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/completion.json)
binds the reports and replay archive. Mathematical authority is
`EC-CANCELLATION-EARLY-GENERATOR-20260914` in
[MATH_STATUS.json](../../MATH_STATUS.json). No fresh-fibre continuation ran.

The [preceding exact audit](SUBGROUP_ACCESSIBILITY_CAUSAL_AUDIT_2026-09-14.md)
found four new-generator anchors that exposed a withheld direction but were
never searched. Removing the relevant generator loses all four finite
witnesses. That motivates changing which compatible anchors are visited first.
It does not turn the [failed comparison](BASIS_AWARE_AMPLIFICATION_2026-09-14.md)
into a successful performance result.

## Charged result

The [final report](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/summary.json)
retains all zero and partial outcomes. Later gains mean gains after the first
positive full cloud, rather than all but one point of that cloud.

| Measurement | Fixed bank | Rebuild and expose new generators first |
|---|---:|---:|
| Three-direction targets | 12/24 | 12/24 |
| Uncapped added directions | 49 | 49 |
| Later-cloud directions | 28 | 28 |
| Complete CPU seconds | 779.756541 | 904.600762 |
| Point calls | 926 | 715 |
| Initial/rebuilt bank component CPU seconds | 93.547619 | 322.120590 |
| Point-call component CPU seconds | 561.806797 | 433.658072 |
| Final replay component CPU seconds | 44.130469 | 56.478093 |

The later-direction/CPU ratio is **0.861990**, with central 97.5% paired
family-bootstrap interval **[0.731502, 1.020819]** and no undefined draws.
The observed rate is 13.80% lower; the interval does not establish a population
disadvantage. The fixed block ratios are 0.832547 and 0.893055. Block 0 has
17 later gains per arm; block 1 has 11 per arm. Both therefore fail the required
strict rate improvement. Neither the aggregate nor repeatability gate passes.

The candidate saves 128.148725 CPU seconds in point calls but spends 228.572971
more on bank construction, plus more final replay and other work. These
components diagnose the cost problem; they are not a counterfactual prediction
of what an optimized engine would recover. All 98 **arm-level** additions have
native independent rank certificates; they are not 98 distinct discoveries.
Across 1641 point calls, 1618 complete and 23 time out. A timeout contributes
its cost and any certified points, without claiming finite-box coverage.

All 24 initial banks agree exactly between arms. The entire pre-first-gain
point-cloud prefix agrees in 23 pairs. The remaining pair,
`4523c54d951c995e467d`, has zero gains in both arms and reaches the same CPU
limit after 46 versus 45 calls; this is not an input or selection mismatch.

## Actual generator-dependent visibility

There are 34 verified rebuilds exporting 631 anchors that use added generators.
The [exposure audit](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/exposure-audit.json)
records 324 such calls among 335 post-refresh calls. All 28 later candidate
directions, on 16 curves, come from these calls. Thus the old failure to reach
the new anchors is repaired in this experiment.

The predeclared [mechanism audit](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/mechanism.json)
binds each actual gained point, pre-call subgroup, required generator, anchor
word, square and exact point map. All 28 points and their negatives are outside
every original-bank box at height 125000. For each row, deleting at least one
participating generator coefficient also loses visibility in the corresponding
factor-free/prime-neighbour dictionary. Block 0 supplies 17 witnesses and
block 1 supplies 11, so the declared finite mechanism gate passes.

A separately frozen, post-protocol
[representative sensitivity](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/representative-protocol.json)
then uses both signs of

\[
P+\sum_i k_iG_i,\qquad k_i\in\{-1,0,1\},
\]

where the `G_i` are **all** generators added above M18 before that call. There
are at most two, giving at most 18 representatives per row. The same dictionary
is checked on all cached original maps, each coefficient-deletion dictionary,
and the executed new chart. No new map preparation, CVP or point search runs.
All 28 old-bank and strict deletion witnesses survive, across **47730 exact
coordinate-square evaluations**; winning maps are independently reconstructed
and mapped back. The [result](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/representative-audit.json)
preserves every representative and minimum. It changes no primary decision rule.

For example, the first recorded witness in each fixed block is:

| Block / case / call | Anchor change | Executed coordinate | Original-bank minimum over translated representatives |
|---|---|---|---:|
| 0 / `fcf6ce8e6a518eb715d2` / 017 | `Q0 + G19` | `[-67767:19469]` | 215855850046008 |
| 1 / `d572443b6b5c9f96a1a8` / 004 | `Q0 + G19` | `[-110855:40413]` | 11333505542 |

Here `Q0` is the recorded old-subgroup word and `G19` is the independently
admitted nineteenth basis vector. The corresponding deleted-anchor minima
are also strictly above 125000; their exact integers, full words and point
coordinates are in the linked packets.

These are finite accessibility changes. They do **not** assert that the same
quotient direction is inaccessible via every representative `P+S`, nor that
the baseline cannot recover other directions. No finite dictionary minimum is
promoted to a whole-coset minimum or arithmetic exclusion. The equal aggregate
gain counts and higher candidate CPU remain the performance result.

The retained `rank_growth.py` diagnostics use numerical canonical-height Gram
data only after execution. Fifteen rows have a one-dimensional added block,
where a nonzero last-increment correlation of one is automatic. Thirteen have
a two-dimensional block. The latter diagnostics depend on the recorded basis
order and are not a target-independent selector or an exact height certificate.

## Exploration beyond the fitted range

The audit records **1054** neighbour calls with exact coordinate witnesses
beyond the old height-125000 range and outside prior completed boxes on the
same anchor; **1042** of those calls complete. They comprise 605 baseline and
449 candidate calls. These are coordinate witnesses, not quartic-square or
elliptic-point existence assertions. A failed boundary probe stays UNKNOWN.
The candidate skips only 21 exactly equivalent completed boxes. Neither arm
uses fitted radius support to discard a distinct unsearched job.

## Fixed scope and comparison

The roster has 24 previously untested whole-j groups, four from each of the six
retained R17 families. Each has an independently checked rank-18 seed and a
retained endpoint of rank at least 21. All preceding scheduler and basis-ablation
CPU groups are excluded. The roster is selected by a fixed hash order from the
retained eligibility audit, without using the new search outcomes.

There are two preassigned validation blocks, each containing two groups per
family. Both blocks are fixed before execution. The second block is not selected
or altered after observing the first. These are known retained-corpus controls,
not fresh rank discoveries or externally pristine population data.

Both arms aim to add **three directions** under the same **40 working CPU-second**
allowance, with 100 hard process CPU and 150 wall seconds per arm. Initial banks,
every rebuild, rational-CVP comparison, native anchor checks, map preparation,
all point calls, complete-cloud reconciliation and rank certification are
charged. Final independent point/policy/rank replay is mandatory and also charged
to complete process CPU. Atomic operations can overrun the working launch cap.

Both arms reconstruct the same initial bank from the same 16 retained generic
anchors, use the same finite places through 500 and retain every independent
direction in each returned cloud. The fixed arm keeps its initial search basis
and bank. The candidate independently rebuilds from the complete enlarged
subgroup after a positive cloud when the target and CPU allowance permit.

After each verified rebuild, the candidate stably partitions its exported
anchors into:

1. anchors with a nonzero coefficient in the most recently admitted block;
2. remaining anchors involving earlier admitted generators;
3. anchors wholly in the original rank-18 subgroup.

The initial order and the relative order within each group remain unchanged.
The bank itself is not edited. Each visited anchor uses its factor-free box,
then its at-most-two retained prime neighbours, all at height 125000 and subject
to the same gain/CPU stops. No residue-weight, radius or target-aware numerical
fit enters either arm. Exact completed-box deduplication survives basis changes;
zero fitted support never suppresses a distinct unsearched box.

## Prerequisites and failure semantics

All 24 starting rank-18 and endpoint rank-at-least-21 certificates pass independent
preflight, which took 28.431219 component CPU seconds. A separate order/fault
regression replays 10 retained point calls without point search, then injects
the actual Sage `AlarmInterrupt` class at the next bank boundary. It preserves
the complete rank-19 cloud, records an UNKNOWN bank and passes final independent
verification. The original failed preparation is not retried or used for search.

One separately labelled development integration uses the already tested
`6f5afe31891f7dc1938b` control. The new order recovers two directions, including
one later gain, and independently verifies the rank-19 replacement bank. It
exhausts the working allowance during the next rank-20 rebuild; that bank remains
UNKNOWN and unused. The complete run retains a certified rank-20 endpoint and
costs 42.246007 internal CPU seconds including final replay. This is software
integration evidence selected from the old mechanism audit, not a validation
outcome. Its original primary-comparison result is untouched.

The declared 40-second control allowance is unchanged despite that development
limit. Failed preparation remains charged and blocks promotion under this
protocol. It is not an arithmetic exclusion and cannot invalidate already
certified point prefixes. Unhandled infrastructure failures stop dispatch;
completed attempts must not be retried or replaced.

## What the generator can change

The existing [pointed-chart identities](CURVE302_V3_ANATOMY_AND_FINITE_ATLAS_2026-09-11.md#the-exact-pointed-quartic)
give, on a short model and away from the separately handled exceptional points,
`t_Q(P)=(y(P)+y(Q))/(x(P)-x(Q))`. After the fixed rational parameter matrix,
the height is the maximum of the absolute primitive numerator and denominator.
Admitting `G` permits an anchor `Q=Q0+cG` that was unavailable in the original
subgroup. Deleting its coefficient returns `Q0`, so the postmortem can compare
the two exact rational representations of the same recovered point.

The [covering-height calculation](POINTED_CHART_HEIGHT_BOUNDS_2026-09-13.md#the-checked-height-inequality)
explains the geometric candidate: changing the anchor changes `2P-Q`, and
bilinearity gives

\[
\widehat h(2P-Q_0-cG)-\widehat h(2P-Q_0)
=c^2\widehat h(G)-2c\langle 2P-Q_0,G\rangle.
\]

This can decrease the leading height term, but chart-dependent real and finite
distortion still controls actual parameter height. The pairing with the missing
`P` is unavailable during selection. Schur reduction alone therefore neither
predicts useful anchors nor proves cheap recovery. The current policy tests the
simple target-independent rule of exposing the latest generators first; exact
coordinate and deletion witnesses are evaluated only after both blocks finish.
The birational recovery chart and the degree-four map `P -> 2P-Q` remain distinct.

## Decision and attribution gates

Across all 24 pairs, later-cloud directions per complete CPU must improve by
at least 10%, with the central 97.5% paired family-bootstrap interval wholly
above 1, using 10000 draws and seed 20260914. Candidate uncapped total directions
and three-direction target completions must both match or exceed the baseline.
Later means after the first positive **full cloud**, so easy first recovery
cannot pass the primary endpoint.

Each preassigned block must separately have a later-direction/CPU ratio above 1
and an actual later gain on an anchor requiring a newly admitted generator.
After both blocks finish, a bounded postmortem must bind the exact generator,
pre-call subgroup, anchor word, gained point, square and coordinate map. For
each such gain it prepares old-bank and single-generator-deletion maps without
target input, then evaluates the discovered point and its negative. Promotion
also requires a strict finite deletion witness in both blocks. These two
representatives do not exhaust a subgroup coset.

The postmortem has a 300-CPU-second hard diagnostic limit and adds no point
searches or CVP-bank searches. It reuses `rank_growth.py` for numerical Schur
and cascade diagnostics only after execution. Missing work or an undefined
comparison cannot pass a gate. All zeros, partial gains, timeouts and failed
preparations remain in the denominators. Shared preflight and attribution CPU
are reported separately and charged entirely to the candidate in conservative
sensitivity calculations.

The two postmortems cost 24.293446 and 23.396765 complete child CPU seconds.
Their component meters are 22.609630 and 22.757197 seconds. Charging all
47.690211 diagnostic CPU seconds and the shared preflight entirely to the
candidate lowers the later-rate ratio to 0.795084. The separately labelled
development integration costs 42.246007 internal CPU seconds and the order/
alarm regression 8.169124 component seconds. These are not a complete account
of development or historical seed/parent creation; packaging overhead is also
outside the arm comparison. The original complete arm charges are unchanged.

Repeatable later-gain advantage and general target-independent prediction remain
unproved. A next implementation gate should address the measured bank cost on
retained transitions, while independently checking the complete subgroup and
every reused or newly selected anchor. Reusing a compatible old map or metric
block requires exact bindings; changing the bank or order requires a distinct
policy, not a retroactive timing correction. A faster reconstruction alone is
still insufficient: a separately frozen disjoint comparison must earn any later-
gain performance claim. The current control is complete and is not restarted.

Broader comparison against an established
implementation also remains open. Cremona's
[eclib/mwrank](https://johncremona.github.io/mwrank/index.html) is available
locally; its [Sage interface](https://doc.sagemath.org/html/en/reference/libs/sage/libs/eclib/interface.html)
supports processing a supplied subgroup and a separate point search, whose
height convention and default saturation differ from this pointed-chart
engine. An external benchmark needs its own matched input, limit and cost
protocol. No state-of-the-art claim follows from an internal comparison.
The local Sage dependency is eclib 20250627; the separate system `mwrank`
executable reports 20231211. The earlier
[R17 canary](../../elkies-k3/R17_PROSPECTIVE_CRT_RANK_JUMP_EXPERIMENT_2026-09-04.md#phase-4-identical-bounded-search)
exhausted 300 seconds during initialization before any search. That retained
failure is an operational boundary, not a point-search baseline or authorization
to repeat the old campaign.

The underlying local-solubility and covering-height theory has established
foundations, including [Fisher–Sills](https://arxiv.org/abs/1103.4944). The present
test concerns exploitation of a changing supplied subgroup within a bounded
search, not novelty of those foundations.

## Entry points

[`cancellation_basis_early.py`](../cas/cancellation_basis_early.py) owns the
frozen roster, preflight, dispatcher and report.
[`verify_cancellation_basis_early.py`](../cas/verify_cancellation_basis_early.py)
independently reconstructs the required visit order and replays full clouds.
[`audit_cancellation_basis_early.py`](../cas/audit_cancellation_basis_early.py)
is the separately metered post-execution attribution stage. A live process
is no longer running for this comparison. The
[replay bundle](../../artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1/replay-bundle.tar.gz)
contains 14490 files and 29432787 compressed bytes, with SHA256
`ce46b023edc15b05571b19e9a2e20176bd0a9bf43c398619c015e480a44733f8`.
All bytes, imports, source seals, 48 arm input bindings and 24 extracted endpoint
proof bindings pass in an empty root. This is a packaging preflight, not another
full arithmetic replay. Failed development work and the original dependency
record preceding the representative supplement are retained verbatim.

For a selected arithmetic replay, extract into a separate workspace with the
matching Sage/PARI runtime and use an output path that does not exist:

```sh
sage -python elliptic-curves/cas/retain_cancellation_basis_early.py replay \
  --case fcf6ce8e6a518eb715d2 --arm basis_refresh \
  --output /tmp/early-generator-arm-replay.json
```

This checks the saved clouds, maps, native ranks and independent bank-selection
receipts without dispatching point searches. It does not recompute the producer
or reference CVP enumeration or supply a new search-completeness proof.
