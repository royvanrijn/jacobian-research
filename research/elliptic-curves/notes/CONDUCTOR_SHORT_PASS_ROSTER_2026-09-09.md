# Conductor-directed native-seed first passes

These bounded searches use exact conductor inequalities only to schedule
certified native seeds. The comparison is pinned at
`artifacts/local/elliptic-curves/future-conductor-payoff-v1/report.json`.
Reaching a listed target requires new independent points and a fresh catalogue
check before claiming a conductor record. Bounded misses give no rank upper
bound. The broader rank29–32 route remains open alongside this roster.

## Inventory181: native M22, conductor target24

The R17 family11952 fibre at `-8/27` has certified lower bound22 in
`compact192_r17_results_v1.json`, entry `11952-019`. Its exact conductor is

```
801759366600408255731113968344306127138129416996405428108788127367893872409350
```

This is below the pinned reported rank24 benchmark, ICARM623. Two new
independent directions would meet that rank threshold; the inequality does
not establish their existence. The native preparation freshly rechecks the
generic17, exact scale, M22 certificate and retained gain history. Its one
productive generic class is51930, with independently replayed exact generic
CVP. Preparation completes in3.760 seconds under120 seconds/1GiB.

The [short-pass runner](../cas/run_short_seed_v3.py) retains V3 selection, both
maps, height125000,10 seconds per invocation, a100-total-invocation budget,
and immediate rebuild after the first standalone certified gain. Search and
replay each have1800 seconds/3GiB and one worker. No public higher point or
retrospective visibility oracle enters selection. Preparation and receipts
are under `artifacts/local/elliptic-curves/curve181-v3-preparation-v1/` and
`curve181-short-v3-discovery-v1/`, with sibling supervisors and frozen sources.
The selected policy exhausts after58 completed invocations, with no gain or
timeout: M22→M22. Search takes33.944 seconds, peak RSS325488640 bytes;
independent replay passes in4.281 seconds, peak302145536 bytes. The
[sealed result](../../artifacts/generated-results/elliptic-curves/curve181_short_v3_v1/result.json)
has terminal reason `FINITE_POLICY_EXHAUSTED_NO_CERTIFIED_GAIN`. There is no
unvisited suffix within this selected policy. No rank24 point subgroup or
conductor improvement is found; the curve's exact rank remains unproved here.

## Inventory54: native M25, conductor target27

The next first pass uses R17 family103b2 at `-538/249`, source
`next24_r17_results_v1.json`, entry `103b2-025`. The pinned conductor payoff
requires two new directions beyond M25 to reach its reported rank27 benchmark.
Preparation freshly verifies the native17, M25 and productive generic class50670,
including independent generic CVP replay, in3.453 seconds under120 seconds/1GiB.
Only the existing certified history enters selection.

The same100-invocation, two-map, height125000 policy runs under
`artifacts/local/elliptic-curves/curve54-short-v3-discovery-v1/`, with inputs in
`curve54-v3-preparation-v1/` and the same1800-second/3-GiB supervisors.
The selected policy exhausts after88 invocations, with no gain or timeout:
M25→M25. Search takes75.874 seconds, peak RSS325197824 bytes; independent
replay passes in8.904 seconds, peak302723072 bytes. The
[sealed result](../../artifacts/generated-results/elliptic-curves/curve54_short_v3_v1/result.json)
has terminal reason `FINITE_POLICY_EXHAUSTED_NO_CERTIFIED_GAIN`; there is no
unvisited suffix within this selected policy. No rank27 or conductor result
is obtained, and no curve-rank upper bound is inferred.

## Inventory43: native M25, conductor target27

The next first pass uses R17 family08f72 at `1361/72`, source
`fresh_r17_paired_results_v1.json`, entry `08f72-070`. Its pinned conductor
target is27, requiring two additional independent directions. Preparation
freshly verifies native17, M25 and the two productive generic classes79534
and126516, with independent generic-CVP replay, in4.808 seconds under
120 seconds/1GiB. No higher point or visibility oracle enters selection.

The same100-total-invocation policy runs under
`artifacts/local/elliptic-curves/curve43-short-v3-discovery-v1/`, with preparation
in `curve43-v3-preparation-v1/` and separate1800-second/3-GiB search/replay
supervisors. Both maps and immediate gain-triggered rebuilding are retained.
The first100 invocations complete without timeout or certified gain: M25→M25.
Search takes76.591 seconds, peak RSS326332416 bytes; independent replay
passes in20.821 seconds, peak302325760 bytes. The
[sealed result](../../artifacts/generated-results/elliptic-curves/curve43_short_v3_v1/result.json)
has terminal reason `CHART_BUDGET_EXHAUSTED`. Its `continuation-queue.json`
retains centre50 /quartic-minimized as the next policy to consider.
The unvisited suffix remains UNKNOWN. No rank27 or conductor result is
obtained; this prefix gives no curve-rank upper bound.

The cached continuation in `curve43-cached-v3-discovery-v1/` completes the
remaining96 invocations,196 cumulative, with no gain or timeout: M25→M25.
Search takes73.969 seconds, peak RSS322695168 bytes; replay passes in3.504
seconds, peak273539072 bytes, inheriting the bound prior full CVP proof and
freshly checking the combined point cloud. The
[final sealed result](../../artifacts/generated-results/elliptic-curves/curve43_cached_v3_v1/result.json)
has terminal reason `FINITE_POLICY_EXHAUSTED_NO_CERTIFIED_GAIN`. The original
queue's `queue-resolution.json` binds this completion. No suffix remains in
this selected policy; no rank27 or conductor result is obtained.
