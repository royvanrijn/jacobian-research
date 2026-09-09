# Curve90: bounded productive-anchor discovery attempt

This is a new search for directions beyond the certified27-point subgroup on
MW16 curve `new-20260906-90`, parameter `-1867/270`. It uses the eight
historically productive generic anchors prepared in the
[performance and transfer study](V3_FUTURE_SEARCH_PERFORMANCE_2026-09-08.md).
The two earlier pointwise masked nulls remain null controls. This attempt
does not claim that the modified selector passed those controls or that a
new direction is guaranteed.

**Completed:** all752 invocations finish without timeout, and independent
replay passes. The adopted subgroup remains27; the925-point retained cloud
also certifies27 modulo2,3,5. This is a bounded no-gain result, not exact rank27.
The search supervisor takes525.047 seconds and independent replay105.913
seconds. The initial landscape takes31.818 seconds and selects376 centres.
Each coordinate policy completes376 invocations; no equal-coordinate pair is
skipped in this run.

Point-backend wall times total220.903 seconds for quartic-minimized maps and
239.594 seconds for factor-free maps. Together they account for87.7% of search
wall time. This describes this run; it is not an end-to-end speed comparison
with a different historical campaign. Individual chart receipts occupy
10,879,954 bytes, and the single complete-cloud snapshot12,237,045 bytes.

## Frozen method and limits

The actual starting basis is the exact generic-first M27 packet, recertified
before search. All eight anchors enter full extension scoring; at M27 this
means16,384 cosets. The existing V3 extension shortlist, exact parity CVP and
centre ordering are reused. Integer-scaled LDL accelerates the search CVP;
independent replay uses the original rational solver.

For each selected centre, the search tries the historical quartic-minimizing
map first, then the factor-free map. Equal projective horizontal maps are
deduplicated exactly. Both use height125000 and a10-second backend cap. The
campaign has at most4096 actual invocations and six epochs, stops at a
certified lower bound32, and rebuilds immediately after a certified gain.
Search and replay each have a7200-second,3-GiB, one-worker supervisor limit.
These are limits, not completion estimates.

Finite quotient tables and point columns persist within the worker. Every
gain requires a fresh standalone complete finite-group certificate. Each
chart is stored once in a hash-linked receipt; one complete cloud is audited
modulo2,3,5 at each epoch end. A larger complete-cloud rank than the adopted
basis stops for reconciliation. Censored boxes and failed finite admissions
remain explicit; no rank upper bound follows.

The final checker regenerates the full landscape with rational CVP, checks
centre order and both maps, replays every retained square/point witness,
binds the complete cloud to its chart receipts, and rechecks gain provenance
and finite rank certificates. A search terminal is not a completed replay.

## Operation and evidence

The search and independent replay are terminal and complete. Their source
files are also retained in `frozen-sources.zip` under the campaign directory.
No new-rank result is asserted here, and no continuation of this exhausted
initial schedule is implied.

All paths below are relative to `research/artifacts/local/elliptic-curves/`:

- `curve90-productive-v3-discovery-v1/protocol.json`: exact scope, limits,
  software versions and narrowly bound inputs/sources;
- `curve90-productive-v3-discovery-v1/progress.json`: last completed chart
  and certified adopted lower bound;
- `curve90-productive-v3-discovery-v1/epoch-*/`: landscapes, immutable chart
  receipts, gain packets and whole-cloud audits;
- `curve90-productive-v3-discovery-v1/terminal.json`: search completion,
  pending separate replay;
- `curve90-productive-v3-discovery-v1/verified.json`: independent replay,
  present only after it completes;
- `curve90-productive-v3-discovery-v1-supervision/`: owned process identity,
  enforced limits and worker log.

Runner: [`run_curve90_productive_v3.py`](../cas/run_curve90_productive_v3.py).
Its `search` action preserves and replays completed chart receipts when
resumed with unchanged inputs and code. Incomplete landscapes are archived
before reconstruction. `replay` is a separate action. Invoke both through
the shared bounded supervisor; do not start a duplicate worker merely because
a progress file is old. A lifetime file lock also prevents overlap.

Existing searches, their frozen source files and all prior proofs remain
unchanged. Rank/conductor novelty must be assessed against the current
catalogue after a new independently verified result, not inferred from the
campaign name or an incremental counter.

## Next candidate, separate from this completed roster

An exact comparison of the retained conductor certificates identifies
`new-20260906-92`, MW16 family02 at3161/432, as a one-direction opportunity:
its certified bound is26 and its exact conductor is about0.063725 times the
reported rank-at-least27 minimum. A new certified27th direction would beat
that reported benchmark by about15.7 times, subject to a fresh catalogue check
when the gain is proved. This is a scheduling payoff, not evidence that the
extra direction exists.

A fresh public catalogue download contains682 curves. Its reported rank27
conductor minimum is still ICARM614, and none of its curves has the rational
`j`-invariant of this local curve92. Missing reported conductors remain missing.
The pinned download, fetch hash, benchmark comparison and exact `j` comparison
are under `artifacts/local/elliptic-curves/future-catalogue-check-v1/`.
The conductor payoff audit is under `future-conductor-payoff-v1/`.
No curve92 search has been launched as part of the curve90 roster.

## Fresh pairwise-parent short pass

After the separate curve92 productive search and four verified25→26 transfers,
a new curve90 experiment uses the100-total-invocation policy. Its seed remains
the certified native M27 at `-1867/270`; only directions beyond those27 count
as new. The [input adapter](../cas/prepare_curve90_parent_input.py) consolidates
the sealed M27, generic metric and eight productive words without modifying
the old preparation or including the masked oracle file.

All28 distinct pairwise XORs of the eight productive masks are nonzero and
outside the original bank. Their generic minima are computed exactly and
independently replayed in1.374 seconds under120 seconds/1GiB. Norm counts are
3.5 (one),5.5 (three),6 (seven),7.5 (three),8 (six),9.5 (eight). All fit within
the existing16-per-shell quota; no additional parent truncation is introduced.
These are untested derived classes, not a claim of inherited productivity.

The [short-pass wrapper](../cas/run_short_seed_v3.py) retains the original
parent-subset landscape and both maps, with height125000,10 seconds per
invocation,100 actual invocations across all adaptive epochs, target32,
immediate rebuild on certified gain, and separate1800-second/3-GiB search and
replay bounds. One worker runs at a time. The prior752-chart search is not
restarted. Inputs and receipts are under
`artifacts/local/elliptic-curves/curve90-parent-input-v1/`,
`curve90-pairwise-preparation-v1/`, and
`curve90-pairwise-short-v3-discovery-v1/`.

Search completes100 invocations in183.010 seconds, peak RSS328642560 bytes,
with no timeouts and no certified gain: M27→M27. The selected1314 centres are
disjoint up to sign from the previous376 centres on the identical ordered
basis; `initial-centre-disjointness.json` binds this exact comparison.
The terminal reason is `CHART_BUDGET_EXHAUSTED`, not policy exhaustion.
Independent replay passes in355.775 seconds, peak RSS303296512 bytes.
The [sealed100-chart result](../../artifacts/generated-results/elliptic-curves/curve90_pairwise_short_v3_v1/result.json)
is exported. The queue recorder preserves centre50, quartic-minimized map,
as the next policy to consider in `continuation-queue.json`. The untested
suffix remains UNKNOWN; this bounded prefix does not rule out rank28 or29.

## Cached second batch

The [cached continuation](../cas/run_cached_seed_v3.py) now extends the verified
prefix in `artifacts/local/elliptic-curves/curve90-pairwise-cached-v3-discovery-v1/`.
Its cumulative cap is200:100 inherited receipts and at most100 new invocations.
The old full-landscape CVP replay remains a hash-bound proof dependency;
all point receipts and the extended cloud are checked again. A gain still
rebuilds immediately, and a new landscape requires full independent CVP replay.
Height125000,10 seconds per invocation, one worker, and separate1800-second /
3-GiB search and replay limits are unchanged. This targets directions beyond
the certified M27 without recomputing the already-verified1314-centre landscape.
The batch completes100 new invocations,200 cumulative, with no timeout or
gain: M27→M27. Search takes74.392 seconds, peak RSS325206016 bytes. Replay
passes in10.000 seconds, peak297091072 bytes, reusing the bound prior full
CVP verification. The
[sealed continuation result](../../artifacts/generated-results/elliptic-curves/curve90_pairwise_cached_v3_v1/result.json)
records the cumulative count and explicit100-chart inheritance.
The original queue's `continuation-head.json` points to this verified
descendant, with centre100 /quartic-minimized as the next policy to consider.
The remaining suffix stays queued while a different native family receives
a first pass; no absence of rank28 or29 is inferred.
