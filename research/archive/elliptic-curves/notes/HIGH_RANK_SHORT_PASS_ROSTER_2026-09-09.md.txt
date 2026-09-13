# Further native M27 first passes

These tests seek new independent directions beyond existing certified M27
seeds, toward the rank29–32 objective. They use the frozen100-total-invocation
V3 policy, both maps at height125000,10 seconds per invocation, and immediate
rebuild after a certified gain. Search and replay each have1800 seconds/3GiB
and one worker. Higher public points and retrospective visibility oracles do
not enter selection. No bounded miss is a curve-rank upper bound.

## Inventory40: family074d9 at2818/1535

Native preparation rechecks the generic17, exact model transport and M27
certificate in `fresh_r17_paired_results_v1.json`, entry `074d9-007`.
All nine historical productive classes are retained:7615,25067,31428,44182,
78315,78981,91675,94715,117749. Their generic CVPs independently replay.
Preparation takes14.423 seconds under120 seconds/1GiB.

A separate [exposure preflight](../cas/prepare_short_pass_exposure.py) freezes
the424-centre landscape without point searches and compares rational centres
up to sign on the identical equation against the retained301-chart adaptive
follow-up and49-chart specialized-parity experiment. Only two centres overlap,
at zero-based indices416 and421. Each centre has at least one invocation, so
neither overlap can enter the initial100-invocation epoch. This comparison is
limited to those two named histories, not every historical search. Preflight
takes41.896 seconds under300 seconds/1GiB; full independent landscape replay
still follows the point search.

Inputs and receipts are under
`artifacts/local/elliptic-curves/curve40-v3-preparation-v1/` and
`curve40-short-v3-discovery-v1/`, with a preflight supervisor and separate
search/replay supervisors. The frozen protocol binds the old histories used
in the comparison. The attempt is intentionally stopped after180.634 seconds
with zero point-search invocations. A live stack inspection identifies the
first `hyperellminimalmodel`/`hyperellred` call in the minimized mapper; the
point-search timeout had not begun. This is a mapping-stage failure, not a
bounded point-search miss. The native seed, preflight, source archive and
`abort.json` are preserved, with no new rank claim.

A separate [bounded map worker](../cas/bounded_map_worker.py) tests the same
first centre in isolated processes with5-second/1-GiB limits. Minimization
hits the strict wall timeout (6.014 seconds including shutdown). The
factor-free map completes in0.698 seconds and passes a separate exact rational
quartic/map-identity check. The diagnostic and both supervisor receipts are
under `artifacts/local/elliptic-curves/curve40-bounded-map-diagnostic-v1/`.
No point search is performed by this diagnostic.

### Replacement with bounded map construction

The [new runner](../cas/run_bounded_maps_seed_v3.py) gives each coordinate map
a separate5-second/1-GiB worker bound, followed by the existing10-second point
search bound. The100-invocation allowance counts actual point searches;
map timeouts are recorded separately. Resource-limited minimized maps do not
prevent trying factor-free coordinates. Both mapper policies are retained.

[Map receipts](../cas/bounded_map_receipts.py) bind inputs, worker/mapper
sources, logs, supervisor outcomes and completed results. Every completed
map receives an independent exact rational quartic/map-identity check.
Replay reads only the attempted map prefix and never restarts minimization,
including for timeout receipts. Any entirely exhausted but map-censored policy
is labelled `CENSORED_FINITE_EXPOSURE`, not complete point coverage. The full
landscape still receives independent rational-CVP replay. No canonical map
minimality or absence of points in a censored map is claimed.

Three [retained-input regressions](../tests/test_bounded_map_receipts.py) pass
in6.952 seconds: timeout replay starts no worker; the exact factor-free map
matches the retained diagnostic; and a deliberately corrupted quartic fails
the independent identity check even after its containing hashes are rebound.

The replacement runs in
`artifacts/local/elliptic-curves/curve40-bounded-maps-v3-discovery-v1/`, with
separate1800-second/3-GiB search and replay supervisors and one sequential map
worker. The interrupted attempt remains unchanged. The replacement completes
100 point-search invocations with no certified gain: M27→M27. No point search
times out; six map constructions hit their resource limit and remain explicitly
incomplete. The terminal reason is `CHART_BUDGET_EXHAUSTED`.

Search takes227.368 seconds, peak RSS691924992 bytes across its process tree;
replay passes in128.764 seconds, peak297160704 bytes. The
[sealed result](../../artifacts/generated-results/elliptic-curves/curve40_bounded_maps_v3_v1/result.json)
binds the100 executed charts, all attempted-map receipts and the full reference
landscape replay. No rank28 direction, conductor improvement or rank upper
bound is obtained. Censored maps and the remaining suffix stay UNKNOWN.

### Preconditioned minimization diagnostic

The installed PARI help documents an optional explicit prime list for
`hyperellminimalmodel`; minimality is then guaranteed only at those primes.
A separate [coordinate policy](../cas/preconditioned_pari_mapping.sage) first
constructs the exact factor-free integral model, then minimizes at2 and3 and
reduces the result. It composes and checks the complete rational transformation;
it makes no global minimality claim.

On the same first centre that stalled the raw minimized mapper, this completes
in0.758 seconds under5 seconds/1GiB. A separate pointed-quartic verifier checks
the exact map identity. The largest coefficient of the discriminant quartic
uses82 bits, versus91 for the factor-free map. This is one diagnostic, not a
point-search result or a calibrated replacement policy. Evidence is under
`artifacts/local/elliptic-curves/curve40-preconditioned-map-diagnostic-v1/`.

The four fixed winning-centre controls expose a limitation of minimization
at2 and3 only: the curve52 gain moves to height313417, outside125000. That
variant is therefore not adopted as a replacement. Its negative calibration
is preserved in `preconditioned-winning-controls-v1/`.

A [full-minimization variant](../cas/preconditioned_full_pari_mapping.sage)
instead applies the default PARI minimization to the factor-free integral
model, with the external timeout retained. On the previously stalled first
centre it completes in0.702 seconds and passes separate exact map verification;
evidence is in `curve40-preconditioned-full-map-diagnostic-v1/`.

The [four retrospective controls](../cas/audit_preconditioned_full_gain_coordinates.py)
complete in5.743 seconds. Positive-point heights are26302 (curve200),175691
(curve113),37261 (curve116), and110235 (curve52). The minimized lane retains
the three in-bound gains; curve113 still needs the unchanged factor-free lane,
which exposed its gain at82957. This is a fixed known-point calibration, not
a prospective discovery or a general equivalence of coordinate policies.
Results are retained in `preconditioned-full-winning-controls-v1/`. This full
preconditioned policy is now integrated as described below.

### Fixed follow-up of the six censored centres

[The targeted runner](../cas/run_preconditioned_censored_followup.py) freezes
exactly the six parent centres whose raw minimization timed out: indices
0,1,2,3,5,7. It inherits their independently verified CVPs, uses the full
preconditioned mapper through [bounded receipts](../cas/preconditioned_map_receipts.py),
and skips any coordinate box already completed in the parent. It searches no
higher point or new parameter and stops at the first certified new direction
for later V3 amplification. The actual supervisors impose120 seconds/2GiB
per search/replay phase, stricter than the inherited outer maxima.

All six new maps and point-search invocations complete without timeout or
gain: M27→M27. None duplicates a completed parent coordinate box. Search
takes16.366 seconds, peak RSS540106752 bytes; replay passes in1.944 seconds,
peak67428352 bytes. It checks inherited proof bindings, the exact alternative
maps and new point witnesses, and fresh whole-cloud certificates. These are
alternative coordinate boxes at the blocked centres, not an assertion that
the unknown old minimized boxes were reproduced.

- [Sealed six-chart follow-up](../../artifacts/generated-results/elliptic-curves/curve40_preconditioned_censored_v1/result.json).
- Raw evidence: `artifacts/local/elliptic-curves/curve40-preconditioned-censored-followup-v1/`.
- [General preconditioned two-map runner](../cas/run_preconditioned_seed_v3.py),
  prepared for future native seeds with the same map and point-search limits.

No rank28 point or conductor improvement is found. The previously censored
centres now have this additional tested coordinate policy; no absence result
or global rank upper bound follows.

## Curve71: prospective preconditioned first pass

Inventory `new-20260906-71`, native family103b2 at `3726/881`, starts from
the certified M27 in `retention24_r17_results_v1.json`, entry `103b2-733`.
Preparation independently re-certifies the seed and its eight productive
parent classes in12.829 seconds under120 seconds/1GiB. The new
[preflight adapter](../cas/prepare_preconditioned_exposure.py) freezes381
centres in36.581 seconds under300 seconds/1GiB, without point searching.
Compared with the earlier301-centre adaptive and49-centre specialized runs,
only indices316 and373 overlap up to sign on the identical equation.
All first50 centres are new relative to those named histories.

The general preconditioned two-map runner completes100 point invocations
over those50 centres: M27→M27, zero map timeouts and zero point timeouts.
Search takes135.395 seconds, peak process-tree RSS477290496 bytes; independent
replay passes in110.443 seconds, peak296165376 bytes. Preparation and preflight
are additional costs. All100 map workers complete, taking65.279 seconds in
aggregate, median0.644 and maximum0.750 seconds. This includes worker startup
and is not a controlled speedup measurement against the raw mapper.

Replay verifies the original rational-CVP landscape, exact map identities,
point witnesses, receipt chain and whole-cloud finite certificates. The mod2,
mod3 and mod5 lower bounds remain27. The100-box cap is exhausted, not the
381-centre policy. A sealed cursor retains centre50, `preconditioned_full`,
with331 centres left; it does not automatically launch another batch.

- [Sealed result](../../artifacts/generated-results/elliptic-curves/curve71_preconditioned_v3_v1/result.json).
- Raw evidence: `artifacts/local/elliptic-curves/curve71-preconditioned-v3-discovery-v1/`.
- Preparation: `artifacts/local/elliptic-curves/curve71-v3-preparation-v1/`.

No new rank28 direction, conductor improvement or rank upper bound is obtained.
The pass validates bounded execution on a fresh seed. Worker startup and
repeated factor-free construction are measured optimization candidates;
their removal has not been benchmarked. Existing historical searches remain
unchanged.

## Curve48: lean-worker prospective pass

The next native seed is inventory `new-20260906-48`, family11952 at
`2828/2015`, certified M27 from `next24_r17_results_v1.json`, entry
`11952-041`. Preparation takes12.968 seconds under120 seconds/1GiB.
Its eight productive parents yield390 selected centres. Exact preflight
takes37.978 seconds under300 seconds/1GiB and finds no overlap up to sign
with the98 centres from its original49-box extended run and subsequent
49-box specialized run. This comparison concerns only those named histories.

The [lean runner](../cas/run_lean_preconditioned_seed_v3.py) initializes PARI
directly through `cypari2` inside the bounded map workers. Before prospective
use, all16 fixed retained maps match exactly and independently verify; three
timeout/replay/corruption regressions pass. Detailed component timings are in
the [performance note](V3_FUTURE_SEARCH_PERFORMANCE_2026-09-08.md).

The prospective pass completes100 point-search invocations at the first50
centres, with zero map or point-search timeouts: M27→M27. Search takes90.348
seconds, peak process-tree RSS159113216 bytes. Its100 map workers take16.304
seconds in aggregate, median0.173 and maximum0.182 seconds. Independent replay
passes in117.099 seconds, peak294400000 bytes. Preparation and preflight are
additional costs; cross-seed timings are not a controlled speedup comparison.

Replay checks the full original rational-CVP landscape, exact map identities,
point witnesses, receipt chain and whole-cloud finite certificates. Mod2,
mod3 and mod5 lower bounds remain27. The
[new cursor recorder](../cas/queue_verified_preconditioned_pass.py) seals the
next action at centre50, `preconditioned_full`, leaving340 centres unvisited.
No continuation is automatically launched. Earlier censored maps, if present,
are explicitly separate from the suffix cursor; this pass has none.

- [Sealed result](../../artifacts/generated-results/elliptic-curves/curve48_lean_preconditioned_v3_v1/result.json).
- Raw evidence: `artifacts/local/elliptic-curves/curve48-lean-preconditioned-v3-discovery-v1/`.
- Preparation: `artifacts/local/elliptic-curves/curve48-v3-preparation-v1/`.

No new rank28 direction, conductor improvement or rank upper bound is obtained.
The fixed-map audit and prospective replay support use of the lean workers in
future passes; the full selected policy is not exhausted.

### Curve48 cached continuation

The [lean continuation](../cas/run_lean_cached_seed_v3.py) preserves the first
pass and inherits its independently verified landscape through exact bindings.
It replays the original100 point receipts without searching them again,
and preserves each inherited map's original worker command/path bindings.
Missing inherited receipts fail closed instead of restarting workers.

The next100 actual point invocations complete at centres50–99, making200
cumulative: M27→M27, zero map or point timeouts. Search takes88.352 seconds,
peak RSS158298112 bytes; replay passes in2.414 seconds, peak71933952 bytes.
Replay independently checks maps, point witnesses, receipt order and the
extended whole cloud; the full CVP proof is inherited from the verified parent.
Mod2, mod3 and mod5 lower bounds remain27. No global upper bound follows.

- [Sealed continuation](../../artifacts/generated-results/elliptic-curves/curve48_lean_cached_v3_v1/result.json).
- Raw evidence: `artifacts/local/elliptic-curves/curve48-lean-cached-v3-discovery-v1/`.
- Its [cursor recorder](../cas/queue_verified_lean_cached_pass.py) retains
  centre100, `preconditioned_full`, with290 centres unvisited. The original
  pass's `continuation-head.json` points to this verified descendant.

The budget is exhausted, not the selected policy. No rank28 direction or
conductor improvement is obtained.

### Curve48 second cached continuation

Chaining the same wrapper from the first cached pass adds100 new invocations
at centres100–149,300 cumulative. M27→M27, zero map or point timeouts.
Search takes97.222 seconds, peak RSS159670272 bytes; replay passes in3.348
seconds, peak80936960 bytes. Original map receipts continue to resolve to
their original source directories across both continuation generations.
Inherited independently verified CVPs are reused; all300 point receipts and
the extended cloud are replayed. The budget is exhausted, not the policy.

- [Sealed second continuation](../../artifacts/generated-results/elliptic-curves/curve48_lean_cached_v3_v2/result.json).
- Raw evidence: `artifacts/local/elliptic-curves/curve48-lean-cached-v3-discovery-v2/`.
- Cursor: centre150, `preconditioned_full`, with240 centres unvisited.
  The first continuation's `continuation-head.json` points here.

No new rank28 direction, conductor improvement or upper bound is obtained.

### Curve48 outside-span parent branch

A separate [302-motivated parent reassessment](PRODUCTIVE_PARENT_SPAN_REASSESSMENT_2026-09-09.md)
prepares16 generic classes outside curve48's old productive span. Their779
centres are disjoint up to sign from the248 centres in the named compared
histories. The first100 invocations complete without timeout or gain,
M27→M27. Search87.931 seconds and full replay225.810 seconds pass. Its729
unvisited centres form a separate queue from the old productive-bank suffix;
neither queue is a completed-policy or rank-upper-bound claim.

[Sealed outside-span pass](../../artifacts/generated-results/elliptic-curves/curve48_complement_v3_v1/result.json).

Its first cached continuation completes100 additional boxes,200 cumulative
within the outside-span branch, with no gain or timeout. Search90.711 seconds
and replay2.513 seconds pass. The679-centre suffix starts at centre100;
the independent earlier full CVP proof is inherited, while the extended
point cloud and exact receipts are checked again.
[Sealed outside-span continuation](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v1/result.json).

Curve48 complementary-parent continuation is now at300 cumulative calls,
still M27, with all100 new calls and their cached replay complete.
The [current result and cursor](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v2/result.json)
advance this branch to centre150; see the [canonical parent reassessment](PRODUCTIVE_PARENT_SPAN_REASSESSMENT_2026-09-09.md).

Curve71 now also has a replayed100-call sampled complementary-parent pass,
still M27 with no timeouts. This is separate from its prior productive-parent
pass. [Result](../../artifacts/generated-results/elliptic-curves/curve71_sampled_shell_v1/result.json);
see the [canonical parent reassessment](PRODUCTIVE_PARENT_SPAN_REASSESSMENT_2026-09-09.md).

Curve48's complementary-parent branch has now completed600 cumulative calls,
still M27, after three further independently replayed cached batches. The
[current cursor](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v5/suffix-queue.json)
is centre300; the canonical parent reassessment records all timings and bounds.
