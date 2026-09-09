# Curve52: certified25→26 in the100-invocation first pass

The native R17 family08f72 fibre at `-164/2855` now has26 independently
certified rational points. Independent replay passes. The new direction was
found on chart73 using quartic-minimized coordinates, generic norm12 class73128
and extension95. On the discovery equation its coordinates are

```
x = 163573931100909861130269627551261840905/355788620652
y = 30351911861295762967751450335528774930611721054492281600/5105240566787269
```

The standalone certificate uses exact finite-group signatures and modular
exclusion of rational2-torsion to prove independence by infinite descent.
This is a lower bound, not an exact rank or a conductor record. The earlier
root-number+1 diagnostic only scheduled this candidate; it did not certify
the new point or the parity of the Mordell–Weil rank.

## Frozen preparation and search

Preparation uses `next24_r17_results_v1.json`, entry `08f72-053`, verifies the
native generic17 and existing M25, and independently replays the generic CVPs
for all seven historical productive masks:
19918,57793,64272,73128,80170,94624,94859. It takes11.247 seconds under a
120-second/1-GiB bound. Neither an external higher-rank point nor a withheld
visibility oracle enters selection.

The new short-pass wrapper retains the parent runner's V3 landscape, both
coordinate maps, cached admission, exact gain proof and immediate rebuild.
The only search-budget change is100 total invocations across all epochs.
Height125000 and10 seconds per invocation remain fixed. Search and replay each
have1800 seconds/3GiB and one worker. Existing runners are unchanged.

The initial landscape selects320 centres. After73 invocations the new direction
triggers a rebuild, which selects325 centres. The remaining27 invocations find
no further certified gain. All100 invocations complete without timeout, and
both retained clouds have finite ranks26 modulo2,3,5. The terminal reason is
`CHART_BUDGET_EXHAUSTED`: the second landscape was not exhausted.

Search takes110.055 seconds, peak RSS339103744 bytes; independent replay takes
117.552 seconds, peak314060800 bytes. Replay regenerates both full landscapes
with the original rational CVP solver and verifies the executed prefix of
maps, point witnesses, gain provenance, and whole-cloud finite certificates.
This provides an actual bounded short-pass result, not a controlled speedup
comparison or a general success rate. Full-landscape replay is now a material
cost relative to this shorter search.

## Evidence and continuation

- [Sealed result, full basis and independent certificates](../../artifacts/generated-results/elliptic-curves/curve52_short_v3_v1/result.json).
- [Short-pass wrapper and replay entry point](../cas/run_short_seed_v3.py).
- [Unchanged parent runner](../cas/run_parent_seed_v3.py).
- [Native R17 preparation](../cas/prepare_r17_productive_seed.py).
- Raw preparation: `artifacts/local/elliptic-curves/curve52-v3-preparation-v1/`.
- Raw search/replay: `artifacts/local/elliptic-curves/curve52-short-v3-discovery-v1/`,
  including frozen sources and a hash-bound `continuation-queue.json`.

The continuation queue records the verified M26 basis, the current landscape,
and the next coordinate policy to consider, including a partially visited
centre. It does not launch or certify the remaining suffix. Any deeper run
must preserve this exposure rather than repeating or silently dropping it.
The reusable [queue recorder](../cas/queue_verified_short_pass.py) reproduces
the same next step in `continuation-queue-v2.json`: centre13, factor-free map.
It also distinguishes a gain at the cutoff, which requires a basis rebuild,
from continuation of an unchanged basis. The original queue record is retained.
The broader rank29–32 and conductor-record objectives remain open.
