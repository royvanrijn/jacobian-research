# Curve200: productive V3 raises the certified lower bound from25 to26

The native MW16 fibre `a1-fibration-02` at `-32999/14074`, inventory
`new-20260906-200`, now has26 exactly certified independent rational points.
Both adaptive epochs pass independent replay. This is a new direction beyond
the retained M25, not an exact-rank proof or conductor-record claim.

## Construction and verification

The seed comes from entry57 of `strata60_mw16_results_v1.json`. Preparation
reconstructs the native sixteen through the compact atlas and exact scale156,
rechecks the M25 certificate, and rechecks every gain in the retained43-chart
attempt. All eight productive generic masks are kept:
15231,20702,22105,25820,29804,32970,33057,34257. Each generic minimum is
independently computed with the integer and original rational CVP solvers.
This is an explicit historical subset, not a complete generic shell bank.

The first epoch scores4096 extension classes and selects361 centres. On its
81st chart (index80), quartic-minimized coordinates supply the26th direction:

```
x = 2241867054295795317370962650905329172441/3992918783052
y = 1454116856212563924640018818514884646454473064646615894240/191939273158077719
```

These coordinates use the discovery equation recorded in the evidence packet,
not its global minimal model. The centre belongs to generic mask22105,
norm10, extension432. That generic class supplied two directions in the old
search. The successful backend invocation takes0.782 seconds. Exact finite
reductions certify independence of the entire26-point basis; the modular
2-torsion exclusion makes the infinite-descent independence argument valid.

The search rebuilds immediately around M26. Its367 selected centres yield734
further chart invocations with no additional certified direction. Both epochs'
complete point clouds have finite ranks26 modulo2,3,5. Overall815 invocations
complete, with zero timeouts. The search takes613.019 seconds, peak RSS
354455552 bytes; independent replay takes139.552 seconds, peak326660096 bytes.
Replay regenerates both full landscapes with the original rational CVP solver,
both coordinate maps and the invocation schedule, replays exact point witnesses,
and checks gain provenance and whole-cloud certificates.

## Reusable implementation and limits

[The prepared-seed runner](../cas/run_productive_seed_v3.py) accepts native MW16
seed preparations instead of requiring another curve-specific runner copy.
It retains cached finite admission, exact integer CVP, immutable chart receipts,
two map policies and immediate rebuilding. The frozen limits are height125000,
10 seconds per chart,4096 actual invocations, target32, at most eight epochs,
one worker, and separate7200-second/3-GiB search and replay limits.

An initial26.951-second launch was stopped after finding an inherited curve92
label in cloud metadata. Its actual equation and seed were curve200. Its frozen
sources and receipts remain under the `discovery-v1` directory with an explicit
abort record. The corrected `discovery-v2` protocol binds family and parameter;
only that run supports this result. Old90/92 runners remain unchanged.

## Evidence and scope

- [Sealed result, full basis, certificates and completion records](../../artifacts/generated-results/elliptic-curves/curve200_productive_v3_v1/result.json).
- [Native preparation](../cas/prepare_curve200_productive.py).
- [Reusable evidence exporter](../cas/export_productive_seed_result.py), which requires successful independent replay and is not itself another arithmetic replay.
- Raw inputs: `artifacts/local/elliptic-curves/curve200-v3-preparation-v1/`.
- Search/replay: `artifacts/local/elliptic-curves/curve200-productive-v3-discovery-v2/`
  and sibling supervision directories; includes the frozen source archive.

Selection used the existing M25 and historical productive classes, without
higher-point oracle input. The previously computed root number+1 was a scheduling
clue only; the new rational point and independence proof establish the increase.
The bounded M26 miss proves no upper bound. This is the first successful transfer
among the optimized [curve90](CURVE90_PRODUCTIVE_V3_SEARCH_2026-09-09.md),
[curve92](CURVE92_PRODUCTIVE_V3_SEARCH_2026-09-09.md), and curve200 trials.
The broader rank29–32 and conductor-record objective remains open.

## Retrospective map check

After completion, an exact four-coordinate audit of the winning centre finds
the new point at height26302 in quartic-minimized coordinates and35687 in
factor-free coordinates. Both lie below125000; the negative point has much
larger height under both policies. This success alone therefore does not show
that running both maps was necessary. The audit performs no point search and
supplies no prospective visibility labels. Its script is
[`audit_productive_gain_coordinates.py`](../cas/audit_productive_gain_coordinates.py),
with the source-bound `winning-coordinate-audit.json` in the completed run directory.
