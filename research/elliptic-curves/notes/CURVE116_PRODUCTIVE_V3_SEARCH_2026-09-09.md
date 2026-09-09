# Curve116: certified lower bound25 to26, independent replay complete

The R17 fibre `103b2` at `-1588/3431`, inventory `new-20260906-116`,
now has26 certified independent rational directions. Both adaptive epochs
pass independent replay. No exact rank or conductor record is claimed.

Preparation uses entry141 (`103b2-014`) of `compact192_r17_results_v1.json`,
reconstructs the native seventeen with exact scale36, and rechecks M25.
All historical gains are rechecked, retaining seven generic masks:
10721,11241,75261,82549,90789,93313,95179. Each generic minimum is computed
with integer CVP and independently with the original rational solver.
The reusable [R17 preparation tool](../cas/prepare_r17_productive_seed.py)
accepts a source certificate and curve ID rather than requiring another copy.

The initial landscape scores1792 extension classes and selects332 centres.
Chart27 (index26), using quartic-minimized coordinates, supplies:

```
x = -14743825871171457660859570975/2028
y = -98658854846281797735707223655959199600175/2197
```

These coordinates use the discovery equation in the evidence packet. Exact
finite signatures certify26 independent columns; the modular exclusion of
rational2-torsion completes the infinite-descent independence argument.
The winning centre has generic mask75261, exact generic norm12, extension251.
An8/10-only parent bank would omit this successful generic class.

The search rebuilds immediately around M26, selects326 centres and completes
652 additional invocations without another gain. All679 invocations complete
without timeout. Both complete clouds have finite ranks26 modulo2,3,5.
Search takes482.586 seconds, peak RSS357109760 bytes; independent replay takes
119.215 seconds, peak331726848 bytes. Replay regenerates both landscapes using
the original rational solver, both maps and the invocation schedule, exact
point witnesses, gain provenance, and whole-cloud certificates.

The frozen limits remain height125000,10 seconds per invocation,4096 actual
invocations, target32, at most eight epochs, one worker, and separate
7200-second/3-GiB search and replay supervisors. The controller runs replay
and export only after successful search completion.

A post-discovery exact coordinate audit puts this point at height37261
quartic-minimized and74522 factor-free at the winning centre. Both are within
125000; the negative point is outside under both maps. This example does not
require both policies for that point, while curve113 supplies a concrete case
where factor-free coordinates alone place the gain within the chosen bound.

## Evidence and continuation

- [Sealed result, full basis, certificates and completion records](../../artifacts/generated-results/elliptic-curves/curve116_productive_v3_v1/result.json).
- [Bounded runner and independent replay](../cas/run_productive_seed_v3_v2.py).
- [Coordinate audit](../cas/audit_productive_gain_coordinates.py).
- Preparation: `artifacts/local/elliptic-curves/curve116-v3-preparation-v1/`.
- Search/replay: `artifacts/local/elliptic-curves/curve116-productive-v3-discovery-v1/`,
  including frozen sources and `winning-coordinate-audit.json`, plus sibling
  supervision directories.

The selected M25 sequence has now produced three verified25→26 transfers:
[curve200](CURVE200_PRODUCTIVE_V3_SEARCH_2026-09-09.md),
[curve113](CURVE113_PRODUCTIVE_V3_SEARCH_2026-09-09.md), andcurve116. The
root-number diagnostic only scheduled them; exact points and independence
proofs establish the gains. These selected outcomes are not a general success
rate or rank-parity theorem.

The next bounded search starts from the already-known native M28 on inventory188
/ ICARM619. Its17→27 history and retained27→28 recovery supply nine productive
classes, prepared by [the native adapter](../cas/prepare_curve188_productive.py).
The old28th direction is seed input, not discovery; only a certified extension
beyond M28 meets that continuation's success threshold. This note reports no
outcome of that separate search. Rank29–32 and conductor-record goals remain open.

## Fixed diagnostic of unadmitted cloud points

Two separately fixed samples test whether finite-span UNKNOWNs conceal another
direction: the eight smallest coordinate-size scores, then eight evenly spaced
order statistics including the extremes. Both select from sign-normalized
nonbasis points in the completed cloud, using maximum absolute numerator times
maximum denominator as the size score. There are15 distinct sampled points.

All16 sample decisions initially have `UNKNOWN_FINITE_COLUMN_IN_SPAN`. The
exact halving-or-cycle classifier proves dependence on M26 in2–4 steps for
every selected point. The small-point diagnostic takes2.148 seconds; the
size-stratified diagnostic takes8.629 seconds, including setup. Both use a
four-step cap and a300-second/1-GiB supervisor. A separate Fraction-arithmetic
replay verifies every final odd-multiple relation without Sage or the halving
classifier. No hidden direction is found in these samples. This does not
classify the whole cloud or bound the curve's rank, and does not justify
silently treating other finite-span UNKNOWNs as dependent.

Scripts: [small sample](../cas/audit_productive_cloud_halving.py),
[stratified sample](../cas/audit_productive_cloud_halving_strata.py),
[independent relation replay](../cas/verify_productive_cloud_relations.py).
Raw protocols, frames, decisions and `relations-verified.json` are retained
under `artifacts/local/elliptic-curves/curve116-cloud-halving-diagnostic-v1/`
and `curve116-cloud-halving-strata-v1/`, with sibling supervision records.
