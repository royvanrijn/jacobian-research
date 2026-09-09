# Curve113: certified lower bound25 to26, independent replay complete

The R17 fibre `08f72` at `-3292/1853`, inventory `new-20260906-113`,
has a new26th independent rational direction. Both adaptive epochs replay
independently. This is a subgroup lower bound, not exact rank or a conductor record.

## Native seed and successful continuation

Preparation starts from `compact192_r17_results_v1.json`, entry96,
`08f72-001`. Exact specialization through the compact six-family atlas and
scale6 reconstructs the native seventeen and verifies the full M25 certificate.
Rechecking all49 historical chart clouds supplies eight productive masks:
17296,18290,19918,57199,60048,81501,101054,107299. Their exact generic minima
are recomputed with integer CVP and independently with the original rational
solver. No full-shell completeness claim is needed for this explicit subset.

The initial landscape scores2048 extension classes and selects377 centres.
The12th chart, index11, supplies the new point using factor-free coordinates:

```
x = -1894014613201859342541817847894344538788871/6337139633045868
y = 69689765089548413761037004541532068865484991424167812874617600/12135789803388143515013
```

These are coordinates on the discovery equation in the evidence packet.
Exact finite signatures certify26 independent columns; modular exclusion of
rational2-torsion validates the infinite-descent independence argument.
The winning centre has generic mask57199, exact generic norm12, extension205.
An8/10-only generic bank would omit this productive class.

Immediate rebuilding around M26 selects375 centres and completes750 further
chart invocations without another gain. Both whole-cloud audits have finite
ranks26 modulo2,3,5. All762 invocations complete without timeout. Search takes
517.262 seconds with peak RSS363544576 bytes; independent replay takes
141.356 seconds with peak336166912 bytes. Replay regenerates both full
landscapes with the original rational CVP solver, both maps and the exact
invocation schedule, point witnesses, gain provenance and cloud certificates.

## What this changes for future searches

The [runner v2](../cas/run_productive_seed_v3_v2.py) reads the generic dimension
from the prepared seed and supports both MW16 and R17. The frozen successful
curve200 runner stays unchanged. Bounds remain height125000,10 seconds per
chart,4096 invocations, target32, at most eight epochs, one worker and separate
7200-second/3-GiB search and replay supervisors.

A post-discovery exact coordinate audit puts the new point at height82957
factor-free and175691 quartic-minimized at the winning centre. Only the first
lies within125000. The negative point lies far outside the bound under both
maps. Thus removing factor-free searches would lose this particular point
at this centre and bound; the audit makes no claim about other points or centres.
For comparison, curve200's gain is visible under both maps. Keep both policies.

## Evidence

- [Sealed result, full basis and certificates](../../artifacts/generated-results/elliptic-curves/curve113_productive_v3_v1/result.json).
- [Native preparation](../cas/prepare_curve113_productive.py).
- Raw preparation: `artifacts/local/elliptic-curves/curve113-v3-preparation-v1/`.
- Search and replay: `artifacts/local/elliptic-curves/curve113-productive-v3-discovery-v1/`,
  including frozen sources and `winning-coordinate-audit.json`, and sibling
  search/replay supervision directories.
- [Coordinate audit](../cas/audit_productive_gain_coordinates.py).

The root-number diagnostic scheduled this candidate; it did not prove the gain.
The completed trials now have two transfers25→26 (curve200 andcurve113), while
curve90 stayed27 andcurve92 stayed26. These selected outcomes are not a success
rate estimate. Rank29–32 and conductor-record objectives remain open.
