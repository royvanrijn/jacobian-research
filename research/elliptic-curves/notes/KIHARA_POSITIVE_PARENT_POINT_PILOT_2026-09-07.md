# Three discovered directions in the positive-parent fibre pilot

All **294 declared point-search boxes** completed on six fixed fibres of
the [three new positive-ratio parents](KIHARA_THREE_POSITIVE_PARENT_EXPANSION_2026-09-07.md).
Independent complete finite-group replays at moduli3 and5 certify three
additional directions on two curves. The strongest lower bound is14.
All six equations are distinct and unmatched in the pinned620-curve
catalogue and201-curve local inventory. This is not a literature-wide
novelty certificate or a near-record result.

Authority: `EC-KIHARA-POSITIVE-PARENT-FIBRE-PILOT-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).
The [report](../../artifacts/generated-results/elliptic-curves/kihara_positive_point_pilot_report_v1.json)
contains exact equations, complete exposure, source bindings and costs.

| Parent p/q | Base T | Exact initial displayed span | Final certified lower bound | Discovered directions | Completed boxes |
|---:|---:|---:|---:|---:|---:|
|1|1|6|7|1|49|
|1|1009/101|7|7|0|49|
|2|1|4|4|0|49|
|2|1009/101|9|9|0|49|
|3|1|12|14|2|49|
|3|1009/101|12|12|0|49|

The rank-at-least14 equation is

```
y^2 = x^3 - 23954333173179822722475*x
            + 1416739547341158145265761630526250.
```

Its [point-cloud certificate](../../artifacts/generated-results/elliptic-curves/kihara_positive_v3_unit_cloud_v1.json)
gives6298 verified rational points, explicit independent column indices
and finite-quotient witnesses. The other gain occurs on

```
y^2 = x^3 - 890938647507*x + 316307846233142894,
```

which has the rational2-torsion point(476338,0). Its lower bound improves
from6 to7. Model transformations remove exact fourth/sixth powers at
primes through997; global minimality is not claimed.

## Intake exposed two restrictions before searching

The six parent/base pairs were fixed before arithmetic: p/q=1,2,3 and
T=1,1009/101. The latter has coordinate height1009; heights depend on the
chosen base coordinate. There was no score, public target, validation-prime
selection or replacement.

The generic section ranks7,9,12 did not all survive specialization. Both
fibres on the first parent have rational2-torsion, so the adaptive worker's
required no-rational2-torsion witness cannot hold. Separately, the mod2
rank of the second parent's unit-fibre seed was only4 instead of9.
A deficient finite quotient alone would establish neither rank loss nor
nonsaturation.

Complete bounded mod3/5 diagnostics gave ranks6,7,4 on the three failed
intakes. To close the two smaller spans, numerical canonical heights
proposed rational coefficient words with denominator at most128. Exact
elliptic group equations then verified every original point as a rational
combination of the selected independent basis, modulo explicit2-torsion.
Thus the two specialization losses7to6 and9to4 are proved. They are not
just weak finite-image lower bounds. An initial vector-indexing error and
its corrected source are preserved.

The [portable seed bundle](../../artifacts/generated-results/elliptic-curves/kihara_positive_fibre_seed_bundle_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/kihara_positive_fibre_seed_replay_v1.json)
check the model scaling, original section specializations, complete finite
quotient groups, torsion exclusion at the chosen modulus and every exact
span relation. The starting ranks are therefore6,7,4,9,12,12. There is no
whole-curve upper-bound claim.

## Torsion-compatible exposure

The existing pointed-quartic backend accepts a rational curve and an
explicit centre without an independent adaptive MWState basis. This
pilot uses that interface with an empty valid rank state. The centres
are independently certified combinations of the known seed basis; an
empty state is not a claim that the known group has rank zero.
No no-two-torsion validator or production admission rule was weakened.

Each fibre receives49 fixed pointed charts, height125000 and at most ten
seconds per chart. All six map files precede the first point attempt.
Canonical heights at384-bit precision, rounding at10^6, unimodular LLL
and numerical CVP propose representatives. Exact parity and rounded-norm
identities, rational group sums and quartic transformations are checked.
No closest-vector optimality or complete parity-coverage claim follows.

The mask domain is all nonzero masks when fewer than2048 exist, otherwise
2048 fixed SHA256 masks. For the rank4 seed there are only15 nonzero
parities. Its remaining34 centres are deterministic distinct representatives
in those same classes, obtained by signed twice-basis translations and
deduplicated up to sign. They are not extra parity classes or rank
directions. This policy was frozen after the arithmetic intake and before
any point search, with the same49-box budget on all six curves.

Every known input point and every returned point is retained. Rank checks
run after the fixed exposure, using moduli3 and5 with independently checked
absence of rational torsion at the proof modulus. The exact initial-span
upper bounds prevent an omitted known dependent point from being counted
as a new direction.

## Independent checks and costs

The new independent checks cover all294 metric/parity/centre choices,
the exact backend transformations and raw square witnesses, and all
**12402 retained points** across the six clouds. Complete finite groups,
quotients and independent columns are reconstructed in a portable verifier.
The two proof moduli agree on every final lower bound.

The [cloud bundle](../../artifacts/generated-results/elliptic-curves/kihara_positive_point_cloud_bundle_v1.json)
and [standalone result](../../artifacts/generated-results/elliptic-curves/kihara_positive_point_cloud_replay_v1.json)
provide the replay inputs:

```sh
sage -python verify_kihara_positive_fibre_seeds.sage \
  --input kihara_positive_fibre_seed_bundle_v1.json --output seed-replay.json
sage -python verify_kihara_positive_point_clouds.sage \
  --input kihara_positive_point_cloud_bundle_v1.json --output cloud-replay.json
python3 elliptic-curves/cas/report_kihara_positive_point_pilot.py --check
```

Total new supervised cost is **218.001065285 seconds**, including intake,
the failed span-proposal attempt, its correction, map preparation, all
point workers and history replays, cloud audits and independent proofs.
Previous parent-construction costs are not counted again. All workers
used a2GiB limit and ran sequentially; per-stage wall limits and the
ten-second per-box limits are retained in the protocols.

The result measures a small fixed portfolio with equal completed point
exposure. It does not establish parent superiority, a score policy,
specialized-rank upper bounds or discovery probabilities. No curve enters
the high-rank inventory, and no larger population or automatic next wave
was launched.
