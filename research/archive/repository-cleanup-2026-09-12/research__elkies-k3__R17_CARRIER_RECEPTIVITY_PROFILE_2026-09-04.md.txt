# Marked-`U` carrier-receptivity profile

Date: 2026-09-04  
Status: exact profile on 34 transported curve/chart cells; untransported coordinates remain `UNKNOWN`

## Purpose

The unit of this audit is a primitive ordered hyperbolic plane

\[
U=\langle F,O\rangle\subset \operatorname{NS}(X)
\]

in the pinned marking, not the abstract Mordell--Weil frame.  Its profile is

\[
\mathcal R(U)=(r_{\rm rigid},r_{8,\rm off},G_{\rm branch},
r_{\rm base},T/\text{minimum classes},C_{\rm equation,parameter}).
\]

Consequently, neither equality of frame classes nor rational `PGL2`
equivalence authorizes copying a quotient, carrier rank, incidence edge, or
base-rank result from one marked `U` to another.  The machine table contains
one record for each of the 43 marked charts and retains literal `UNKNOWN`
values where the marked transport has not been performed.

## Anchor quotient-labelled profile

The original twelve exact displayed quotients live on four anchor charts.
Ranks in the third column are ranks in the displayed exceptional quotient,
not ranks of the full specialized Mordell--Weil group.  Table labels abbreviate
the canonical `norm12-orbit-*` marked-`U` labels.

| marked `U` | curve: displayed quotient rank | rigid transfer rank | fitted norm-eight off-diagonal audit | low-genus base audit | marked minimum/Tate data | equation rank / support / input bits |
|---|---|---|---|---|---|---|
| `074d9` | `351:8, 356:12, 376:5, 377:6, 385:12` | `351:6, 356:2, 376:5, 377:5, 385:2` | 96 direction tests, all eight off-diagonal ranks `0` | complete rigid cross-fibre atlas | exact marked norm-eight/norm-ten counts only; Tate quotient `UNKNOWN` | `1 / 7 / 2544` |
| `11952` | `12:12, 395:11` | `12:0, 395:2` | 23 direction tests, both off-diagonal ranks `0` | exact bounded 1,143-class laboratory | complete marked minimum histogram; conditional zero-Tate-class exclusion on 17 product characters; quotient otherwise `UNKNOWN` | `11 / 9 / 3482` |
| `08f72` | `363:10, 364:11, 378:7` | `363:2, 364:1, 378:6` | 56 direction tests, all six off-diagonal ranks `0` | `UNKNOWN` | exact marked norm-eight/norm-ten counts only; Tate quotient `UNKNOWN` | `15 / 10 / 3468` |
| `103b2` | `393:9, 404:10` | `393:2, 404:1` | `UNKNOWN` | complete native rigid-cover atlas | complete marked minimum histogram; Tate quotient `UNKNOWN` | `30 / 14 / 4687` |

The exact specialization-parameter height ranges in these four charts are,
respectively, `93--104`, `98--99`, `122--125`, and `64--70` bits.  The JSON
stores each rational parameter itself, its numerator and denominator sizes,
and the full chart equation-complexity vector.

After the same-curve panel and the completed rank-28 tranche, the six
coordinates have marked data on `12/43` charts for complete rigid transfer
ranks, `3/43` for fitted off-diagonal norm-eight ranks, `12/43` for at least
one exact branch-character/incidence layer, `3/43` for low-genus base ranks,
`12/43` for at least one marked minimum-class row, and `43/43` for equation
complexity.  No chart has a complete Tate quotient.  The fifteen exact
quotient-labelled curves occur in 112 curve/chart parameter matches across
all 43 charts; 34 of those cells on fourteen charts now have saturated
quotient transports, and 31 cells on twelve charts have complete fixed-cover
span audits.

This already separates arithmetic families carried by the same abstract
frame.  Within alternate Q80, for example, the exact rigid rank vectors on
`11952` and `08f72` are `(0,2)` and `(2,1,6)`, while only `11952` currently
has the bounded low-genus base audit.  These are marked-chart statements;
unknown coordinates on either side are not negative results.

## Same-curve controlled panel

The balanced alternate-Q80 subpanel now transports five fixed public curves
through every chart in their rational-`PGL2` class: curves 12 and 395 on six
charts, and curves 363, 364, and 378 on four.  This adds nineteen non-native
transports, for 24 exact cells in total.  Every cell uses a saturated generic
MW17 basis and the complete 39,147-carrier inventory.

For each fixed curve, the generic subgroup is the same primitive sublattice
of the ordered displayed public group: every chart basis differs from the
native anchor by an integral determinant-`+/-1` matrix.  The complete split
count, quotient-span rank, and carrier-visible extension lattice are also
identical across presentations.  What varies is the deterministic carrier
cost order.  The first splitting carrier moves by factors 19.13 for curve 395
and 22.89 for curve 363; curve 12 remains the exact zero-split control on all
six charts.  This makes prefix visibility and first-hit position measured
search affordances, not intrinsic rank mechanisms.  See
[`R17_SAME_CURVE_MARKED_U_PANEL_2026-09-04.md`](R17_SAME_CURVE_MARKED_U_PANEL_2026-09-04.md).

## Branch-incidence graph

The stored graph is bipartite: carrier nodes on one side and exact
quotient-labelled specialization nodes on the other.

- The complete fixed inventories on all twelve targets give 70 rational
  rigid-cover incidence edges: 46 on `074d9`, 2 on `11952`, 18 on `08f72`,
  and 4 on `103b2`.  No rigid carrier in these inventories is incident to two
  current targets on the same marked chart.  In particular, the two covers
  at rank-29 curve 356 and the two at rank-29 curve 385 are disjoint, so this
  rigid layer supplies no common rank-at-least-18 corridor through those two
  records.
- The 75 fitted norm-eight quotient-basis directions give 75 diagonal edges.
  Evaluating every such source direction at every other exact target on the
  same marked chart gives 175 exact off-diagonal square tests and no split
  direction: 96 tests on `074d9`, 23 on `11952`, and 56 on `08f72`.  The
  `074d9` experiment selects one target-blind canonical norm-eight trace,
  then fits twelve members separately at curves 356 and 385 before freezing
  them.  Curve 12 is recorded as not applicable to this matrix because it is
  a `11952` fibre; comparing it requires a separate common-K3 birational
  transport, not substitution into the `074d9` base.
- The zero off-diagonal ranks apply only to these fitted carriers and current
  targets.  They are not an exhaustion of the norm-eight frontier and are not
  a held-out prediction score, because each pencil was fitted at its diagonal
  target.
- The inherited branch-character layer is exact on all ten alternate-Q80
  marked charts (1,198 covers summed over distinct base coordinates).  All
  ten now have complete 39,147-class smooth norm-ten layers, each with no
  equal-cover collision, internal three-character relation, or match to the
  twelve committed characters.  On `135b7`, `10f72`, `09952`, and `0ae21`,
  the old degree-one section span first had index two; an exact rational
  height-four half now saturates each marking.
- Marked `103b2` separately has its complete 39,120-class smooth norm-ten
  layer, with distinct characters, no internal three-character relation, and
  no formal committed-catalogue product match.

## Low-genus and Tate boundary

On marked `103b2`, all 39,120 rigid quadratic carriers have rational conic
bases and hence surface rank at least 18.  Their 765,167,640 distinct pairs
have genus-one `V4` base and surface rank at least 19.  Among the 5,566 pairs
with the immediate rational-point gate, the certified base-Jacobian rank
lower bounds have maximum 9.  A stored lower bound zero means that the
declared pass found no point; it is not an exact rank-zero result.

On marked `11952`, the exact bounded laboratory contains 1,143 distinct
rational conics, 652,653 genus-one pairs of surface rank at least 19, and
248,225,691 genus-five triples of surface rank at least 20.  Of the 64
selected base-rank screens, 62 completed and 17 have exact rank one.  This
does not classify the remaining native carrier classes.

The exact frame minimum-norm histograms are retained only as reference data.
They become marked minimum-class data on `103b2` and `11952`, and only the
norm-eight/norm-ten counts have marked evidence on `074d9` and `08f72`.  The
complete marked norm-ten count is exact on all ten alternate-Q80 charts;
their other minimum rows remain `UNKNOWN`.  Every full Tate quotient remains
`UNKNOWN`.  On `11952`, the separate product-character
reduction excludes the zero Tate class for height-eight sections under its
stated gates; it does not determine a quotient dimension, a nonzero class,
or existence.

## Next quotient-only transport tranche

The former rank-28 tranche is complete: curves 11, 391, and 423 each have an
exact displayed quotient `Z^11`.  Among the remaining 54 rows without a
displayed quotient transport, the maximum current rank lower bound is 27.
The next quotient-only tranche is:

| curve | marked `U` | parameter height | rank jump over generic 17 |
|---:|---|---:|---:|
| 67 | `08234` | 99 bits | at least 10 |
| 416 | `08234` | 104 bits | at least 10 |

The other 52 transports are deferred.  Selection uses only the maximum rank
lower bound; public family hit counts are deliberately absent because the
historical search exposure and its denominator are unknown.  The two selected
quotient fields remain `UNKNOWN` until their saturated
chart-specific generic bases and exact specialization relations have been
compiled.

## Artifacts and replay

- Full 43-chart profile: [`../artifacts/generated-results/elkies-k3-r17-carrier-receptivity-profiles-v1.json`](../artifacts/generated-results/elkies-k3-r17-carrier-receptivity-profiles-v1.json)
- Compact comparison table: [`../artifacts/generated-results/elkies-k3-r17-carrier-receptivity-profiles-v1.tsv`](../artifacts/generated-results/elkies-k3-r17-carrier-receptivity-profiles-v1.tsv)
- Same-curve exact panel: [`../artifacts/generated-results/elkies-k3-r17-same-curve-marked-u-panel-v1.json`](../artifacts/generated-results/elkies-k3-r17-same-curve-marked-u-panel-v1.json)
- Builder/checker: [`scripts/build_r17_carrier_receptivity_profiles.py`](scripts/build_r17_carrier_receptivity_profiles.py)

```bash
.venv/bin/python elkies-k3/scripts/build_r17_carrier_receptivity_profiles.py
.venv/bin/python elkies-k3/scripts/build_r17_carrier_receptivity_profiles.py --check
```

This profile is a middle-layer audit.  It changes no theorem status in
`MATH_STATUS.json` and makes no inference from bounded search counts.
