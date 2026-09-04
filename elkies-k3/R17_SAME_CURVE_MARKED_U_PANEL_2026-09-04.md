# Same-curve / different-marked-`U` carrier panel

Date: 2026-09-04

Status: exact 24-cell balanced subpanel; 19 new non-native transports

## Question and design

The public-curve atlas now has 112 exact curve/chart parameter matches across
all 43 marked charts for the fifteen curves with exact displayed quotient
labels.  Repeated appearances of one curve offer a within-curve control: the
elliptic curve over `QQ`, its ordered displayed independent-point subgroup,
and its rank lower bound are held fixed while the marked presentation and its
search coordinates change.

The exact balanced subpanel uses the two alternate-Q80 rational-`PGL2`
classes.  Every one of their ten charts has a saturated generic MW17 basis
and a complete 39,147-class smooth norm-ten carrier inventory.  Thus no chart
is compared using a shorter or different-size census.

| rational `PGL2` class | fixed curves | charts | cells |
|---|---|---:|---:|
| `11952` | 12, 395 | 6 | 12 |
| `08f72` | 363, 364, 378 | 4 | 12 |

For each of the 24 cells the replay:

1. specializes the saturated generic MW17 basis;
2. expresses it in the same ordered public-point coordinates and verifies all
   relations by exact elliptic-curve group law;
3. proves the displayed quotient by Smith form;
4. evaluates all 39,147 quadratic carriers over `QQ`;
5. expresses every splitting point in the fixed public group and computes its
   exact quotient span.

Canonical heights propose integral relations only.  They do not enter any
retained proof.

## Exact result

| curve | rank lower bound | presentations | displayed quotient | complete splits | exact split span | first-split priority range | spread |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 29 | 6 | 12 | 0 | 0 | none | -- |
| 395 | 28 | 6 | 11 | 2 | 2 | 1,513--28,940 | 19.13x |
| 363 | 27 | 4 | 10 | 3 | 2 | 850--19,453 | 22.89x |
| 364 | 28 | 4 | 11 | 2 | 1 | 5,739--30,681 | 5.35x |
| 378 | 24 | 4 | 7 | 13 | 6 | 788--4,486 | 5.69x |

The complete count and the quotient-span rank are constant for every fixed
curve.  The stronger lattice comparison is constant as well:

- the 17-dimensional generic bases on all presentations of a fixed curve are
  related to the native anchor by integral matrices of determinant `+/-1`;
- their primitive generic sublattices inside the displayed public group are
  therefore literally equal, not merely equal in rank;
- adjoining every splitting carrier point gives the same primitive visible
  extension lattice on every presentation of that curve.  Its ranks for
  curves 12, 395, 363, 364, and 378 are respectively 17, 19, 19, 18, and 23.

This corrects an important premise of the proposed experiment.  These ten
rational-`PGL2`-equivalent charts do change the marked presentation and basis
words, but after the exact `QQ` identification they do **not** change the
specialized generic subgroup or the complete carrier-visible extension
lattice.  They are therefore an especially clean test of presentation cost,
not ten independent intrinsic subgroup treatments.

## What changes: search exposure

`priority_rank` is the deterministic lexicographic order by generic-section
addition-chain upper bound, support, maximum coefficient, `L1` norm, and then
the coefficient vector.  It is a reproducible search-cost proxy, not a
measured runtime.

That proxy is strongly presentation-dependent even though the exact output
lattice is fixed:

- curve 363 has a splitting carrier in the first 1,000 rows only on `0ae21`;
  its first hit is 850 there and 19,453 on `135b7`;
- curve 395 has a splitting carrier in the first 5,000 rows only on `098fc`;
  its first hit ranges from 1,513 on `098fc` to 28,940 on `091e4`;
- curve 364 has no hit before 30,681 on the native `08f72` order, but hits at
  5,739 on `09952`;
- curve 378 has the same thirteen splitting carriers and six-dimensional
  quotient span everywhere, while the first 10,000 rows expose two of them on
  `08f72`/`09952` and four on `135b7`/`0ae21`.

Equation cost and carrier exposure are distinct affordances.  For example,
`08f72` has equation-rank 15 but does not expose curve 364 until priority
30,681, whereas equation-rank-38 `09952` exposes it at 5,739.  A globally
"cheap equation" need not give a cheap carrier basis for the target at hand.

## Operational consequence

For PGL2-equivalent marked presentations, the search should operate on the
union of chart prefixes and deduplicate exact carrier outputs in the fixed
public group.  A target-blind scheduler should minimize the declared carrier
cost across presentations, rather than selecting one canonical chart and
assuming its prefix is representative.  Once a complete transport is made,
the quotient span must be reported separately from the priority at which it
was found.

Within this panel, first-hit position, prefix hit count, carrier label, and
equation complexity are search affordances.  The quantities that survive the
within-curve presentation changes are the complete split count, the exact
quotient-span rank, the saturated generic subgroup, and the full visible
extension lattice.  Their stability is exact for these five curves; it is not
a theorem for the other 88 parameter matches.

## Claim boundary

The certificate proves 24 exact transports and complete carrier censuses.  It
does not prove that a displayed public subgroup is the full Mordell--Weil
group, does not give an exact rank upper bound, does not infer historical
search exposure, and does not generalize presentation invariance beyond the
ten complete alternate-Q80 inventories tested.

## Artifacts and replay

- Full certificate:
  [`../artifacts/generated-results/elkies-k3-r17-same-curve-marked-u-panel-v1.json`](../artifacts/generated-results/elkies-k3-r17-same-curve-marked-u-panel-v1.json)
- Compact 24-cell table:
  [`../artifacts/generated-results/elkies-k3-r17-same-curve-marked-u-panel-v1.tsv`](../artifacts/generated-results/elkies-k3-r17-same-curve-marked-u-panel-v1.tsv)
- Exact replay:
  [`scripts/certify_r17_same_curve_marked_u_panel.sage`](scripts/certify_r17_same_curve_marked_u_panel.sage)

```bash
PYTHONPATH=elliptic-curves/cas sage -python \
  elkies-k3/scripts/certify_r17_same_curve_marked_u_panel.sage --check
```

This is a carrier/search audit.  It changes no theorem status in
`MATH_STATUS.json`.
