# Chart coverage on the public28 control

Follow-up: the [fixed own27 point-search control](INVENTORY188_EXCEPTIONAL_DIRECTION_RECOVERY_2026-09-07.md) recovers the known28 direction on chart5. The representative misses below remain correct, but do not imply failure to recover their quotient direction.

The [current catalogue audit](CURRENT_CATALOGUE_AND_PUBLIC28_2026-09-07.md)
independently reproduces28 public points on our earlier local27 curve188,
the11952 fibre110314/102227. These are retrospective geometry diagnostics,
with no new point search or changes to completed candidate experiments.

## Fixed one-basis translations

Let P be certified extra public point26 (zero-based) and B_0,...,B_26 the
previous local basis. The frozen roster is P and P±B_i:55 distinct points,
each outside the old subgroup. Both signs are located in all49 original
completed125000-height charts. All5,390 observations are outside. None
improves the unshifted minimum21,632,242,813,382.

[The exact audit and replay](../../artifacts/generated-results/elliptic-curves/inventory188_public28_translates_v1.json)
and [independent Sage check](../../artifacts/generated-results/elliptic-curves/inventory188_public28_translates_sage_v1.json)
pass. Sage independently checks every addition and inverse rational coordinate
map. Its first wrapper passed a point dictionary to a constructor and failed
before producing evidence; the source and log remain preserved. The V2
wrapper explicitly reads x and y. Longer words and other representatives
remain untested.

## Geometry from the existing27-point subgroup

A separate geometry-only control uses the existing specialized-parity policy
on the old27-point subgroup, without public-point inputs. It retains the
sample domain `full11952-specialized-followup-v1`:2,048 distinct27-bit masks
with nonzero quotient above the17-point generic prefix. Canonical heights
at384 bits, rounding at10^6, unimodular LLL and CVP determine the sample.
The49 largest computed norms are selected, with parity ties. All maps are
frozen before the public witness is read. This is not a complete enumeration
of2^27 masks or a CVP-optimality theorem.

The exact audit checks the positive rounded metric, unimodular change,
all2,048 parity/norm transports, the fixed selection, every rational centre
and every quartic identity. Geometry construction takes2.7852 seconds; its
exact map/witness audit takes1.7087 seconds. There is no point enumeration.
[The map and visibility certificate](../../artifacts/generated-results/elliptic-curves/inventory188_own27_geometry_visibility_v1.json)
passes, but its98 signed observations of public point26 are all outside125000.
The minimum height is104,137,837,718,785,499,398.

## Both certified extra representatives

The retained finite signatures separately certify rank28 for the old27 basis
plus public point26 and plus public point27. They are the complete two public
points individually detected as extra by those modulo2 signatures. The
rank27 finite results for the other26 public points do not classify those
points as lying in the old rational span.

[The two-witness comparison](../../artifacts/generated-results/elliptic-curves/inventory188_two_witness_chart_comparison_v1.json)
rechecks both28-point certificates and all392 coordinates: two points, both
signs,49 original generic17 charts and49 own27 charts.

| Public point | Best height in generic17 charts | Best height in own27 charts |
|---|---:|---:|
| 26 | 21,632,242,813,382 | 104,137,837,718,785,499,398 |
| 27 | 173,970,974,235,857 | 2,186,776,913,276,331 |

Neither chart set contains either representative at125000. Own27 gives worse
coordinates for both. This is one curve and two witnesses, not a general
policy or sensitivity theorem. Other representatives may be cheaper. No
point-search budget is increased and no further search is launched.

The remaining question is whether a different centre selection or a better
representative in the same quotient supplies affordable coverage. More known
points and larger computed parity norms alone do not establish that coverage.
Any proposed policy change needs its own fixed test; public control points
stay outside prospective selection inputs.

## Reproduction

```sh
python3 elliptic-curves/cas/audit_inventory188_public28_translates.py check
python3 elliptic-curves/cas/audit_inventory188_two_witness_charts.py --check
```

Independent Sage sources are `verify_inventory188_public28_translates_v2.sage`
and `audit_inventory188_own27_geometry.sage`. Use an isolated copy rather than
overwrite frozen outputs. Local protocols and supervision records are under
`artifacts/local/elliptic-curves/inventory188-public28-translates-v1/` and
`inventory188-own27-geometry-control-v1/`.
