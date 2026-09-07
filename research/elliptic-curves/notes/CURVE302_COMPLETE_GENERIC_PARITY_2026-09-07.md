# Complete generic17 centre-class geometry on302

The [exact audit](../../artifacts/generated-results/elliptic-curves/curve302_complete_generic_parity_geometry_v1.json)
checks all131071 nonzero parities of the original17-section subgroup on302.
It produces49 factor-free charts,48 absent from the previous2048-class sample.
This is a geometry result; no point search or additional rank gain is claimed.
The active [48-fibre campaign](DET1092_RECORD_SCALE_CAMPAIGN_2026-09-07.md)
keeps its frozen sampled policy.

The [completed recovery calibration](CURVE302_RECOVERED_SUBGROUP_CALIBRATION_2026-09-07.md)
uses49 centres selected from2048 sampled parities, only about1.56percent of
the nonzero generic17 classes. Its success recovers seven exceptional
directions through subsequent subgroup enlargement. The present audit tests
the earlier centre-selection stage without another parameter population.

Before execution, the protocol fixes every nonzero mask1..131071, the same
384-bit rounded height metric, deterministic numerical CVP representative,
49 largest computed norms with smaller-mask tie-breaking, and the existing
factor-free map. It accepts only the original generic17 seed as an external
artifact input. The artifact-access guard precedes arithmetic imports and its
recorded reads pass inspection. Neither public exceptional points nor the
previous sample enter geometry selection. The previous sample enters only
the subsequent comparison, after the49 new maps have been frozen.

Every parity representative, rounded norm and integral basis transport is
checked exactly, followed by49 exact rational centre sums and expanded binary
quartic substitution identities with square scaling. The131071 parities occur
once each in immutable chunks. All2048 old sample entries reproduce exactly,
so the difference comes from completing the candidate classes rather than
changing the metric or CVP procedure.

| Selected centres | Minimum computed norm | Maximum computed norm |
|---|---:|---:|
|49 from2048 sampled classes|149668095|169049461|
|49 from131071 classes|161236614|170482149|

Only one chosen centre overlaps. These are rounded-form integer norms, not
certified covering radii. Numerical CVP choices remain unproved optimal;
complete parity enumeration is only in the supplied subgroup, not the unknown
full Mordell–Weil group. Higher computed centre norms do not prove lower point
height, improved recovery or a new independent direction. A subsequent bounded
recovery experiment must measure those endpoints before prospective adoption.

Geometry completes in22.800151072 supervised seconds under a180-second/2GiB
limit. The separate exact verification takes2.282146919seconds under a120-second
limit. All intermediate chunks, maps, access logs and supervision remain under
`artifacts/local/elliptic-curves/curve302-complete-generic-parity-v1/`.
Sources are `curve302_complete_generic_parity_geometry.sage` and
`verify_curve302_complete_generic_parity.py`. There is no automatic enlargement
of the active point campaign or repeat of a failed/censored stage.
