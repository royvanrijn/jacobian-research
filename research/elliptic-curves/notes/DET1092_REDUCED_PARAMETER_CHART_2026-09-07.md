# A smaller arithmetic chart for the determinant1092 MW17 parent

The [certified determinant1092 parent](CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md)
admits an exact rational parameter change whose integral short equation has
largest coefficient sizes **90 bits in A and 135 bits in B**. All seventeen
generic sections transport exactly. This is the same elliptic fibration on the
same K3 surface, not another parent or an increase of generic rank.

The portable inputs, equation, sections, and independent proof are in
[`det1092_reduced_parameter_chart_v1`](../../artifacts/generated-results/elliptic-curves/det1092_reduced_parameter_chart_v1/manifest.json).
The standalone Sage verifier proves the function-field identities and all
section transports without project arithmetic imports. A replay from copied
inputs produces byte-identical equation/section and proof files.

## Why change the parameter chart?

The [previous invariant-height audit](BOUNDED_PRIME_EXPOSURE_AND_HEIGHT_FLOOR_2026-09-07.md)
proved that the raw fibres at T=1 and 1009/101 require completed-square integral
quartics of at least 255 and 438 coefficient bits. More effort reducing those
same fibres cannot remove that invariant obstruction. A fixed projective
height-two diagnostic, including the excluded anchor only as a diagnostic,
revealed unusually large removable arithmetic scalings. This motivated an
equation-only parameter change, before selecting any further fibres.

Write the original short equation as y²=x³+A(T)x+B(T), of degrees 8 and 12.
First clear coefficient denominators with an exact Weierstrass scaling.
For the resulting binary forms F and G, factor the common coefficient content
of the two Hessians and their cross-Jacobian:

```
56 F F'' - 49(F')²
132 G G'' - 121(G')²
12 F'G - 8 FG'
```

If both forms reduce to powers of a common linear form, these covariants
vanish. Their content therefore supplies necessary prime support for this
complete-collapse condition. It is not a proof that all partial-collapse
routes have been considered.

At this frozen support, retain only index-p parameter substitutions followed
by a removable Weierstrass scaling p^k with k≥2. Nineteen strict contractions
were found. Each reduces the homogeneous resultant by p^(96(k−1)), with the
full rational matrix/scaling identity checked exactly. Then compare seventeen
weighted two-dimensional LLL bases, with exponents −8 through 8, by the exact
coefficient score max(max|Aᵢ|³,max|Bᵢ|²). Its selected value has 270 bits.
This finite rule does not prove global minimality, an optimal coefficient
height, or completeness of routes with equal-resultant intermediate moves.

The selected matrix and scaling are

```
a = -36201755177827284938903000
b = -41418390914593414769282000
c = -1016449687002447546724417437
d = 1324257248229139113399816775
u = 127912441776881731337678240068563891280789386291887300898201991681/6
```

Set T=(as+b)/(cs+d), h=cs+d. The new coefficients are h⁸A(T)/u⁴
and h¹²B(T)/u⁶. The old short coordinates are
x=u²X/h⁴ and y=u³Y/h⁶. The certificate also checks the original long-to-short
coordinate change, equality of j-invariants, and every transported section.
Two generic sections have rational-function coordinates; no polynomial-section
assumption is needed.

## Fixed point-exposure experiment

Retain the same two fixed-coordinate choices as the earlier pilot, now in the
new chart: s=1 and 1009/101. No arithmetic score, held-out validation prime,
public exceptional point, record rank, or catalogue match enters selection.
Both remain in the experiment; there are no outcome-based replacements.

| New coordinate | Original T | Reduced j numerator/denominator bits | Quartic coefficient floor |
|---|---|---:|---:|
| 1 | −5544296149458621407727500/21986254373335111905385667 | 217 / 194 | 32 |
| 1009/101 | 496473517765873968232324500/10876192098955201514531167069 | 486 / 466 | 77 |

These original parameter addresses are far outside the compact-height region,
while their invariant arithmetic heights are much smaller than the earlier
raw-chart pilot. This is a comparison of different fibres, not a violation of
the invariant floor for a fixed elliptic curve.

At s=1 the seventeen section values span exactly fifteen directions. With
zero-based section indices, exact Sage group arithmetic proves

```
P14 = P2 + P5 - P6 + P9 - P13
P15 = -P3 - P5 + P7 + P8 - P12 + 2 P13.
```

The subset P0,…,P13,P16 has an independent finite certificate. At the outer
fibre all seventeen section values are certified independent. These are
statements about the specialized section spans, not upper bounds on the
whole curves. Every discovered gain must therefore be measured from 15 and
17 respectively.

The new protocol retains the full-blind-regression factor-free mapper,
2048 deterministic nonzero parity samples, and 49 largest computed-norm
centres per fibre. All 98 maps freeze before points. Each chart receives
height 125000 and ten seconds, without rank stopping or an adaptive wave.
Point membership, histories, and complete-cloud independence are separate
exact checks. Catalogue comparison is post-search only.

The [completed exposure certificate](../../artifacts/generated-results/elliptic-curves/det1092_reduced_chart_point_pilot_v2.json)
records all 98 completed boxes, both exact histories, and twelve independent
certificate stages. The complete retained point clouds certify the same lower
bounds modulo 2, 3, and 5, with standalone Sage replay:

| Coordinate | Initial certified span | Final certified lower bound | Retained point witnesses | Actual quartic coefficient bits |
|---|---:|---:|---:|---:|
| 1 | 15 | 15 | 11,578 | 45–58 |
| 1009/101 | 17 | 17 | 535 | 87–97 |

Neither curve matches the pinned 626-equation ICARM snapshot or the 201-curve
local inventory. This post-search comparison does not establish literature-wide
novelty, and these bounds do not qualify as new near-record results. The high-rank
inventory is unchanged. The new chart returns substantial point witnesses where
the previous raw-chart fibres returned none, but adds **zero certified directions**.
Different fibres and different specialized section ranks prevent a causal
claim about detector efficiency from this two-fibre comparison.

The exposure and certification stages consumed **607.078143972 seconds**,
including the failed first intake. Parameter reduction, section verification,
and the preserved diagnostic attempts consumed another **12.683273788 seconds**;
the independent portable chart replay used **0.743798093 seconds**. The combined
supervised total is **620.505215853 seconds**. These are supervised stage wall
times, not a claim to include editing, unsupervised diagnostic commands, or
report-generation overhead. No larger parameter population or following point
wave was launched by this protocol.

## Preserved failures and replay

The first generic verifier wrongly required polynomial coordinates for every
section. The corrected version permits rational functions and verifies their
exact equations. The first reduced-chart intake wrongly required all seventeen
specializations to be independent; it stopped before any point search. The
second protocol preserves both fibres and uses the certified fifteen-point
unit subset, while recording and replaying all original section values. Its
centre-policy prose retains the legacy phrase “17-point”; the operative row
dimensions, parity moduli, seeds, and all checks use 15 and 17 respectively.
The frozen protocol is preserved, with this descriptive erratum.

The first finite-kernel diagnostic used Fraction(Sage QQ), which constructed
invalid numerator/denominator values in this runtime. Broad ValueError handling
then discarded every mod2 prime. Its empty mod2 report is invalid and is not
evidence. Version two converts through exact rational strings and requires
nonempty seventeen-column rows. The original mod3/mod5 diagnostic used serialized
strings and was unaffected. Signed-kernel relation proposals did not close the
span; bounded numerical-height proposals did, after exact group verification.
All attempts and their supervision costs are retained.

Replay the portable chart with

```
python3 elliptic-curves/cas/package_det1092_reduced_parameter_chart.py --check
python3 elliptic-curves/cas/report_det1092_reduced_chart_point_pilot.py --check
```

To re-execute its standalone arithmetic, copy the bundle to a fresh directory
and run `verify.sage` with its `proof-protocol.json`, `source-parent.json`, and
`chart-search.json`, choosing fresh output and certificate paths. The input
protocol binds the verifier and both input hashes.
