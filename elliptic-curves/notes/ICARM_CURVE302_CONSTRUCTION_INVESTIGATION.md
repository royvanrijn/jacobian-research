# ICARM curve 302: construction and H3/R17 provenance investigation

Status: **exactly excluded from the published R17 chart; broader public
provenance incomplete**.

## Bottom line

Curve 302 and 31 independent points are exactly reproducible; see
[`ICARM_CURVE302_RANK31.md`](ICARM_CURVE302_RANK31.md). The public material
currently inspected gives attribution and a conditional exact-rank statement,
but no family equation, search parameter, specialization map, or decomposition
of the 31 points into generic and exceptional sections.

The published Elkies--Klagsbrun rank-17 equation is now available locally and
has been tested exactly.  The primitive degree-24 equation
`j_R17(t)=j_302` has an irreducible degree-24 reduction modulo `397`, hence is
irreducible over `QQ` and has no rational root.  Curve 302 is therefore not a
direct rational specialization of that published chart.  This does not rule
out another K3 fibration, an isogenous construction, or a different family.
See
[`ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md`](ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md).

## Exact recognition fingerprint

For exact family matching, use the global minimal invariants from the rank
certificate. The reduced `j`-invariant is

```text
375212263011874190418465904591842149883143397000379116040077721209939307800933802063374036449207052554068812677819360841029222896769044613030052189320099029335648176287279769780720048237450063678337025
/
91178667460631761509802129711614912708907877152843688131591561513957139867216918982080782638846592571825918112525868294328087917747020295905203091292094298537560118649108877196489360560362979328
```

with SHA-256

```text
5939208330113d89ae063d62053f0c8383e18b3a564919b86f86a02a4d13a550.
```

This exact equality test has now been performed against the published `R17`
`j`-map.  There is no rational match.  For any other proposed family, a
positive match would still need a Q-isomorphism and section transport;
equality of `j` alone is insufficient because of twists.

The bad-reduction fingerprint is

```text
I15@2, I4@3, IV@5, I6@7, I4@11, I5@13,
I2@19, I2@23, I3@29, I2@37, I2@41, I2@73,
I2@131, I2@167,
I1@7547, I1@632881, I1@966509,
I1@18145679437533309132469,
I1@767028866604834801397681553,
I1@30580600452196904409276223329355584892025407195996968868775951126238056443210297.
```

This profile is useful for rejecting proposed specializations whose minimal
invariants cannot have the same local valuations.

## Relation to the reconstructed H3 route

The reconstructed H3 programme starts from the level-474 `E7+E8/MW2` Kumar
source and reaches the recovered rootless `MW17` lattice.  Independently, the
published compact rootless equation and all seventeen generic sections are now
available locally.  Its exact `j`-map supplies the exclusion above: curve 302
does not occur at a rational parameter in this published chart.

Had curve 302 been a nonsingular specialization of that generic rank-17
family, its 31-point subgroup would have contained fourteen directions beyond
the generic rank:

```text
31 - 17 = 14 exceptional specialization directions.
```

The exact `j`-exclusion shows that no such `17+14` decomposition exists through
the published rational parameter chart.

## Exact specialization exclusion

The replay constructs

```text
c4_R17(t)^3*Delta_302-c4_302^3*Delta_R17(t).
```

After primitive normalization it has degree 24, retains degree 24 modulo 397,
and is irreducible over that finite field.  Gauss's lemma proves irreducibility
over `QQ`; the nonzero leading coefficient also excludes the point at
infinity.  Thus there is no rational specialization parameter to which an
isomorphism or section transport could be attached.

## Best next calculations

The current high-leverage route is not a blind 32nd-point search on curve 302:
the public BSD+GRH calculation already predicts exact rank 31. Instead:

- obtain the discoverers' construction record and test any different proposed
  family by the same exact invariant gate;
- keep curve 302 as an external regression for finite-quotient independence,
  height-lattice diagnostics and residual descent;
- pursue rank 32 in the published `R17` family through its certified rational
  parameters and quotient-targeted specialization pipeline, not through a
  nonexistent curve-302 parameter.

The existing bounded six-root Mestre census has not been rerun for curve 302 in
this update. No negative Mestre-family conclusion is claimed.
