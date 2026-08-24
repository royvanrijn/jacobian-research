# ICARM curve 302: construction and H3/R17 provenance investigation

Status: **public provenance incomplete; no K3 specialization claim**.

## Bottom line

Curve 302 and 31 independent points are exactly reproducible; see
[`ICARM_CURVE302_RANK31.md`](ICARM_CURVE302_RANK31.md). The public material
currently inspected gives attribution and a conditional exact-rank statement,
but no family equation, search parameter, specialization map, or decomposition
of the 31 points into generic and exceptional sections.

It is therefore not justified to identify curve 302 with the Elkies--Klagsbrun
rank-17 K3 family reconstructed under [`../../elkies-k3/`](../../elkies-k3/).
That family remains an obvious object to test because it produced the preceding
rank-record programme, but this is a research hypothesis, not provenance.

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

This is the first cheap exact equality test against any proposed one-parameter
specialization. A positive match still needs a Q-isomorphism and section
transport; equality of `j` alone is insufficient because of twists.

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

The current H3 programme starts from the level-474 `E7+E8/MW2` Kumar source
and has an exact degree-two neighbour corridor to the recovered rootless
`MW17` lattice. At repository head
`2ebb1f612bb04f74d25786e35c8a30eab5a7bedf`, the selected q24 horizontal
section over `QQ` is exact, but its resolved Riemann--Roch pencil, binary-quartic
Jacobian, and `D12/MW5` child fibre certificate are still open. Equation-level
reconstruction is therefore not complete through the full generic rootless
family. Consequently there is currently no exact function `j_R17(t)` or full
list of 17 generic section coordinates to evaluate at curve 302.

If curve 302 is a nonsingular specialization of that generic rank-17 family,
then its 31-point subgroup must contain fourteen directions beyond the generic
rank. This count is only conditional on the family identification:

```text
31 - 17 = 14 exceptional specialization directions.
```

No such `17+14` decomposition has been produced.

## Required exact specialization certificate

A genuine H3/R17 identification should contain all of the following.

1. An explicit rational point on the H3 base and every neighbour parameter
   needed to reach the rootless fibration.
2. A specialization of the generic rootless Weierstrass equation whose
   `c4,c6` are related to curve 302 by one rational scaling/twist-compatible
   change of variables.
3. An exact Q-isomorphism to the public global minimal model.
4. Exact transport of the seventeen generic sections.
5. A finite-reduction rank check separating the transported rank-17 subgroup
   from fourteen additional directions.

Until this exists, notes should say "candidate specialization" at most.

## Best next calculations

The current high-leverage route is not a blind 32nd-point search on curve 302:
the public BSD+GRH calculation already predicts exact rank 31. Instead:

- finish the selected H3 equation route far enough to expose the generic
  rootless `MW17` model and its `j`-map;
- solve the exact equation `j_R17(t)=j_302`, including points at infinity and
  every rational chart;
- for any rational solution, test the Q-isomorphism and all seventeen section
  transports before running exceptional-point searches;
- use the recovered parameter, if any, to search nearby rational
  specializations for rank 32 or a substantially smaller rank-31 model;
- obtain the discoverers' construction record, which could immediately rule
  the H3 hypothesis in or out.

The existing bounded six-root Mestre census has not been rerun for curve 302 in
this update. No negative Mestre-family conclusion is claimed.
