# Surface--fibration--character identity graph

Date: 2026-09-04

Status: design contract plus exact determinant-948 seed; no theorem-status
change

## Outcome

The construction state is no longer a flat sequence of gates.  It has three
layers:

```text
surface identity
  S=(X,T,NS,A_NS,O_Cliff,Gamma*,period point),

fibration presentation
  F=(U,W,R,MW,marking,k(P1_U) -> k(X)),

cover/character construction
  C=(C -> P1_U,G,chi,branch data,X_C,M_chi,Tate data).
```

Changing a primitive `U` on one fixed K3 changes the second line, not the
first.  A base extension creates a new surface above `X`; the surface-level
objects pull back to it but need not remain primitive, exhaustive, or even
isometric with their old intersection form.  This distinction is part of the
data model rather than a prose convention.

The machine-readable contract is
[`data/object-graph/surface-fibration-character-graph-v1.schema.json`](data/object-graph/surface-fibration-character-graph-v1.schema.json).
The first seed is
[`data/object-graph/determinant-948-surface-graph-v1.json`](data/object-graph/determinant-948-surface-graph-v1.json).
It registers published R17, alternate Q80, hidden `103b2`, and the noncyclic
`4A1/MW13` fibration under one surface identity.  It also deliberately gives
the `103b2` fibration and the `103b2` quadratic character different IDs.

## What is persistent

For two marked Jacobian fibrations `U_1,U_2` on the same fixed marked K3
surface, the following records are inherited through a `HOP_TO` edge:

| layer | persistent data | presentation-dependent data |
|---|---|---|
| K3 surface | `X`, `T`, `NS`, `A_NS = A_T(-1)`, Clifford order, stable marking group, marked period point | none |
| elliptic fibration | common ambient `NS` and its Galois action | `U`, `W=U^perp(-1)`, fibre roots `R`, `MW=W/R`, zero, base coordinate, equation |
| search | exact divisor/Galois objects already identified in a common marked lattice | equation size, section words, carrier prefix order, first-hit cost |

This is consistent with the exact same-curve panel: after transport, its
complete visible extension lattice is fixed while first-hit priority changes
by as much as `22.89x`.  Presentation is therefore a search affordance, not a
new surface identity.

## The character correction

A square class

```text
[q(t_U)] in k(P1_U)^*/k(P1_U)^{*2}
```

belongs first to the base function field of the chosen fibration.  There is
no automatic map from `k(P1_{U_1})` to `k(P1_{U_2})`.  Its canonical
surface-level image is instead the restriction

```text
res_U([q]) in k(X)^*/k(X)^{*2}
```

along `k(P1_U) -> k(X)`.  A character for `U_1` transports to a base
character for `U_2` only after proving that this surface-field class descends
to `k(P1_{U_2})`.  Literal comparison of two polynomials is insufficient;
the certificate must include the field embeddings or a birational map and an
exact square-class identity.

The same caution applies to rational base change.  Under `t=phi(s)`, the
recorded operation is restriction

```text
[q(t)] -> [q(phi(s))].
```

It is not recorded as equality unless `phi` induces an isomorphism of the
named function fields.  A pullback can split a previously nontrivial class.

## Conditional commuting square

Let `f_i:X -> B_i` be two fibrations and let

```text
alpha_i in H^1(k(B_i),mu_2).
```

If their restrictions have been proved equal in `H^1(k(X),mu_2)`, then the
normalizations of

```text
X times_{B_i} C_i
```

are birational over `X`.  On a common smooth resolution `Y`, the deck action
on divisors is one common Galois object.  The two pulled-back fibrations give
different quotients of that divisor module by their respective trivial and
fibre-root lattices.  The rational representation statement is then the
Galois-equivariant Shioda--Tate identity already proved as Theorem A2 in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).

This gives the desired commuting diagram only under the displayed
surface-square-class equality.  Integral exactness additionally requires the
primitive closures, torsion, reducible components, height scaling, and
`|G|`-primary graph glue on `Y`.  Base change can enlarge `NS(Y)`, and the
pullback intersection form scales by the degree.  Neither integral MW
exactness nor saturation is inherited merely from the rational diagram.

## Identity keys and certified edges

The registry uses stable mathematical keys:

- a surface key names a marked arithmetic surface/period point, not a
  Weierstrass equation;
- a fibration key includes its primitive marked `U` inside that surface's
  `NS`;
- a base-character key includes its ambient function-field ID and normalized
  class;
- a surface-character key is the pulled-back class in `k(X)^*/k(X)^{*2}`;
- a divisor key uses coordinates in a common marked divisor lattice together
  with its field and Galois orbit;
- equation charts, polynomial normal forms, and section words are aliases or
  presentations, never identity keys.

Every `CERTIFIED_EXACT` node or edge has hash-pinned evidence.  Important
edge meanings are:

```text
surface --HAS_FIBRATION--> fibration
fibration --HAS_FRAME--> frame
fibration --HOP_TO--> fibration
fibration --HAS_BASE_FUNCTION_FIELD--> field
character --HAS_AMBIENT_FIELD--> field
base character --PULLS_BACK_TO--> surface character
surface character --DESCENDS_TO_BASE--> fibration
fibration --SUPPORTS_CARRIER--> carrier
carrier --BASE_CHANGES_TO--> surface cover
surface cover --HAS_FIBRATION--> pulled-back fibration
pulled-back fibration --HAS_MW_BLOCK--> character block
divisor --PRESENTED_BY--> pulled-back fibration
divisor --SAME_DIVISOR--> divisor
```

An exact `HOP_TO` edge is accepted only when both endpoints belong to the same
surface node.  An exact `DESCENDS_TO_BASE`, `SAME_SURFACE_CHARACTER`, or
`SAME_DIVISOR` edge requires a transport witness in a common ambient field or
lattice.  Missing transport is `UNKNOWN`, never inferred from matching labels.

## Determinant-948 seed

The seed contains the common surface fingerprint

```text
surface id:       K3-8188cdcda8c57b2d
T:                [[-2,0,1],[0,4,0],[1,0,118]], det -948
A_NS:             Z/948, with 2-primary part Z/4
Clifford algebra: (2,237/4), ramified at 2 and 3
integral order:   reduced discriminant 474, local level 79
stable group:     O+(T)^* = ker(O+(T) -> O(A_T))
```

The exact rational marking proves a rational point on the abstract stable
period curve, but an independent explicit model and coordinates on that
stable cover remain unavailable.  The period node is therefore
`CERTIFIED_PARTIAL`, not silently promoted from the explicit genus-two
Atkin--Lehner quotient point.

The registered same-surface fibrations are:

| fibration | root lattice | MW rank | frame class |
|---|---:|---:|---|
| published R17 | rootless | 17 | published R17 |
| alternate `11952` | rootless | 17 | alternate Q80 |
| hidden `103b2` | rootless | 17 | published R17 |
| noncyclic `4A1` | `4A1` | 13 | noncyclic `4A1/MW13` |

The quadratic `103b2` carrier is currently attached to the published-R17
base.  Its visible cover lattice has invariant block `R17(2)`, one
anti-invariant line `<16>`, and one order-two graph class.  The complete
anti-invariant MW group is not known.  Descent of its surface character to
the alternate-Q80, hidden-`103b2`, or `4A1` base is registered as `UNKNOWN`;
the shared label does not supply such a proof.

## Future-utility profile

Arithmetic admission and construction utility remain separate orders.  An
`ARITHMETICALLY_EXCLUDED` source cannot be rescued by an attractive frame or
character score, and an `UNKNOWN` feature is not scored as zero.

After admission, attach a versioned profile

```text
U(S) = (
  A_NS,p and exact local forms,
  certified primitive-nef U classes and rootless classes,
  certified U-graph connectivity and equation costs,
  allowed bridge/character glue modules,
  certified low-genus carriers,
  intersections of base-character images inside H^1(k(X),mu_n),
  exact character-block and saturation data
).
```

Before enough calibration exists, ranking is Pareto or lexicographic:

```text
arithmetic admission
  -> fibration richness
  -> character/carrier richness
  -> presentation-specific compilation cost.
```

Writing a numerical product of the three factors is allowed only after the
normalization, treatment of unknown values, and benchmark objective are
declared.  The utility profile is a scheduling diagnostic and does not change
mathematical status.

For every new exact construction, ingestion should perform four checks:

1. resolve or create the surface ID before creating a fibration ID;
2. store the fibration's literal marked `U` and base-field embedding;
3. normalize every character in its named ambient field and pull it back to a
   surface-character node;
4. compare divisors only in a common marked lattice on `X` or a certified
   common cover.

This turns cross-presentation transport into an explicit work queue rather
than an implicit hope.

## Fail-closed boundaries

- Same determinant or discriminant form does not identify a surface.
- An isometry class of frames does not identify a marked `U`.
- Matching polynomial text across unnamed coordinates does not identify a
  function-field extension.
- A character on one fibration base is not automatically a character on
  another base.
- A section on a base-changed surface is not a divisor on the original K3.
- Rational character decomposition does not determine the integral MW
  lattice or its saturation.
- A bounded carrier census supplies lower bounds and search affordances, not
  global nonexistence.
- Utility fields never change `UNKNOWN` or an arithmetic admission decision.

## Replay

```bash
python3 elkies-k3/scripts/validate_surface_fibration_character_graph.py
```

The validator checks the schema-level contract, node and edge types, common
surface inheritance for fibration hops, evidence hashes and JSON assertions,
and the exact lower-bound counts in the determinant-948 utility profile.

