# Different-NS foundry objective: arithmetic MW17 on NS0024

Date: 2026-09-04.

Status: **OPEN**.

## Milestone

The determinant-948 equation programme has reached the maximal noncyclic
bridge

```text
published R17
  -> noncyclic 4A1/MW13
  -> published R17,
```

with an explicit characteristic-zero equation, thirteen saturated rational
sections, and target-free reverse selection.  Another equation on that same
Neron--Severi lattice is therefore no longer the principal foundry milestone.

The next milestone is

```text
the first arithmetic MW17 fibration on a different Neron--Severi lattice
that is found by the planner rather than supplied as its target.
```

The preferred test class is `NS0024`, of determinant `950`.  Its determinant
already proves that a successful surface is not isometric to the
determinant-948 R17 surface.  The existence of suitable geometric frames and
positive-lattice routes is proved, but the arithmetic fibration in the
displayed milestone remains unknown.

## First missing gate

Before another planner or compiler campaign, supply one equation-facing
characteristic-zero source over `QQ` with a **rational rank-19 marking**.
Concretely, the source package must contain:

1. an explicit smooth elliptic K3 equation over `QQ` and its rational
   fibre and effective zero;
2. nineteen independent divisor classes defined over `QQ`, including the
   source `U`, the reducible-fibre components, and the source sections;
3. their exact intersection matrix, a determinant-one identification with
   `NS0024`, and a proof that the geometric Picard rank is exactly `19`;
4. the complete component and section incidence data needed by the marked-`U`
   planner and the equation compiler.

Here “rational” means fixed by `Gal(Qbar/QQ)` and represented in the explicit
equation, not merely a basis of `NS(X_Qbar) tensor QQ`.

No modular model, Frobenius-fixed rank estimate, abstract frame Gram, or
characteristic-zero equation without the full marking passes this gate.

## Existing NS0024 inputs and their boundaries

Three existing routes are useful, but none supplies the missing source
package.

- `A3+A4+A6/MW4` has a certified thirteen-edge degree-two marked lattice
  corridor to the catalogue frame `NS0024-F005`, and the first equation
  compiler adapter is prepared.  It has no characteristic-zero equation with
  the required rational rank-19 marking.
- `D5+E8/MW4` is attached to a short completed-core sequence ending in a new
  determinant-950 rootless frame.  The primes `17,13,7` are core Kneser
  neighbour primes, not elliptic-neighbour degrees.  Moreover, the requested
  first marked-`U` edge is absent in the certified degree-two through
  degree-four search boxes.  This route cannot be substituted for a marked
  elliptic corridor.
- `2E8/MW1` has an Inose source formula on the degree-475 isogeny locus, but
  the direct `QQ` specialization is impossible: `475` is absent from the
  complete Mazur--Kenku list of rational cyclic-isogeny degrees.  An
  Atkin--Lehner or quadratic descent is not excluded, but would still have to
  descend all nineteen divisor classes individually.  Its minimum
  source-section pole order is `473`.

These are alternative source strategies on the same determinant-950 class.
They must not be spliced together unless an exact unimodular change of marking
on one explicit surface is supplied.

The active source strategy is therefore `A3+A4+A6/MW4`: recover a common
characteristic-zero producer for the semistable surface and its four resolved
sections, then certify the rational marking and geometric Picard rank.  The
`D5+E8/MW4` frame remains a geometric control, and quotient descents of the
Inose source remain a separate speculative route rather than a substitute for
that reconstruction.

## Planner and compiler protocol

Once the source gate closes, run the marked-`U` realization planner with:

- the explicit source `(NS,U,W)` marking;
- a rootless rank-17 endpoint predicate and determinant `950`;
- a declared low-degree search box, beginning with old-fibre degree two;
- no target frame Gram, historical neighbour line, or target-overlap
  fingerprint.

A retained planner result must give a literal primitive `U'` in the same
Neron--Severi lattice and independently pass nefness and effective-zero gates.
The determinant-950 rootless frames already stored in the repository are
post-selection landing tests only.  Using one of them as planner input does
not meet the “found by the planner” milestone.

If the selected edge has old-fibre degree two, feed its exact divisor and
physical incidence data to the universal degree-two compiler.  A
characteristic-zero Weierstrass equation is required; a modular compilation
is a discovery certificate only.  If the planner selects another degree, that
result remains a marked route and requires a separately certified compiler.

## Arithmetic acceptance certificate

The open problem closes only when one certificate chain proves all of the
following on the compiled endpoint:

1. the equation defines the same characteristic-zero K3 with
   `NS(X_Kbar) = NS0024`;
2. the selected fibration has no reducible fibres, so its geometric
   Mordell--Weil rank is `17`;
3. seventeen displayed sections lie in `E(QQ(t))`, are independent, and have
   saturated Shioda height lattice of determinant `950`;
4. torsion is excluded and the arithmetic Mordell--Weil group is exactly
   `Z^17`;
5. the planner transcript records target-free selection before endpoint
   construction and equation comparison.

The rank-19 rational source marking may also close the arithmetic rank step by
equivariant rank transfer, but it does not remove the requirement for the
endpoint equation and saturated seventeen-section basis in this milestone.

## Canonical supporting records

- [`R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md`](R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md)
  — determinant-948 closure that triggers this pivot.
- [`LATTICE_FOUNDRY_REPORT_2026-09-01.md`](LATTICE_FOUNDRY_REPORT_2026-09-01.md)
  and [`NS0024_EDGE1_COMPILER_PREPARATION_2026-09-01.md`](NS0024_EDGE1_COMPILER_PREPARATION_2026-09-01.md)
  — thirteen-edge route and prepared first degree-two compiler input.
- [`NS0024_NEW_ROOTLESS_SOURCE_ROUTE_2026-09-03.md`](NS0024_NEW_ROOTLESS_SOURCE_ROUTE_2026-09-03.md)
  and [`RELATIVE_U_BRIDGE_LIFTING_2026-09-03.md`](RELATIVE_U_BRIDGE_LIFTING_2026-09-03.md)
  — completed-core route, Inose source option, and marked-`U` bounded
  obstruction.
- [`NS0024_DIRECT_QQ_INOSE_OBSTRUCTION_2026-09-04.md`](NS0024_DIRECT_QQ_INOSE_OBSTRUCTION_2026-09-04.md)
  — exact removal of the direct `X0(475)(QQ)` source route and the remaining
  quotient-descent boundary.
- [`ARITHMETIC_RANK_TRANSFER_2026-09-03.md`](ARITHMETIC_RANK_TRANSFER_2026-09-03.md)
  — exact arithmetic promotion theorem and current NS0024 fail-closed result.
- [`MARKED_U_REALIZATION_PLANNER_2026-09-03.md`](MARKED_U_REALIZATION_PLANNER_2026-09-03.md)
  and [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
  — planner contract, degree-two compiler, and theorem layer.
