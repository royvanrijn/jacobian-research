# Different-NS foundry objective: first planner-found arithmetic MW17

Date: 2026-09-04.

Status: **OPEN**.

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-SOURCE-PRECURSOR 2e115b35c30a8cea -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-MARKED-RATIONAL-PARAMETER-SCAN ca678e520745dd3c -->
<!-- status-consumer: OP-K3-DIFFERENT-NS-ARITHMETIC-MW17 c8566fbe8f4dc838 -->

## Milestone

The determinant-948 equation programme is closed through the maximal
noncyclic bridge

```text
published R17 -> 4A1/MW13 -> published R17,
```

including a characteristic-zero equation, thirteen saturated rational
sections, and target-free reverse selection. Another equation on that same
Neron--Severi lattice is therefore not the next foundry milestone.

The live milestone is

```text
the first arithmetic rootless MW17 fibration on a different rank-19
Neron--Severi lattice that is selected by the marked-U planner without a
supplied target frame.
```

Determinant-950 `NS0024` is no longer a candidate. A full rational marking
would make its primitive `2E8/MW1` Inose fibration descend and hence give a
noncuspidal, non-CM rational point on `X0+(475)`. Momose's theorem excludes
such a point. The exact argument is in
[`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md).

## Preferred working class: NS0031

The preferred replacement is determinant-1184 `NS0031`. This priority is
not inferred from determinant alone. It is the strongest current candidate
that simultaneously has a positive equation-level source precursor, a
low-pole complete source basis, rootless MW17 frames, and an exact physical
same-NS corridor.

The stable source key is the pair

```text
(artifacts/generated-results/
 elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json,
 NS0031-S001).
```

The artifact qualification is essential: `source_id` values are shard-local
and are not globally unique. This source has

```text
root type       A1+2A7
source MW rank  2
height Gram     [[2,1],[1,41/8]]
basis poles     [0,1]
```

Its normalized square-twist chart over `GF(7)` contains two complete marked
MW2 pairs on model 157. At one pair, the 59-equation system in 52 variables
has Jacobian rank 51 with a unit maximal minor, and explicit coordinates
solve every equation through `7^8`. The eight rows outside that minor are
now forced on the localized marked scheme by an exact discriminant/node
identity and the fibre/component orders. Consequently the full model-157
germ is a one-parameter formally smooth `ZZ_7` branch. This remains a formal
local result: it does not algebraize the branch or produce a rational
characteristic-zero point.

The exact marking-level corridor is

```text
A1+2A7 -> A1+2A3+A5 -> 5A1+A3 -> 4A1 -> A1
       -> rootless NS0031-F017.
```

All five edges have old-fibre degree two and exact unimodular transport. The
corridor proves physical same-lattice reachability, but `F017` is only a
post-selection control for this objective. Feeding `F017`, its Gram matrix,
or this route to the planner as a target would not meet the milestone.

`NS0005` is not the replacement lead. Its similarly normalized `A1+2A7`
charts at 5 and 7 contain individual generator types but no pair with the
required mutual height. That bounded negative result is not a global
obstruction, but it is strictly weaker source evidence than the positive
`NS0031` chart.

## First missing gate

Algebraize the `NS0031` marked formal branch and produce one equation-facing
characteristic-zero source over `QQ` with a rational rank-19 marking. The
package must contain:

1. an explicit smooth elliptic K3 equation over `QQ`, with rational fibre and
   effective zero;
2. the two source sections and all reducible-fibre components as nineteen
   individual `QQ`-rational divisor classes;
3. their exact intersection matrix, a determinant-one identification with
   determinant-1184 `NS0031`, and a proof that the geometric Picard rank is
   exactly 19;
4. complete component and section incidence data for the marked-`U` planner
   and equation compiler.

The next narrow calculation is to algebraize the certified model-157 formal
branch, for example as a characteristic-zero section-first or
rational-surface-base-change curve containing that residue disk, and then
find and certify a rational point on it. A longer finite `7`-adic lift, or a
purely formal family without algebraization, does not pass this gate.

The first exact rational-coordinate scan on the formal branch is negative in
its declared box: all 247 reduced values `m9=n/d` with `|n|,d<=40` in the
model-157 residue disk lift through `7^40`, but none simultaneously rationally
reconstructs all 52 coordinates. This is a bounded miss, not a rational-point
obstruction. It redirects the next calculation toward an explicit algebraic
model of the marked curve rather than a blind enlargement of the same box.

## Target-free planner and compiler protocol

Once the source gate closes, run the marked-`U` planner with only:

- the explicit source `(NS,U,W)` marking;
- a rootless rank-17 endpoint predicate of determinant 1184;
- a declared low-degree search box beginning with old-fibre degree two;
- no target Gram, target frame identifier, historical route, or endpoint
  overlap fingerprint.

The selected primitive `U'` must independently pass nefness and
effective-zero gates. Existing `NS0031` rootless frames and the certified
`F017` corridor may be used only after selection to identify or compare the
landing class.

Compile the selected moves over `QQ`. The universal degree-two compiler is
the preferred backend when applicable; another selected degree requires its
own exact compiler certificate. Modular equations remain discovery data,
not characteristic-zero endpoints.

## Arithmetic acceptance certificate

The problem closes only when one certificate chain proves:

1. the compiled equation is the same characteristic-zero K3 with geometric
   Neron--Severi lattice `NS0031`;
2. the selected fibration has no reducible fibres;
3. seventeen displayed sections lie in `E(QQ(t))`, are independent, and
   generate a saturated Shioda lattice of determinant 1184;
4. torsion is excluded, so the arithmetic Mordell--Weil group is `ZZ^17`;
5. the planner transcript proves target-free selection occurred before
   endpoint construction and comparison.

A finite-field marked point, finite-precision or purely formal lift, abstract
frame Gram, preselected target, unsaturated section subgroup, or endpoint
equation without the rational rank-19 source marking does not close the
problem.

## Canonical supporting records

- [`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
  — theorem excluding the former determinant-950 candidate over `QQ`.
- [`LATTICE_FOUNDRY_EQUATION_FIRST_SHORTLIST_2026-09-02.md`](LATTICE_FOUNDRY_EQUATION_FIRST_SHORTLIST_2026-09-02.md)
  — exact NS0031 source, marked finite-field pair, tangent/lift certificate,
  and complete degree-three comparisons.
- [`NS0031_MARKED_FORMAL_BRANCH_2026-09-04.md`](NS0031_MARKED_FORMAL_BRANCH_2026-09-04.md)
  — exact dependence of the eight omitted residual rows and the resulting
  one-parameter formally smooth `ZZ_7` marked branch.
- [`SECTION_FIRST_NORMAL_FORM_COMPILER_2026-09-02.md`](SECTION_FIRST_NORMAL_FORM_COMPILER_2026-09-02.md)
  and [`RATIONAL_SURFACE_BASE_CHANGE_AUDIT_2026-09-02.md`](RATIONAL_SURFACE_BASE_CHANGE_AUDIT_2026-09-02.md)
  — equation chart and model-157 quadratic-base-change structure.
- [`SAME_NS_COMPILER_ROUTES_2026-09-02.md`](SAME_NS_COMPILER_ROUTES_2026-09-02.md)
  — exact five-edge physical control corridor to `NS0031-F017`.
- [`MARKED_U_REALIZATION_PLANNER_2026-09-03.md`](MARKED_U_REALIZATION_PLANNER_2026-09-03.md)
  and [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
  — planner and compiler contracts.
