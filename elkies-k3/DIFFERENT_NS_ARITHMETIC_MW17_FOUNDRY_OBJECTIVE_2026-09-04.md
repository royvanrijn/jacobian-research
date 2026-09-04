# Different-NS foundry objective: first planner-found arithmetic MW17

Date: 2026-09-04.

Status: **OPEN**.

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-SOURCE-PRECURSOR 2e115b35c30a8cea -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-MARKED-RATIONAL-PARAMETER-SCAN ca678e520745dd3c -->
<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->
<!-- status-consumer: EC-K3-DIFFERENT-NS-ARITHMETIC-GATE-RERANK f968ac0d6fa311fb -->
<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER 40745008c2fe2a80 -->
<!-- status-consumer: OP-K3-DIFFERENT-NS-ARITHMETIC-MW17 e93bdd3228be30d0 -->

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

## Two arithmetic exclusions

Determinant-950 `NS0024` is not a candidate. A full rational marking would
make its primitive `2E8/MW1` Inose fibration descend and hence give a
noncuspidal, non-CM rational point on `X0+(475)`. Momose's theorem excludes
such a point. The exact argument is in
[`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md).

Determinant-1184 `NS0031` is also not a candidate. Its exact split Clifford
order gives the modular curve

```text
X_ns(4) x_{X(1)} X_0(37).
```

Vélu's determination of `X_0(37)(QQ)`, followed by an exact mod-4 Frobenius
test at `19`, excludes a lift of either noncuspidal rational point. Thus a
full rational `NS0031` marking is impossible. The proof and replay are in
[`NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md).

The one-parameter formally smooth `ZZ_7` model-157 branch and the physical
corridor

```text
A1+2A7 -> A1+2A3+A5 -> 5A1+A3 -> 4A1 -> A1
       -> rootless NS0031-F017
```

remain exact local and geometric controls. They cannot be promoted to the
required arithmetic source over `QQ`. Longer lifting, a larger rational
parameter scan, or algebraization of that branch is no longer a live attack
on this milestone.

## First missing gate

Rerank the remaining different-NS frames with the rational-marking arithmetic
gate before equation search. Determinant `720` is the strongest existing
lattice/corridor control, but its known rational `3A5/MW2` equation is not a
source for the determinant-720 lattice: the displayed determinant-720
sublattice saturates with index `6` to determinant `20`, with rational
3-torsion and a rational half-section.

The reranked planner artifact makes this fail-closed: 66 surfaces pass the
exact lattice filters, the two proved arithmetic obstructions are removed,
and 64 candidates remain. After the same-NS determinant-948 control, the
determinant-720 surface is the first unresolved row.

The separate arithmetic-marking classifier sharpens the type of those rows:
one surface is `ARITHMETICALLY_POSSIBLE` (the already-realized
determinant-948 control), two are `ARITHMETICALLY_EXCLUDED`, and the remaining
63 are `UNKNOWN`. Its new different-NS equation-agent handoff is therefore
empty. `UNKNOWN` rows stay in the arithmetic-curve queue and may not be sent
to an equation agent. See
[`RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md`](RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md).

The next narrow calculation is therefore an exact arithmetic moduli decision
for the determinant-720 marking: identify the stable marked curve and decide
whether any rational noncuspidal point has saturated determinant `720` rather
than landing on a proper overlattice. Only after that decision should the
coefficient search resume. If determinant `720` is excluded, repeat the same
pre-screen down the determinant-aware queue.

The first positive source package must contain:

1. an explicit smooth elliptic K3 equation over `QQ`, with rational fibre and
   effective zero;
2. source sections and all reducible-fibre components giving nineteen
   individual `QQ`-rational divisor classes;
3. their exact intersection matrix, a determinant-one identification with
   the selected new Neron--Severi lattice, and a proof that the geometric
   Picard rank is exactly 19;
4. complete component and section incidence data for the marked-`U` planner
   and equation compiler.

The old NS0031 rational-coordinate scan remains a bounded negative record:
all 247 reduced values `m9=n/d` with `|n|,d<=40` in the model-157 residue disk
lift through `7^40`, but none simultaneously rationally reconstructs all 52
coordinates. The modular obstruction supersedes that search as an arithmetic
route.

## Target-free planner and compiler protocol

Once a source passes the arithmetic gate, run the marked-`U` planner with
only:

- the explicit source `(NS,U,W)` marking;
- a rootless rank-17 endpoint predicate of the selected NS determinant;
- a declared low-degree search box beginning with old-fibre degree two;
- no target Gram, target frame identifier, historical route, or endpoint
  overlap fingerprint.

The selected primitive `U'` must independently pass nefness and
effective-zero gates. Existing rootless frames and certified corridors may be
used only after selection to identify or compare the landing class.

Compile the selected moves over `QQ`. The universal degree-two compiler is
the preferred backend when applicable; another selected degree requires its
own exact compiler certificate. Modular equations remain discovery data, not
characteristic-zero endpoints.

## Arithmetic acceptance certificate

The problem closes only when one certificate chain proves:

1. the compiled equation is the same characteristic-zero K3 with the selected
   geometric Neron--Severi lattice;
2. the selected fibration has no reducible fibres;
3. seventeen displayed sections lie in `E(QQ(t))`, are independent, and
   generate a saturated Shioda lattice of the selected determinant;
4. torsion is excluded, so the arithmetic Mordell--Weil group is `ZZ^17`;
5. the planner transcript proves target-free selection occurred before
   endpoint construction and comparison.

A finite-field marked point, finite-precision or purely formal lift, abstract
frame Gram, preselected target, unsaturated section subgroup, or endpoint
equation without the rational rank-19 source marking does not close the
problem.

## Canonical supporting records

- [`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
  — theorem excluding determinant `950` over `QQ`.
- [`NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
  — split-Clifford/fibre-product theorem excluding determinant `1184` over
  `QQ`.
- [`RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md`](RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md)
  — batch classifier, marking-level proof boundary, and empty fail-closed
  different-NS equation handoff.
- [`LATTICE_FOUNDRY_EQUATION_FIRST_SHORTLIST_2026-09-02.md`](LATTICE_FOUNDRY_EQUATION_FIRST_SHORTLIST_2026-09-02.md)
  and [`NS0031_MARKED_FORMAL_BRANCH_2026-09-04.md`](NS0031_MARKED_FORMAL_BRANCH_2026-09-04.md)
  — retained NS0031 finite-field, finite-lift, and formal-local controls.
- [`SAME_NS_COMPILER_ROUTES_2026-09-02.md`](SAME_NS_COMPILER_ROUTES_2026-09-02.md)
  — exact determinant-720 and NS0031 physical control corridors.
- [`GOLAY_OCTAD_LATTICE_DESIGN_2026-09-01.md`](GOLAY_OCTAD_LATTICE_DESIGN_2026-09-01.md)
  — exact determinant-720 lattice and the saturation rejection of the known
  rational `3A5` point.
- [`MARKED_U_REALIZATION_PLANNER_2026-09-03.md`](MARKED_U_REALIZATION_PLANNER_2026-09-03.md)
  and [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
  — planner and compiler contracts.
