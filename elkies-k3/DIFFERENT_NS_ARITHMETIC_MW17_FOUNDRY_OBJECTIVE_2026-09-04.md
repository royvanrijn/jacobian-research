# Lane B arithmetic-first foundry objective

Date: 2026-09-04.

Status: **OPEN**.

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-SOURCE-PRECURSOR 2e115b35c30a8cea -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-MARKED-RATIONAL-PARAMETER-SCAN ca678e520745dd3c -->
<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->
<!-- status-consumer: EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION 972f591d2885f9ba -->
<!-- status-consumer: EC-K3-DIFFERENT-NS-ARITHMETIC-GATE-RERANK 252991e141c42e55 -->
<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER 93e6c5626d369572 -->
<!-- status-consumer: EC-K3-DET378-QQ-MARKING-OBSTRUCTION 1e910f72f54ac228 -->
<!-- status-consumer: EC-K3-ARITHMETIC-FIRST-MARKED-T-FOUNDRY d94b3dbddf5cb529 -->
<!-- status-consumer: OP-K3-DIFFERENT-NS-ARITHMETIC-MW17 9f40eebe50b66ea4 -->

## Milestone

The determinant-948 equation programme is closed through the maximal
noncyclic bridge

```text
published R17 -> 4A1/MW13 -> published R17.
```

Another equation on that Néron--Severi lattice is not the next Lane B
milestone. The construction programme now starts from the rank-three
transcendental lattice and admits equation work only after the arithmetic
marking gate.

The primary success target is stronger than a plain new MW17 equation:

```text
a different-NS arithmetic MW17 fibration over QQ(t), together with a
certified positive-rank low-genus carrier base and an exact new section
after base change.
```

The stretch target is an integral `V4`-stable Mordell--Weil lattice with
character ranks

```text
17 + 1 + 1 + 1,
```

including its complete 2-primary graph glue. This would target generic rank
20 over the connected `V4` carrier rather than finding two characters and
leaving the product character undecided.

## Determinant 720 is closed first

The determinant-720 Golay surface is no longer an unresolved source. Its
literal transcendental lattice has determinant `-720` and

```text
A_NS = Z/2 + Z/6 + Z/60.
```

The exact split-Clifford calculation identifies the primitive-similarity
norm-one curve as `X_0(15)`. Computing the literal discriminant action gives
image `S3`; imposing the stable kernel adds the full mod-2 identity
condition. After rational conjugation the exact marked curve is

```text
X_0(60),       genus 7.
```

Its twelve rational points are cusps: the Mazur--Kenku rational cyclic
isogeny classification excludes degree `60`. Hence there is no full rational
determinant-720 marking and no positive point to hand to an equation agent.

The four noncuspidal rational points on the useful quotient `X_0(15)` do not
lift. The known rational `3A5/MW2` equation explains the boundary: its
displayed determinant-720 subgroup has index six, with rational `3`-torsion
and a rational half-section, and saturates to determinant `20`.

The theorem and exact replay are in
[`GOLAY_DET720_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](GOLAY_DET720_QQ_MARKING_OBSTRUCTION_2026-09-04.md).

## Global arithmetic-first order

The foundry order is now

```text
rank-three T
  -> exact integral Clifford order and O+(T)^*
  -> full marked Shimura/modular curve
  -> rational noncuspidal non-CM point
  -> saturated rational marking of NS = T^perp
  -> rootless-frame test+marked-U planning
  -> equation compilation.
```

The generated `T`-first planner enforces this order on all 827 catalogue
rows without using rootless-frame data in its arithmetic priority. It
propagates four exact exclusions and the already-realized determinant-948
positive control. The remaining arithmetic research queue has 822 rows: 62
from the old rootless-MW17 subcatalogue and 760 not yet screened at the NS
stage. Twenty-three currently have coarse genus at most two. That coarse genus
is a prioritization diagnostic only; it is not the genus of the stable
marking curve and proves neither a rational point nor a rational marking.

The split determinant-`378` row is now closed: its coarse `X_0(7)` becomes
`X_0(63)` after imposing the literal stable kernel, and all its rational
points are cusps. The next exact-coarse calculations are determinant `256`
and `512`, with coarse curves `X_0(2)` and `X_0(4)`; their literal stable
kernels are still unknown.
The priority tuple is

```text
(full marked genus 0/1, occasionally 2;
 rational non-CM point;
 small stable-kernel index;
 small similarity/marking gap).
```

Only a certified positive row may enter `NS=T^perp`, saturation, rootless,
or equation work. See the generated
[`../artifacts/generated-results/elkies-k3-arithmetic-first-marked-t-foundry-v1.json`](../artifacts/generated-results/elkies-k3-arithmetic-first-marked-t-foundry-v1.json).

## Arithmetic MW17 plus carrier certificate

A primary construction closes only when one certificate chain proves:

1. a characteristic-zero elliptic K3 over `QQ` with geometric Picard rank
   exactly 19 and nineteen individually rational divisor classes;
2. an exact saturated identification of its Néron--Severi lattice with the
   selected `T^perp`, including source `U` and incidence data;
3. target-free selection of a primitive nef rootless `U'` and an exact
   characteristic-zero endpoint equation;
4. seventeen displayed independent sections that generate the saturated
   Mordell--Weil lattice over `QQ(t)`, with torsion excluded;
5. a smooth geometrically integral carrier `C/QQ` of genus one, or another
   explicitly justified low genus, a nonconstant map `C -> P1`, and a
   certified positive-rank Jacobian;
6. an exact section on the pullback surface that is independent of the
   invariant MW17 lattice, with its height and character block certified.

If the last section is not saturated, the result remains a rank-lower-bound
certificate rather than the complete carrier objective.

## Prescribed integral V4 target

For the stretch target, begin with the group action rather than with a list
of unrelated quadratic covers. Specify a connected `V4` cover of the base
and four rational representation blocks

```text
M_1 of rank 17,
M_chi1, M_chi2, M_chi3 each of rank 1.
```

The construction data must also prescribe every allowed half-sum across the
four blocks. Exact local discriminant forms and isotropic graph subgroups at
`2` must determine an integral even `V4`-stable overlattice before equation
search. A completion certificate must then display one section in each
nontrivial character, prove independence and the full height pairing, prove
the claimed graph glue and saturation, and exclude torsion. Those assertions
give generic rank 20 over the `V4` carrier.

The current alternate-Q80 laboratory is a control, not this completion. It
has two nontrivial rational characters, two exact half-sum glues of total
index four, 64 rational genus-one `V4` bases, and seventeen base Jacobians of
rank one. Complete zero-Tate-class product inversion and bounded deeper
searches do not produce the third section. The nonzero product class, its
section, and full `V4` lattice saturation remain `UNKNOWN`. The relevant
boundaries are recorded in
[`INTEGRAL_RANK_TRANSFER_GLUE_CALCULUS_2026-09-02.md`](INTEGRAL_RANK_TRANSFER_GLUE_CALCULUS_2026-09-02.md),
[`R17_PRODUCT_TATE_COHOMOLOGY_REDUCTION_2026-09-04.md`](R17_PRODUCT_TATE_COHOMOLOGY_REDUCTION_2026-09-04.md), and
[`R17_RATIONAL_V4_DEEP_TRACE_EXHAUSTION_2026-09-04.md`](R17_RATIONAL_V4_DEEP_TRACE_EXHAUSTION_2026-09-04.md).

## Fail-closed boundaries

- A coarse Clifford curve is not a full marked curve.
- A rational point without a non-CM check and saturated rational marking is
  not an arithmetic source.
- A finite-field point, finite-precision lift, or formal local branch is not
  a rational marking.
- A rootless frame is not inspected merely because its lattice score is
  attractive; its `T` row must pass first.
- A supplied target, an unsaturated section subgroup, or a carrier with only
  a heuristic rank does not close the objective.
- Two quadratic characters do not imply a third product-character section.

## Replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_qq_marking_obstruction.sage --check
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det378_qq_marking_obstruction.sage --check
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage --check
python3 elkies-k3/scripts/build_rank7_determinant_aware_ranking.py --check
python3 elkies-k3/scripts/build_arithmetic_first_marked_t_foundry.py --check
```

## Canonical supporting records

- [`GOLAY_DET720_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](GOLAY_DET720_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
  — exact stable curve and determinant-720 exclusion.
- [`DET378_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](DET378_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
  — literal `A4` spin image, stable `X_0(63)`, and split determinant-378 exclusion.
- [`RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md`](RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md)
  — marking decisions and empty equation-agent handoff.
- [`DETERMINANT_AWARE_FOUNDRY_RANKING_2026-09-02.md`](DETERMINANT_AWARE_FOUNDRY_RANKING_2026-09-02.md)
  — retained rootless-subcatalogue accounting after arithmetic rejection.
- [`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
  and [`NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
  — the other exact rational-marking exclusions.
- [`MARKED_U_REALIZATION_PLANNER_2026-09-03.md`](MARKED_U_REALIZATION_PLANNER_2026-09-03.md)
  and [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
  — downstream planner and compiler contracts, used only after arithmetic
  admission.
