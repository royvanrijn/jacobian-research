# Lane B arithmetic-first foundry objective

Date: 2026-09-04.

Status: **OPEN**.

<!-- status-consumer: EC-K3-DET1236-GENUS2-RATIONAL-POINTS 5a3c84eb9f7f0604 -->
<!-- status-consumer: EC-K3-DET1236-MARKED-SHIMURA-CURVE 185d31609e7702fc -->

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-SOURCE-PRECURSOR 2e115b35c30a8cea -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-MARKED-RATIONAL-PARAMETER-SCAN ca678e520745dd3c -->
<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->
<!-- status-consumer: EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION 972f591d2885f9ba -->
<!-- status-consumer: EC-K3-DIFFERENT-NS-ARITHMETIC-GATE-RERANK d569364c553007a2 -->
<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER 6043be45b20f8241 -->
<!-- status-consumer: EC-K3-DET378-QQ-MARKING-OBSTRUCTION 1e910f72f54ac228 -->
<!-- status-consumer: EC-K3-ARITHMETIC-FIRST-MARKED-T-FOUNDRY 9e9e0a1a8ac7c088 -->
<!-- status-consumer: EC-K3-DET500-DET750-QQ-MARKING-OBSTRUCTIONS 14498ad134ffa60e -->
<!-- status-consumer: OP-K3-DIFFERENT-NS-ARITHMETIC-MW17 b8ef1932e51636fa -->

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

## Two rootless arithmetic rows are closed

The determinant-500 and determinant-750 rootless-MW17 rows are now decided
at the exact marked-curve gate. Their literal lattices are

```text
T_N = U(5) + <10N> = 5(U+<2N>),       N=2,3.
```

The coarse curves `X_0(2)` and `X_0(3)` forget projective level five. Exact
Clifford and discriminant-action calculations give stable groups

```text
Gamma_0(N) intersection +/-Gamma(5).
```

After conjugation these are `Gamma_H(50)` and `Gamma_H(75)`, where `H` is the
inverse image of `{+/-1}` modulo five. The exact marked curves have genera
four and nine and degree-two maps to `X_0(50)` and `X_0(75)`. Each has exactly
four rational cusps and no rational noncuspidal point by Mazur--Kenku.
Therefore both rows are `ARITHMETICALLY_EXCLUDED`, satisfying the first
two-row arithmetic-classification gate without authorizing equation work.
See
[`DET500_DET750_QQ_MARKING_OBSTRUCTIONS_2026-09-04.md`](DET500_DET750_QQ_MARKING_OBSTRUCTIONS_2026-09-04.md).

## Determinant 1236 is the positive-directed marked-curve target

The content-one cyclic row `K3-6d288cfad55e0d15` now has an exact Phase-1
certificate. Its literal Clifford order has Eichler pair `(D,N)=(6,103)`,
and its Atkin--Lehner action exhausts `O(A_T)` for `A_T=Z/1236Z`. The exact
projective stable curve is

```text
C_1236 = X_0^6(103)/<w_618>,       genus 6.
```

It maps with degree two to the explicit genus-two curve

```text
y^2 = 1944*x^6 + 441*x^4 - 90*x^2 + 9,
```

and with degree four to the rank-one elliptic curve `618f1`. The exact marked
curve itself has two rational discriminant-`-3` CM points, but both are
Picard-rank-20 specializations. No non-CM rational lift from the genus-two
quotient has been certified.

An exact bielliptic quadratic-Chabauty/Mordell--Weil-sieve certificate proves
that these are all fourteen rational points on the genus-two quotient,
mapping to `+/-G`, `+/-3G`, `+/-4G`, and `+/-10G` on `618f1`. It also resolves
both rational fixed fibers: one comes from the rational discriminant-`-3` CM
pair, while the other comes from a quadratic discriminant-`-24` CM pair and
has no rational point upstairs. Thus the remaining covering descent starts
with twelve explicit non-fixed rational quotient points.

Exact Jacquet--Langlands accounting gives
`Jac(C_1236) ~ 618a1 x ... x 618f1`, all six factors of rank one. The
genus-two quotient contributes `618e1 x 618f1`, while the missing marking
cover has Prym `618a1 x ... x 618d1`. Thus classical Chabauty is exactly at
rank equal to genus; the split-Jacobian quadratic-Chabauty dimension screen
passes, but it still requires the explicit degree-two cover.

The row therefore has Phase-2 status `UNRESOLVED_FOR_EXPLICIT_REASON` and
remains `UNKNOWN` in the three-way equation-dispatch classifier. Its next
task is the degree-two covering descent, not a K3 equation. See
[`DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md`](DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md).

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
propagates six exact exclusions and the already-realized determinant-948
positive control. The remaining arithmetic research queue has 820 rows: 60
from the old rootless-MW17 subcatalogue and 760 not yet screened at the NS
stage. Twenty-one currently have coarse genus at most two. That coarse genus
is a prioritization diagnostic only; it is not the genus of the stable
marking curve and proves neither a rational point nor a rational marking.

The split determinant-`378` row is now closed: its coarse `X_0(7)` becomes
`X_0(63)` after imposing the literal stable kernel, and all its rational
points are cusps. Among actual rootless-MW17 candidates, determinants `500`
and `750` are also closed by their exact `X_H(50)` and `X_H(75)` curves. The
next global exact-coarse calculations are determinant `256`
and `512`, with coarse curves `X_0(2)` and `X_0(4)`; their literal stable
kernels are still unknown.
Within the already-rootless queue, determinant `1236` is now the strongest
positive-directed row because its literal stable kernel, genus-six marked
curve, genus-two quotient model, fourteen rational quotient points, and the
two fixed CM fibers are exact. Its six-factor Jacobian and four-factor Prym
are also exact. The remaining unknown is specifically the rational lift
through one degree-two cover at the complete list of twelve non-fixed points.
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
  elkies-k3/scripts/certify_det500_det750_qq_marking_obstructions.sage --check
sage -- elkies-k3/scripts/certify_det1236_genus2_rational_points.sage --check
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage --check
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
- [`DET500_DET750_QQ_MARKING_OBSTRUCTIONS_2026-09-04.md`](DET500_DET750_QQ_MARKING_OBSTRUCTIONS_2026-09-04.md)
  — literal mod-five stable kernels, exact `X_H(50)`/`X_H(75)` curves, and two
  rootless-MW17 arithmetic exclusions.
- [`DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md`](DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md)
  — exact genus-six stable marked curve, low-genus quotient tower, two rational
  CM controls, complete fourteen-point genus-two quotient, and the degree-two
  rational-lift obstruction.
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
