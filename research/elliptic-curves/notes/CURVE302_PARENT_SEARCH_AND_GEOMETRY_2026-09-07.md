# Curve302 parent: faster screening and sharper construction constraints

The [explicit determinant1092 parent](CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md)
now has a reusable parameter-screening adapter and two further construction
results. Its geometric Picard rank is **exactly19**, its full geometric MW
group is the displayed rational MW17, and its MW lattice has automorphism
group `{I,-I}`. These are properties of the recovered parent; they do not
identify the discoverers' unpublished procedure.

## Reuse the successful search infrastructure

The [shared runtime](SHARED_RESEARCH_RUNTIME.md) already separates finite-field
facts, parameter scores, curve normalization and point search. The
[factor-free control](BLIND_FACTOR_FREE_CONTROL_AND_PROSPECTIVE_EXPOSURE_2026-09-07.md)
recovers a known28th direction from an original27 subgroup, using49 frozen
charts. The [earlier reassessment](RANK_JUMP_REASSESSMENT_2026-09-05.md)
distinguishes candidate incidence, coordinate visibility and global solubility.
Those lessons transfer to1092; a successful parameter score alone does not
show that extra points exist or are visible in the current chart budget.

The new adapter reuses the existing `family_traces` implementation rather
than introducing another point-count engine. It stores exact traces for all
points of `P1(Fp)` at the first16 admissible primes through281: eight discovery
primes and eight held-out primes. There are **3,506 entries**, including all
sixteen infinity fibres. Every smooth trace was independently verified with
PARI elliptic cardinality. The nine singular entries remain `null`.
Primes with unsuitable reduction in this presentation are skipped, with the
reason recorded; they do not exclude the surface or any rational fibre.

For a reduced rational parameter `m/n`, the lookup is `m/n mod p` when
`p` does not divide `n`, and the infinity entry otherwise. The degree `(8,12)`
homogeneous Weierstrass scaling proves that these are the correct local
traces. This avoids evaluating large rational coefficients or recounting
elliptic curves separately for each parameter.

The fixed529-parameter engineering benchmark gave:

| Stage | Observed time |
|---|---:|
| Cold table preparation | 0.1434 seconds |
| Recount reductions for529 supplied parameters | 0.1088 seconds |
| Warm table lookup for the same parameters | 0.00147 seconds |

This is about **74 times faster for warm trace lookup** in the retained run
(68 times in the initial run). Cold preparation amortizes after roughly700
parameters against that comparator. These timings exclude process startup
and file loading. They do **not** measure point-search speed, predict rank,
or establish that a larger population improves discovery odds.

The [trace tables](../../artifacts/generated-results/elliptic-curves/curve302_parent_trace_tables_v1.json)
and [benchmark record](../../artifacts/generated-results/elliptic-curves/curve302_parent_trace_benchmark_v1.json)
are retained. Given a text file with one rational parameter per line:

```sh
python3 elliptic-curves/cas/curve302_parent_trace_cache.py --parameters parameters.txt
```

Output separates the two bands and records missing primes. Its score is the
finite sum `-sum(a_p log(p)/p)` over available primes. It is a scheduling
heuristic; `rank_bound` is always `null`. The held-out band can only serve as
a holdout if it is kept out of candidate selection.

The [text export](../../artifacts/generated-results/elliptic-curves/curve302_parent_trace_tables_v1.txt)
also works with the existing `newfamily/scan_rational_nagao_tables.cpp`.
All seven coprime positive parameters in the3-by3 adapter control were retained
and checked, including both bands' exact integer score-unit sums. Thus the
existing compiled rational scanner can consume this family without changes.
No production sweep or parameter promotion was performed. The separately
prepared two-fibre pilot and its frozen inputs were left untouched.

Rebuild and independently verify the tables with:

```sh
sage -python elliptic-curves/cas/prepare_curve302_parent_trace_cache.sage
```

## Geometric Picard rank is exactly19

Authority: `EC-CURVE302-PARENT-GEOMETRIC-PICARD19`.

The new checker first verifies all seventeen section identities and the full
height Gram directly from the equation. Together with the fibre and zero
section these give19 independent rational divisor classes. A good reduction
therefore leaves only a **three-dimensional** orthogonal complement in `H2`.
Its eigenvalues have the form

\[
 \epsilon p,\alpha,\beta,\qquad \epsilon\in\{1,-1\},\quad \alpha\beta=p^2.
\]

If `s1,s2` are the first two traces after removing the19 known divisor
directions, then

\[
 s_2=s_1^2-2\epsilon p s_1.
\]

The two counts determine a unique sign at both tested primes. This
reconstructs the **entire** Frobenius characteristic polynomial from exact
fibre counts and the known divisor lattice, without controlled reduction:

| Prime | Surface count over `Fp` | Surface count over `Fp²` | Remaining quadratic |
|---|---:|---:|---|
|149|24636|493345524|`T²+248T+22201`|
|151|25276|520355556|`T²+244T+22801`|

At either prime the full polynomial is
`(T-p)^19 (T+p)` times the displayed quadratic. The quadratic roots divided
by `p` are not roots of unity: their rational trace is not an integer.
The calculation at149 independently recovers the earlier full
controlled-reduction polynomial, not merely two otherwise insufficient
moments. The new argument uses the19 already verified divisor directions.

The Tate theorem for K3 surfaces in characteristic at least5 gives geometric
Picard rank20 in both reductions. Over `Fp²` the Artin–Tate formula determines
their NS discriminants modulo rational squares as

\[
 248^2-4\cdot149^2=-27300,\qquad
 244^2-4\cdot151^2=-31668.
\]

Their ratio is **29/25**, which is not a square. If the characteristic-zero
geometric NS rank were20, its specialization would have full rank in both
reductions and force the discriminants to have the same squareclass.
Consequently the geometric Picard rank is19. This is the
[van Luijk two-prime argument](https://pub.math.leidenuniv.nl/~luijkrmvan/ps/picone.pdf),
using the [proved Tate theorem](https://arxiv.org/abs/1301.6326), not an
unproved conjectural assumption.

Because there are no reducible fibres, Shioda–Tate now gives

\[
 \operatorname{rank}E(\overline{\mathbb Q}(t))
 =\operatorname{rank}E(\mathbb Q(t))=17.
\]

The existing no-overlattice obstruction proves that the displayed rational
basis is also the **full geometric basis**. The full geometric NS lattice is
`U + (-MW17)`, with determinant1092, and all its classes are rational.
The transcendental rank is3. In the rank19 Shimura/K3 moduli interpretation,
this surface is not a Picard20/CM point. Locating its explicit coordinate on
`X(546)/<w546>` remains a different problem.

This closes a search direction: extending the constant field cannot supply
an eighteenth independent section on this fibration. Every elliptic fibration
with section on this same K3 has geometric MW rank at most17. A nonconstant
base change creates another surface and is not excluded by that statement.

The [certificate](../../artifacts/generated-results/elliptic-curves/curve302_parent_geometric_picard19_v1.json)
replays in a bounded single-worker process:

```sh
sage -python elliptic-curves/cas/verify_curve302_parent_geometric_picard.sage
```

Finite fibre cardinalities use Sage/PARI. No full Frobenius backend is needed.
The original parent certificate and raw149 output remain unchanged.

## A quadratic construction obstruction

Authority: `EC-CURVE302-PARENT-QUADRATIC-DESCENT-GATE`.

PARI's complete lattice automorphism calculation gives `Aut(MW17)={I,-I}`.
An independent graph calculation supplies the same upper bound: enumerate
all2,436 norm-four vectors, join pairs of inner product two, and compute
the automorphisms of this158,496-edge graph. Its group has order two.
The norm-four shell spans rank17, so every lattice isometry acts faithfully
on this graph. Both signs are actual lattice isometries.

Suppose this pointed elliptic fibration descended along a quadratic rational
base map. Its deck involution would act by an integral height isometry.
The invariant subgroup is precisely the MW group of the source. Since the
only possible actions are `I` and `-I`, its rank must be17 or0. Therefore
**no source of arithmetic generic rank1–16 can produce this fibration by
quadratic base change**. In particular, the standard MW8 rational-source
construction is excluded. The same reasoning applies geometrically now that
the full geometric group is known.

This does not exclude a rank-zero source, another fibration on the K3,
higher-degree covers or multiparameter constructions. It does not establish
the discoverers' provenance. It does make the direct quartic/section-incidence
reconstruction more informative: there is no hidden nontrivial `8+9`
invariant/anti-invariant splitting of this MW lattice.

The [automorphism certificate](../../artifacts/generated-results/elliptic-curves/curve302_parent_quadratic_descent_gate_v1.json)
has two independent group-order checks:

```sh
sage -python elliptic-curves/cas/verify_curve302_parent_quadratic_descent_gate.sage
```

## Subsequent calibration and deployment

The [generic17-only adaptive calibration](CURVE302_RECOVERED_SUBGROUP_CALIBRATION_2026-09-07.md)
now completes196 boxes, reaches24 and identifies seven of the fourteen known
exceptional directions by exact identities checked only after the searches.
The [record-scale campaign](DET1092_RECORD_SCALE_CAMPAIGN_2026-09-07.md)
has completed its million-parameter intake and deployed the calibrated policy
on48 fixed prospective fibres. The following motivation predates those runs.

The most informative point-search control is302 with **only the17 generic
images** supplied to centre selection and execution. Keep the other14 public
directions outside the worker and use them only to evaluate recovery after
the fixed run. This tests the1092 lattice and its actual coordinate visibility;
the earlier27-to28 control tested a different surface and subgroup.

Such a control should compare frozen equal-budget centre policies and retain
every completed or censored box. Its outcome would guide whether effort should
go into point charts, adaptive subgroup enlargement or a larger parameter
population. The subsequent runs above address that control; no rank-success
calibration is claimed for the screening scores themselves.
