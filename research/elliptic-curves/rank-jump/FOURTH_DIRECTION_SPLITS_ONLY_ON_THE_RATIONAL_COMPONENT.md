# The fourth lift splits on the rational component of the relation scheme

The successful +7 quartet has a complete finite simultaneous-solubility
certificate under its retained three-cover relation. The degree-twelve
relation scheme has parameter-field degrees 1+11. Adding the fourth
root changes its component degrees to

\[
\boxed{1+1+22.}
\]

The two rational components lie over the **same** parameter t=−288/65
and differ only by the sign of the fourth root. The degree-eleven
component does not acquire that root in its residue field: its fourth
value has nonsquare norm, with an obstruction at 29.

The two matched negative relation schemes instead give single degree-24
fields. This is an exact comparison of finite arithmetic constructions,
not a bounded failure to find their points.

The result closes the conditional chain

\[
\text{a rational realization of the fixed successful triple relation}
\ \Longrightarrow\ t=-288/65
\ \Longrightarrow\ \text{four rational native lifts}
\ \Longrightarrow\ \text{three certified quotient directions}.
\]

It does not explain all seven retained quotient directions on that
original fibre. Nor does it explain why the selected relation has a
rational component in the first place.

## Fixed inputs and scope

Use A=`01333`, B=`0b2d0`, D=`19e45`, and the fourth cover C=`13109`.
The three schemes impose

\[
P_A-P_B+P_D=S_i(t)
\]

on the first three native lifts. Their generic words were selected from
the frozen lattice enumeration before comparing with the oracle relation.
The [matched-translate certificate](MATCHED_TRANSLATES_SEPARATE_RATIONAL_SPLITTING_FROM_LATTICE_CAPACITY.md)
proves that each scheme has proper intersection degree twelve and
verifies all its native roots and original elliptic-curve group relations.
The positive index 1 is the conjugate of the earlier successful relation.

For index 1,

\[
S_1=P_7-P_8+P_9-P_{12}+P_{14}-P_{16}+P_{17}
\]

in the published generic basis. Its parameter polynomial is

\[
F_1(t)=(t+288/65)H_{11}(t),
\]

with H11 irreducible over Q. The controls F0 and F2 are irreducible of
degree twelve. These coefficients and their previous proofs are imported
unchanged.

The additional primitive square polynomial is

\[
q_C(t)=459069811561+92081791006t-2589226631t^2.
\]

This is a previously retained successful cover label. Its selection is
retrospective; the calculation is not a new prospective selector.

## The complete fourth-root algebra

For each irreducible factor F of a parameter polynomial, let

\[
K=\mathbf Q[t]/(F),\qquad b=q_C(\bar t).
\]

The corresponding fourth-lift algebra is K[z]/(z²−b). Exact gcds
show that b is nonzero on every component. The projection is therefore
étale of degree two over the entire finite relation scheme. There is no
branch component or excluded denominator to add to this accounting.

We compute

\[
N_{K/\mathbf Q}(b)=\operatorname{Res}(F,q_C)
\]

for the monic factors F. The independent verifier recomputes the norm
as the determinant of multiplication by qC in the power basis of K.
If b were square in K, its rational norm would be square. This one-way
implication suffices for all three nonrational components here.

| Relation | Parameter component | Fourth-value norm | Fourth-root algebra |
|---|---:|---|---|
| 0 | Degree 12 | Negative | A degree-24 field |
| 1 | Rational, t=−288/65 | (1018829/65)² | Q × Q |
| 1 | Degree 11 | A 29-adic unit congruent to 15 modulo 29, a nonsquare | A degree-22 field |
| 2 | Degree 12 | Negative | A degree-24 field |

The negative norms certify nonsquareness over the control residue fields;
they do not describe a real obstruction to the original elliptic fibre.
For H11, the norm obstruction means the fourth value cannot be square
at every completion above 29. No individual prime ideal above 29 has
been designated without further factorization.

To verify complete component degrees independently, the checker uses
the primitive element w=t+sqrt(qC(t)) and forms

\[
\operatorname{Res}_t(F(t),(w-t)^2-q_C(t)).
\]

Each resulting polynomial has degree twice deg F and is squarefree.
Thus w has the full number of distinct conjugates and generates the
corresponding field in the nonsquare cases. Exact factorization confirms
the complete degree-24 patterns

\[
\begin{array}{c|ccc}
\text{relation index}&0&1&2\\\hline
\text{fourth-lift component degrees}&24&1+1+22&24.
\end{array}
\]

All polynomial coefficients and norm values are retained in the two
certificates. No number-field class group computation or point search
was necessary.

## What this says about the successful block

There are two arithmetic splitting steps in the displayed construction:
the triple relation scheme acquires a rational component, and the
fourth-root quadratic algebra splits over that component. The fourth
root is **not** supplied by a rational function on the entire finite
triple-relation scheme, since it is absent over H11.

This does not disprove a rational-point implication for the particular
successful relation. On the contrary, F1 has only one rational root,
and qC is square there, so every rational realization of that relation
does extend to the fourth cover. The distinction is between a uniform
algebraic lift across all components and the arithmetic splitting of
the selected rational component.

The first three displayed points have quotient rank two. The two fourth
root signs are conjugate native points whose sum is a generic trace.
They supply the same additional quotient line, not two extra directions.
The [exact quartet certificate](PAIRED_SOLUBILITY_AND_SPECIALIZATION_COLLAPSE.md)
proves that the full quartet has quotient rank three over the marked
rank-17 subgroup. Its generated specialized subgroup has rank twenty,
inside the retained rank-24 subgroup of the original fibre.

Thus the evidence supports a **two-direction construction extended by
one additional soluble cover**. It does not establish statistical
independence of the two splitting events, or show that the fourth
direction is an unrelated accident. In particular, whether the
degree-22 field is a compositum of the degree-eleven field with one
constant quadratic field remains a separate arithmetic question.

## Mechanism and remaining implication

The strongest concrete mechanism now has a complete finite expression:
the successful marked construction has a rational component whose
remaining quadratic lift splits. Its matched controls have no rational
component even after adjoining the fourth root. This is information
about **global solubility**, not point-search visibility.

The missing prospective implication is a condition on generic cover
labels and translates that forces, or reliably identifies, such a
rational component with enough further split lifts. The matched controls
already show that the same trace norms and intersection degree do not
force it. A larger Selmer space or an auxiliary product point alone
does not provide that implication.

This finite scheme is cut out by an additional generic group relation.
Its complete rational-point calculation does not enumerate all rational
points on the genus-17 carrier: other points can obey other relations,
or no relation of the prescribed type. The other four retained quotient
directions on the +7 fibre remain outside this construction.

## Replay

The bounded [protocol](FOURTH_LIFT_ON_RELATION_PROTOCOL.json) fixes the
three existing relation schemes and one fourth cover. The square-norm
fallback was not needed: all nonrational components had nonsquare norms.

```sh
sage -python elliptic-curves/rank-jump/fourth_lift_on_relation.py check
sage -python elliptic-curves/rank-jump/verify_fourth_lift_on_relation.py check
```

The immutable outputs are
`rank_jump_fourth_lift_on_relation_norms_v1.json` and
`rank_jump_fourth_lift_on_relation_verification_v1.json` under
`artifacts/generated-results/elliptic-curves/`. No active search file,
candidate population, original parameter search or online computation
was involved.
