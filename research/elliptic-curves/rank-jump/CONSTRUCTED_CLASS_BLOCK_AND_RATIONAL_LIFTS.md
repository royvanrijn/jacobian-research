# A constructed class-group block with two rational lifts

The subsequent [blind recovery and V3 comparison](BLIND_CONSTRUCTED_CLASS_RECOVERY_2026-09-12.md)
recovers rational points on both frozen covers without reading the known point
oracle. This report retains the original class-block proof and its explicitly
retrospective first solubility evaluation.

The two equation/generic-only strict classes on MW16-05 at `t=3/17` now
have a certified arithmetic interpretation and rational lifts. Their
square-root ideals add **two ordinary ideal-class2-torsion directions**
beyond the six generic strict images. Both explicit genus-one covers have
rational points, verified retrospectively using the already certified
point group and then checked independently with rational arithmetic.

The source classes were constructed before exceptional points were
admitted. The final solubility evaluation is explicitly oracle-assisted.
This distinction preserves the result's meaning: we constructed extra
classes without points, then tested the fixed output against known points.

This report supersedes the unresolved solubility statement in the
[first positive construction report](TWO_CONSTRUCTED_STRICT_CLASSES_AND_302.md).
The class-creation condition on t and the additional-class construction
on302 remain open.

## An actual elementary class-group factor

For each strict class beta, let `J_beta` be the unique fractional ideal
with `(beta)=J_beta^2`. The ordinary half-ideal map is well-defined modulo
squares. Its kernel consists of appropriate unit squareclasses; its image
lies in `Cl(K)[2]`. This is the same map used in the
[earlier half-ideal analysis](STRICT_HALF_IDEALS_AND_UNIT_KERNEL.md), but the
two new inputs here came from equation-only principal dependencies.

We constructed reduced representatives for the six generic strict half
ideals and the two additions. Generic factors use `(sqrt(Norm(gamma)),gamma)`
away from the discriminant, with exact corrections at every bad prime.
The constructed factors use the already certified complete ideal valuations.
Every ideal reduction retains and verifies its principal multiplier.

The [direct Artin calculation](../../artifacts/generated-results/elliptic-curves/rank_jump_constructed_class_direct_artin_v1.json)
evaluates ten generic ordinary unramified characters and the two new strict
characters on those eight ideals. At each prime-ideal factor P it uses
`chi_beta(P)=(beta,pi_P)_P`, with a verified local uniformizer. Since the
whole character is unramified, this is its Frobenius value, even when
individual factors in the expression for beta are not local units.

The character matrix has rank8; its six generic columns have rank6. The
[separate block verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_constructed_half_block_verification_v2.json)
checks an exact principal-square identity for every displayed ideal, proves
all27 rational primes in the Artin factorizations, and replays the29 prime
ideal records and eight complete factorizations. It constructs eight
combinations of the twelve characters whose evaluation
matrix is the identity. Consequently

\[
 \boxed{\operatorname{Cl}(K)\cong(\mathbf Z/2)^8\oplus C}
\]

for some uncomputed finite group C, with the first six displayed factors
coming from the generic strict classes and the last two from the newly
constructed classes. This does not compute the whole class group. It proves
that the two additions have new, non-2-divisible ideal-class information;
they cannot be explained solely by additional unit squareclasses.

On the two new ideals, the two new characters themselves give the matrix

\[
 \begin{pmatrix}0&1\\1&1\end{pmatrix}.
\]

Its determinant is1 over F2. **This is not the Cassels--Tate pairing.** In
particular its nonzero diagonal is compatible with both covers being
rationally soluble, as certified below.

The first Artin attempt required a cyclic odd-norm quotient with all
individual factors invertible, and completed only one column. Its optional
direction list also included negative entries, outside PARI's documented
nonnegative convention; no coverage or bounded-miss claim is based on those
trials. Half-ideal construction used default-direction reductions with
exact identities. The completed replacement factors the eight fixed
reduced ideals and uses local Hilbert symbols. It adds no new ideal targets
or reduction directions. The partial result is preserved.

## Explicit covers and exact rationality

The half-ideal reductions give `J_original=A*J_reduced`. Thus division by
`A^2` preserves the squareclass. Expanding the resulting fixed expression
gives representatives with maximum coefficient sizes2464 and1766bits,
respectively. Both satisfy the exact identity
`(beta_reduced)=J_reduced^2`.

The [compaction certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_constructed_class_compaction_v1.json)
records those representatives and the two quadrics of each cover. The
square-equivalence circuits are retained in the checkpoint archive.
Expansion uses PARI's documented
[factored-element and extended-ideal interfaces](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html#nffactorback).

As before, write

\[
 \beta(u+v\theta+w\theta^2)^2=Q_0+Q_1\theta+Q_2\theta^2.
\]

The cover is `Q2=0, Q1+s^2=0` in P3. The retrospective evaluator first
matches good-prime characters against the frozen22-point certificate.
It then adds the indicated existing points and verifies an actual square
root in K:

\[
 \boxed{\beta\xi^2=4x(P)-\theta.}
\]

Character agreement is only a candidate-matching step; this identity and
the point equations are the rational-solubility certificate.

The [oracle-labelled evaluation](../../artifacts/generated-results/elliptic-curves/rank_jump_constructed_class_oracle_solubility_v1.json)
contains each rational point P, the coefficients of xi, and a primitive
rational point on the corresponding cover. A
[Sage-free verifier](verify_constructed_cover_points.py) independently
checks the elliptic equation, the cubic-field identity using polynomial
arithmetic over exact fractions, and both cover quadrics. Both pass.

Thus the two classes are in the rational Kummer image, not nontrivial Sha
classes. This conclusion is unconditional and uses no complete Selmer
upper bound. It does not claim a newly discovered curve or point, nor
increase the reference's already known rank lower bound.

No smaller common simultaneous-solubility carrier for these two covers has
been constructed. Their rational points prove simultaneous solubility of
this fixed block; they do not supply a point-free condition on t.

## What the construction explains, and what it still does not

On this reference we now have the verified chain

\[
 \text{equation-only principal dependencies}
 \Longrightarrow \text{two new strict classes and half-ideal factors}
 \xRightarrow{\text{retrospective exact witnesses}}
 \text{two rational lifts}.
\]

The first arrow is **incidence**. The second is **solubility**, certified
with known-point evaluation inputs admitted after the constructor froze.
Neither statement is a visibility feature or a rank predictor.

The missing leading arrow is still

\[
 \text{an interpretable specialization condition on }t
 \Longrightarrow\text{the new principal dependencies / class factors}.
\]

The observed closing relation atoms are not a proposed cause of the rank
jump. Smooth norms made an existing arithmetic dependency computable.
The explicit nondivisible2-torsion block is the new mathematical object;
we have not proved which change in ramification or class-group structure
forces it to appear at this parameter.

## Transfer to302

The [maximal norm-form transfer](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_maximal_norm_form_v1.json)
first reproduces the reference's exact frozen norm form, then applies the
same maximal-order multiplication and Hessian reduction to302. The identity
`Norm(a*u+v*w)=a^2*F(m,n)`, including the integral basis and SL2 transport,
is checked coefficientwise. No principal-relation or point search is run.

| Constructor arithmetic | +6 reference | 302 |
|---|---:|---:|
| Cubic-field discriminant bits |237|496|
| Maximum reduced binary-cubic coefficient bits |65|156|
| Generic ordinary unramified character dimension |10|6|
| Generic strict dimension |6|0|

The first two rows describe computational scale only. They are not rank
features. The exact302 parent and its rank14 displayed quotient remain the
primary mechanism case. The larger norm-form coefficients do not justify
copying the reference's bounded relation budget and treating a miss as an
incidence result.

The next constructive requirement on302 is one new global principal
dependency yielding a nonzero strict class. Its generic strict dimension
zero makes the final independence target especially clear. Once such a
class exists explicitly, its half ideal and cover can be checked by the
same pipeline used here.

## Current ranking and replay

1. **Supported arithmetic mechanism:** additional elementary ordinary
   ideal-class2-torsion blocks carrying strict classes. Two new directions
   have now been constructed and rationally lifted on the reference.
2. **Unresolved creation mechanism:** the parameter-dependent ramification,
   principalization or integral-cover event that forces those new factors.
   The explicit block makes this a concrete comparison target.
3. **Separate solubility mechanism:** a global cover point, or an equivalent
   complete argument. Artin rank alone does not decide this; here exact
   retrospective points do.
4. **Computational lesson only:** small unresolved ideal representatives
   and appropriate character primes improve class construction. Their
   success is not evidence of a prospective rank selector.

The [new evidence archive](../../artifacts/generated-results/elliptic-curves/rank_jump_class_block_lift_evidence_v1.zip)
and [manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_class_block_lift_evidence_v1.json)
retain the bounded workers and square-equivalence circuits. Restore only
in an empty replay checkout. Cheap independent checks are:

```sh
sage -python elliptic-curves/rank-jump/verify_constructed_half_block_v2.py check
python3 elliptic-curves/rank-jump/verify_constructed_cover_points.py check
```

All changes are rank-jump-specific. Active search work and its mathematical
status entries remain untouched. Oracle-labelled outputs are forbidden as
prospective selector inputs.
