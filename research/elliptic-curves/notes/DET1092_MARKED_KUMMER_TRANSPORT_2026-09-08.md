# Marked elliptic Kummer transport and two norm obstructions

## Result and scope

**Verified application.** The first302 unlock has a nonzero elliptic
2-Kummer class modulo the inherited MW17 image. All nine source-selected
generic-point controls have zero relative elliptic class: their marked
points are exactly generic section0. Nevertheless the corresponding
Jacobian classes are independent of their rank17 Picard images in every
case, by the completed
[first-witness](DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md) and
[control-panel](DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md) proofs.

**New deduction.** Two natural attempts to transport a class without first
resolving the marked point fail exactly:

1. Symmetrizing the two branches of an RR fibre returns the known centre's
   elliptic Kummer class, hence zero modulo the inherited image.
2. Each of the ten sextic descent fields has no proper intermediate field.
   In particular it has no elliptic cubic subfield on which to define a
   relative field norm. The elementary tensor-norm replacement is trivial
   on norm-square representatives.

The subsequent [torsion-module obstruction](#broader-obstruction-irreducible-2-torsion)
also rules out nonzero Jacobian-to-elliptic homomorphisms and nonzero
torsion-module-induced Selmer transfer on these ten members. Other auxiliary
algebras, nonlinear predicates and more general correspondences are not
excluded. No map from an unevaluated genus-two Selmer class to an exceptional
elliptic direction is constructed here.

## Exact marked map

**Verified application.** The literal parent model is

\[
E_t:y^2+xy+y=x^3+x^2+a_4(t)x+a_6(t).
\]

Set `X=4x`, `Y=8y+4x+4`. Then

\[
Y^2=F_t(X)=X^3+5X^2+(16a_4(t)+8)X+64a_6(t)+16.
\]

The generic coefficient functions, exact centre, all ten marked points,
and all polynomial coefficients are retained in the
[construction certificate](../../artifacts/generated-results/elliptic-curves/det1092_marked_kummer_transport_v1.json).
There is no new point search. The first point is reconstructed by the
completed frozen-chart replay; the other nine use the already frozen
source-only selector `u=r0(t), v=0`.

**Established literature.** For the smooth odd-degree model, the usual
2-Kummer representative is `beta_P=X_P-theta` in
`K_t=Q[theta]/F_t(theta)`, with norm `Y_P^2`. Even-degree sextic Jacobian
descent instead has a scalar quotient and may have a kernel; these two
objects must not be identified. See
[Poonen--Schaefer, introduction and sections5,11--13](https://math.mit.edu/~poonen/papers/descent.pdf).

**Verified application.** Write the generic centre as `Z=P_w` to distinguish
it from the genus-two curve. The RR line `f0+f1*x+f2*y=0` becomes

\[
L(X)=\left(1-\frac{2f_1}{f_2}\right)X+4-\frac{8f_0}{f_2}.
\]

It passes through `-Z`, `P`, and `Q=Z-P`. Let

\[
G(X)=\frac{F_t(X)-L(X)^2}{X-X_Z}
     =(X-X_P)(X-X_Q).
\]

All point identities, line incidences and norms are checked exactly. The
independent replay computes norms as determinants of multiplication matrices,
not by the constructor's resultants. Nonzero denominators, nonzero ordinates,
smooth elliptic cubics and squarefree sextics are explicitly verified on all
ten cases. Outside this patch a different representative is required; no
undefined value is treated as a zero class.

## Why the unordered pair loses the direction

**New deduction, symbolic verification.** Evaluation at `theta` gives

\[
\boxed{(X_Z-\theta)G(\theta)=L(\theta)^2},\qquad
\beta_P\beta_Q\beta_Z=L(\theta)^2.
\]

Thus `delta(P)+delta(Q)=delta(Z)`, already inherited. The pair's relative
class is zero. The two individual branches have the *same* relative class,
since `Q=Z-P` and signs disappear modulo2. Therefore a rational split branch
does give an unambiguous relative elliptic class, but the symmetric product
does not compute that class before splitting.

The first relative class is separated from MW17 by the previously certified
finite-group characters (matrix ranks17 and18). Its branch orientation is
checked to agree with the original witness or its centre-complement.
The old strict-filtration character annihilates `MW17+K_strict` and evaluates
to1 on this class. It is therefore **not strict modulo MW17 for that frozen
elliptic filtration**. This is not an identification with an MW16 ideal class
or an independently defined Jacobian strict kernel.

## Why a direct sextic-to-cubic norm is unavailable

**Verified application.** Put `A=Q[T]/q(T;u,0)`. At a squarefree good
reduction, factor degrees give a Frobenius cycle type. The retained witnesses
are:

| Marked case | Sextic six-cycle prime | Sextic 1+5 prime | Irreducible cubic prime |
|---|---:|---:|---:|
|302 first unlock|83|179|31|
|scale-0131232|17|109|23|
|scale-0257585|53|43|23|
|scale-0487239|71|23|23|
|scale-0177036|17|131|23|
|scale-0043332|17|53|29|
|scale-0590501|23|113|17|
|scale-0290097|17|107|29|
|scale-0748009|67|23|67|
|302 generic-section control|79|109|31|

**New deduction.** The irreducible sextic reduction proves `A` is a field
and the Galois action on its six embeddings is transitive. A nontrivial
block system must have two blocks of size3 or three blocks of size2. A
5-cycle would act trivially on the set of blocks, but its length5 orbit
cannot fit in a block of size2 or3. Hence the action is primitive. By the
subgroup/block correspondence, `A/Q` has no proper intermediate fields.
The checker also exhausts the15 pair partitions and10 triple partitions
of six labels and verifies that none is fixed by a5-cycle.

Every displayed elliptic cubic is independently irreducible, so `K_t` is a
degree3 field. It cannot embed in `A`. Thus `Norm_(A/K_t)` is **not defined**;
this is not a failed numerical class-group calculation.

**New deduction.** Base-extending the multiplication matrix gives

\[
\operatorname{Norm}_{A\otimes K_t/K_t}(\beta)
 =\operatorname{Norm}_{A/\mathbf Q}(\beta).
\]

For the norm-square representatives of rational Jacobian classes the right
side is a rational square, so this route produces a trivial cubic
squareclass. The claim concerns this elementary tensor norm, not all maps
that might use a larger algebra or additional marked data.

## Broader obstruction: irreducible 2-torsion

**New deduction, independently verified application.** All ten Jacobians in
the table are simple over `Q`. For every elliptic curve `E/Q`,

\[
\operatorname{Hom}_{G_{\mathbf Q}}(J[2],E[2])=0,
\qquad \operatorname{Hom}_{\mathbf Q}(J,E)=0.
\]

Both statements also hold in the reverse direction. In particular, none
of these ten RR curves admits a nonconstant map over `Q` to an elliptic
curve, of any degree. This strengthens the earlier degree-two-quotient
obstruction for the first witness; it is not an assertion of absolute
simplicity for these ten members.

**Established literature.** For a squarefree sextic, `J[2]` is the space of
even subsets of its six branch points modulo complementation, compatibly
with Galois action. This is the `p=2,d=6` case of
[Poonen--Schaefer, section6, Proposition6.2 and its concluding identification](https://math.mit.edu/~poonen/papers/descent.pdf).
Its dimension over `F2` is four.

**New deduction.** The previously certified `1+5` reduction supplies an
element with branch permutation `(01234)(5)`. Its action on `J[2]` has
characteristic polynomial

\[
\Phi_5(X)=X^4+X^3+X^2+X+1,
\]

which is irreducible over `F2`. Thus `J[2]` is an irreducible Galois
module. An independent combinatorial verification finds orbit lengths
`1,5,5,5` on its16 elements. A proper nonzero invariant subspace would
have1,3 or7 nonzero vectors, none a union of the length5 nonzero orbits.

Any nonzero equivariant map from this irreducible four-dimensional module
to a two-dimensional elliptic2-torsion module would be injective, which is
impossible. In the other direction its nonzero image would be a proper
invariant subspace. The checker independently tests all six possible
`GL2(F2)` actions and all256 linear maps in each direction: in each case
only the zero map intertwines the actions.

If `J` contained an elliptic subvariety over `Q`, its2-torsion would be a
Galois-stable two-dimensional subspace. Hence `J` is `Q`-simple. A nonzero
homomorphism `J -> E` would have an elliptic connected kernel, and a nonzero
homomorphism `E -> J` would have an elliptic image, both impossible. With
the inherited rational basepoint, a map from the RR curve to an elliptic
curve factors through its Jacobian after subtracting its basepoint image.
This proves the asserted curve-map obstruction.

**New deduction: a necessary constant-field degree.** Let `L/Q` be any
finite number field with degree not divisible by5. The image of `G_L` in
the finite `J[2]` Galois image has index dividing `[L:Q]`. Its order is
therefore still divisible by5, so it contains an element of order5. Any
nonidentity order5 action in dimension4 over `F2` has minimal polynomial
`Phi5` and is irreducible. Thus the same simplicity and Hom obstructions
hold over `L` as well. An extension allowing an elliptic factor must have
degree divisible by5; that condition is **necessary, not sufficient**.
Quadratic, cubic and quartic constant extensions cannot repair this route.

**Scope.** The zero torsion-module Hom group means there is no nonzero
map on `H1` or Selmer groups *induced by a coefficient-module homomorphism*
`J[2] -> E[2]`. It does not say that no abstract linear map can be written
between finite Selmer vector spaces. It does not compute either Selmer
group or exclude nonlinear incidence conditions. The marked RR point
still maps to the varying elliptic surface; that is not a morphism of the
whole RR curve to one fixed elliptic fibre.

The [construction](../../artifacts/generated-results/elliptic-curves/det1092_rr_torsion_transport_gate_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_torsion_transport_gate_replay_v1.json)
reuse exactly the ten saved `1+5` witnesses. The replay verifies the
degree-five factors using Frobenius-power irreducibility tests rather than
the constructor's factorization routine, then uses bit subsets and an
exhaustive3072-map linear check rather than the constructor's constraint
matrices. Both Sage10.9 jobs completed in under a second, with separate
25-second caps and no new prime or point search.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rr_torsion_transport_gate.sage
```

## Rational classes, Selmer classes, and what remains unknown

**Established interpretation and verified application.** The inherited
basepoint has Abel--Jacobi class zero. The known first witness and the nine
generic-point controls give genuine rational Jacobian classes outside the
Picard-image rational span, not merely locally surviving candidates. Their
true Kummer classes lie in the rational image inside `Sel_2(J)` and map to
zero in `Sha(J)[2]`. Full Selmer groups and other Sha classes remain
**NOT_COMPUTED**. A nonzero relative Jacobian class is not synonymous with
a nonzero Sha class or an extra marked elliptic point.

The field and pair calculations explain why these particular class-transfer
shortcuts do not repair the failed specificity test. A predictive invariant
available before rational branch recovery remains **unconstructed**; no
positive incidence conjecture is promoted to a theorem.

## Reproduction and limits

```sh
sage -python research/elliptic-curves/cas/verify_det1092_marked_kummer_transport.sage
```

The [independent replay certificate](../../artifacts/generated-results/elliptic-curves/det1092_marked_kummer_transport_replay_v1.json)
binds all inputs and the checker hash. Sage10.9 completed construction and
independent replay in about0.6seconds each, under separate25-second caps.
There are ten already certified marked cases, one symbolic identity, and
44 retained cubic-prime attempts within the existing64-prime list, stopping
at the first irreducible cubic in each case. All20 sextic gates reuse earlier
recorded witnesses and are recomputed exactly. No full Galois-group or
class-group computation is needed. Passing output is immutable; a failed
assertion or timeout cannot generate a passing certificate.

No point search, global Selmer calculation, parameter sweep, or V3/V4 pilot
mutation was performed. The old independence and strict certificates are
consumed as proof dependencies, not rediscovered by another point campaign.
