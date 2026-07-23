# Controlled-divisor equivariant Keller suspensions

## Research assessment

As of 23 July 2026, the proposed lead is promising after a scope correction.
Its most defensible form is a **one-boundary incidence-extraction
conjecture**, followed by a separate **two-gauge rigidity problem**.

The main correction is that three notions currently compressed into
"one irreducible controlled divisor" must be separated:

1. the irreducible critical divisor of the two-dimensional invariant
   quotient;
2. the nonproperness hypersurface in the target;
3. the boundary prime in the canonical finite normalization
   \[
   \partial_F=
   \left(\operatorname{Norm}_{\mathbb A^3}k(\mathbb A^3)
   \setminus\mathbb A^3\right)_{\mathrm{red}}.
   \]

They coincide in the foundational cubic only after applying the marked-root
and reconstruction identifications.  A classification statement should say
which one is assumed irreducible and should include completeness of the
canonical boundary list.

There is also an important evidence correction.  The noncubic
root-engineered quadratic-gauge maps are not equivariant for the fixed
grading `(1,-1,-2)`: their higher seed decorations mix weights.  Only the
cubic seed is the foundational equivariant map.  The standard cancellation
family *is* equivariant, with source weights

\[
(m,-1,-m-1)
\]

and target weights `(-m-1,-1,m)`, but every noncubic standard cancellation
map has a second canonical boundary image `P=0`.  Thus an intrinsic
**exactly-one-boundary-prime** hypothesis excludes the obvious equivariant
non-incidence examples.  This is the real reason the one-boundary clause has
classification force.

## Verified literature input

Shaska's paper
[Graded Keller maps and the Jacobian Conjecture](https://arxiv.org/abs/2607.20210)
was submitted on 22 July 2026.  It proves over `C` that:

1. a `G_m`-equivariant Keller map with all source weights of one sign is an
   automorphism; the target degrees are automatically a permutation of the
   source weights;
2. every `G_m`-equivariant Keller map of `C^2`, for every nontrivial weight
   signature, is an automorphism;
3. for source weights `(1,-1,-2)` and target degrees `(-2,-1,1)`, put
   \[
   u=xy,\qquad v=x^2z.
   \]
   If the positive-weight target component is `G_3=x\Lambda(u,v)` and the
   quotient invariants are
   \[
   P=G_2G_3,\qquad Q=G_1G_3^2,
   \]
   then
   \[
   \det JG=\pm\Lambda^{-2}\operatorname{Jac}_{u,v}(P,Q).
   \]
   Hence the Keller condition is equivalent to
   \[
   \operatorname{Jac}_{u,v}(P,Q)=\kappa\Lambda^2.
   \tag{1}
   \]

Equation (1) is exactly the proposed quotient starting point.  It is a
reduction to a plane map with prescribed divisorial Jacobian, not yet a
classification of that plane map.

## The incidence lemma

For a parametrized plane curve `(X(S),Y(S))`, the marked-line incidence

\[
(S,B)\longmapsto\bigl(B,Y(S)-BX(S)\bigr)
\tag{2}
\]

has determinant

\[
-(Y'(S)-BX'(S)).
\tag{3}
\]

Conversely, a coordinate-preserving plane map is of the form (2), up to a
target shear, precisely when its second coordinate is affine-linear in the
preserved coordinate `B`.  Thus the proposed normal form is equivalent to a
nontrivial assertion:

> the intrinsic marked-root and one-boundary data force an
> affine-linear coefficient pencil on the quotient.

If a source chart has reciprocal Jacobian `D^{-1}` and

\[
D=D_0(P,S)-\kappa(P,S)Q,
\]

then choosing `X_S=lambda*kappa` and

\[
\beta=\frac{Y_S-\lambda D_0}{X_S},\qquad
B=Q+\beta,\qquad C=Y-BX
\]

gives

\[
Y_S-BX_S=\lambda D
\]

and a constant total determinant.  This proves the construction criterion.
It does **not** prove that an arbitrary one-boundary equivariant Keller map
admits such coordinates.

The two displayed gauges are:

\[
\begin{array}{c|c|c}
\text{gauge}&X(S)&\text{critical normalization}\\ \hline
\text{weighted tangent incidence}&S&\mathbb A^1\\
\text{quadratic cancellation gauge}&S^2&\mathbb G_m.
\end{array}
\]

The first is genuinely equivariant in all degrees.  The second is an
all-degree incidence construction but is equivariant for the fixed
three-weight action only at its foundational cubic seed.

## Recommended conjecture

The statement becomes falsifiable in the following form.

> **Equivariant one-boundary incidence conjecture.**  
> Let `F:C^3 -> C^3` be a nonproper `G_m`-equivariant Keller map for a
> linear mixed-sign action whose invariant quotient is a smooth affine
> plane.  Assume:
>
> 1. the generic inverse extension has a primitive marked root `S` whose
>    reconstruction open is the distinguished affine source;
> 2. the canonical finite-normalization boundary is complete and has one
>    geometrically integral prime, mapping birationally to an irreducible
>    target nonproperness divisor;
> 3. the boundary link is height-one saturated, its valuation is primitive,
>    and the order-two quotient Jacobian in (1) splits into two simple
>    ledger factors rather than one double factor;
> 4. the residue coordinate and transverse conormal class are saturated,
>    and there is no nontrivial target divisor ledger.
>
> Then, after equivariant triangular source and target changes, the quotient
> core is a marked-line incidence (2), and the orbit coordinate or source
> chart contributes the reciprocal factor to (3).  Ordinary stable
> left-right equivalence may be applied only after this equivariant normal
> form has been extracted.

The last sentence matters.  Arbitrary stable left-right equivalence need not
preserve the displayed torus action, so "equivariant classification up to
stable equivalence" should not use the grading as though it remained fixed.

A stronger second conjecture can then ask when `X` is equivalent to `S` or
`S^2`.  That is not a formal consequence of the incidence conjecture:
`X_S` is the coefficient of the marked variable in the controlled divisor,
so `deg X>=3` is exactly a higher-coefficient reciprocal chart.

## What is already proved internally

The repository contains two nearby scoped theorems.

1. The
   [one-boundary chart theorem](ONE_BOUNDARY_CHART_CLASSIFICATION.md)
   proves the polynomial/reciprocal sign dichotomy after height-one
   saturation and boundary monotonicity.  It also shows by countermodels
   that divisor minimality alone does not imply polynomial extension or
   determine tangential valuations.
2. The
   [strict quartic dichotomy](../extended-geometry/DEGREE_FOUR_MARKED_ROOT_CLASSIFICATION.md)
   proves that every straightened elementary quartic marked-root
   presentation is either a weighted quartic or standard cancellation type
   `(m,r)=(2,1)`.
3. The
   [degree-four incidence-suspension classification](../verified/INCIDENCE_SUSPENSION_DEGREE_FOUR_CLASSIFICATION.md)
   exhausts every root-preserving, `P`-fibration-preserving affine rechart of
   the quadratic reciprocal chart through `deg X<=4`.  The only controlled
   divisors are the `S,S^2,S^3` coefficients giving `X=S^2,S^3,S^4`;
   the latter two have unavoidable `Q/S` and `Q/S^2` pullback poles.  Hence
   only the quadratic gauge polynomializes in this bounded class.

None of these theorems extracts the suspension square from the bare canonical
finite normalization.  The missing implication is still:

\[
\text{intrinsic one-boundary marked cover}
\Longrightarrow
\text{coordinate-preserving, affine-linear quotient pencil}.
\tag{4}
\]

For the foundational geometric-degree-three case, the internal
[minimal-boundary program](MINIMAL_BOUNDARY_CLASSIFICATION.md) proves that
both known branches collapse to one polynomial left-right class once either
the suspension gateway or the cubic finite-normalization gateway is
available.

## Smallest test: weights `(1,-1,-2)`, coordinate degree at most seven

For target weights `(-2,-1,1)`, every equivariant map can be written on the
invariant chart as

\[
F=(x^{-2}A(u,v),x^{-1}B(u,v),xC(u,v)).
\]

Keeping **every** polynomial source monomial of total degree at most seven
gives support sizes

\[
\#A=10,\qquad \#B=9,\qquad \#C=7.
\]

After fixing the linear part `(z,y,2x)` and the two effective diagonal
gauges, 21 coefficients remain.  The invariant determinant equation gives
37 coefficient equations.

The exact checker
[`explore_equivariant_degree7_tangent.py`](../scripts/explore_equivariant_degree7_tangent.py)
finds at the foundational point:

\[
\operatorname{rank}(dI)=20,\qquad
\dim T=1.
\tag{5}
\]

All ten coefficients added beyond the published sixteen-monomial,
`z`-linear slice vanish in that unique tangent vector.  The tangent is the
same known dual-number direction and has a nonzero quadratic obstruction on
its literal affine line.  Thus the enlarged support produces no new
first-order branch at the foundational map.

This is a local necessary test, not component exhaustion.  A direct
21-variable Gröbner basis did not finish in the short exploratory run.
The published narrower support has the stronger exact result

\[
\mathcal O_{\mathrm{Keller,gauge}}
\simeq\mathbb Q[\epsilon]/(\epsilon^2),
\]

so its unique closed nonautomorphic point is the foundational map.

No known geometric-degree-four family lies in the coordinate-degree-seven
box:

- the displayed weighted quartic has profile `(12,11,4)`;
- the root-engineered quadratic quartic has profile `(7,26,24)`;
- standard cancellation type `(2,1)` already has
  `deg P=10` and `deg Q=8`.

Consequently a degree-four equivariant Keller map in this box would be a
genuinely new class, not a small representative of a known quartic family.

## Sharp failure modes

The first counterexample search should target the following in order.

1. **Higher-power one-place core.**  A linear section with quotient
   Jacobian `D^3` has inverse degree four and ramification index four.  Its
   plane core exists; no polynomial one-boundary Keller suspension is known.
2. **Higher incidence coordinate.**  A reciprocal chart with coefficient
   `kappa(S)` of degree at least two gives `deg X>=3`.
3. **Unsaturated two-place marking.**  On a `G_m` critical normalization,
   a nonmonomial Laurent residue can evade the primitive `Y=Q-sP` lift.
4. **Nontrivial target ledger.**  A target divisor may absorb part of the
   quotient Jacobian and invalidate the simple reciprocal reconstruction.
5. **Additional reconstruction divisor.**  This is the proposed
   two-independent-divisor failure, but it lies outside an exact
   one-boundary theorem and should be treated as the next classification
   stratum rather than its first internal counterexample.

The first two are the smallest genuine internal failures of the corrected
conjecture.

## Proof and computation program

1. Linearize the source and target torus actions and compute the invariant
   quotient.  Restrict initially to weights for which both invariant rings
   are polynomial.
2. Apply the quotient determinant identity (1), factor its divisor, and
   compare the two copies of the controlled prime with the canonical
   boundary valuation.
3. Prove height-one saturation, boundary monotonicity, and residue/conormal
   saturation from the marked reconstruction open.
4. Extract a preserved quotient coordinate.  This is the substantive
   two-variable classification step; without it the line-incidence form is
   not justified.
5. Prove affine-linearity in that coordinate and apply the incidence lemma.
6. Reconstruct the orbit coordinate from the primitive boundary valuation
   and eliminate any target ledger.
7. Only then quotient by triangular gauge and compare stable boundary
   invariants.
8. Computationally, finish the 21-variable degree-seven open chart by
   modular Gröbner bases and rational reconstruction.  Then test its
   remaining boundary charts and compute geometric degree from the quotient
   function-field extension.

The degree-three target is realistic: the full-support tangent calculation,
the existing cubic marking extraction, and the collapse of both known cubic
branches all point to rigidity.  Degree four is the correct discovery test,
but the coordinate-degree-seven cutoff currently contains no known positive
example.
