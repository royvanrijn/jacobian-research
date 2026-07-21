# Foundational facts

## Exact facts checked in this repository

Put `u=1+xy` and write `F=(a,b,c)`, with

\[
a=u^3z+y^2u(4+3xy),\quad
b=y+3xu^2z+3xy^2(4+3xy),\quad
c=2x-3x^2y-x^3z.
\]

Exact symbolic expansion gives `det DF = -2`. Thus `F` is étale at every
finite complex point. Exact rational substitution also gives

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

Consequently `F` is a polynomial Keller map but is not injective and hence not
a polynomial automorphism. Subject only to checking the displayed arithmetic,
this is a counterexample to the usual Jacobian conjecture in dimension 3.
Appending identity coordinates gives counterexamples in every dimension
`n >= 3`. Multiplying one output by `-1/2` normalizes the determinant to 1.

Ordinary total degrees of the three coordinates are `(7,6,4)`.

## Three-minimum polynomial and escape at infinity

For the collision value `q=(-1/4,0,0)`, define

\[
\mathcal L(x,y,z)=(4F_1+1)^2+16F_2^2+16F_3^2
                  =16\lVert F-q\rVert^2.
\]

This is a degree-14 polynomial with integer coefficients, and

\[
\nabla\mathcal L=32\,DF^T(F-q).
\]

Since `DF` is invertible everywhere, its critical set is exactly `F^{-1}(q)`.
That fiber consists of the three points displayed above.  At each of them,

\[
\nabla^2\mathcal L=32\,DF^TDF>0.
\]

Thus the three points are nondegenerate global minima and there are no other
finite critical points.

The failure of compactness has the following explicit certificate.  For
`R != 0`, put

\[
X_R=\left(R,-\frac1R,\frac5{R^2}\right).
\]

Exact substitution gives

\[
F(X_R)=\left(0,\frac2R,0\right),\qquad
\mathcal L(X_R)=1+\frac{64}{R^2},
\]

and

\[
\nabla\mathcal L(X_R)
=\left(-\frac{392}{R^3},\frac{264}{R},0\right).
\]

Consequently, as `|R| -> infinity`, the points `X_R` escape to infinity while
`\mathcal L(X_R) -> 1` and `\nabla\mathcal L(X_R) -> 0`.  This is an explicit
Palais--Smale sequence at level 1.

## Function-field model

Let `K=C(a,b,c)` and `L=C(x,y,z)`. On `x != 0`, set

\[
t=y+1/x,\qquad P(T)=cT^3-2T^2+bT-2a.
\]

Direct identities give

\[
b=4t+2/x-3ct^2,\qquad 2a=ct^3-2t^2+bt,
\qquad P'(t)=3ct^2-4t+b=2/x.
\]

If `r=P'(t)`, the source is recovered rationally by

\[
x=2/r,\quad y=t-r/2,\quad
z=5r^2/4-3tr/2-cr^3/8.
\]

Thus `L=K(t)` and the extension has degree at most 3. The explicit three-point
fiber plus étaleness implies nearby targets have at least three distinct
preimages, so the generic degree is exactly 3. Therefore `P` is the minimal
polynomial after monic normalization.

The discriminant is

\[
\operatorname{Disc}_T(P)=-4Q,
\quad Q=27a^2c^2-18abc+16a+b^3c-b^2.
\]

The available same-day structural analysis argues that `Q` is nonsquare and
irreducible, giving Galois closure group `S_3` and a non-normal cubic extension.
The repository verifies the discriminant identity; irreducibility over the
function field is recorded as a proof argument, not delegated to a finite test.

At the collision target, `c=0` and the cubic degenerates to
`-2T^2+1/2`. Two roots are `t=+/-1/2`; the third is at infinity and corresponds
to the preimage with `x=0`. This is not finite ramification: `P'(t)=2/x` shows
that a repeated root forces `x` to escape to infinity.

The full calculation is recorded in `IMAGE_AND_NONPROPERNESS.md`.  In
particular, with
`Gamma=V(3bc-4,12a-b^2)`, the exact image is `C^3 minus Gamma`, the fiber
cardinalities are respectively `3,1,0` on `Q != 0`, `Q=0` off `Gamma`, and
`Gamma`, and the nonproperness set is exactly `V(Q)`.  Moreover `Gamma` is the
singular locus of `V(Q)`.

## Explicit normal-form reductions

The standard reductions have been executed explicitly:
`CUBIC_HOMOGENEOUS_REDUCTION.md` gives a 95-dimensional map `I+H` with `H`
cubic homogeneous, and `CUBIC_LINEAR_REDUCTION.md` gives a 451-dimensional
Druzkowski map `X-(AX)^{*3}` with `rank(A)=95`. Both have determinant one and
stored exact rational collisions. These are explicit normal forms of the same
counterexample mechanism, not logically independent counterexamples.

## Audited downstream conjectures

By contraposition of the cited implication/equivalence theorems, the verified
3D counterexample also implies: the Mathieu conjecture fails for `SU(3)`; the
Gaussian Moments Conjecture fails in at least one finite dimension; Zhao's
quartic homogeneous Vanishing Conjecture fails in at least one finite
dimension; and the all-dimensional Image Conjecture fails.  C15 now supplies
an explicit 190-variable quartic HN witness, an explicit Special
Image-Conjecture pair, and a third-Weyl-algebra Dixmier witness.  Least
dimensions and the least exceptional Vanishing exponent remain unknown. See
`C15_INDEPENDENT_AUDIT.md` and `DIRECT_CONSEQUENCES.md`.
