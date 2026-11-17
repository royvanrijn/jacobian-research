# Coordinate obstructions for low degree Keller collisions

Three simple necessary conditions explain the small sparse censuses. The
first two work in arbitrary degree. The third combines an elementary
degree-controlled reduction with Moh's plane theorem.

## A univariate coordinate

If Fi is a nonlinear polynomial in xi alone, then its nonconstant derivative
divides det JF. A nonconstant polynomial cannot multiply another polynomial
to a nonzero constant. Thus such an Fi is impossible in a Keller map in
characteristic zero.

This is a necessary condition on a completed support. A partial univariate
row can still acquire mixed terms; discarding that partial branch before
accounting for possible additions would be unsound.

## A hyperplane normal derivative

Suppose Fi is divisible by xi. On xi=0 every off-diagonal derivative in row i
vanishes, while its normal derivative is

\[
N=(F_i/x_i)|_{x_i=0}.
\]

Restricting the determinant identity to xi=0 gives N times a polynomial
cofactor equal to a nonzero constant. Consequently N must be a nonzero
constant.

For a normalized exact support, if all nonlinear terms of Fi contain xi
and at least one has xi-exponent exactly one, N is nonconstant: the remaining
source monomials are distinct and their coefficients are nonzero. This
closes that support in every degree. The sole non-shear eight-term residual
has F2=y(1+dx) and violates this condition.

The condition does not exclude xi+xi^2 H merely because it is reducible.
There its normal derivative on xi=0 is still one.

## A constant nonzero Jacobian pairing

Let F be a three-variable Keller map with maximum ordinary coordinate
degree D. Suppose constant vectors v,w satisfy

\[
v^T JF(X)w=c\ne0
\]

identically in X. Then F has no collision whenever D<=10.

Choose linear source coordinates whose first direction is w and linear
target coordinates whose first output is v^T F/c. In those coordinates
the first output has the form

\[
x_1+H(x_2,x_3),\qquad \deg H\le D.
\]

The triangular source substitution x1=t-H(x2,x3) makes that output equal
to t. The substitution has a polynomial inverse and determinant one.
The transformed map is

\[
(t,P(t,x_2,x_3),Q(t,x_2,x_3)),
\]

up to coordinate order, with nonzero constant determinant. For each fixed
t, the last two outputs form a plane Keller map. If delta=max(1,deg H),
their degrees are at most D*delta<=D^2.

Two colliding source points must lie in the same t-fibre after these
invertible transformations. But D<=10 gives D^2<=100, and
[Moh's degree-100 plane theorem](https://www.math.purdue.edu/~ttm/jacobian)
makes that plane map injective. This is a contradiction. Moh's primary
author account explicitly states the maximum-degree-at-most-100 range
for both plane components; the publication is
T.T. Moh, Journal fuer die reine und angewandte Mathematik 340 (1983),
140-212.

The argument applies over every characteristic-zero field: all map,
linear-transformation and collision coefficients lie in a finitely
generated subfield, which embeds in C. Invertible source changes preserve
distinctness of the colliding points. No global polynomial inverse for a
family of plane fibres is needed.

This yields a support-independent exclusion: in raw degree at most ten,
a three-variable collision cannot have even a linear target combination
that is affine with constant nonzero slope along a constant source
direction. In particular a coordinate xi+H(other variables), including
a fixed identity coordinate, is excluded.

The pairing property is invariant under invertible affine source and
target changes. This gate is therefore stronger than checking three
displayed diagonal derivative entries.

## Exact direction incidence for normalized maps

When JF(0)=I, write

\[
JF(X)=I+\sum_{\alpha\ne0}M_\alpha X^\alpha.
\]

A constant nonzero pairing exists exactly when the finite bilinear system

\[
v^TM_\alpha w=0\quad(\alpha\ne0),\qquad v^Tw=1
\]

has a solution. The last normalization is legitimate by rescaling w.
For a degree-at-most-ten Keller collision this incidence system must be
empty over the algebraic closure.

The [control checker](verify_coordinate_gates.py) builds this complete
incidence for the known degree-seven counterexample after tangent
normalization G=(F3/2,F2,F1). Both Singular and SymPy return the unit ideal
for its 27 equations, as required. It also verifies a nonlinear shear
control by an exact inverse substitution and Jacobian chain rule.

A point in the direction incidence would exclude a low-degree Keller
candidate. Emptiness alone would not construct a collision.

## Sharp dimension four for this bounded special form

Among Keller counterexamples of maximum coordinate degree at most ten
that admit a constant nonzero Jacobian pairing, the smallest dimension
is **four**.

Dimension one is affine. In dimension two, straightening the affine-shear
coordinate leaves a one-variable map with nonzero constant derivative,
so it is injective without any degree restriction. Dimension three is
excluded by the preceding D^2 bound and Moh's plane theorem.

For existence in dimension four, take the known
[foundational map](../../foundations/verified/FOUNDATIONAL_GEOMETRY.md),
normalize it as G=(F3/2,F2,F1), and append an identity coordinate:

\[
\widehat G(x,y,z,t)=(G(x,y,z),t).
\]

This has determinant one, degree seven, and the inherited distinct
collision points with t=0. The pairing v=w=e4 is constantly one. The
[exact control receipt](coordinate-gates-validation.json) verifies its
determinant, collision and coordinate degrees (4,6,7,1).

Thus dimension four is attained in this explicitly bounded special class.
For unbounded degree, the three-variable obstruction would require an
unrestricted plane theorem; that question remains outside this result.

