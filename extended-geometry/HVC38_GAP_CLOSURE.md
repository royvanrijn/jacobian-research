# Nonlinear-pivot and coordinated left-right closure at dimension 38

## Scope

The [dimension-38 cross-construction audit](HVC38_CROSS_CONSTRUCTION_FRONTIER.md)
left three nearby routes toward a dimension-36 quartic HN example:

1. nonlinear completion of the public eleven-variable \(d\)-pivot;
2. higher target-degree completion of the local twelve-variable \(z_8\)-pivot;
3. coordinated source and target automorphisms of the local map.

This note closes bounded versions of all three routes.  It does **not**
produce a dimension-36 example or prove a global lower bound.

## A common square obstruction

Let \(R\) denote the pullback of the target completion polynomial.

For the public lift, complete

\[
G_7=d-xy
\]

by a target polynomial in the six outputs whose pullbacks are independent
of \(d\).  On the completed pivot slice, \(d=xy-R\), and the fourth
component has the exact form

\[
G_4=a+(xy)^2-R^2.                                                \tag{1}
\]

For the local twelve-variable map, complete

\[
K_8=z_8+z_1z_2
\]

by a target polynomial in the eight outputs whose pullbacks are independent
of \(z_8\).  On that completed pivot slice,

\[
K_4=z_4+(z_1z_2)^2-R^2.                                        \tag{2}
\]

In either case, if the restricted map has degree at most three, the highest
homogeneous part of \(R^2\) first forces

\[
\deg R\le2.
\]

Unique factorization then forces the quadratic part to be

\[
R_2=\pm xy
\quad\text{or}\quad
R_2=\pm z_1z_2,                                                  \tag{3}
\]

respectively.  Thus the nonlinear completion problem becomes a linear
filtered-subalgebra calculation.

## Target degree through eight

For each target degree \(e=1,\ldots,8\), take every target monomial of
degrees one through \(e\) in the allowed outputs and map it to its source
terms of degree at least three.  Exact low-degree output coordinates give a
visible characteristic-zero kernel.  Reduction modulo the good prime

\[
p=1000003
\]

has the complementary full rank, so the modular and rational kernels have
the same dimension.

For the public \(d\)-pivot, the cumulative column counts and ranks are

\[
\begin{array}{c|rrrrrrrr}
e&1&2&3&4&5&6&7&8\\ \hline
\#\text{columns}&6&27&83&209&461&923&1715&3002\\
\operatorname{rank}&3&24&80&206&458&920&1712&2999.
\end{array}
\]

The kernel is always the three-dimensional span of outputs
\(G_8,G_{10},G_{11}\).  Their quadratic parts have rank three; adjoining
\(xy\) raises the rank to four.

For the local \(z_8\)-pivot,

\[
\begin{array}{c|rrrrrrrr}
e&1&2&3&4&5&6&7&8\\ \hline
\#\text{columns}&8&44&164&494&1286&3002&6434&12869\\
\operatorname{rank}&3&39&159&489&1281&2997&6429&12864.
\end{array}
\]

The kernel is always the five-dimensional span of outputs
\(K_7,K_9,K_{10},K_{11},K_{12}\).  Their quadratic parts have rank five;
adjoining \(z_1z_2\) raises the rank to six.

Together with (1)--(3), this excludes both nonlinear completed pivots for
every target polynomial of degree at most eight.

## Coordinated quadratic source and target directions

Now use the local map \(K\).  Simultaneously shear the jointly affine source
block

\[
(z_4,z_7,z_9,z_{10},z_{11})
\]

by arbitrary quadratic polynomials in

\[
(z_1,z_2,z_3,z_5,z_6,z_8,z_{12}).
\]

This gives \(5\binom{8}{2}=140\) source columns.  Add all elementary
quadratic target directions

\[
Y_i\longmapsto Y_i+\epsilon Y_jY_k,
\qquad j,k\ne i,
\]

giving another \(12\binom{12}{2}=792\) columns.  The combined
degree-at-least-four coefficient matrix therefore has 932 columns.

Modulo \(1000003\), its rank is \(896\), so its kernel has dimension
thirty-six.  Sparse echelon reduction gives thirty-six relations with
support at most five.  Rational reconstruction lifts every relation to
\(\mathbb Q\), and direct substitution verifies that all terms of degree at
least four cancel exactly.

On the cubic coefficient matrix of this thirty-six-parameter linearized
family, a fixed \(5\times5\) minor has determinant \(12\).  One exact Schur
entry is the constant

\[
-2.
\]

Consequently every degree-preserving linearized source-target combination
still has cubic-output rank at least six.

## Finite triangular integration

Twenty-two of the thirty-six sparse relations contain one elementary target
direction.  Integrate each by:

1. applying the corresponding simultaneous triangular source shear; then
2. applying the elementary triangular target shear to the precomposed map.

Seventeen relations integrate to degree-three one-parameter families.
For the remaining five, every high-degree coefficient has common parameter
factor \(\lambda^2\), so only \(\lambda=0\) preserves degree three.

The seventeen integrable directions can be combined.  Their target changes
modify only outputs \(1,2,3\), using products of outputs
\(7,8,9,11,12\); hence the combined target map is triangular.  The exact
degree-three locus of the resulting seventeen-parameter source-target
family is defined by 48 equations of parameter degrees two and three.

A selected \(6\times6\) cubic minor is

\[
-(t_3t_9-1)
\left(2t_0t_{13}+9t_3t_7-6t_3t_9+6\right).                     \tag{4}
\]

Over \(\mathbb Q\), the Gröbner basis of the 48 degree-three equations
together with (4) is

\[
(1).
\]

Therefore the degree-three locus never meets the zero locus of this minor.
Every member of the combined genuine triangular source-target family has
cubic-output rank at least six.

## Consequence and remaining gaps

No \(n+r\le17\) construction is obtained.  The exact new conclusions are:

- the public \(d\)-pivot is obstructed through target degree eight;
- the local \(z_8\)-pivot is obstructed through target degree eight;
- the full thirty-six-dimensional linearized quadratic left-right kernel
  has cubic rank at least six;
- the exact degree-three locus of the derived seventeen-parameter finite
  triangular family also has cubic rank at least six.

The nearest remaining routes are now materially broader:

- target completion of degree at least nine;
- a different jointly affine source block or quadratic source basis;
- nonlinear target generators rather than elementary quadratics;
- a non-nested state realization with tail score at most seventeen;
- a symmetric Hessian-nilpotent lift that avoids full dimension doubling.

## Verification

Run

```bash
make verify-hvc38-gap-closure
```

The checker reconstructs both maps, all bounded target monomials, the
good-prime ranks, all 932 source-target columns, the thirty-six sparse
relations, their rational lifts, the finite triangular integrations, and
the final characteristic-zero Gröbner basis.  The generated record is
[`hvc38_gap_closure.json`](../artifacts/generated-results/hvc38_gap_closure.json).

No Lean formalization, external review, global minimality statement, or
claim beyond the displayed families is made.
