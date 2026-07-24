# Jelonek's coefficient-space corollary and the deformation project

Work over `C`.  For integers `n,d>=1`, let

\[
 {\cal V}(n,d)=
 \bigl(\mathbb C[x_1,\ldots,x_n]_{\leq d}\bigr)^n
\]

be the affine coefficient space of polynomial self-maps of degree at most
`d`, and let

\[
 X(n,d)=\{F\in{\cal V}(n,d):\det DF=1\}                 \tag{1}
\]

be the coefficient scheme cut out by the coefficients of `det DF-1`.
Statements about irreducible components below concern the irreducible
components of the underlying reduced variety.  The scheme structure becomes
essential in the tangent-space project in Section 3.

## 1. Jelonek's component theorem

Jelonek proves that the subset

\[
 {\cal A}(n,d)=
 \{F\in X(n,d):F\text{ is a polynomial automorphism}\}
\]

is Zariski closed in `X(n,d)`.  Consequently every irreducible component
`Omega` of `X(n,d)` satisfies exactly one of the following:

1. `Omega` is contained in `A(n,d)`; or
2. `Omega\setminus A(n,d)` is a nonempty dense open subset, so a general
   point of `Omega` is a counterexample to the Jacobian conjecture.

He also proves

\[
 \dim\Omega\geq n^2-1.                                 \tag{2}
\]

These statements are Jelonek's
[Lemma 2.1 and Theorem 2.2](https://arxiv.org/abs/2607.20597),
not new results of this repository.

## 2. Immediate corollary for the explicit maps

Let

\[
 F_N:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q},
 \qquad N\geq3,
\]

be any of the explicit determinant-one weighted maps constructed in
[the all-degree rational-fiber theorem](../verified/ALL_DEGREE_RATIONAL_FIBERS.md).
It has a complete fiber of `N` distinct points, hence its complexification is
not injective and is not a polynomial automorphism.

Fix any coefficient bound `d>=deg(F_N)`, and regard `F_N` as a point of
`X(3,d)`.  If `Omega` is any irreducible component of `X(3,d)` containing
that point, then `Omega` cannot be contained in the closed automorphism
locus.  Jelonek's theorem therefore gives:

> **Coefficient-space corollary (Jelonek).** Every irreducible component of
> `X(3,d)` containing `F_N` is generically noninvertible and has dimension at
> least
> \[
> 3^2-1=8.
> \]

More precisely, the nonautomorphisms form a dense open subset of every such
component.  The corollary supplies neither the number of components through
`F_N` nor their dimensions beyond the lower bound.  It also does not imply
that `F_N` is a smooth or reduced point.

## 3. The ordinary tangent space

The incremental coefficient-space project starts with the full scheme (1),
not with a weighted-support ansatz.  Put `F=F_N` and write a first-order
deformation as

\[
 F_\epsilon=F+\epsilon G,\qquad
 G\in{\cal V}(3,d),\qquad \epsilon^2=0.
\]

Linearizing the determinant gives

\[
 \det D(F+\epsilon G)
 =1+\epsilon\,\operatorname{tr}
   \bigl(\operatorname{adj}(DF)\,DG\bigr).              \tag{3}
\]

Thus the linearized Jacobian operator is

\[
 L_F:{\cal V}(3,d)\longrightarrow
 \mathbb C[x,y,z]_{\leq3(d-1)},\qquad
 G\longmapsto
 \operatorname{tr}\bigl(\operatorname{adj}(DF)\,DG\bigr), \tag{4}
\]

and

\[
 \boxed{T_FX(3,d)=\ker L_F,\qquad
 \dim T_FX(3,d)=3\binom{d+3}{3}-\operatorname{rank}L_F.} \tag{5}
\]

Equation (4), expanded monomial by monomial, is the requested system of
linearized Jacobian equations.  Since `det DF=1`, one may equivalently replace
`adj(DF)` by `(DF)^{-1}`; the adjugate form is preferable for exact sparse
coefficient arithmetic.

Jelonek's bound already implies

\[
 \dim T_FX(3,d)\geq8,                                  \tag{6}
\]

but equality in (6) would not by itself prove that there is a unique smooth
eight-dimensional component through `F`.

## 4. The unrestricted infinitesimal quotient collapses

Let `U` and `V` be polynomial vector fields on the target and source,
respectively.  The dual-number target and source automorphisms
`id+epsilon U` and `id+epsilon V` induce the first-order left--right
variation

\[
 \delta_{U,V}F=U\circ F+DF\cdot V.                     \tag{7}
\]

It belongs to the fixed-Jacobian tangent space precisely when

\[
 (\operatorname{div}U)\circ F+\operatorname{div}V=0,   \tag{8}
\]

and it belongs to the chosen coefficient box when its degree is at most
`d`.  Define the unrestricted infinitesimal left--right subspace by

\[
 {\cal O}^{LR}_{F,d}
 =
 \left\{
 U\circ F+DF\cdot V:
 \begin{array}{l}
 \deg(U\circ F+DF\cdot V)\leq d,\\
 (\operatorname{div}U)\circ F+\operatorname{div}V=0
 \end{array}
 \right\}
 \subseteq T_FX(3,d),                                  \tag{9}
\]

where `U,V` range over all polynomial vector fields.  No rank computation is
needed: the source part alone exhausts the tangent space.

Indeed, for any `G in T_FX(3,d)`, put

\[
 V=(DF)^{-1}G=\operatorname{adj}(DF)G.                 \tag{10}
\]

This is a polynomial vector field, and

\[
 F\circ(\operatorname{id}+\epsilon V)
 =F+\epsilon\,DF\cdot V
 =F+\epsilon G.                                        \tag{11}
\]

The chain rule, or the linearized determinant identity, gives
`div(V)=L_F(G)=0`.  Thus (11) is a determinant-one source automorphism over
the dual numbers.  Consequently

\[
 \boxed{{\cal O}^{LR}_{F,d}=T_FX(3,d),\qquad
 {\cal N}^{LR}_{F,d}=0.}                               \tag{12}
\]

This is the first-order form of the repository's
[formal orbit-triviality theorem](FORMAL_ORBIT_TRIVIALITY.md).  In
particular, the proposed unrestricted quotient

\[
 {\cal N}^{LR}_{F,d}
 =T_FX(3,d)/{\cal O}^{LR}_{F,d}
\]

cannot detect any of the stable parameters.

## 5. What ordinary coefficient deformation theory can detect

The collapse (12) sharply separates two questions.

The raw tangent space (5) remains useful.  Together with local equations it
can show that `F_N` is singular or nonreduced in `X(3,d)`, and primary
decomposition of the completed local ring can show that several components
pass through it.  Jelonek then says that every one of those components is
generically counterexample-producing and has dimension at least eight.

By contrast, the normalized weighted seed space has dimension `N-3`, and
decorated normalization detects an `(N-3)`-dimensional family of stable
left--right classes on a nonempty clean open.  Differentiating the family at
the chosen seed gives

\[
 T_{H_N}{\cal S}_N\longrightarrow T_{F_N}X(3,d),
 \qquad \dim T_{H_N}{\cal S}_N=N-3.                    \tag{13}
\]

Composing (13) with the unrestricted LR quotient is identically zero by
(12).  This is not a contradiction.  The dual-number source automorphism in
(11) need not algebraize to one polynomial automorphism family that
trivializes a positive-dimensional seed family, and its required complexity
can grow without a uniform bound.  Stable moduli live in precisely that
formal-to-algebraic and bounded-complexity gap.

Thus the meaningful comparisons are:

- the rank of (13) inside the raw coefficient tangent space;
- `dim T_{F_N}X(3,d)` versus the dimensions and tangent spaces of the local
  irreducible components;
- quadratic and higher local obstructions to raw coefficient directions;
- bounded-degree or Rees-filtered contact spaces before the cutoff is
  allowed to grow, followed by the algebraization/no-escape problem for an
  entire parameter family.

For a vector-field cutoff `b`, one may define

\[
 {\cal O}^{LR,\leq b}_{F,d}\subseteq T_FX(3,d)          \tag{14}
\]

by imposing `deg U,deg V<=b` in (9).  These spaces measure the complexity
needed to kill a tangent direction, but they are not stable-moduli tangent
spaces.  In fact (10) gives a finite cutoff that kills every fixed
first-order direction.  Detecting stable moduli therefore requires uniform
higher-order or family-level algebraization, not merely a larger
first-order matrix.

For the foundational case `N=3`, the stable-parameter count is zero.  The
existing
[sixteen-monomial coefficient calculation](FOUNDATIONAL_WEIGHTED_COEFFICIENT_SCHEME.md)
finds a nonreduced tangent direction in a much smaller normalized ansatz and
then proves that direction is polynomial source-orbit tangent.  It does not
compute (5) for the full bounded-degree scheme.

## 6. Revised exact computational deliverables

For each selected `F_N`, beginning with the foundational map and the minimal
bound `d=deg(F_N)`:

1. enumerate the coefficient variables of `G` and build the sparse matrix of
   `L_{F_N}` from (4);
2. compute its rank and a certified basis of `T_{F_N}X(3,d)`;
3. differentiate the `N-3` seed parameters inside that raw tangent space and
   record their canonical divergence-free source trivializers (10);
4. compute quadratic Kuranishi equations and, where feasible, the minimal
   primes of the completed local coefficient ideal;
5. only for the separate complexity question, compute the filtered spaces
   (14) and higher-order contact growth, then test whether the formal
   trivializers algebraize uniformly along the seed family.

Sparse ranks should first be audited modulo several good primes and then
certified over `Q`.  The output should distinguish ambient degree `d`,
geometric degree `N`, raw tangent dimension, local component dimensions,
filtered orbit dimensions, and higher-order obstruction ranks.  The
unrestricted first-order LR-normal dimension is always zero and should not
be reported as a candidate stable-moduli dimension.
