# Affine-plane / Schubert reduction of the final HC4 bridge

## Status

This note continues `HC4-MR` after the Frobenius closure of the regular `[4]`
stratum.  It does **not** assume that the moving Jordan frame is affine
parallel.

> **Theorem HC4RSD77 — affine-plane middle foliation.**  Let
> \[
> S=\operatorname{Hess}\psi,\qquad
> T=\operatorname{Hess}A,\qquad
> N=S^{-1}T
> \]
> be a polynomial regular nilpotent `[4]` pencil with `det S` a nonzero
> constant, and assume the complete Jordan flag is Frobenius.  On the generic
> constant-rank locus put
> \[
> E_2=\ker N^2=\operatorname{im}N^2.
> \]
> Then `E_2` is autoparallel for the ambient flat affine connection:
> \[
> \nabla_XY\in E_2\qquad(X,Y\in E_2).
> \tag{0.1}
> \]
> Hence every connected leaf of `E_2` is an open subset of an affine
> two-plane in `A^4`.

> **Theorem HC4RSD78 — Schubert dichotomy.**  The direction map of these
> affine two-planes
> \[
> \sigma:B\dashrightarrow \operatorname{Gr}(2,4)
> \]
> from the two-dimensional leaf space has rank two in the genuinely moving
> branch.  Its differential is of common-image/upper-triangular type.  Thus
> its image is locally contained in one of the two standard two-dimensional
> Schubert leaves of `Gr(2,4)`:
>
> 1. all plane directions contain one fixed line; or
> 2. all plane directions lie in one fixed three-space.
>
> Case 1 gives a constant affine invariant and returns to the already-closed
> fixed-direction branch.  In Case 2, on any affine Grassmann chart the
> moving plane direction is represented by a row `(P,Q)` and the HC4
> transverse determinant identity gives
> \[
> J(P,Q)=c\ne0
> \tag{0.2}
> \]
> in local quotient coordinates.  Thus the only nonparallel middle-foliation
> geometry is a **rational plane-Keller quotient**.

The remaining globalization problem has therefore become very precise:
prove that the rational pair `(P,Q)` supplied by Case 2 extends polynomially
across the Grassmann-chart boundary, or show that failure of such extension
forces Case 1.  No local differential twisting remains beyond this pole
problem.

## 1. Why `E_2` is affine

Choose at a generic point an `S`-adapted Jordan frame

\[
Ne_1=0,\quad Ne_2=e_1,\quad Ne_3=e_2,\quad Ne_4=e_3
\]

with anti-diagonal `S`.  The exact linear Codazzi system for both `S` and
`T=SN`, together with Frobenius of `E_2=<e_1,e_2>`, forces

\[
\Gamma^3_{i j}=\Gamma^4_{i j}=0
\qquad(i,j\in\{1,2\}).
\tag{1.1}
\]

Equation (1.1) is tensorial as the vanishing of the affine second fundamental
form of the leaf.  Hence (0.1) holds on the generic locus.  Integral leaves
of an autoparallel distribution in affine space are affine subspaces.

The same first-order system also sharpens the previous twist count: the three
entries that were individually allowed are not independent.  In fact

\[
\Gamma^3_{4,1}=\Gamma^2_{3,1}.
\tag{1.2}
\]

Thus there are only two independent projective kernel-motion modes before the
Grassmann reduction.

## 2. Grassmann tangent geometry

For a two-plane `L` in a four-space,

\[
T_L\operatorname{Gr}(2,4)=\operatorname{Hom}(L,V/L).
\]

The HC4 Krylov calculation says that the two independent first derivatives of
`L=E_2` have a common one-dimensional image (equivalently, after dualizing,
a common one-dimensional kernel).  The rank-two integral manifolds of this
rank-one tangent cone are exactly the classical alpha/beta Schubert planes:

- planes containing a fixed line;
- planes contained in a fixed three-space.

This gives the two cases in HC4RSD78.

## 3. Keller pair in the moving Schubert chart

In the second Schubert case choose constant ambient coordinates so the fixed
three-space is `W`, and choose an affine Grassmann chart.  The direction plane
is represented by

\[
L(u)=\operatorname{graph} B(u),
\]

where only one row of the `2 x 2` matrix `B` moves.  Write that row as

\[
(P(u_1,u_2),Q(u_1,u_2)).
\]

The transverse projective derivative computed in `HC4RSD72` has determinant a
nonzero constant, while `HC4RSD74` kills its lower-left entry.  In these
coordinates this derivative is exactly

\[
D(P,Q).
\]

Therefore

\[
\det D(P,Q)=c\in K^*.
\]

So a nonparallel affine-plane foliation does not leave an arbitrary surface
problem: it produces a plane Keller pair on every Grassmann chart.

## 4. Exact remaining bridge

The previous master theorem asked for an abstract `affine-or-Keller` bridge.
HC4RSD77--78 identify it concretely:

\[
\boxed{
\text{parallel affine }2\text{-planes}
\quad\text{or}\quad
\text{a rational plane Keller pair from their direction map}.}
\]

The only unresolved issue is **polynomialization across Schubert-chart
boundaries**.  This is the global algebraic problem to attack next; trying to
force the whole Jordan frame affine-parallel is unnecessarily strong and is
formally false.
