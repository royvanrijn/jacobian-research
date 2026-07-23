# Natural torus Hamiltonian reduction of the foundational lift

The foundational Keller map has a canonical algebraic torus symmetry.  This
note determines exactly what its Hamiltonian reduction does and isolates the
remaining rank-two Dixmier problem.  The result improves both the naive
positive expectation and the naive singular-quotient obstruction:

* the full affine moment reduction is not affine four-space;
* a primitive cotangent chart removes the stabilizer and retains the complete
  collision;
* the target quotient on that chart is standard affine four-space;
* the source quotient is an explicit smooth affine fourfold whose polynomial
  triviality is open.

## 1. Equivariance and moment maps

Write the foundational map as

\[
 F=(A,B,C):\mathbb A^3_{x,y,z}\longrightarrow\mathbb A^3.
\]

It is equivariant for

\[
 t\cdot(x,y,z)=(t^{-1}x,ty,t^2z)
\]

on the source and

\[
 t\cdot(A,B,C)=(t^2A,tB,t^{-1}C)
\]

on the target.  On the two cotangent bundles the moment maps are

\[
\begin{aligned}
 \mu_{\rm src}&=-xp_x+yp_y+2zp_z,\\
 \mu_{\rm tgt}&=2A\alpha+B\beta-C\gamma.
\end{aligned}
\tag{1}
\]

If \(\widehat F(x,p)=(F(x),DF(x)^{-T}p)\), equivariance gives the exact
identity

\[
 \widehat F^*\mu_{\rm tgt}=\mu_{\rm src}.             \tag{2}
\]

Thus every level admits a classical Hamiltonian reduction.

## 2. Why the full affine quotient is not \(\mathbb A^4\)

For every value \(c\), the target level contains

\[
 (A,B,C;\alpha,\beta,\gamma)
 =(-1/4,0,0;-2c,0,0).                                \tag{3}
\]

Only coordinates of weights \(2\) and \(-2\) are nonzero, so its stabilizer
is \(\mu_2\).  The transverse symplectic slice consists of the four
weight-odd coordinates \(B,C,\beta,\gamma\), on which \(\mu_2\) acts by
minus one.  The coarse quotient is locally

\[
 \mathbb A^4/\{\pm1\}.                                \tag{4}
\]

Its maximal ideal has ten independent quadratic generators modulo its
square, whereas its dimension is four.  Hence it is singular.  In
particular, no full affine moment reduction is polynomial affine four-space.

This obstruction alone is not enough to discard the torus route: (3) lies on
the nonprimitive momentum stratum \(\beta=0\).

## 3. The primitive momentum chart

The target cotangent coordinate \(\beta\) has weight \(-1\).  On
\(\beta\ne0\), every orbit has a unique representative with \(\beta=1\).
The moment equation becomes

\[
 B=c-2A\alpha+C\gamma.                                \tag{5}
\]

Consequently the target reduction is exactly

\[
 \mathbb A^4_{A,C,\alpha,\gamma}.
\]

Before fixing the gauge, invariant coordinates are

\[
 A\beta^2,\qquad C/\beta,\qquad
 \alpha/\beta^2,\qquad\gamma\beta.                    \tag{6}
\]

Their reduced brackets give the two standard canonical pairs
\((A,\alpha)\) and \((C,\gamma)\).

For the source, use the polynomial momentum coordinates

\[
 p=DF(x)^T(\alpha,\beta,\gamma)^T.
\]

The inverse-Jacobian matrix has constant determinant, so this is a polynomial
coordinate change in the momentum variables.  Equations (1)--(2), followed
by the gauge \(\beta=1\), identify the source quotient with the smooth affine
hypersurface

\[
\boxed{
 X_c=\left\{
 B(x,y,z)+2\alpha A(x,y,z)-\gamma C(x,y,z)=c
 \right\}
 \subset\mathbb A^5_{x,y,z,\alpha,\gamma}.}
\tag{7}
\]

The reduced map is

\[
 \pi_c:X_c\longrightarrow\mathbb A^4,\qquad
 (x,y,z,\alpha,\gamma)\longmapsto
 (A(x,y,z),C(x,y,z),\alpha,\gamma).                  \tag{8}
\]

It is the free Hamiltonian reduction of the cotangent lift on this chart,
so it is an etale symplectic morphism between smooth fourfolds.

## 4. The collision survives without orbit identification

The foundational fiber over \((-1/4,0,0)\) is

\[
\begin{aligned}
 b_0&=(0,0,-1/4),\\
 b_+&=(1,-3/2,13/2),\\
 b_-&=(-1,3/2,13/2).
\end{aligned}
\tag{9}
\]

Set

\[
 (\alpha,\beta,\gamma)=(-2c,1,0).
\]

All three points satisfy (7), and (8) sends them to

\[
 (-1/4,0,-2c,0).                                     \tag{10}
\]

They are distinct quotient points.  Indeed, \(\beta\) has primitive weight
\(-1\), so the slice \(\beta=1\) meets every orbit exactly once.  The
transformation \(t=-1\), which exchanges \(b_+\) and \(b_-\) on the zero-odd-
momentum stratum, changes \(\beta=1\) to \(\beta=-1\) and therefore does not
identify the two displayed gauge representatives.

Thus the natural torus simultaneously supplies a free reduction chart and a
complete three-point reduced collision.

## 5. The localized target quantum reduction is \(A_2\)

Let the target Weyl algebra have

\[
 [\alpha,A]=[\beta,B]=[\gamma,C]=1
\]

and all mixed commutators zero.  Use the ordered quantum moment element

\[
 M=2A\alpha+B\beta-C\gamma.                           \tag{11}
\]

Localize at \(\beta\).  The following four elements commute with \(M\):

\[
\begin{aligned}
 Q_1&=A\beta^2,&P_1&=\alpha\beta^{-2},\\
 Q_2&=C\beta^{-1},&P_2&=\gamma\beta.
\end{aligned}
\tag{12}
\]

There are no ordering corrections because \(\beta\) commutes with the
\((A,\alpha)\) and \((C,\gamma)\) pairs.  Directly,

\[
 [P_i,Q_j]=\delta_{ij},\qquad
 [P_1,P_2]=[Q_1,Q_2]=0.                               \tag{13}
\]

For exhaustion, Fourier-transform the remaining pair by putting

\[
 U=\beta,\qquad V=-B,\qquad [V,U]=1.
\]

Then

\[
 M=2A\alpha-VU-C\gamma.
\]

In the reduction \(M=c\), invertibility of \(U\) eliminates \(V\):

\[
 V=(2A\alpha-C\gamma-c)U^{-1}.                        \tag{14}
\]

Every weight-zero normally ordered Laurent monomial is generated by (12),
and (13) has the PBW basis of the second Weyl algebra.  Hence the localized
target quantum Hamiltonian reduction is exactly

\[
 \boxed{A_2.}                                         \tag{15}
\]

The unresolved algebra is the corresponding reduction on the source side,
where the inverted primitive element is the exotic semi-invariant
\(\Phi_F(\beta)=\delta_B\).

## 6. The remaining polynomial and quantum questions

The torus search has therefore reached a sharper endpoint than a generic
symmetry question:

> Is the smooth symplectic hypersurface \(X_c\) polynomially symplectic-
> isomorphic to \(\mathbb A^4\), compatibly with (8)?

A positive answer would turn (8) into another four-dimensional polynomial
symplectic Keller counterexample.  The repository already has such a
counterexample by cotangent-graph restriction, so the main value would be
quantum: the construction is now first-class rather than second-class.

The corresponding fixed-rank question is now one-sided:

> Is quantum Hamiltonian reduction of the source \(A_3\), localized at
> \(\delta_B=\Phi_F(\beta)\), also \(A_2\), compatibly with the explicit
> target identification (12)--(15)?

This must be checked at the normalizer/Ore-localization level.  The global
graph-centralizer theorem does not answer it, because it concerns the
second-class graph pair and the unlocalized one-moment constraint
\(\delta_S=0\), whereas the present chart inverts a primitive momentum
semi-invariant.

There are now three finite next tests.

1. Compute invariants of \(X_c\) that would obstruct \(X_c\simeq\mathbb A^4\):
   units, locally nilpotent derivations, Makar--Limanov/Derksen invariants,
   and factoriality.
2. Search for a triangular elimination of (7), first at \(c=0\), using the
   adapted coordinates \((X,Q,Z)\) from the rank-two symplectic descent.
3. Compute the source normalizer after adjoining \(\delta_B^{-1}\), using
   \((F_i,\delta_i)\) for the image subalgebra but exhausting the ambient
   \(A_3[\delta_B^{-1}]\), not only \(\Phi_F(A_3)[\delta_B^{-1}]\).

The exact polynomial certificate is
[`verify_natural_torus_reduction.py`](../scripts/verify_natural_torus_reduction.py).
