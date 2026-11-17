# A twelve-variable degree-three Keller counterexample

## The result

Over every characteristic-zero field, the following polynomial map
\(K:\mathbb A^{12}\to\mathbb A^{12}\) has degree three, identity linear
part, Jacobian determinant one, and a rational collision:

\[
\begin{aligned}
K_1={}&z_1+\tfrac12z_1^2z_{11}-\tfrac32z_1^2z_2-z_1z_{12}z_3-z_{11}z_{12},\\
K_2={}&z_2+12z_1z_2^2-z_1z_2z_9-6z_1z_3z_8+3z_1z_3
       +3z_1z_5z_8-9z_2^2z_6\\
     &\quad+3z_3z_6z_8-3z_5z_6-z_8z_9,\\
K_3={}&z_3-z_1z_{10}z_2+3z_1z_2z_3+z_1z_7z_8-z_{10}z_8
       -3z_2^2z_4-7z_2^2z_8+4z_2^2\\
     &\quad-3z_2z_3z_6+z_2z_5z_6+z_3z_4z_8-z_4z_5-z_6z_7,\\
K_4={}&z_4-2z_1z_2z_8-z_8^2,\\
K_5={}&z_5+z_1z_2z_3+3z_2^2,\\
K_6={}&z_6+z_1^2z_2,\\
K_7={}&z_7+3z_2z_3-z_2z_5,\\
K_8={}&z_8+z_1z_2,\\
K_9={}&z_9+6z_1z_3-3z_1z_5-3z_3z_6,\\
K_{10}={}&z_{10}-z_1z_7+7z_2^2-z_3z_4,\\
K_{11}={}&z_{11}+z_1z_3,\\
K_{12}={}&z_{12}-\tfrac12z_1^2.
\end{aligned}                                                       \tag{1}
\]

The colliding points and common image are

\[
\begin{aligned}
p={}&(0,0,-1/4,0,0,0,0,0,0,0,0,0),\\
q={}&(1,-3/2,13/2,-9/4,3,3/2,99/4,3/2,-3/4,-45/8,-13/2,1/2),\\
K(p)={}&K(q)=p.
\end{aligned}                                                       \tag{2}
\]

This is an exact upper-bound theorem, not a minimality or priority claim.
The derivation uses A. MacFarlane's \(F_{13}\), whose provenance is pinned
in the [twenty-variable audit](../extended-geometry/MACFARLANE_G20_DIMENSION_REDUCTION_AUDIT.md).

## Coordinate-pair restriction

The mechanism is broader than a fixed invariant polynomial.  Let
\(F:\mathbb A^{n+1}\to\mathbb A^{n+1}\) be Keller.  Suppose determinant-one
source and target automorphisms put it in the relative form

\[
\Phi(z,s)=\bigl(K(z)+sL(z,s),s\bigr).                              \tag{3}
\]

Then the source level \(s=0\) maps to the target level \(s=0\), and

\[
D\Phi(z,0)=
\begin{pmatrix}
DK(z)&L(z,0)\\
0&1
\end{pmatrix}.
\]

Therefore \(\det DK=\det D\Phi=1\).  Any collision of \(\Phi\) contained in
the common level \(s=0\) descends to a collision of \(K\).  Unlike a
pullback-fixed invariant \(P\circ F=P\), the source coordinate and target
coordinate may be different.  This is why the earlier invariant search did
not see the reduction.

For MacFarlane's map, use

\[
T(x)=(x_1,\ldots,x_{12},\,x_{13}+x_2^2)
     =(x_1,\ldots,x_{12},F_{13,13}(x))                            \tag{4}
\]

on the source and

\[
A(y)=(y_1,y_2,y_3,\,y_4-y_8^2,y_5,\ldots,y_{13})                 \tag{5}
\]

on the target.  Both are triangular automorphisms with determinant one.
Exact substitution gives

\[
A\circ F_{13}\circ T^{-1}(z,s)
=
\bigl(
K_1,K_2,K_3,K_4+s(2z_{12}-z_1^2),K_5,\ldots,K_{12},s
\bigr).                                                          \tag{6}
\]

The source graph \(x_{13}=-x_2^2\) by itself would create the quartic
\(z_1^2z_2^2\) in component four.  The square completion \(y_4-y_8^2\)
removes it because \(F_{13,8}=x_8+x_1x_2\).  Equation (6) is the exact
factorization certificate for (1).

## What the construction adds to the search language

There are three distinct reduction tests:

1. a fixed invariant asks for \(P\circ F=P\);
2. a coordinate-pair restriction asks for a source coordinate \(h\) and a
   target coordinate \(g\) with \(g\circ F=h\);
3. a degree-preserving coordinate-pair restriction may first require target
   completion to cancel the high-degree defect created by solving \(h=0\).

The second condition is strictly weaker than the first because \(g\) and
\(h\) need not be the same polynomial.  In the present construction,
\(g=y_{13}\) while \(h=x_{13}+x_2^2\).  The earlier fixed-invariant audit was
therefore correct but searched the wrong equality.

This also turns the degree cleanup into exact linear algebra.  After a
candidate source graph is substituted, choose a bounded monomial basis
\(P_\alpha\) for triangular target shears.  The coefficients of all terms
above the desired degree in

\[
F_i-\sum_\alpha c_\alpha P_\alpha(F)
\]

are linear in the unknown \(c_\alpha\).  For the successful graph
\(x_{13}=-x_2^2\), the defect is \(z_1^2z_2^2\), and the quadratic target
basis finds \(P=y_8^2\).  This suggests a backward compiler transition
consisting of “graph deletion plus target completion,” scored before
declaring a degree increase fatal.

## Homogeneous consequence

Write \(K=z+Q+C\), with \(Q\) quadratic and \(C\) cubic.  Exact coefficient
row reduction gives

\[
\dim\operatorname{span}\{C_1,\ldots,C_{12}\}=6,
\qquad C_7=\cdots=C_{12}=0.                                     \tag{7}
\]

Using the first six coordinate directions as \(B\), rank-compressed
homogenization therefore gives the nineteen-variable map

\[
G_{19}(z,w,\tau)=
\left(
z+\tau Q(z)+\tau^2B w,\;
w-(C_1(z),\ldots,C_6(z)),\;
\tau
\right).                                                         \tag{8}
\]

Every nonlinear term in (8) is homogeneous cubic.  Companion cancellation
at \(\tau=1\), together with

\[
z+\tau Q(z)+\tau^2C(z)=\tau^{-1}K(\tau z),
\]
proves \(\det DG_{19}=1\).  The points
\((p,C(p),1)\) and \((q,C(q),1)\) collide.  Thus the certified upper bounds
become

\[
n_{\mathrm{degree}\le3}\le12,\qquad
n_{\mathrm{cub}}\le19,\qquad
n_{\mathrm{HN},4}\le38,                                         \tag{9}
\]

where the last inequality is the homogeneous cotangent lift of (8).

## Verification

Run

```bash
make verify-macfarlane-f12
```

The primary checker reconstructs (6) from the pinned \(F_{13}\) formulas,
computes the determinant of (1) as an exact sparse polynomial, verifies
(2), proves the cubic-output rank in (7), and constructs and checks (8).
The independent checker hard-codes (1) in a separate standard-library
sparse-polynomial implementation and recomputes its degree, collision, and
exact determinant without SymPy or shared repository algebra code.

The generated record is
[`macfarlane_f12_coordinate_pair_reduction.json`](../artifacts/generated-results/macfarlane_f12_coordinate_pair_reduction.json).
No Lean formalization or external review is claimed.
