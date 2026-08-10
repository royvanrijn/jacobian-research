# The degree-108 JC2 no-vertical branch as a finite Belyi deformation problem

## Scope

Guccione--Guccione--Horruitiner--Valqui reduced the then-unresolved plane
Jacobian degree pair `(72,108)` to two combinatorial families.
Their `(9,27),(2,3)` family is eliminated in Section 5 of arXiv:2204.14178;
the remaining family is `(8,28),(3,2)`.  Proposition 4.3 gives, after a
birational Laurent coordinate reduction, a pair `P,Q` with

\[
[P,Q]_{x,y}=x^2
\]

and two possible Newton polygons.  This note treats the first/no-vertical-edge
polygon

\[
N(P)=\{(0,0),(1,0),(8,14),(8,16)\},
\]
\[
N(Q)=\{(0,0),(2,1),(12,21),(12,24)\}.
\tag{0.1}
\]

It converts its top equation into a finite Belyi problem.

## 1. A symplectic Laurent coordinate change

Set

\[
X=xy^2,\qquad s=y^{-1}.
\tag{1.1}
\]

Then

\[
\det\frac{\partial(X,s)}{\partial(x,y)}=-1,
\]

so

\[
[P,Q]_{X,s}=-X^2s^4.
\tag{1.2}
\]

Every monomial transforms by

\[
x^iy^j=X^is^{2i-j}.
\]

The polygons (0.1) therefore imply the exact forms

\[
P=A(X)s^2+B(X)s+C(X),
\tag{1.3}
\]

\[
Q=D(X)s^3+E(X)s^2+F(X)s+G(X),
\tag{1.4}
\]

with

\[
\begin{array}{c|ccc}
 &A&B&C\\ \hline
\deg_X&8&8&8\\
\operatorname{ord}_X&1&1&0
\end{array}
\]

and

\[
\begin{array}{c|cccc}
 &D&E&F&G\\ \hline
\deg_X&12&12&12&12\\
\operatorname{ord}_X&2&2&1&0.
\end{array}
\]

Thus the original bivariate system is a five-equation system in seven
**univariate** polynomials.

## 2. The five Wronskian equations

Expanding (1.2) by powers of `s` gives

\[
\begin{aligned}
s^4:&\quad -2AD'+3DA'=-X^2,\\
s^3:&\quad -2AE'-BD'+3DB'+2EA'=0,\\
s^2:&\quad -2AF'-BE'+3DC'+2EB'+FA'=0,\\
s^1:&\quad -2AG'-BF'+2EC'+FB'=0,\\
s^0:&\quad -BG'+FC'=0.
\end{aligned}
\tag{2.1}
\]

The first equation is equivalently

\[
\boxed{2AD'-3DA'=X^2.}
\tag{2.2}
\]

The remaining four equations are successive deformation equations over a
solution `(A,D)` of (2.2).

## 3. The top equation is Belyi

Equation (2.2), together with

\[
\deg A=8,\quad \operatorname{ord}_0A=1,
\qquad
\deg D=12,\quad \operatorname{ord}_0D=2,
\tag{3.1}
\]

has a useful global interpretation.  Put

\[
\beta(X)=\frac{D(X)^2}{A(X)^3}.
\tag{3.2}
\]

At every nonzero root `alpha` of `A`, (2.2) gives

\[
D(\alpha)A'(\alpha)\ne0,
\]

so every such root is simple and is not a root of `D`.  Similarly every
nonzero root of `D` is simple and disjoint from `A`.  Hence

\[
A=X\bar A,\qquad D=X^2\bar D,
\]

with `bar A` squarefree of degree 7, `bar D` squarefree of degree 10 and
`gcd(bar A,bar D)=1`.

After cancelling the common power at zero,

\[
\beta=X\frac{\bar D^2}{\bar A^3}
\]

has degree 21.  Moreover

\[
\frac{\beta'}{\beta}
=\frac{2AD'-3DA'}{AD}
=\frac{X^2}{AD}.
\tag{3.3}
\]

Therefore its ramification is completely determined:

* over `0`: ten points of ramification index 2 and the simple point `X=0`;
* over `infinity`: seven points of ramification index 3;
* at `X=infinity`: since `beta'=Theta(X^{-18})`, the local ramification index
  is 17; the remaining four points in that fiber are simple.

The Riemann--Hurwitz contributions are

\[
10+7(3-1)+(17-1)=40=2\cdot21-2,
\]

so there is no other ramification.  Thus `beta` is a degree-21 Belyi map with
passport

\[
\boxed{
(2^{10},1),\qquad(3^7),\qquad(17,1,1,1,1).
}
\tag{3.4}
\]

Consequently the top `(A,D)` layer is a finite Hurwitz/Belyi problem rather
than a positive-dimensional polynomial ansatz.

## 4. Dessin enumeration

Fix the permutation of the third branch value to be a 17-cycle with four fixed
points.  If `sigma_0` has type `(2^10,1)` and `sigma_1` has type `(3^7)`, the
relation

\[
sigma_0 sigma_1 sigma_\infty=1
\]

forces a simple matching problem on the 17-cycle.  Exact enumeration, modulo
the rotations centralizing the 17-cycle and permutations of its four fixed
points, leaves five center-set types:

\[
\begin{aligned}
&(0,3,7,11),\\
&(0,3,7,12),\\
&(0,3,8,11),\\
&(0,3,8,13),\\
&(0,3,9,13),
\end{aligned}
\tag{4.1}
\]

where indices are read modulo 17.  The exact branch-cycle replay verifies that
these are five transitive dessin types with trivial automorphism groups;
orientation reversal is part of their arithmetic Galois action, not an extra
quotient in this enumeration.

Thus the no-vertical-edge `(8,28)` residue has only finitely many top dessins.
The four lower equations in (2.1) are linear successively in the deformation
pairs `(B,E)`, `(C,F)`, and `G` once `(A,D)` is fixed.  This is the natural
next computation: solve those deformation spaces for each Belyi type instead
of the original 72-coefficient system.

## 5. Connection with the curvature approach

The Laurent transformation does not preserve a constant Jacobian: the reduced
pair has bracket `-X^2s^4`.  Nevertheless the exact target-shear curvature law
from `JC2_DUAL_SCHUR_CURVATURE_OBSTRUCTION.md` remains valid with the actual
Jacobian polynomial in place of the constant `c`; the correction term is the
square of the bracket.

More importantly, (2.2) is itself the `2:3` approximate-root resonance which
appears in the generalized Magnus expansion of a `(72,108)` pair.  The Belyi
reduction therefore packages the same `2:3` obstruction globally rather than
as a large coefficient system.

## 6. Exact implementation and proof boundary

The finite program proposed above is now implemented and closes this
no-vertical-edge Laurent branch.  The intrinsic quintic graph reconstructs
the five conjugate `(A,D)` pairs exactly; its irreducible quintic has Galois
closure group `S_5`.  The distinguished `X=0` point and the polynomial
normalization are compatible with every conjugate, so the arithmetic shortcut
does not remove a dessin before deformation.

Over the quintic field the successive `(B,E)`, `(C,F)`, and derivative-`G`
maps have ranks `17`, `18`, and `12`, with kernel dimensions `2`, `3`, and
`0`; the constant of `G` is a separate target-translation kernel.  After
solving these maps, the last two levels give only 25 sparse degree-three/four
equations in five parameters.  On the required `deg(B)=8` open their exact
Singular ideal is `(1)`.

The canonical proof and reproduction commands are in
[`plane-jc/JC2_72_108_BELYI_DEFORMATION_CLOSURE.md`](plane-jc/JC2_72_108_BELYI_DEFORMATION_CLOSURE.md).
This conclusion is conditional on the audited Proposition-4.3 Laurent
reduction and the complete intrinsic first-block graph; it is not a
stand-alone proof of the general degree reduction or of `JC(2)`.  The other
Laurent polygon is closed by the separate certified `(72,108)` calculations
recorded in `MATH_STATUS.json`.
