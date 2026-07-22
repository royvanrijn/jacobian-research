# The normalized linear–quadratic factorization bridge

This note connects the intrinsic marked-root construction directly to the
expanded foundational Keller map.  It isolates three independent facts:

1. multiplication together with the resultant is étale on the coprime locus;
2. the normalized factorization slice is polynomially isomorphic to `A^3`;
3. multiplication on that slice is linearly equivalent to the displayed map.

Work over a field `k` of characteristic zero.  Write

\[
L=aT+bS,
\qquad
Q=cT^2+dTS+eS^2.
\]

Their product is

\[
LQ=(ac)T^3+(ad+bc)T^2S+(ae+bd)TS^2+(be)S^3,
\]

and their resultant is

\[
R(L,Q)=a^2e-abd+b^2c.
\]

## 1. The ambient coefficient–resultant map

Define

\[
\Theta:\mathbb A^2\times\mathbb A^3\longrightarrow\mathbb A^4\times\mathbb A^1,
\qquad
(L,Q)\longmapsto(LQ,R(L,Q)).
\]

### Proposition

The morphism `Theta` is étale on the open set `R(L,Q) != 0`.

### Proof

Suppose a tangent vector `(dot L,dot Q)` is killed by the differential of
multiplication.  Then

\[
\dot L\,Q+L\,\dot Q=0.
\]

When `R(L,Q) != 0`, the factors are coprime.  Hence `L` divides `dot L`; since
both are linear, there is a scalar `lambda` such that

\[
(\dot L,\dot Q)=\lambda(L,-Q).
\]

This is exactly the infinitesimal relative-scaling direction.  The resultant
has bidegree `(2,1)`, so along this direction

\[
dR_{(L,Q)}(L,-Q)=(2-1)R(L,Q)=R(L,Q),
\]

which is nonzero.  Thus the differential of `Theta` has trivial kernel, and
source and target both have dimension five.  Therefore it is an isomorphism.

In coordinates this is the determinant identity

\[
\det D\Theta=-R(L,Q)^2.
\]

The conceptual content is that multiplication forgets exactly the relative
scaling `(L,Q) -> (lambda L,lambda^{-1}Q)`, while the resultant detects that
missing direction.

## 2. The normalized slice

Put

\[
X_{\mathrm{fac}}=
\left\{(a,b,c,d,e)\in\mathbb A^5:
R(L,Q)=1,\ ad+bc=1
\right\}.
\]

Equivalently, `X_fac` is the base change of `Theta` to the affine target slice

\[
\{[LQ]_{T^2S}=1\}\times\{R=1\}.
\]

Consequently the multiplication morphism

\[
\mu:X_{\mathrm{fac}}\longrightarrow\mathbb A^3,
\qquad
(a,b,c,d,e)\longmapsto(ac,ae+bd,be)
\]

is étale.  No separate calculation near `a=0` is required: étaleness is
preserved by base change.

## 3. Polynomial coordinates on the slice

Define

\[
\Phi:\mathbb A^3_{a,y,z}\longrightarrow X_{\mathrm{fac}}
\]

by

\[
\begin{aligned}
b&=1+ay,\\
c&=1-\frac32ay+a^2z,\\
d&=\frac12y-az+\frac32ay^2-a^2yz,\\
e&=-2z+4y^2-4ayz+3ay^3-2a^2y^2z.
\end{aligned}
\]

Direct substitution gives

\[
ad+bc=1,
\qquad
a^2e-abd+b^2c=1.
\]

The inverse is polynomial on the whole slice, including `a=0`:

\[
\boxed{
 y=2bd-ae,
 \qquad
 z=2d^2+ce+6bd^2+3bce-\frac92e.
}
\]

Substitution in both directions proves

\[
\boxed{X_{\mathrm{fac}}\cong\mathbb A^3.}
\]

The residual torus action

\[
\lambda\cdot(a,b,c,d,e)
=(\lambda a,b,c,\lambda^{-1}d,\lambda^{-2}e)
\]

becomes the linear action

\[
\lambda\cdot(a,y,z)=(\lambda a,\lambda^{-1}y,\lambda^{-2}z).
\]

Thus the formerly apparent Laurent coordinates are globally regular
homogeneous coordinates of weights `1,-1,-2`.

## 4. The Keller map obtained by multiplication

Compose the slice isomorphism with multiplication and omit the normalized
constant coordinate `ad+bc=1`:

\[
G=\mu\circ\Phi:\mathbb A^3_{a,y,z}\longrightarrow\mathbb A^3,
\qquad
G=(ac,ae+bd,be).
\]

The ambient étaleness argument already proves that `G` is étale.  Exact
expansion gives the sharper normalization

\[
\boxed{\det DG=-1.}
\]

Let `F_orig` be the foundational polynomial in the repository.  Define linear
maps

\[
A(z_1,z_2,z_3)=\left(z_1,z_2,-\frac12z_3\right),
\qquad
B(u_1,u_2,u_3)=(u_3,2u_2,2u_1).
\]

Then

\[
\boxed{F_{\mathrm{orig}}=B\circ G\circ A.}
\]

Indeed, after `a=z_1`, `y=z_2`, and `z=-z_3/2`, one has

\[
F_{\mathrm{orig}}=(G_3,2G_2,2G_1).
\]

Since `det A=-1/2`, `det B=-4`, and `det DG=-1`, this also recovers

\[
\det DF_{\mathrm{orig}}=-2.
\]

## 5. Geometric interpretation

A generic cubic has three distinct linear factors.  For each choice of linear
factor, the condition `R(L,Q)=1` fixes its relative scaling uniquely.  Hence
`mu` is generically three-to-one.  The three ingredients of the counterexample
are therefore cleanly separated:

\[
\begin{array}{rcl}
\text{étaleness}
&\Longleftarrow&
\text{coprimality plus the resultant weight},\\
\text{generic degree three}
&\Longleftarrow&
\text{three choices of marked root},\\
X_{\mathrm{fac}}\cong\mathbb A^3
&\Longleftarrow&
\text{the tangent non-osculating normalized slice}.
\end{array}
\]

The exact symbolic certificate is
[`scripts/verify_normalized_factorization_slice.py`](../scripts/verify_normalized_factorization_slice.py).
