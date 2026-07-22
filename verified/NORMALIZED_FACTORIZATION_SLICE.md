# The normalized linear–quadratic factorization bridge

This note connects the intrinsic marked-root construction directly to the
expanded foundational Keller map.  It isolates three independent facts:

1. multiplication together with the resultant is étale on the coprime locus;
2. the normalized factorization slice is polynomially isomorphic to `A^3`;
3. multiplication on that slice is linearly equivalent to the displayed map.

Work over a field `k` of characteristic zero, and put

\[
V_i=\operatorname{Sym}^i(k^2).
\]

Write

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
source and target both have dimension five.  Therefore
\(d\Theta_{(L,Q)}\) is an isomorphism, so \(\Theta\) is etale at `(L,Q)`.

In coordinates this is the determinant identity

\[
\det D\Theta=-R(L,Q)^2.
\]

The conceptual content is that multiplication forgets exactly the relative
scaling `(L,Q) -> (lambda L,lambda^{-1}Q)`, while the resultant detects that
missing direction.

## 2. The normalized slice and the projective open

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

Let

\[
U=\bigl(\mathbb P(V_1)\times\mathbb P(V_2)\bigr)
\setminus\bigl(\{R=0\}\cup\{[LQ]_{T^2S}=0\}\bigr).
\]

This is exactly the projective source open in the marked-root construction.
Projectivization gives an isomorphism

\[
\boxed{X_{\mathrm{fac}}\xrightarrow{\sim}U.}
\]

Indeed, take a projective pair `([L],[Q])` in `U`, choose arbitrary
representatives, and put

\[
m=[LQ]_{T^2S},\qquad r=R(L,Q).
\]

Both are nonzero.  Under independent rescaling

\[
(L,Q)\longmapsto(\lambda L,\mu Q),
\]

they transform as

\[
m\longmapsto\lambda\mu m,
\qquad
r\longmapsto\lambda^2\mu r.
\]

The unique rescaling for which both normalized values equal one is

\[
\boxed{
\lambda={m\over r},
\qquad
\mu={r\over m^2}.
}
\]

This construction is independent of the chosen representatives.  If they are
replaced by `(alpha L,beta Q)`, then `m` and `r` become
`alpha beta m` and `alpha^2 beta r`, while the new normalization scalars become
`lambda/alpha` and `mu/beta`.  Thus the normalized affine pair is unchanged.
The formulas are regular on `U`, so they give the inverse to projectivization.

This identifies the normalized complete intersection literally, not merely
birationally, with the geometric source open.

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

## 6. Unequal factor degrees

The étaleness mechanism is not special to degrees one and two.  For positive
integers `p,q`, define

\[
\Theta_{p,q}:V_p\times V_q\longrightarrow
V_{p+q}\times\mathbb A^1,
\qquad
(A,B)\longmapsto(AB,\operatorname{Res}(A,B)).
\]

### Unequal-degree coefficient–resultant theorem

If `p != q`, then `Theta_(p,q)` is étale on the coprime locus.

### Proof

If

\[
\dot A B+A\dot B=0
\]

and `A,B` are coprime, then the same divisibility argument gives

\[
(\dot A,\dot B)=\lambda(A,-B).
\]

The resultant has bidegree `(q,p)`.  Therefore

\[
d\operatorname{Res}_{(A,B)}(A,-B)
=(q-p)\operatorname{Res}(A,B).
\]

In characteristic zero this is nonzero when `p != q`.  It kills the unique
kernel direction of multiplication, and the source and target dimensions are
both `p+q+2`.

The normalization step reveals why **consecutive** degrees are especially
natural.  Suppose `q>p`, let `m` be any nonzero linear functional of the
product `AB`, and put `r=Res(A,B)`.  Under factor rescaling,

\[
m\longmapsto\lambda\mu m,
\qquad
r\longmapsto\lambda^q\mu^p r.
\]

After imposing `lambda mu m=1`, the resultant condition becomes

\[
\lambda^{q-p}={m^p\over r},
\qquad
\mu={1\over\lambda m}.
\]

Thus:

- if `q-p=1`, there is a unique algebraic normalization,
  \[
  \lambda={m^p\over r},\qquad
  \mu={r\over m^{p+1}};
  \]
- if `q-p>1`, a residual `mu_(q-p)` ambiguity remains;
- if `q=p`, the resultant has weight zero along relative scaling and cannot
  repair the differential kernel of multiplication.

So unequal degrees explain étaleness, while consecutive degrees explain the
canonical normalized affine slice used here.

The exact symbolic certificate is
[`scripts/verify_normalized_factorization_slice.py`](../scripts/verify_normalized_factorization_slice.py).
