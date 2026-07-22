# The normalized factorization model

This note packages the foundational counterexample into three propositions:
an explicit affine-space source, a coordinate-free étaleness mechanism, and
an exact linear comparison with the announced polynomial.  The longer
[factorization bridge](NORMALIZED_FACTORIZATION_SLICE.md) supplies the
projective interpretation and the extension to arbitrary unequal degrees.

Work over a field `k` of characteristic zero.  Put

\[
 L=aT+bS,\qquad Q=cT^2+dTS+eS^2.
\]

Then

\[
 \operatorname{Res}(L,Q)=a^2e-abd+b^2c,\qquad
 [LQ]_{T^2S}=ad+bc.
\]

## Proposition 1: the normalized factorization variety

Define

\[
 X=\left\{(a,b,c,d,e)\in\mathbb A^5:
 a^2e-abd+b^2c=1,\quad ad+bc=1\right\}.
\]

There is a polynomial isomorphism

\[
 \Phi:\mathbb A^3_{a,y,z}\xrightarrow{\sim}X
\]

given by

\[
\begin{aligned}
 b&=1+ay,\\
 c&=1-\frac32ay+a^2z,\\
 d&=\frac12y-az+\frac32ay^2-a^2yz,\\
 e&=-2z+4y^2-4ayz+3ay^3-2a^2y^2z.
\end{aligned}
\]

Its inverse is the polynomial map

\[
 \boxed{y=2bd-ae,\qquad
 z=2d^2+ce+6bd^2+3bce-\frac92e.}
\]

### Proof

Substitution of the four displayed forward formulas gives

\[
 ad+bc=1,\qquad a^2e-abd+b^2c=1.
\]

Substituting them into the two inverse formulas returns `y` and `z` exactly.
Conversely, substitute the inverse formulas into the forward formulas and
reduce modulo

\[
 (a^2e-abd+b^2c-1,\ ad+bc-1).
\]

The four remainders are zero.  Thus both compositions are identities in the
coordinate rings.  In
particular this is a global polynomial certificate, including the divisor
`a=0`; no localization at `a` occurs.

### Vertical derivation and polynomial slice

There is a more structural way to read the same formulas.  Hold `a,b` fixed
and take the saturated primitive extension of the generic common tangent
vector to the two equations in the `c,d,e` variables.  The two gradients are

\[
 (b,a,0),\qquad (b^2,-ab,a^2),
\]

and their cross product has the common factor `a`.  Removing that factor gives
the globally defined derivation

\[
 \boxed{D=a^2\partial_c-ab\partial_d-2b^2\partial_e.}
\]

Equivalently, `D_0=D/2` is induced by the elementary factor operation

\[
 Q\longmapsto Q+\frac{t}{2}L(aT-2bS).                \tag{1}
\]

Indeed, the added quadratic has coefficient vector
`(a^2/2,-ab/2,-b^2)`.  Adding a multiple of `L` to `Q` does not change the
resultant modulo `L`, and a direct multiplication shows that (1) also
preserves `[LQ]_(T^2S)`.  Thus the vertical `G_a`-action has an intrinsic
factorization interpretation.

It preserves both defining equations:

\[
 D(ad+bc-1)=0,\qquad D(a^2e-abd+b^2c-1)=0.
\]

Moreover `D(a)=D(b)=0` and `D^2` kills every ambient coordinate, so `D` is
locally nilpotent.  Its additive action is

\[
 (a,b,c,d,e)\longmapsto
 (a,b,c+t a^2,d-tab,e-2tb^2).
\]

The inverse coordinates above satisfy

\[
 D(y)=0,\qquad D(z)=1.                              \tag{2}
\]

There is also a mechanical Bézout construction of the slice.  Since
`b-ay=1` in `k[X]` and

\[
\begin{aligned}
 D_0(-e)&=b^2,\\
 D_0(4yd)&=-2aby,\\
 D_0(2cy^2)&=a^2y^2,
\end{aligned}
\]

the polynomial

\[
 \widetilde z=4yd+2cy^2-e
\]

satisfies

\[
 D_0(\widetilde z)=(b-ay)^2=1.
\]

Reduction by the defining ideal gives `tilde z=2z`, so this is exactly the
same affine coordinate with the generator normalized by `D_0=D/2`.

Thus `z` is a global polynomial slice for the vertical `G_a`-action.  The
slice theorem gives

\[
 k[X]=\ker(D)[z].
\]

On the section `z=0`, the two invariant coordinates `a,y` give

\[
\begin{aligned}
 b&=1+ay,&
 c&=1-\frac32ay,\\
 d&=\frac12y+\frac32ay^2,&
 e&=4y^2+3ay^3.
\end{aligned}
\]

Hence `ker(D)=k[a,y]` and

\[
 \boxed{k[X]=k[a,y,z].}
\]

This recasts the explicit elimination as a “vertical locally nilpotent
derivation plus polynomial slice” argument.  The field is obtained
canonically by saturating the Jacobian cross-product by its common factor
`a`.  At `a=0` the fibers of `X -> A^2_(a,b)` jump in dimension, so no claim
is made that `D` spans their entire tangent space there.

The primitive saturated generator is unique up to a nonzero scalar after the
projection and orientation are fixed.  The `G_a`-action is not unique in the
broader sense: every `f(a,y)D` with `f in ker(D)` is again locally nilpotent
and has the same invariant ring when `f` is nonzero.  A nonconstant `f` may
introduce additional fixed fibers, so only invariant **units** preserve the
everywhere-free slice action.

## Proposition 2: the coefficient--resultant map is étale

Let `V_i=\operatorname{Sym}^i(k^2)` and define

\[
 \Theta:V_1\times V_2\longrightarrow V_3\times\mathbb A^1,
 \qquad (L,Q)\longmapsto(LQ,\operatorname{Res}(L,Q)).
\]

Then `Theta` is étale wherever `L` and `Q` are coprime.

### Proof

Suppose `(\dot L,\dot Q)` is killed by the differential of multiplication:

\[
 \dot L\,Q+L\,\dot Q=0.
\]

Coprimality implies that `L` divides `\dot L` and `Q` divides `\dot Q`.
The degree bounds therefore give

\[
 (\dot L,\dot Q)=\lambda(L,-Q).
\]

This is the relative-scaling direction forgotten by multiplication.  The
resultant has bidegree `(2,1)`, so

\[
 d\operatorname{Res}_{(L,Q)}(L,-Q)
 =(2-1)\operatorname{Res}(L,Q)
 =\operatorname{Res}(L,Q).
\]

It is nonzero on the coprime locus.  Thus `d\Theta_(L,Q)` is injective.
Source and target both have dimension five, so `d\Theta_(L,Q)` is an
isomorphism and `Theta` is étale at `(L,Q)`.  In coordinates the same fact is
the identity

\[
 \boxed{\det D\Theta=-\operatorname{Res}(L,Q)^2.}
\]

Let

\[
 H=\{[C]_{T^2S}=1,\ \rho=1\}
 \subset V_3\times\mathbb A^1,
\]

where `C` is the cubic coordinate and `rho` the resultant coordinate.  Then
`X=\Theta^{-1}(H)`, and the map

\[
 X\longrightarrow H\simeq\mathbb A^3,\qquad
 (L,Q)\longmapsto(ac,ae+bd,be)
\]

is the base change of `Theta` along `H`.  It is therefore étale.  Conceptually,

\[
\boxed{\text{multiplication forgets relative scaling;\quad
the resultant detects relative scaling}.}
\]

## Proposition 3: comparison with the foundational polynomial

In the coordinates of Proposition 1, define

\[
 G(a,y,z)=(ac,\ ae+bd,\ be),
\]

where `b,c,d,e` are the displayed polynomials in `a,y,z`.  Then

\[
 \boxed{\det DG=-1.}
\]

Let `F_{\mathrm{orig}}` denote the foundational polynomial map and set

\[
 A(z_1,z_2,z_3)=\left(z_1,z_2,-\frac12z_3\right),\qquad
 B(u_1,u_2,u_3)=(u_3,2u_2,2u_1).
\]

Then

\[
 \boxed{F_{\mathrm{orig}}=B\circ G\circ A.}
\]

Equivalently,

\[
 F_{\mathrm{orig}}(z_1,z_2,z_3)
 =(G_3,2G_2,2G_1)
 \left(z_1,z_2,-\frac12z_3\right).
\]

Since `det A=-1/2` and `det B=-4`, this also gives

\[
 \det DF_{\mathrm{orig}}
 =(\det B)(\det DG)(\det A)=-2.
\]

The construction is therefore summarized by the entirely explicit chain

\[
 \boxed{\text{normalized coprime factorization}
 \longrightarrow X\simeq\mathbb A^3
 \xrightarrow{\,G\,}\mathbb A^3
 \xrightarrow{\text{linear target change}}F_{\mathrm{orig}}.}
\]

## Exact certificate

The canonical symbolic verifier is
[`scripts/verify_normalized_factorization_slice.py`](../scripts/verify_normalized_factorization_slice.py).
It checks:

- both defining equations after the forward substitution;
- both compositions, with the reverse composition reduced by a Gröbner basis
  of the defining ideal;
- the primitive vertical locally nilpotent derivation, `D(y)=0`, and
  `D(z)=1`;
- `det DG=-1`;
- the exact linear equivalence with `F_orig`;
- the optional five-dimensional identity
  `det DTheta=-Res(L,Q)^2`.

It is part of `make verify-core`.  Keeping one executable certificate avoids
two copies of the same formulas drifting apart.
