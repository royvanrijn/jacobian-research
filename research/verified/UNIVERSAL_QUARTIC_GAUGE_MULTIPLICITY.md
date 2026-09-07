# Universal quartic multiplicity by power-shifted gauges

> **All-degree extension.**  The
> [common power-shift theorem](UNIVERSAL_POWER_SHIFTED_GAUGE_MULTIPLICITY.md)
> applies the same extra `P`-power to every decoration of degree at least
> four and proves the corresponding result uniformly for all `N>=4`.  The
> quartic theorem below is its triangular first case.

The trace-chord construction proves universal quartic multiplicity over
number fields, but its defining quadric can be anisotropic over a general
characteristic-zero field.  This note gives a second mechanism which has no
quadratic-form input.

The new parameter does not move the primitive generator.  It changes the
power with which the first target coordinate enters the quartic decoration.
That change is invisible on the selected fiber, where the coordinate equals
one, but it changes an intrinsic lattice index on the normalized ramified
boundary.

## The theorem

> **Universal quartic gauge-multiplicity theorem.**
> Let `K` be any characteristic-zero field and let `A` be any rank-four
> finite etale `K`-algebra.  Then
> \[
>   \boxed{|\mathcal R_K(A)|=\infty.}
> \]
> All classes may be represented by determinant-one polynomial maps
> `A^3_K -> A^3_K` of geometric degree four.  They are indexed by
> `m>=0` and are separated by the intrinsic odd lattice indices
> \[
>   2m+5.
> \]

Thus the anisotropic biquadratic algebra over
`\mathbb Q((a))((b))` from the
[low-rank boundary note](LOW_RANK_MULTIPLICITY_BOUNDARIES.md) still occurs
in infinitely many stable classes.  Anisotropy obstructs the weighted
trace-chord presentation, not quartic multiplicity itself.

## 1. A quartic root presentation

Every finite etale algebra over the infinite field `K` is monogenic.  Choose
a monic squarefree quartic `f(T)` with

\[
 A\simeq K[T]/(f).
\]

Choose `a in K` such that

\[
 f'(a)f'''(a)\ne0.                                    \tag{1.1}
\]

This is possible because `f'` and `f'''` are nonzero polynomials and `K` is
infinite.  Write

\[
 G(S)=f(a+S)-f(a)
     =g_1S+g_2S^2+g_3S^3+g_4S^4.                     \tag{1.2}
\]

Then `g_1g_3g_4!=0`.

For each integer `m>=0`, replace the usual diagonal quartic lift by

\[
 \boxed{
 G_{P,m}(S)
 =g_1S+P(g_2S^2+g_3S^3)+g_4P^{m+4}S^4.
 }                                                       \tag{1.3}
\]

At `P=1`, every lift in (1.3) restricts to the same polynomial `G(S)`.
The exponent `m` is therefore invisible to the selected inverse quotient.

## 2. The polynomial Keller maps

Put

\[
 t=1+xy,\qquad
 q=t^2z+\frac{g_1}{g_3}y^2(1+3t),\qquad P=tq.
                                                               \tag{2.1}
\]

Define `F_m=(P,B_m,C_m)` by

\[
\begin{aligned}
B_m={}&y+3\frac{g_3}{g_1}xq
          +2\frac{g_2}{g_1}tq
          +4\frac{g_4}{g_1}
             t^{m+2}x^2q^{m+4},\\
C_m={}&x(5-3t)-\frac{g_3}{g_1}x^3z
          -2\frac{g_4}{g_1}
             t^mx^4q^{m+4}.
                                                               \tag{2.2}
\end{aligned}
\]

These are polynomials.  Their coordinate degrees are

\[
 \deg(F_m)=(7,\,7m+26,\,7m+24).                       \tag{2.3}
\]

On `t!=0`, set

\[
 S=\frac{x}{t},\qquad Q=y+xq,\qquad
 D=1-S(Q-PS)=\frac1t.                                  \tag{2.4}
\]

Direct substitution gives the inverse equation

\[
 \boxed{
 E_m(P,B,C;S)
 =G_{P,m}(S)-\frac{g_1}{2}(BS^2+C)=0
 }                                                       \tag{2.5}
\]

and, on the source incidence,

\[
 \boxed{\partial_SE_m=g_1D.}                           \tag{2.6}
\]

Equivalently, at fixed `P` the map is the marked-line incidence for

\[
 \left(S^2,\frac{2G_{P,m}(S)}{g_1}\right).
\]

The plane Jacobian is `-2D`; the reciprocal chart (2.4) has Jacobian
`D^{-1}`.  Hence

\[
 \boxed{\det DF_m=-2.}                                 \tag{2.7}
\]

Scaling the second target coordinate by `-1/2` gives determinant one and
does not change the stable class.

The polynomial (2.5) is irreducible over `K(P,B,C)`: it is primitive and
linear in `C`, whose coefficient is a nonzero constant.  It has degree four
in `S`, and its generic roots are simple and reconstruct through (2.4).
Thus every `F_m` has geometric degree four.

## 3. The fixed complete fiber

At

\[
 y_m=\left(1,0,-\frac{2f(a)}{g_1}\right),              \tag{3.1}
\]

the inverse equation is independent of `m`:

\[
\begin{aligned}
 E_m(y_m;S)
 &=G(S)+f(a)\\
 &=f(a+S).                                             \tag{3.2}
\end{aligned}
\]

It is squarefree.  For every root `r`, equation (2.6) gives

\[
 D_r=\frac{f'(a+r)}{g_1}\ne0,
\]

so every root lies in the reconstruction chart and produces exactly one
regular source point.  Since the geometric degree is four, no further
point or inverse sheet remains.  Scheme-theoretically,

\[
 \boxed{
 F_m^{-1}(y_m)
 \simeq\operatorname{Spec}K[S]/(f(a+S))
 \simeq\operatorname{Spec}A.
 }                                                       \tag{3.3}
\]

Thus (3.3) is a complete fiber for every `m`.

## 4. The canonical boundary is still two-vertex

Away from `P=0` and the repeated-root discriminant, (2.5)--(2.6) reconstruct
every inverse root regularly.  Over the generic point of `P=0`, the lower
Newton polygon has vertices

\[
 (0,0),\quad(2,0),\quad(3,1),\quad(4,m+4),             \tag{4.1}
\]

and slopes

\[
 0,\quad1,\quad m+3.                                   \tag{4.2}
\]

The first block is the degree-two affine branch, the second is the
degree-one affine branch, and the last is one unramified boundary branch.
On the last branch,

\[
 v(P)=1,\qquad v(S)=-(m+3),\qquad
 v(D)=-2m-5,\qquad v(q)=-2m-4<0,                       \tag{4.3}
\]

so it is genuinely outside affine source space.  The degree sum is

\[
 2+1+1=4.
\]

The other boundary image is the discriminant.  It receives one generic
ramified boundary prime of index two.  Hence the two vertices are
intrinsically ordered: the discriminant vertex has ramification index two,
while `P=0` has one unramified boundary prime.

Delete `P=0` and normalize the ramified target stratum.  A repeated root
`r` satisfies

\[
 B(r)=\frac{G_{P,m}'(r)}{g_1r},\qquad
 C(r)=\frac{2G_{P,m}(r)-rG_{P,m}'(r)}{g_1}.             \tag{4.4}
\]

As in the ordinary quadratic gauge, `r!=0`,
`dC+r^2dB=0`, and the pole of `B` at `r=0` has odd order one.  Therefore
(4.4) is birational and

\[
 \operatorname{Norm}(Z_\Delta\setminus Z_0)
 \simeq\operatorname{Spec}K[P^{\pm1},r^{\pm1}]
 \simeq\mathbb G_m^2.                                  \tag{4.5}
\]

## 5. The stable lattice index

Divide by `g_1` and write `a_j=g_j/g_1`.  Relative differentiation in
(4.4) gives

\[
 \operatorname{Fitt}_0\Omega_{\widetilde Z_\Delta^\circ/
                                  Z_\Delta^\circ}
 =(J_m),
\]

where

\[
 \boxed{
 J_m(P,r)=r^2B_r
 =-1+3a_3Pr^2+8a_4P^{m+4}r^3.
 }                                                       \tag{5.1}
\]

All three coefficients are nonzero.  Its Laurent support is

\[
 \mathcal S_m=\{(0,0),(1,2),(m+4,3)\}\subset\mathbb Z^2.
                                                               \tag{5.2}
\]

The affine lattice generated by the differences in (5.2) has index

\[
 \boxed{
 \iota_m=
 \left|
 \det\begin{pmatrix}1&2\\m+4&3\end{pmatrix}
 \right|
 =2m+5.
 }                                                       \tag{5.3}
\]

Stable normalization functoriality transports the ordered ramified stratum
and its relative Fitting divisor.  Stabilization only takes (4.5) and
(5.1) times affine space, so its units remain Laurent monomials in `P,r`.
An isomorphism acts on exponent vectors by an affine
`GL_2(\mathbb Z)` transformation.  Such a transformation, and multiplication
of `J_m` by a Laurent unit, preserve the index (5.3).

Therefore a stable polynomial left--right equivalence
`F_m~F_n` would imply

\[
 2m+5=2n+5,
\]

and hence `m=n`.  The maps `F_0,F_1,F_2,...` are pairwise stably
inequivalent.  Together with (3.3), this proves the theorem.

## 6. Scope and local types

The argument uses only:

1. monogenicity of finite etale algebras over an infinite field;
2. the quartic root-engineered incidence;
3. the intrinsic normalized-boundary Fitting divisor.

It does not use a trace form, a rational point on a quadric, a Galois-group
classification, or a restriction on ramification.  Consequently the same
construction applies uniformly to split, `1+3`, `2+2`, cyclic, `D_4`,
`A_4`, and `S_4` quartic algebras, as well as wildly ramified dyadic
quartics after base change to their characteristic-zero local field.

The old trace-chord family remains useful: it gives much smaller coordinate
degrees and a continuous weighted modulus when its quadric is isotropic.
The power-shifted family supplies universality by a discrete boundary
modification invisible on `P=1`.

## 7. Exact regression

Run

```bash
.venv/bin/python scripts/verify_universal_quartic_gauge_multiplicity.py
```

The checker verifies the polynomial formulas, determinant, inverse and
derivative identities for the first power shifts, the fixed quartic fiber,
Newton slopes, Fitting polynomial, and lattice indices.  Monogenicity and
stable normalization functoriality are written inputs.
