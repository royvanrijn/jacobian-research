# Triple-root affine completion and its determinant obstruction

The triple-root coefficient slice has an affine-space completion preserving
two distinct torus boundary divisors.  This is the first normalized source
for which the source-side affine gap disappears.  However, the minimal
polynomial clearing of the oriented coefficient map creates three additional
affine Jacobian divisors.  Its target is already affine three-space and has
no nonconstant units, so target reparametrization cannot absorb them.

This isolates the next obstruction on the target side.

## 1. The toric source

Use the triple-root slice

\[
 r_{13}=r_{23}=p_0=1.                               \tag{1}
\]

Write

\[
 L_1=aX+bY,\qquad
 L_3=cX+dY,\qquad ad-bc=1,
\]

and `L_2=L_1+tL_3`.  The condition `p_0=1` is

\[
 ac(a+tc)=1.
\]

It forces `a,c` to be nonzero and gives

\[
 d=\frac{1+bc}{a},\qquad
 t=\frac{1-a^2c}{ac^2}.                             \tag{2}
\]

Therefore

\[
\boxed{
 Y_{p_0}\simeq\mathbb G_m^2\times\mathbb A^1
}
                                                            \tag{3}
\]

with coordinates `(a,c,b)`.  Its evident completion is

\[
 \overline Y_{p_0}=\mathbb A^3_{a,b,c},              \tag{4}
\]

and the missing primes `(a=0)` and `(c=0)` remain distinct.

## 2. Laurent oriented coefficients

Substitution of (2) gives

\[
 p_1=\frac{-a^2c+3bc+2}{ac},                         \tag{5}
\]

\[
 p_2=
 \frac{-2a^2bc^2-a^2c+3b^2c^2+4bc+1}{a^2c^2},       \tag{6}
\]

\[
 p_3=
 \frac{-b(bc+1)(a^2c-bc-1)}{a^3c^2},                \tag{7}
\]

and the oriented discriminant coordinate is

\[
 D=r_{12}=\frac{1-a^2c}{ac^2}.                      \tag{8}
\]

These are regular on the torus source (3), but have poles on its affine
completion.

## 3. Minimal polynomial clearing

Clear exactly the displayed denominator monomials:

\[
\begin{aligned}
 A&=ac\,p_1=-a^2c+3bc+2,\\
 B&=a^2c^2p_2
   =-2a^2bc^2-a^2c+3b^2c^2+4bc+1,\\
 C&=a^3c^2p_3
   =-b(bc+1)(a^2c-bc-1),\\
 E&=ac^2D=1-a^2c.
\end{aligned}                                       \tag{9}
\]

The four polynomial outputs satisfy exactly

\[
 A^2-3B-E^2+E-1=0.                                  \tag{10}
\]

In characteristic zero, (10) solves for `B`, so the cleared target is

\[
 \overline T_{\rm clr}\simeq\mathbb A^3_{A,C,E}.     \tag{11}
\]

The induced polynomial map on the affine source completion is

\[
 \Phi_{\rm clr}:\mathbb A^3_{a,b,c}
 \longrightarrow\mathbb A^3_{A,C,E}.                \tag{12}
\]

## 4. Exact determinant

Direct differentiation gives

\[
\boxed{
 \det\frac{\partial(A,C,E)}{\partial(a,b,c)}
 =
 -6abc\,(bc+1)(a^2c-bc-1).
}                                                     \tag{13}
\]

The desired completion divisors `a=0` and `c=0` occur, but so do the three
additional affine divisors

\[
 b=0,\qquad bc+1=0,\qquad a^2c-bc-1=0.              \tag{14}
\]

They are not target units: the target (11) is affine three-space, whose
unit group consists only of constants.  A polynomial automorphism of the
target has constant nonzero Jacobian, so postcomposing (12) cannot remove
any factor in (13).

### Theorem 4.1

The minimal monomial clearing of the triple-root oriented chart extends to
a polynomial map between affine three-spaces, but it is not Keller.
Moreover, no polynomial reparametrization of its cleared affine target can
make it Keller.

Thus filling the two source torus boundaries alone is insufficient.  A
successful target absorption must mix the boundary coordinates
non-minimally—before pole clearing—or introduce a nonlinear complete
intersection whose residue form cancels all five factors in (13).

The latter repair is explicit.  The
[affine-source oriented construction](AFFINE_SOURCE_TRIPLE_ROOT_COX_MAP.md)
adjoins the cleared determinant as a target orientation coordinate.  It
restores constant residue Jacobian and turns four of the five target branch
components into distinct dicritical images.

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/verify_triple_root_affine_completion.py
```

The checker verifies the source parametrization, every Laurent coefficient,
the cleared target equation, the affine-three target identification, and
the five-factor determinant.
