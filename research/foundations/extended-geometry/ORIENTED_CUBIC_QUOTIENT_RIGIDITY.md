# Quotient rigidity of the oriented cubic Cox chart

The tangent-hyperplane oriented source is not affine three-space, and no
stabilization by ordinary affine coordinates makes it affine space.  Its
natural affine quotient is therefore a genuine descent rather than a hidden
change of coordinates.  Moreover, factor-permutation symmetry forces every
nontrivial quotient which preserves the pair of dicritical boundary
valuations to exchange and identify them.

This closes the natural symmetric-quotient route to an affine Keller map
with two independently marked dicritical divisors.  It does not rule out a
nonsymmetric affine modification.

Work over a characteristic-zero field `k`.

## 1. The tangent source

Let

\[
 L_i=u_iX+v_iY,\qquad
 r_{ij}=u_iv_j-v_iu_j,
\]

and let `p_1` be the coefficient of `X^2Y` in `L_1L_2L_3`.  The oriented
tangent source is

\[
 Y_{\rm tan}^{\rm or}
 =
 \{r_{13}=r_{23}=p_1=1\}\subset\mathbb A^6.          \tag{1}
\]

Write

\[
 L_1=aX+bY,\qquad L_3=cX+dY.
\]

The first normalization says

\[
 ad-bc=1,                                            \tag{2}
\]

and every `L_2` satisfying `r_(23)=1` is uniquely

\[
 L_2=L_1+tL_3.                                      \tag{3}
\]

On (2), direct coefficient extraction gives

\[
 p_1=B_0+tB_1,
\]

\[
 B_0=a^2d+2abc=a(3ad-2),\qquad
 B_1=bc^2+2acd=c(3ad-1).                            \tag{4}
\]

Thus (1) is the hypersurface

\[
 B_0+tB_1=1
\]

over `SL_2`.

## 2. Exact motivic class

Put

\[
 H=(B_1=0)\subset SL_2,\qquad
 Z=(B_1=0,\ B_0=1).
\]

Projection to `SL_2` has one point over `SL_2-H`, no point over `H-Z`, and
an affine line over `Z`.  Hence

\[
 [Y_{\rm tan}^{\rm or}]=[SL_2]-[H]+\mathbb L[Z].    \tag{5}
\]

Equation (4) splits `H` into two disjoint pieces:

\[
 H_1=(c=0)\simeq\mathbb G_m\times\mathbb A^1,
\]

\[
 H_2=(3ad=1)\simeq\mathbb G_m^2.
\]

On `H_1`, the equation `B_0=1` gives `a=d=1`, so
`Z_1` is `A^1`.  On `H_2`, it gives `a=-1`, `d=-1/3`, and
`bc=-2/3`, so `Z_2` is `G_m`.  Consequently

\[
 [H]=(\mathbb L-1)\mathbb L+(\mathbb L-1)^2,
\qquad
 [Z]=\mathbb L+(\mathbb L-1).
\]

Since `[SL_2]=L^3-L`, (5) becomes

\[
 \boxed{
 [Y_{\rm tan}^{\rm or}]=\mathbb L^3+\mathbb L-1.
 }                                                   \tag{6}
\]

Over `C`, its compactly supported Hodge--Deligne polynomial is therefore

\[
 E_c(Y_{\rm tan}^{\rm or};u,v)=(uv)^3+uv-1.         \tag{7}
\]

This differs from `(uv)^3`; more strongly, after multiplying by `(uv)^n`
it still differs from `(uv)^(n+3)` for every `n>=0`.  Therefore

\[
 \boxed{
 Y_{\rm tan,\mathbb C}^{\rm or}\times\mathbb A^n
 \not\simeq\mathbb A^{n+3}_{\mathbb C}
 \quad\text{for every }n\geq0.
 }                                                   \tag{8}
\]

Ordinary affine stabilization cannot straighten the oriented source.

## 3. Rigidity of factor-permutation quotients

The two dicritical boundary valuations are labeled by the collision edges

\[
 e_{13}=\{1,3\},\qquad e_{23}=\{2,3\}.               \tag{9}
\]

The symmetric group acts faithfully on the three edges of the complete
graph on `{1,2,3}`.  The pointwise stabilizer of both displayed edges is

\[
 \operatorname{Stab}(e_{13})
 \cap\operatorname{Stab}(e_{23})
 =
 \{1\}.                                             \tag{10}
\]

The setwise stabilizer of the two-element set `{e_(13),e_(23)}` is

\[
 \{1,(12)\}.                                        \tag{11}
\]

It follows that a factor-permutation quotient which retains both boundary
valuations as individually labeled divisors must be trivial.  The only
nontrivial factor-permutation symmetry preserving their union is `(12)`,
and it exchanges them.

For the tangent chart, `(12)` is precisely the regular involution

\[
 (L_1,L_2,L_3)\longmapsto(L_2,L_1,L_3).
\]

Its quotient is the normalized linear--quadratic slice `A^3`, while the
target quotient forgets the sign of the oriented discriminant and is the
coefficient space `A^3`.  The two dicritical primes have one image in that
quotient.  Thus the affine descent proved in
[the affine-descent theorem](ORIENTED_CUBIC_AFFINE_DESCENT.md) is the unique
nontrivial factor-permutation quotient preserving the dicritical pair as a
set, and it necessarily merges the pair.

## 4. Equivariant affine stabilization does not repair the quotient

One might try to retain the lost sign by adjoining an affine representation
`V` before taking the involution quotient.  This also has a local
obstruction.

At the generic ramification point of

\[
 Y_{\rm tan}^{\rm or}\longrightarrow S\simeq\mathbb A^3,
\]

choose a transverse parameter `delta`.  The involution acts by
`delta -> -delta`.  Decompose a linear affine suspension as

\[
 V=V^+\oplus V^-,
\]

where the involution is `+1` on `V^+` and `-1` on `V^-`.  On the normal
tangent space to the fixed locus, the number of negative eigenvalues is

\[
 1+\dim V^-.                                        \tag{12}
\]

In characteristic zero, the local quotient by an involution is smooth only
when the involution is a pseudoreflection, hence only when (12) equals one.
Therefore smoothness forces

\[
 V^-=0.                                             \tag{13}
\]

For example, with one attempted sign coordinate `z`, the local invariant
ring contains

\[
 A=\delta^2,\qquad B=\delta z,\qquad C=z^2,\qquad
 B^2=AC,                                            \tag{14}
\]

the singular quadratic cone.

If (13) holds, the involution is trivial on every new affine coordinate and

\[
 (Y_{\rm tan}^{\rm or}\times V)/\langle\iota\rangle
 \simeq S\times V\simeq\mathbb A^{3+\dim V}.
\]

But then the quotient still places the two exchanged dicritical valuations
in one orbit and one prime.  Hence a linear equivariant affine suspension
cannot simultaneously keep the quotient smooth and retain nontrivial
orientation data.

## 5. Consequence and remaining route

The following three operations cannot produce an affine-space Keller map
retaining both marked dicriticals:

1. a polynomial change of coordinates on the tangent source;
2. stabilization by any number of affine coordinates; or
3. a nontrivial quotient induced by permutation of the three factors; or
4. a linear equivariant affine stabilization followed by that quotient.

The first two are excluded by (7)--(8), the third by (10)--(11), and the
fourth by (12)--(14).

The surviving geometric route must therefore change the boundary model:
for example, a nonsymmetric affine modification, an enlargement carrying
additional nontrivial boundary characters, or a different nonsymmetric
monodromy cover.  This result does not claim that such a modification is
impossible.

Changing to another linear coefficient hyperplane does not help.  The
[linear-hyperplane classification](LINEAR_HYPERPLANE_COX_CLASSIFICATION.md)
proves that every nonzero linear slice belongs to one of three binary-cubic
orbit types and that all three normalized sources are stably non-affine.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_oriented_cubic_quotient_rigidity.py
```

The checker verifies the tangent coefficient formula, the two strata in
`H` and `Z`, the class `L^3+L-1`, exact point counts over three prime
fields, the pointwise and setwise stabilizers of the two dicritical edges,
and the singular local invariant ring produced by an added sign coordinate.
