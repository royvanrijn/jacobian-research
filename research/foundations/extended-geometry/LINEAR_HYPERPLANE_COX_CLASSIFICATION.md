# Classification of all linear coefficient hyperplanes

No linear coefficient hyperplane can straighten the normalized
three-factor oriented source to affine three-space.  Over an algebraically
closed characteristic-zero field there are exactly three source types,
classified by the root multiplicities of the dual binary cubic, and none is
even stably affine.

This closes the entire alternative-linear-slice route left open by the
first oriented cubic construction.

## 1. General linear slice

Let

\[
 L_i=u_iX+v_iY,\qquad
 P=L_1L_2L_3=\sum_{j=0}^3p_jX^{3-j}Y^j,
\]

and fix a nonzero coefficient functional

\[
 \ell(P)=\lambda_0p_0+\lambda_1p_1+
          \lambda_2p_2+\lambda_3p_3.                 \tag{1}
\]

Consider the normalized source

\[
 Y_\ell=
 \{r_{13}=r_{23}=\ell(P)=1\}\subset\mathbb A^6.      \tag{2}
\]

As in the earlier charts, write

\[
 L_1=aX+bY,\qquad L_3=cX+dY,\qquad ad-bc=1,
\]

\[
 L_2=L_1+tL_3.
\]

Then

\[
 \ell(L_1L_2L_3)=B_{0,\ell}+tB_{1,\ell},             \tag{3}
\]

where

\[
 B_{0,\ell}=\ell(L_1^2L_3),\qquad
 B_{1,\ell}=\ell(L_1L_3^2).
\]

Putting

\[
 H_\ell=(B_{1,\ell}=0)\subset SL_2,\qquad
 Z_\ell=(B_{1,\ell}=0,\ B_{0,\ell}=1),
\]

gives the universal class formula

\[
 [Y_\ell]=[SL_2]-[H_\ell]+\mathbb L[Z_\ell].         \tag{4}
\]

## 2. Only three geometric orbit types

Projectively, a nonzero element of `(Sym^3 k^2)^*` is a binary cubic.
Over an algebraically closed field its roots have one of exactly three
multiplicity patterns:

1. `(3)`, a triple root;
2. `(2,1)`, a double and a simple root;
3. `(1,1,1)`, three distinct roots.

The group `PGL_2` is triply transitive on `P^1`, so each pattern is one
orbit.  A `GL_2` representative can be rescaled to determinant one, and an
overall scalar in `ell` is absorbed by the factor rescaling

\[
 (L_1,L_2,L_3)\longmapsto
 (sL_1,sL_2,s^{-1}L_3),
\]

which preserves `r_(13)=r_(23)=1` and scales `P` by `s`.  Consequently the
three geometric isomorphism types of (2) are represented by

\[
 \ell=p_0,\qquad \ell=p_1,\qquad \ell=p_0+p_3.       \tag{5}
\]

## 3. Triple-root class

For `ell=p_0`,

\[
 B_0=a^2c,\qquad B_1=ac^2.
\]

The locus `H` is the disjoint union `(a=0)` and `(c=0)`.  Each is
`G_m times A^1`, and `B_0` vanishes on both, so `Z` is empty.  Formula (4)
gives

\[
\boxed{
 [Y_{p_0}]
 =\mathbb L^3-2\mathbb L^2+\mathbb L
 =\mathbb L(\mathbb L-1)^2.
}                                                     \tag{6}
\]

## 4. Double-root class

For `ell=p_1`, the exact stratification in the
[quotient-rigidity theorem](ORIENTED_CUBIC_QUOTIENT_RIGIDITY.md) gives

\[
\boxed{
 [Y_{p_1}]=\mathbb L^3+\mathbb L-1.
}                                                     \tag{7}
\]

Here `H` is the disjoint union
`G_m times A^1` and `G_m^2`, while `Z` is the disjoint union
`A^1` and `G_m`.

## 5. Squarefree class

For `ell=p_0+p_3`, the exact calculation in the
[oriented cubic chart](ORIENTED_CUBIC_COX_CHART.md) gives over the original
ground field

\[
 [Y_{p_0+p_3}]
 =
 \mathbb L^3-2\mathbb L-
 [\operatorname{Spec}k[\xi]/(\xi^2-\xi+1)].          \tag{8}
\]

After base change to an algebraically closed field, the final degree-two
scheme splits into two points:

\[
\boxed{
 [Y_{p_0+p_3}]=\mathbb L^3-2\mathbb L-2.
}                                                     \tag{9}
\]

## 6. Stable-affine obstruction

Over `C`, the compactly supported Hodge--Deligne polynomials of the three
types are, with `q=uv`,

\[
 q^3-2q^2+q,\qquad q^3+q-1,\qquad q^3-2q-2.          \tag{10}
\]

None equals `q^3`.  Multiplication by `q^n` cannot turn any of them into
`q^(n+3)`.  Therefore:

### Theorem 6.1

For every nonzero linear functional `ell` on the binary-cubic coefficient
space and every `n>=0`,

\[
\boxed{
 (Y_\ell)_\mathbb C\times\mathbb A^n
 \not\simeq\mathbb A^{n+3}_\mathbb C.
}                                                     \tag{11}
\]

In particular, changing the coefficient hyperplane—including a
nonsymmetric linear choice—cannot produce the desired affine source.  Any
successful modification must be nonlinear in the coefficient variables or
must change the normalization/boundary model itself.

In fact, a homogeneous nonlinear coefficient hypersurface does not help
either.  The
[coefficient-hypersurface obstruction](COEFFICIENT_HYPERSURFACE_COX_OBSTRUCTION.md)
shows that degree `d>1` leaves a free residual `mu_d` Cox stabilizer, whose
Euler-characteristic divisibility is incompatible with affine space.

## 7. Search audit and reproduction

Before extracting the orbit proof, all 6,928 primitive projective
coefficient vectors of height at most five were enumerated.  None had
`A^3` point counts simultaneously over `F_5`, `F_7`, and `F_11`.  This
bounded search is only a regression; Theorem 6.1 follows from the complete
orbit classification and Hodge realization, not from the finite samples.

Run

```bash
.venv/bin/python scripts/verify_linear_hyperplane_cox_classification.py
```

The checker verifies the three `SL_2` stratifications, their motivic
classes, finite-field realizations, the three Hodge obstructions, and the
height-five search.
