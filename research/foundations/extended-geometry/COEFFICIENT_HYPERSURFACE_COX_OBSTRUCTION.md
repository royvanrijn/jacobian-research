# Obstruction for every homogeneous coefficient hypersurface

No homogeneous hypersurface in the binary-cubic coefficient space can be
used as the third normalization boundary to produce affine three-space.
Degree one is excluded by the complete binary-cubic orbit classification.
Every higher degree has a residual finite Cox stabilizer acting freely, and
compactly supported Euler characteristic excludes affine space even after
ordinary stabilization.

## 1. Higher-degree coefficient boundary

Let

\[
 P=L_1L_2L_3=p_0X^3+p_1X^2Y+p_2XY^2+p_3Y^3
\]

and let

\[
 F(p_0,p_1,p_2,p_3)
\]

be a nonzero homogeneous polynomial of degree `d>=1`.  Consider

\[
 Y_F=
 \{r_{13}=r_{23}=F(P)=1\}\subset\mathbb A^6.         \tag{1}
\]

On `(P^1)^3`, the two resultant boundaries and the pullback of `F=0`
have class rows

\[
 (1,0,1),\qquad(0,1,1),\qquad(d,d,d).               \tag{2}
\]

Their class matrix has Smith normal form

\[
 \operatorname{diag}(1,1,d),                        \tag{3}
\]

and determinant `-d`.  Thus the normalization is integrally unimodular
only for `d=1`.

## 2. The residual Cox stabilizer

For `zeta in mu_d`, act by

\[
 (L_1,L_2,L_3)
 \longmapsto
 (\zeta L_1,\zeta L_2,\zeta^{-1}L_3).               \tag{4}
\]

Both selected resultants remain unchanged:

\[
 r_{13}\longmapsto r_{13},\qquad
 r_{23}\longmapsto r_{23}.
\]

The product cubic scales by

\[
 P\longmapsto\zeta P,
\]

so homogeneity gives

\[
 F(P)\longmapsto\zeta^dF(P)=F(P).
\]

Hence (4) is the residual `mu_d` action predicted by (3).

The action is free on (1).  Indeed, `r_(13)=1` forces `L_1` and `L_3` to be
nonzero.  A point fixed by (4) would satisfy `zeta L_1=L_1`, and therefore
`zeta=1`.

## 3. Euler-characteristic obstruction

Assume `d>1` and base change to `C`.  If `Y_F` were affine three-space, its
free finite action would give a degree-`d` finite etale quotient

\[
 \mathbb A^3\longrightarrow\mathbb A^3/\mu_d.
\]

Compactly supported Euler characteristic is multiplicative for a finite
etale cover:

\[
 \chi_c(\mathbb A^3)
 =
 d\,\chi_c(\mathbb A^3/\mu_d).                       \tag{5}
\]

But `chi_c(A^3)=1`, so (5) would express `1` as `d` times an integer, which
is impossible.

The same argument applies after adjoining any number of affine coordinates
with trivial `mu_d` action, because

\[
 \chi_c(\mathbb A^{3+n})=1.
\]

Thus, for `d>1`,

\[
\boxed{
 (Y_F)_\mathbb C\times\mathbb A^n
 \not\simeq\mathbb A^{3+n}_\mathbb C
 \quad(n\geq0).
}                                                     \tag{6}
\]

## 4. Completion with the degree-one theorem

When `d=1`, the residual group in (4) is trivial, but the
[linear-hyperplane classification](LINEAR_HYPERPLANE_COX_CLASSIFICATION.md)
applies.  Its triple-root, double-root, and squarefree source classes are

\[
 \mathbb L^3-2\mathbb L^2+\mathbb L,\qquad
 \mathbb L^3+\mathbb L-1,\qquad
 \mathbb L^3-2\mathbb L-2,
\]

and none becomes affine after stabilization.

Combining the two cases proves:

### Theorem 4.1

For every nonconstant homogeneous coefficient polynomial `F` and every
`n>=0`, the normalized source (1) does not become affine space after
ordinary stabilization.

Therefore the remaining construction cannot come from merely replacing the
linear coefficient hyperplane by a nonlinear projective coefficient
hypersurface.  It must instead change the boundary-class presentation—for
example by adding primitive Cox coordinates before taking the level,
altering the selected resultant tree, or using a boundary modification not
pulled back from one coefficient hypersurface.

Adding such coordinates as an affine vector bundle is also insufficient.
The
[Cox vector-bundle obstruction](COX_VECTOR_BUNDLE_STABILIZATION_OBSTRUCTION.md)
shows that the degree-one Hodge defect survives every affine-bundle fiber
and that arbitrary linear fiber weights retain the higher-degree free
`mu_d` action.

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/verify_coefficient_hypersurface_cox_obstruction.py
```

The checker verifies the Smith form (3), the weights in (4), freeness on
the normalized source, and the Euler-characteristic divisibility
obstruction.
