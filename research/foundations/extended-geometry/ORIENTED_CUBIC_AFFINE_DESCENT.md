# Affine descent of the oriented cubic Cox chart

The oriented three-linear-factor chart admits an affine-space descent for a
special hyperplane: normalize the tangent coefficient `p_1`, then forget the
ordering of two factors.  The quotient source is the normalized
linear--quadratic factorization slice, hence affine three-space, and the
descended multiplication map is the foundational cubic Keller map.

This solves the affine-transfer problem only after identifying the two
dicritical divisors.  It does not produce an affine-space Keller map with
two independently marked dicritical primes.

The [quotient-rigidity theorem](ORIENTED_CUBIC_QUOTIENT_RIGIDITY.md)
further proves that the tangent oriented source has class
`L^3+L-1`, is not stably affine, and that this involution is the unique
nontrivial factor-permutation symmetry preserving the dicritical pair as a
set.

Work over a characteristic-zero field `k`.

## 1. Tangent-hyperplane oriented chart

Let

\[
 L_i=u_iX+v_iY,\qquad
 P=L_1L_2L_3=p_0X^3+p_1X^2Y+p_2XY^2+p_3Y^3,
\]

and `r_(ij)=u_iv_j-v_iu_j`.  Replace the symmetric hyperplane used in the
first oriented chart by

\[
 m_{\mathrm{tan}}=p_1.
\]

This trilinear polynomial is irreducible.  Viewed as a polynomial linear in
`v_3`, its coefficient is `u_1u_2` and its constant term is
`u_3(v_1u_2+u_1v_2)`; these are coprime.

Define

\[
 Y^{\mathrm{or}}=
 \{r_{13}=r_{23}=p_1=1\}\subset\mathbb A^6,           \tag{1}
\]

and the oriented coefficient target

\[
 T^{\mathrm{or}}=
 \{p_1=1,\ D^2=\operatorname{Disc}(P)\}.              \tag{2}
\]

As before,

\[
 \Psi^{\mathrm{or}}:Y^{\mathrm{or}}\longrightarrow
 T^{\mathrm{or}},\qquad
 (L_1,L_2,L_3)\longmapsto(P,D=r_{12})                 \tag{3}
\]

is polynomial and has constant hypersurface-residue Jacobian `-1/2` on the
smooth target locus.  It has generic degree three and two dicritical
divisors corresponding to the missing `R_(13)` and `R_(23)` branches.

## 2. Linear--quadratic quotient

Put

\[
 L=L_3,\qquad Q=L_1L_2.
\]

Then

\[
 \operatorname{Res}(L,Q)
 =\operatorname{Res}(L_3,L_1)
  \operatorname{Res}(L_3,L_2)
 =r_{13}r_{23}=1,                                    \tag{4}
\]

and

\[
 [LQ]_{X^2Y}=p_1=1.                                  \tag{5}
\]

Consequently the involution

\[
 \iota:(L_1,L_2,L_3)\longmapsto(L_2,L_1,L_3)         \tag{6}
\]

has quotient

\[
 S=
 \left\{
 (L,Q):
 \operatorname{Res}(L,Q)=1,\quad
 [LQ]_{X^2Y}=1
 \right\}.                                           \tag{7}
\]

The normalized linear--quadratic factorization theorem gives an explicit
polynomial isomorphism

\[
 \boxed{S\simeq\mathbb A^3.}                         \tag{8}
\]

It is the slice used by the foundational map.

On the target, `iota` changes

\[
 D=r_{12}\longmapsto-r_{12}.
\]

Thus

\[
 T^{\mathrm{or}}/(D\mapsto-D)
 =
 \{p_1=1\}\simeq\mathbb A^3.                         \tag{9}
\]

## 3. Cartesian discriminant diagram

For a linear form and a quadratic form,

\[
 \operatorname{Disc}(LQ)
 =\operatorname{Res}(L,Q)^2\operatorname{Disc}(Q).
\]

On (7), this becomes

\[
 \operatorname{Disc}(P)=\operatorname{Disc}(Q).      \tag{10}
\]

Ordering the two linear factors of `Q` is exactly the double cover obtained
by adjoining a square root of `Disc(Q)`.  Hence, after normalization, the
square

```text
 Y^or  ---------------- Psi^or ---------------->  T^or
  |                                                |
  | quotient by L1 <-> L2                         | D -> -D
  v                                                v
 S ~= A^3  -------------- G ------------------->  A^3_(p0,p2,p3)
```

is Cartesian over the squarefree locus and extends by normalization across
the generic double-root divisor.

The bottom map is multiplication

\[
 G(L,Q)=LQ.                                          \tag{11}
\]

The ambient coefficient-resultant determinant is

\[
 \det D(\operatorname{coeff}(LQ),\operatorname{Res}(L,Q))
 =-\operatorname{Res}(L,Q)^2.
\]

Taking the residue along (7) gives

\[
 \boxed{\det DG=-1.}                                 \tag{12}
\]

Under the explicit coordinates of the normalized factorization slice, `G`
is linearly equivalent to the foundational noninjective cubic Keller map.

## 4. What happens to the dicritical divisors

Upstairs, the two escaping primes are

\[
 \mathcal D_{13}\quad\text{and}\quad\mathcal D_{23}.
\]

The involution (6) exchanges them:

\[
 \iota(\mathcal D_{13})=\mathcal D_{23}.             \tag{13}
\]

Therefore their images in the affine quotient form one prime.  The finite
collision divisor `R_(12)` is the ramification divisor of the vertical
double cover but is not a ramification divisor of either horizontal Keller
map on its smooth chart.

This proves:

> The natural affine descent of the oriented cubic Cox chart is the
> foundational Keller map, and it necessarily merges the two Cox
> dicritical divisors.

Thus the affine-space gap is not “find any descent”; that descent is now
explicit.  The remaining problem is to descend or modify the chart while
retaining two independently distinguishable boundary primes.

## 5. Arithmetic consequence

The oriented top map separates the sign of the cubic discriminant but has
cyclic degree-three fiber behavior.  The bottom foundational map forgets
that sign and has full symmetric cubic monodromy.  Neither is a permutation
on the finite fields covered by their split-fiber profiles.

The affine quotient therefore does not combine the two desired features:
dicritical independence is lost, and exceptional permutation monodromy is
not gained.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_oriented_cubic_affine_descent.py
```

The checker verifies the resultant and discriminant identities, the two
quotient involutions, the ambient determinant, the normalized-slice
relations, and the exact linear equivalence with the foundational map.
