# Three new Airy matrix directions at cap26, and their extension obstruction

Fix the canonical full polynomial Long symbols \(f=(R,T,D,S)\) in
source \((X,Q,P,Z)\), with \(\{P,X\}=\{Z,Q\}=1\).
Use the matrix Moyal product and
\[
 t=XQ,\quad M(t)=\begin{pmatrix}0&1\\t&0\end{pmatrix},\quad
 V=Q\partial_P+X\partial_Z,\qquad A_i=M(t)Vf_i .
\]
Assign source weights \((-1,1,4,2)\), parameter weight 3, image
weights \((-1,2,4,1)\), and matrix-index weight zero.

**Theorem.** No graded \(2\)-by-\(2\) matrix jet
\[
 F_i=f_iI_2+hA_i+h^2B_i+h^3C_i+h^4E_i+O(h^5),\qquad
            \deg B_i,\deg C_i,\deg E_i\leq26,             \tag{1}
\]
satisfies every canonical Weyl commutator through \(h^5\).
The first matrix coefficient \(A\) is fixed, but every entry in all
later images, including the first image, is allowed.

Cap26 admits THREE traceless second-correction freedoms absent at
cap25. The proof retains these directions and the six scalar
freedoms until exact compatibility equations eliminate them.
It uses no cap-preserving gauge assumption.

Consequently a graded full lift with this prescribed first coefficient
requires common spatial bound at least 27. No positive cap27 point,
other first coefficient, scaled connection, nongraded obstruction,
bounded formal-gauge theorem, or DC1/DC2 conclusion is asserted.

## Complete spaces and the positive lower point

Let \(\delta_i=\{f_i,-\}\) and
\(\mathcal D_i=\delta_i+\operatorname{ad}A_i\).
The unknown coefficient at \(h^m\), \(m\geq2\), contributes at
commutator order \(h^{m+1}\) by
\[
                 (\mathcal LC)_{ij}
                         =\mathcal D_iC_j-\mathcal D_jC_i. \tag{2}
\]
The entire scalar basis at order \(m\), for image weight \(w_i\), is
\[
 X^\alpha Q^\beta P^\gamma Z^\nu,\qquad
 \alpha=\beta+4\gamma+2\nu+3m-w_i\geq0,\quad
 \alpha+\beta+\gamma+\nu\leq26.                           \tag{3}
\]
Every monomial and every matrix entry are included.

Use the admitted cap25 positive point \(B^*,C^*\) in
[the frozen matrix-jet source](AIRY_LONG_CAP25_MATRIX_JET_AND_FIXED_EXTENSION_OBSTRUCTION.md).
It consists of actual original-source polynomials, is scalar at
order two and traceless at order three, and satisfies every relation
through \(h^4\). In its notation,
\[
 B^*=\left(g_0-\frac2{27}K_{R^4a}+\frac t2V^2f\right)I_2,
 \qquad K_G=(\{G,f_i\})_i .
\]
Its full coefficient lists are replay inputs, not an assumed
classification of the larger cap26 space.

Define the original global polynomials
\[
 \begin{split}
 B_0&=3X^2P+2(1-3XQ)Z,\\
 Y&=R(Q+3XQ^2-XB_0/3)-3XQ,\\
 J_4&=9Y^4+66Y^3+36(XQ)Y^2+150(XQ)Y
                          +171(XQ)^2-380XQ .
 \end{split}
\]
In particular \(R^4a=R^3(Y+2)\) is polynomial; no rational phase-chart
admission is assumed.

## The full second coefficient has nine parameters

At \(m=2\), (3) has component counts \((80,114,141,102)\),
totaling 437. The scalar map has rank \(431/437\) and the exact
six-dimensional kernel
\[
 K_{R^3},\ K_{R^4S},\ K_{X^3},\
                         K_{R^4a},\ K_{X^4Q},\ K_{R^3J_4}. \tag{4}
\]
The scalar rank minor and exact kernel classification come from
[the complete scalar cap26 certificate](GRADED_LONG_CAP26_EVEN_JET_OBSTRUCTION.md);
the new replay checks them again.

The traceless map has 1311 columns and 8917 rows. Its retained
1308-square minor has determinant \(310392\) modulo \(1000003\).
It has the following three exact independent kernel tuples:
\[
 T_j=(\mathcal D_iG_j)_i,\qquad
 G_1=X^3E_{21},\quad G_2=X^3M(t),\quad
 G_3=X^3(H+2tM(t)),\quad H=\operatorname{diag}(1,-1).      \tag{5}
\]
All their components have degree at most 26; their maximum component
degrees \((R,T,D,S)\) are \((-,18,26,14)\).
The replay checks every coefficient in these tuples. On three
specified coefficient coordinates their values are respectively
\[
                         (-18,0,0),\quad(0,0,36),\quad
                         (0,108,72).
\]
The resulting determinant \(69984\ne0\) proves independence.
Thus the rank is exactly 1308 and (5) exhausts the traceless kernel.

These are polynomial covariant primitives at one order.
They are not automatically removable by a gauge preserving cap26.
For example \(G_3^2=X^6(1+4t^3)I_2\); a formal inverse for
\(I_2+hG_3\) is not a supplied globally polynomial inverse.

The entire second coefficient is therefore
\[
 \begin{split}
 B={}&B^*+
   (uK_{R^3}+vK_{R^4S}+\gamma K_{X^3}
              +eK_{R^4a}+\varepsilon K_{X^4Q}
              +\zeta K_{R^3J_4})I_2\\
   &+s_1T_1+s_2T_2+s_3T_3.                              \tag{6}
 \end{split}
\]
All nine parameters are retained. Relative to a scalar first
particular \(g_0\), the coefficient of \(K_{R^4a}\) is
\(\delta=e-2/27\). The witness records and the checker verifies the
full rational transition from the compact scalar basis (4) to the
normalized scalar cap26 basis.

## Exact h4 equations eliminate the new matrix parameters

At \(m=3\), (3) has scalar counts \((54,80,102,71)\), totaling
307. Its scalar map has rank \(306/307\), with exact kernel
\[
                             K_{R^6}=(0,0,-6R^5,0).
\]
A retained 306-square minor has determinant \(474270\) modulo
\(1000003\). The traceless map has 921 columns and 7230 rows.
Its 921-square minor has determinant \(75001\) modulo this prime,
so it is injective.

The complete known \(h^4\) right side after (6) has fifteen nonzero
parameter-coefficient vectors and 7448 matrix-polynomial terms.
It is linear in the six scalar parameters. The traceless second
coefficients contribute both linear terms and all six quadratic
monomials \(s_as_b\), \(a\leq b\). These quadratic terms are retained.

First project to the ordinary scalar coefficient trace. An exact
rational functional on 234 specified rows annihilates ALL 307
scalar columns and evaluates the right side to
\[
                                  s_3.                  \tag{7}
\]
Thus \(s_3=0\).

After precisely that forced substitution, five retained rational
functionals annihilate ALL 921 traceless columns. Their respective
right-side values are
\[
                         e,\quad\varepsilon,\quad\zeta,
                                  \quad s_1^2,\quad s_2^2.       \tag{8}
\]
Their supports have sizes \(711,711,670,714,714\).
These are exact identities on every parameter coefficient, not
conclusions inferred from a sampled point or a Groebner output.

Over any characteristic-zero field, (7)--(8) force
\[
                   s_1=s_2=s_3=0,\qquad
                   e=\varepsilon=\zeta=0.               \tag{9}
\]
Individual square identities are used; there is no inference from
a sum of squares and no exclusion of complex parameter branches.

The bounded scalar and traceless kernels now show that the ENTIRE
family through \(h^4\) is
\[
 \begin{split}
 B&=B^*+(uK_{R^3}+vK_{R^4S}+\gamma K_{X^3})I_2,\\
 C&=C^*+vMV(K_{R^4S})+\gamma MV(K_{X^3})+wK_{R^6}I_2.    \tag{10}
 \end{split}
\]
All \(u,v,\gamma,w\) remain. Unlike cap25, the \(\gamma\)-term fits
the third-coefficient cap26 space. Direct coefficientwise substitution
verifies (10), including all matrix entries. Thus this is both a
necessary and a sufficient description at order \(h^4\), over
fields; no equality of nonreduced schemes is claimed.

## The entire surviving family fails at h5

The next scalar counts are \((34,54,71,47)\), totaling 206.
Allowing four matrix entries gives 824 unknown columns, including
all 136 first-image entries.

Use the full matrix Moyal coefficient formula, including even
bidifferential orders; no scalar odd-only formula is substituted.
With (10), the known \(h^5\) right side is scalar and has precisely
the five monomials
\[
                            1,\quad\gamma,\quad v,\quad
                                    \gamma^2,\quad v\gamma.     \tag{11}
\]
Their scalar support sizes are \(201,273,96,301,138\).
Every \(u,w\) coefficient cancels identically, without removing
those parameters by a gauge.

A retained 158-row rational scalar functional annihilates every
one of the 206 derivative columns and evaluates (11) to
\[
                                (1,0,0,0,0).            \tag{12}
\]
Lift it by ordinary coefficient trace divided by two. This gives
a 316-row matrix functional. The checker rebuilds ALL 824 full
covariant columns (2), including off-diagonal and first-image
columns, and verifies that the lifted functional annihilates each.
Ordinary commutators of the classical polynomial matrix coefficients
have zero trace here; cyclic quantum matrix trace is not assumed.

Thus every member of (10), with every possible next matrix
coefficient in (3), would require \(0=1\). This proves the theorem.

## Replay and exact boundary

The [witness](airy_first_coefficient_cap26_witness.json) retains the
two traceless rank minors, scalar third-coefficient minor, staged
functionals (7)--(8), and final functional (12), with row locators.
It identifies the positive point and scalar cap26 certificate.
The compact kernels (4)--(5) are regenerated from their original
polynomial formulas.

The [checker](check_airy_first_coefficient_cap26.py) reconstructs
all column maps, every parameter coefficient, and all 824 final
columns. Rational identities use Fraction arithmetic. The installed
python-flint package is required only for exact modular determinant
replay; missing it is an error, never a skipped rank check.
Run from the repository root:

    .venv/bin/python research/dixmier/dixmier-formal-lifting/check_airy_first_coefficient_cap26.py

The [receipt](airy_first_coefficient_cap26_controls.json) pins the
proof, checker, witness, reused arithmetic engine and exact inputs.
Modular nonzero minors supply rational rank lower bounds; the
opposite bounds use the exact independent kernels. Field-valued
parameter exclusions use the rational functionals above.

The theorem fixes precisely \(F_1=MVf\). It does not extend the
cap25 displayed scaled-family corollary to cap26, classify other
connections or first coefficients, or assert anything at larger
caps. Formal verification, independent whole-proof replay, external
review, bounded formal-gauge equivalence, scalar realization and
DC1/DC2 conclusions remain unclaimed.
