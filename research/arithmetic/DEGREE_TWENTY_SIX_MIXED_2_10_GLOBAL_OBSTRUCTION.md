# A global obstruction to the degree-26 reciprocal split (2,10)

The normalized reciprocal degree-26 mixed case with sign degrees
\((2,10)\), in either orientation, has no model over \(\mathbb Q\).
More strongly, its entire four-equation parameter closure has no point
over \(\mathbb Q_{61}\).

This conclusion covers arbitrary denominators and bad reduction. It
uses two complete finite-field censuses and a weighted valuation
argument. No saturation, integrality, or good-reduction hypothesis is
imposed on a hypothetical model.

The square-parameter cases \((4,8)\), \((6,6)\), and the nonsquare
conjugate \((6,6)\) case remain open here. General nonreciprocal
degree-26 pencils are outside this obstruction.

## The exact four-variable system

Use the normalization and sign convention of the
[mixed coupled-ODE note](MIXED_RECIPROCAL_COUPLED_ODE.md). In degree 26,

\[
 h=6,\quad L=25,\quad N(1)=N'(1)=0,\quad \rho^2=\lambda,\quad r=1/\rho.
\]

The residual critical polynomial has paired roots under
\(u\leftrightarrow\lambda/u\). Over the coefficient field with
\(\rho\) adjoined, write

\[
 \Sigma(v)=v^6P(w)Q(w),\quad w=v+v^{-1},\quad
 X+2(w+2)X'=25(D_6-rD_5)P Q,\quad
 Y+2(w-2)Y'=25(C_6-rC_5)P Q.
\]

In the raw sign split \((2,10)\), unequal degrees force
\(\lambda\) to be a square in the coefficient field. Indeed, a
nontrivial quadratic conjugation would interchange the two aggregate
sign factors. Choose the sign of \(\rho\) so that

\[
 P=w-c,\qquad \deg Q=5,\qquad
 Y=Q^2V,\qquad V=w^2+aw+b.
\]

The equation for \(Y\) is

\[
 4(w-2)VQ'+\bigl(V+2(w-2)V'\bigr)Q
             =25(C_6-rC_5)(w-c).                 \tag{1}
\]

The coefficients of a monic \(Q\) are determined from the top down,
with constant divisors \(21,17,13,9,5\). Two compatibility equations
remain, from the coefficients of \(w^0,w^1\). Adjoin the equations

\[
                     \frac{X(c)}{c+2}=0,\qquad
                     \frac{F(r)}{2r}=0.          \tag{2}
\]

The removed factors are nonzero on the good locus. The mixed note
proves that these four equations are necessary and sufficient on its
explicit nondegenerate open set, including the full reconstruction of
\(F\) and \(N\).

Let \(E_0,E_1,E_2,E_3\) be their primitive integer numerators in the
fixed order (constant compatibility, linear compatibility, first
equation in (2), second equation in (2)). Their exact coefficients
are retained in the authoritative
[system file](degree26_mixed_2_10_system.json).
Its companion checker and proof were independently replayed before
the present argument. The new checker uses those exact integer
polynomials, without dropping any variable factor.

For the variable order \((a,b,c,r)\), use weights

\[
                       (1,2,1,1).
\]

Their ordinary total degrees and weighted degrees are

| Equation | Ordinary degree | Weighted degree |
|---|---:|---:|
| \(E_0\) | 6 | 7 |
| \(E_1\) | 6 | 6 |
| \(E_2\) | 11 | 11 |
| \(E_3\) | 12 | 12 |

Write \(H_i\) for the weighted-leading part of \(E_i\). The whole
four-equation closure, including degenerate parameters, is used below.

## Lemma 1. A weighted residue obstruction

Let \(f_1,\ldots,f_s\in\mathbb Z_p[x_1,\ldots,x_d]\) and assign
positive integer weights \(w_1,\ldots,w_d\). Suppose the common zeros
of their weighted-leading forms over \(\mathbb F_p\) consist only
of the origin. Then every common zero of the original equations
over \(\mathbb Q_p\) has integral coordinates.

To prove this, suppose \(x\) is a common zero and

\[
                   m=\min_i\frac{v_p(x_i)}{w_i}<0.
\]

Choose a finite totally ramified extension containing a number \(q\)
with \(v_p(q)=-m>0\). For example, adjoining a root of
\(T^\ell-p\), where \(\ell\) is the least common multiple of the
weights, suffices. This polynomial is Eisenstein, and the residue
field remains \(\mathbb F_p\).

Set \(X_i=q^{w_i}x_i\). All \(X_i\) are integral, and at least one
is a unit. If \(d_j\) is the weighted degree of \(f_j\), then

\[
 q^{d_j}f_j(q^{-w_1}X_1,\ldots,q^{-w_d}X_d)
 =\sum_\alpha c_{j,\alpha}\,
              q^{d_j-\sum_iw_i\alpha_i}X^\alpha.
\]

After reduction, every term of lower weight vanishes. Thus the
nonzero residue vector of the \(X_i\) is a common zero of the
weighted-leading forms over \(\mathbb F_p\), a contradiction.

Consequently, if the original equations also have no affine point
over \(\mathbb F_p\), they have no point over \(\mathbb Q_p\), and
in particular no rational point.

In the present weights, the only possible denominator in \(m\) is
two, so the extension \(\mathbb Q_{61}(\sqrt{61})\) suffices.
It is totally ramified, so no enumeration over \(\mathbb F_{61^2}\)
is required. When the minimum is half-integral, only \(b\) can
remain nonzero after reduction; the explicit term
\(H_1(0,b,0,0)=124215b^3\equiv19b^3\pmod{61}\) also directly excludes
that case.

## Theorem 2. No point over \(\mathbb Q_{61}\)

Complete exact enumeration over \(\mathbb F_{61}\) gives:

| System | Points satisfying its first two equations | Points satisfying all four |
|---|---:|---|
| \(E_0=E_1=E_2=E_3=0\) | 3758 | none |
| \(H_0=H_1=H_2=H_3=0\) | 10921 | only \((0,0,0,0)\) |

These are counts in all of \(\mathbb F_{61}^4\), including every
degenerate parameter and every vanishing coefficient. Applying
Lemma 1 first excludes nonintegral points over \(\mathbb Q_{61}\).
Reducing an integral point would then contradict the first row.
Thus the entire affine closure has no \(\mathbb Q_{61}\) point.

The finite assertions are retained as two independent exact algorithms:

* The Python enumeration loops over all \(a,b,r\), then solves the
  first two equations, both linear in \(c\). It explicitly retains
  the cases where one or both leading coefficients vanish, including
  the case in which every \(c\) must be tested.
* The separate C++ enumeration loops over all four coordinates,
  testing every \(c\) literally. It selects the weighted-leading
  terms independently from the full integer-polynomial input.

Both implementations return the complete point sets and the two
candidate counts in the table. They use exact arithmetic modulo 61,
with no numerical root finding or probability estimate.

## Corollary 3. Rational reciprocal degree-26 models must avoid (2,10)

Every good model in this normalized raw sign class would supply a
point of the four-variable system over its coefficient field.
For a rational model, \(\lambda\) lies in \(\mathbb Q\): its
unique equal-value pairs are Galois invariant, so their common product
is fixed. Unequal aggregate sign degrees make \(\lambda\) a rational
square, as above. The choice of \(\rho\) then makes
\(a,b,c,r\in\mathbb Q\). This contradicts Theorem 2.

Changing the sign of \(\rho\) interchanges the sign factors, so
both orientations are covered. The same reasoning excludes this
normalized good sign class over \(\mathbb Q_{61}\) itself.

This is stronger than the prior real sign-balance exclusion of
\((2,10)\), because no real-root assumption is used. It is also
independent of the prior normalized good-reduction obstruction
at 3: every denominator pattern and bad reduction at 61 is
accounted for by the weighted argument.

## Replay and boundaries

~~~text
.venv/bin/python research/arithmetic/verify_degree_twenty_six_mixed_2_10_global.py --write-receipt research/arithmetic/degree_twenty_six_mixed_2_10_global_verification.json
~~~

Default replay requires SymPy and g++. It reconstructs the
weighted-leading forms from the retained integer equations, runs
both complete enumerations, compares their exact point sets and
candidate counts, and retains source hashes. The C++ source is
[separate and reviewable](verify_degree_twenty_six_mixed_2_10_finite_field.cpp).
The written valuation proof supplies the passage from the finite
censuses to all \(\mathbb Q_{61}\) points.

No complete rational classification of the remaining mixed cases
is established. In particular this does not prove the absence of
rational degree-26 pencils, or of reciprocal models in the remaining
\((4,8)\), \((6,6)\) classes. No external novelty, formal verification
or external human review is asserted.
