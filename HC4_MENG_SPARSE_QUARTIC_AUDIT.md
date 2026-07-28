# Sparse quartic nonlinear descent toward `HC_4`

## Status

After the unit-pivot nonlinear descent in
[the toric audit](HC5_NONLINEAR_TORIC_DESCENT.md), no homogeneous quartic
vertical-Hamiltonian correction supported on at most four monomials both
retains the normalized Meng--Yang collision and has constant nonzero Hessian
determinant.  The same remains true after adjoining an arbitrary scalar
multiple of one cubic monomial or an arbitrary linear combination of two
cubic monomials.

An exact characteristic-zero continuation also excludes all corrections
supported on three cubic monomials.  It reconstructs the 234 quartics over
\(\mathbb Q\), repeats the odd-layer rank calculation over \(\mathbb Q\),
and proves that every remaining rational determinant ideal is a unit ideal.
An additional exact finite-field continuation excludes support four over
\(\mathbb F_{1000003}\); its characteristic-zero promotion remains open.

This is a bounded-support obstruction, not a proof of `HC_4`.  It does not
exclude dense quartics, cubic supports of size at least four, quadratic
renormalizations, degree at least six, or non-coordinate coisotropic
embeddings.

## 1. Collision-first normalization

At the common critical level \(\sigma=-19/2\), the corrected nonlinear
descent has the form

\[
 f(x,y)+2yr+4xs.
\]

The base term \(f(x,y)\) may be removed by a vertical Hamiltonian correction
without changing the unit \(t\)-pivot.  Thus use

\[
 \psi_0=2yr+4xs,\qquad
 H_0=\operatorname{Hess}\psi_0,\qquad \det H_0=64.
\]

The two transported Meng--Yang points are

\[
 p=\left(1,-\frac32,6,\frac{81}{8}\right),\qquad -p.
\]

Let \(h\) be homogeneous quartic.  Since both \(\nabla\psi_0\) and
\(\nabla h\) are odd, the two gradients collide exactly when

\[
 \nabla h(p)=-H_0p.
\]

For a quartic monomial \(c\,w^e\), multiply the \(i\)-th equation by \(p_i\)
and put \(d_e=c\,p^e\).  The entire collision constraint becomes

\[
 \sum_e d_e e
 =
 \left(-\frac{81}{2},18,18,-\frac{81}{2}\right).
\]

It is therefore imposed before any Hessian determinant is expanded.

## 2. Exhaustive supports of size at most four

There are 35 monomials of degree four in four variables.  Exact linear
algebra on every support of size at most four gives:

\[
\begin{array}{c|r|r}
\text{support size}&\text{isolated collision solutions}
&\text{one-parameter collision families}\\ \hline
1&0&0\\
2&3&0\\
3&328&1\\
4&42622&514
\end{array}
\]

Thus only 42,953 isolated coefficient choices and 515 one-parameter
families reach the determinant test.

The isolated choices are reduced modulo \(1000003\).  The exponent minors
and collision denominators are nonzero modulo this prime, so a rational
constant-determinant solution would survive reduction.  Exact determinant
evaluation at nine points rejects every choice.

Each rank-deficient support is then reconstructed over \(\mathbb Q\).
All 515 are genuinely one-parameter families.  For every family, take the
gcd in \(\mathbb Q[\tau]\) of the equations

\[
 \det\bigl(H_0+\operatorname{Hess}h_\tau(w_j)\bigr)-64=0
\]

at the same nine points.  Every gcd is a unit.  Hence no member of any
family has constant Hessian determinant.

No generic 35-parameter determinant is formed: collision linear algebra and
pointwise principal-part rejection precede every exact family calculation.

## 3. One cubic monomial

Let \(h=h_4+\lambda h_3\), where \(h_4\) is one of the preceding collision
quartics and \(h_3\) is a cubic monomial.  Since \(\nabla h_3\) is even, it
does not change the gradient difference between \(p\) and \(-p\).

The degree-eight part of the Hessian determinant is

\[
 \det\operatorname{Hess}h_4.
\]

It cannot be cancelled by \(h_3\).  Modular principal-part evaluation leaves
only 232 of the 42,953 isolated quartics.  Among the 515 families, exact
rational gcds leave only

\[
\begin{array}{c|c}
\text{quartic support}&\text{family parameter}\\ \hline
\{yrs^2,y^2r^2,xyrs,x^2yr\}&-81/8\\
\{xr^2s,xyrs,xy^2s,x^2s^2\}&9/2.
\end{array}
\]

Both resulting quartics have identically zero Hessian determinant.  For each
of the 234 surviving quartics and each of the 20 cubic monomials, determinant
evaluation produces equations in \(\lambda\).  Their modular gcds for the
isolated quartics and rational gcds for the two family members are all units.
Thus no one-monomial cubic correction completes this sparse quartic layer.

## 4. Two cubic monomials

For each of the 234 quartic principal-part survivors, the checker considers
all \(\binom{20}{2}=190\) cubic supports:

\[
 h=h_4+\lambda m_a+\mu m_b.
\]

Thus 44,460 two-cubic coefficient planes are tested.  The determinant grading
avoids a generic bivariate elimination.

1. Degrees seven and one are linear in \((\lambda,\mu)\).  Rank-two systems
   are inconsistent; rank-one systems fix a cubic direction and their full
   univariate determinant gcds are units.
2. In rank zero, degree six is affine-linear in
   \[
   L=\lambda^2,\quad M=\lambda\mu,\quad N=\mu^2,
   \]
   subject to \(M^2=LN\).  Isolated points, transverse conic intersections,
   and conic rulings all fail exact modular full-determinant tests.
3. The higher-nullity cases pass to degree four on the conic.  Modulo
   \(M^2=LN\), this coefficient is linear in
   \[
   L^2,LM,LN,MN,N^2,L,M,N,1.
   \]
   Degree two adds further affine-linear equations in \(L,M,N\).

Only four bivariate coefficient families survive these graded gates.  For
each, the complete determinant evaluation ideal in
\(\mathbb F_{1000003}[\lambda,\mu]\) has Gröbner basis \(\{1\}\).  Because all
collision denominators and exponent minors are nonzero at this prime, no
rational two-cubic correction can have constant Hessian determinant.

## 5. Three cubic monomials over \(\mathbb Q\)

For each of the 234 surviving quartics, the continuation checks all
\(\binom{20}{3}=1140\) triples, hence 266,760 quartic/triple pairs.  The
degree-seven and degree-one determinant coefficients give a linear system in
the three cubic coefficients.  The finite-field discovery calculation and
the characteristic-zero reconstruction have the same rank census:

\[
\begin{array}{c|r}
\text{rank}&\text{quartic/triple pairs}\\ \hline
0&5480\\
1&53364\\
2&130508\\
3&77408
\end{array}
\]

Rank three is inconsistent.  Of the rank-two null lines, 129,588 lie on a
support-at-most-two boundary already handled by the preceding checker.  Each
of the remaining 920 genuine three-support lines has unit full-determinant
gcd.

Of the rank-one null planes, 50,412 likewise reduce to a two-support
boundary.  The full determinant evaluation ideal of every one of the 2,952
genuine planes is the unit ideal.  Finally, the determinant on each of the
5,480 rank-zero coefficient spaces is interpolated in the 35 monomials of
total degree at most four in three parameters; every resulting
three-parameter evaluation ideal is the unit ideal.

In the characteristic-zero pass, the line gcds need at most five evaluation
points, the plane ideals at most seven, and the three-space ideals at most
nine.  All gcds and Gröbner bases are computed over \(\mathbb Q\), so there
is no denominator-lifting caveat.

## 6. Four cubic monomials over the certificate field

There are
\[
234\binom{20}{4}=1{,}133{,}730
\]
quartic/quadruple pairs.  The determinant-degree-seven and degree-one linear
system has rank census
\[
\begin{array}{c|r}
\text{rank}&\text{quartic/quadruple pairs}\\ \hline
0&5430\\
1&79396\\
2&353740\\
3&504818\\
4&190346.
\end{array}
\]

Rank four is inconsistent.  Among the rank-three lines, 504,352 are
support-at-most-three boundaries and 466 are genuinely four-support; every
genuine line has unit full-determinant gcd.  Rank two contains 347,658
boundary planes and 6,082 genuine planes, whose full determinant ideals are
units using at most seven evaluation points.  Rank one contains 71,440
boundary three-spaces and 7,956 genuine three-spaces, whose ideals are units
using at most eight points.

For rank zero, the zero-signature cubic set has size six for 230 quartics and
size twelve for four exceptional quartics.  Their four-subsets give exactly
5,430 four-parameter spaces.  Every full determinant ideal is a unit using
at most twelve evaluation points.

This exhausts cubic support four over \(\mathbb F_{1000003}\).  Unlike
Section 5, this calculation has not yet been repeated over \(\mathbb Q\), so
it is recorded as a finite-field obstruction rather than a
characteristic-zero theorem.

## 7. Remaining search space

A homogeneous correction of odd degree cannot change the gradient
difference between \(p\) and \(-p\), because its gradient is even.  The next
homogeneous collision-carrying layer is therefore degree six.  Within the
quartic chart, the next characteristic-zero extension has four cubic
monomials, while the next finite-field extension has at least five: the
quartic part carries the collision, while the cubic Hessian may cancel
determinant terms.

Dense quartic supports and quadratic renormalizations also remain outside
the present bounded-support theorem.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_hc4_meng_sparse_quartic_obstruction.py
```

Continue through three cubic monomials over the certificate field with:

```bash
.venv/bin/python scripts/verify_hc4_meng_three_cubic_rank_gate.py
```

Promote that calculation to characteristic zero with:

```bash
.venv/bin/python scripts/verify_hc4_meng_three_cubic_characteristic_zero.py
```

Continue through four cubic monomials over the certificate field with:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic.py
```

The two targeted stages are:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic_rank_gate.py
.venv/bin/python scripts/verify_hc4_meng_four_cubic_rank_zero.py
```
