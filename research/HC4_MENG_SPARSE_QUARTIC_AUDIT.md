# Sparse quartic nonlinear descent toward `HC_4`

## Status

The support-bounded theorem chain is

\[
 \texttt{HC5T1}\longrightarrow\texttt{HC4MQ1}
 \longrightarrow\texttt{HC4MCK}.
\]

It is now subsumed by the dense mixed-degree theorem `HC4CQ1` in
[`HC4_MENG_DENSE_CUBIC_QUARTIC.md`](HC4_MENG_DENSE_CUBIC_QUARTIC.md).

After the unit-pivot nonlinear descent in
[the toric audit](HC5_NONLINEAR_TORIC_DESCENT.md), no homogeneous quartic
vertical-Hamiltonian correction supported on at most four monomials both
retains the normalized Meng--Yang collision and has constant nonzero Hessian
determinant.  The same remains true after adjoining an arbitrary scalar
multiple of one cubic monomial or an arbitrary linear combination of two
cubic monomials.

Exact characteristic-zero continuations first exclude cubic supports three
and four, then parameterize the complete odd-layer kernel inside the full
20-dimensional cubic space.  Descending spatial-coefficient ideals prove
that no arbitrary homogeneous cubic correction works for any of the 234
quartic principal parts.

The support-three and support-four theorems `HC4MC3` and `HC4MC4` remain
proved exact checkpoints.  `HC4MCK` subsumes them, and `HC4CQ1` in turn
subsumes `HC4MCK`.  Sections 3--7 record useful independent coefficient
regressions; they are not prerequisites for the dense structural theorem.

This bounded-support chain is not a proof of `HC_4`.  Section 8 separately
closes the complete 35-dimensional homogeneous quartic chart, including
arbitrary nondegenerate quadratic renormalizations, by reducing it to the
de Bondt--van den Essen theorem for homogeneous symmetric-Jacobian maps.
Mixed homogeneous degrees, corrections beyond the sextic support bounds
recorded in the companion audit, and non-coordinate coisotropic embeddings
remain open.

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

## 6. Four cubic monomials over \(\mathbb Q\)

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

The characteristic-zero checker reconstructs the same rank and boundary
strata over \(\mathbb Q\).  Exact rational line gcds and Gröbner ideals are
units in every genuine nullspace.  The maximum evaluation-prefix lengths
for ranks three, two, one, and zero are respectively five, seven, eight, and
eleven.  Thus cubic support four is excluded in characteristic zero, without
a denominator-lifting caveat.

## 7. Full cubic kernel over \(\mathbb Q\)

Support-by-support enumeration is unnecessary after Section 6.  For each of
the 234 quartics, form the complete ten-row odd determinant signature on the
20-dimensional cubic coefficient space.  The exact ranks and kernel
dimensions are

\[
\begin{array}{c|r|r}
\text{quartics}&\text{odd rank}&\text{kernel dimension}\\ \hline
229&8&12\\
1&7&13\\
4&4&16.
\end{array}
\]

The four rank-four quartics factor as \(u^3L\), where \(u\) is one of the
four coordinates and \(L\) is linear.  Their odd kernel is the complete
\(H_0\)-harmonic cubic space.  In coordinates paired by
\[
\partial_y\partial_r+\tfrac12\partial_x\partial_s,
\]
the checker uses the triangular chart
\[
h_3=a+u b+u^2c+u^3d
\]
and obtains \(b,c\) by integrating the harmonic equation along the variable
paired with \(u\).

For every kernel, the checker computes
\[
\det\operatorname{Hess}(\psi_0+h_4+h_3)-64
\]
symbolically, extracts all spatial coefficient polynomials, clears rational
denominators, and sends the coefficient ideals to Singular in descending
spatial degree.  The first unit layers are

\[
\begin{array}{c|r}
\text{first unit layer}&\text{quartics}\\ \hline
6&62\\
5&16\\
4&156.
\end{array}
\]

Thus every full cubic-kernel determinant ideal is the unit ideal over
\(\mathbb Q\).  No homogeneous cubic correction, sparse or dense, completes
the Meng descent when the collision quartic is supported on at most four
monomials.

## 8. Complete homogeneous quartic chart

This correction space is also a nonlinear mixed Schur ansatz.  At the
nonzero critical level \(\sigma=-19/2\), make the triangular pivot
translation

\[
 T=t+H(x,y,r,s).
\]

The partial Legendre reduction in \(T\) differs from the old reduction by
\(-\sigma H\).  Thus every quartic \(h_4\) below is realized by taking
\(H=-h_4/\sigma\); it is not merely an unrelated deformation of the
four-variable potential.

Decompose a homogeneous quartic according to its degree in the base
variables \(x,y\) and dual variables \(r,s\):

\[
 h_4=\sum_{a=0}^4 h_{a,4-a}.
\]

All five summands are now allowed simultaneously, giving all 35 homogeneous
quartic monomials.  Put

\[
 \psi=\psi_0+h_4=2yr+4xs+h_4.
\]

Over \(\mathbb C\), the explicit linear substitution

\[
\begin{aligned}
x&=(z_1+iz_2)/\sqrt8,&s&=(z_1-iz_2)/\sqrt8,\\
y&=(z_3+iz_4)/2,&r&=(z_3-iz_4)/2
\end{aligned}
\]

sends the quadratic term to

\[
 \frac12(z_1^2+z_2^2+z_3^2+z_4^2).
\]

Consequently the transformed gradient has the form

\[
 F(z)=z+\nabla\widetilde h_4(z),                       \tag{8.1}
\]

where \(\nabla\widetilde h_4\) is homogeneous cubic and \(JF\) is symmetric.
If \(\det\operatorname{Hess}\psi\) is a nonzero constant, then \(F\) is a
Keller map.  The dimension-four theorem of
[de Bondt and van den Essen](https://doi.org/10.1016/S0022-4049(03)00223-8)
therefore makes \(F\) a polynomial automorphism.  An invertible linear
change transports equality of gradients and distinctness of points, so the
Meng collision cannot survive.

The same argument starts from any nondegenerate quadratic form: over an
algebraic closure its Hessian is congruent to the identity.  Thus arbitrary
quadratic renormalizations do not reopen a homogeneous quartic chart.

> **Complete homogeneous-quartic obstruction.**  In four variables, no
> nondegenerate quadratic potential plus an arbitrary homogeneous quartic
> can have constant nonzero Hessian determinant and a gradient collision.

This is theorem `HC4HQ1`.  It uses the cited low-dimensional theorem as an
external structural input; it is not a new proof of that theorem.

### Independent dense mixed regression

Before invoking the complete theorem, an exact coefficient calculation
closed the genuinely mixed space

\[
 h_{1,3}+h_{2,2}+h_{3,1}
\]

has dimension \(8+9+8=25\).  This is a dense coefficient space, not a
support-bounded search.  Two one-sided enlargements have dimension 30:
adjoin either the arbitrary pure-base term \(h_{4,0}\) or the arbitrary
pure-dual term \(h_{0,4}\).

For each of these three spaces, impose the collision equation

\[
 \nabla h_4(p)=-H_0p
\]

and the spatial-degree-two part of the determinant identity,

\[
 \operatorname{tr}\!\left(
   \operatorname{adj}(H_0)\operatorname{Hess}(h_4)
 \right)=0.                                           \tag{8.2}
\]

The combined linear system has rank 14.  It leaves respectively 11, 16, and
16 parameters.  After this exact substitution, the remaining spatial
coefficient ideals of

\[
 \det\operatorname{Hess}(\psi_0+h_4)-64
\]

have 262, 273, and 273 generators.  Singular computes the unit ideal over
\(\mathbb Q\) in all three cases.  Therefore:

> **Dense mixed-quartic regression.**  No homogeneous quartic correction
> containing arbitrary mixed bidegrees \((1,3),(2,2),(3,1)\), with either
> pure sector also allowed arbitrarily, retains the normalized Meng
> collision and has constant nonzero Hessian determinant.

That coefficient theorem is `HC4MDQ`.  Its calculation does not allow
\(h_{4,0}\) and \(h_{0,4}\) simultaneously, but `HC4HQ1` subsumes it and
closes their interaction.

## 9. Remaining search space

A homogeneous correction of odd degree cannot change the gradient
difference between \(p\) and \(-p\), because its gradient is even.  It can,
however, interact with a quartic in the Hessian determinant.  The complete
cubic--quartic interaction is now excluded without a support restriction by
`HC4CQ1`.  The older `HC4MCK` coefficient theorem remains an independent
exact regression over the 234 sparse quartic principal parts.

The next homogeneous collision-carrying layer is degree six.

The parallel
[sparse sextic audit](HC4_MENG_SPARSE_SEXTIC_AUDIT.md) excludes every
sextic-only collision carrier supported on at most four monomials.  Its
mixed theorem `HC4MQS6` also excludes every zero-gradient sextic of that
support size over the 234 quartic principal parts treated here.  The joint
theorem `HC4JQS4` allows sextic terms alongside every quartic component
within combined support four and excludes every such genuinely mixed
correction.  The later theorem `HC4E46` in
[`HC4_SOURCE_DUAL_BIGRADING.md`](HC4_SOURCE_DUAL_BIGRADING.md) removes
this support bound and closes the full homogeneous degrees \(4+6\).

Thus the pure homogeneous quartic chart, its quadratic renormalizations,
the full mixed degrees \(3+4\), and the full mixed degrees \(4+6\) are
closed.  The remaining coordinate-chart problem begins with simultaneous
cubic and sextic interaction, including degrees \(3+4+6\), followed by
non-coordinate coisotropic embeddings.

## Reproduction

The shortest exact replay of the canonical chain is:

```bash
.venv/bin/python scripts/verify_hc5_nonlinear_toric_descent.py
.venv/bin/python scripts/verify_hc4_meng_sparse_quartic_obstruction.py
.venv/bin/python scripts/verify_hc4_meng_full_cubic_kernel.py
```

The last command requires `Singular` on `PATH`.

Verify the exact reduction of the complete homogeneous quartic chart to the
de Bondt--van den Essen theorem with:

```bash
.venv/bin/python scripts/verify_hc4_meng_full_quartic_reduction.py
```

Check the complete dense cubic--quartic reduction with:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_dense_cubic_quartic_reduction.py
```

Check the dense mixed-quartic theorem, including its two one-sided pure
enlargements, with:

```bash
.venv/bin/python scripts/verify_hc4_meng_dense_mixed_quartic.py
```

This command also requires `Singular` on `PATH`.

The following support-three and support-four commands are retained as
targeted historical regressions; they are not prerequisites for the
full-kernel checker.  Check three cubic monomials over the certificate field
with:

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

Promote the four-cubic calculation to characteristic zero with:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic_characteristic_zero.py
```
