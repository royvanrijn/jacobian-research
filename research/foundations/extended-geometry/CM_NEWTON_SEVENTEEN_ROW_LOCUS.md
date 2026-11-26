# Projections of the conductor-18 seventeen-row blind locus

Use the [three-pole reduction](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#three-extra-pole-inputs-for-odd-degree-five-and-even-degree-one)
with the minus-top convention. On \(b\ne0\) the early rows
solve for the three correction coefficients and for the
coefficients of angles 3, 4 and 5. The remaining equations are
\(S_7,\ldots,S_{16}\) in
\[
              (r,\ell,h,v,w)=(b^2,b_1/b,b_2/b,c,e).
\]
The solved input has coherent coordinate one. A common zero
with \(r\ne0\) therefore gives a formal coherent direction
invisible to all seventeen rows. The question whether this
locus is empty over \(\overline{\mathbb Q}\) remains open.

## A single norm for the common-root condition

Suppose \(A/R\) is free of rank \(d\) and \(\rho_0,\ldots,\rho_{s-1}\)
are residuals in \(A\). In any basis write \(M_\rho\) for multiplication
by \(\rho\), and put
\[
 q(T)=\sum_{j=0}^{s-1}T^j\rho_j,\qquad
 N(T)=\det M_{q(T)} .
\]
On a geometric fibre its finite algebra is a product of local
Artinian algebras. On one such factor multiplication by \(q\)
acts as its residue scalar on each piece of the maximal-ideal
filtration. Thus the norm is the product of each scalar
\(q(T)(x)\) to the length of the corresponding factor. It is
identically zero exactly when a single fibre point annihilates
every \(\rho_j\). This includes nonreduced fibres. Applying the
adjugate identity to the vector of the unit also gives
\(N(T)\in(\rho_0,\ldots,\rho_{s-1})A[T]\). The coefficients
and any \(d(s-1)+1\) distinct rational evaluations generate
the same ideal by interpolation.

The rank-twelve chart has nine remaining rows and needs at
most 97 evaluations. In the normalized coordinates \(S_7\)
is cubic in \(r\), with leading coefficient
\(-23055471104/584486340375375\). Dividing by this constant
gives a rank-three chart over \(\mathbb Q[\ell,h,v,w]\)
and reduces the evaluation bound to 25. The exact
\(r=0\) slice certificate below means this projection can use
the whole cubic fibre.

## A quadratic chart and two proper component cuts

**Proposition (two projected cuts for the full residual system).**
The normalized blind locus projects under elimination of \(w,j\)
into the union of the two proper cuts specified below.
Over their admissible locus the projection has finite fibres.

Set \(j=h-17w\). In the variables \(r,\ell,j,v,w\), the
first residual is quadratic in \(w\). Write it as \(Kw^2+Lw+C\),
with \(K=-12829696r/194828780125125\). This gives a rank-two
algebra on \(r\ne0\), over \(\mathbb Q[r,r^{-1},\ell,j,v]\).
Write the residual remainders as \(a_k+d_k w\).
For a remainder \(a+d w\), \(K\) times its norm is
\(Ka^2-Lad+Cd^2\).
For all the remainders to have a common root, it is necessary
and sufficient that every individual norm and every bracket
\(\Delta_{i,k}=a_i d_k-a_k d_i\) vanish. Indeed the coefficient
pairs then have span at most one. When there is a nonzero
pair, its norm test finds the shared root. Both a bracket
and the cleared norm belong to the residual ideal on this chart.

The distinction from separate norms occurs in this family.
Fix \(r=\ell=1\). At \(j=37,v=13\bmod107\) the two
cleared norms for rows 8 and 9 are zero. Their Jacobian in
\(j,v\) is \(69\), while their bracket is \(33\) and the
cleared norm of the sum is \(55\). The simple point lifts
to algebraic characteristic-zero parameters with the same
two norms zero and no shared root.

Use \(N_8\) for the cleared norm of row 8, \(D=\Delta_{8,9}\),
and \(E=\mathrm{Res}_j(a_8,d_8)\). The checker certifies
irreducibility of \(N_8,D,E\) over \(\mathbb Q\) by complete
degree-preserving reductions at \(65521\). In particular, the
first prefix projection \(N_8=0\) is irreducible, and \(D=0\)
gives a proper next cut. On the chart with
\(\ell(\ell-8v)\ne0\), the \(j\)-degrees of \(N_8,D\)
are 7 and 6, with leading coefficients constant multiples
of \(\ell(\ell-8v)\) and \(\ell\).

Let \(P=\mathrm{Res}_j(N_8,D)\), computed at these fixed degrees.
It has the factorization \(P=\ell E^2 F\) up to a nonzero
rational factor. Both \(j\)-degrees drop at \(\ell=0\),
so the Sylvester determinant there is zero. To see the
double factor along \(E\), work at its generic point.
The common root of \(a_8,d_8\) is unique and is simple for \(D\).
These properties follow from the full-degree linear
Euclidean remainder and the nonzero derivative of \(D\) at the
shared root \(j=19\) over \((r,\ell,v)=(1,1,85)\bmod107\).
The checker verifies that the specialization is on a nonzero
reduction of \(E\), and that the needed leading coefficients
and the slope of that linear remainder are units there.
Consequently those quantities are nonzero at the generic
point over \(\mathbb Q\). In the completed discrete valuation
ring at \(E\), lift the simple root of \(D\). At its lift both
\(a_8,d_8\) are divisible by \(E\), so \(N_8\) is divisible by
\(E^2\) there. The leading coefficient of \(D\) is a unit
at \(E\), hence \(E^2\) divides the resultant.

For weights \((2,1,2,1)\) on \((r,\ell,j,v)\),
\(N_8,D,E\) have weights \(16,13,17\).
The resultant weight is at most \(6\cdot16+7\cdot13-2\cdot6\cdot7=103\).
The exact modular quotient \(F\) has weight 68 and is
irreducible at \(65521\). The weight bound and
\(\mathbb Z_{(65521)}\)-primitivity of the factors give
irreducibility over \(\mathbb Q\).

Every remaining point projects, after eliminating \(w,j\),
to one of these two proper cuts in \((r,\ell,v)\):
\[
\begin{array}{ll}
 E=0,& \mathrm{Res}_j(D,N_9)=0,\\
 F=0,& \mathrm{Res}_j(N_8,\Delta_{8,10})=0 .
\end{array}
\]
For the second equation the row-10 remainder may be scaled
by a power of \(K\) to clear its denominators. In each
case the indicated resultant is nonzero modulo the
component equation, certified by reduction at \(65521\).
The pairs of equations therefore cut the two projected
surfaces to subsets of dimension at most one. The highest
weighted part of \(E\) at \(r=0\) is a \(331\)-unit times
\(\ell^2(\ell-8v)^{15}\). Thus the two off-\(L_3\) boundary
points in the [weighted boundary analysis](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#the-projective-boundary-at-331)
cannot be the reduction of a point on \(E=0\). The \(j\)
and \(w\) fibres are finite, using the nonzero \(j\)-leading
coefficient of \(N_8\) and the quadratic chart. Both
\(\ell=0\) and \(\ell=8v\) are excluded for the full system
by the slice certificates. The displayed two cuts are
necessary equations and may contain points that are
removed by the remaining residuals.

## Finite projection to the even coefficients

**Proposition (finite projection and a two-coefficient integrality criterion).**
The normalized residual locus is finite over a proper closed
subset of the \((v,w)\)-plane. At a place above \(331\) its
points cannot have both \(v\) and \(w\) integral. In particular,
the seventeen-shift theorem holds for every fixed algebraic
phase with \(c,e\) integral at one place above \(331\), regardless
of the coefficients of \(B\).

Filter \(r,\ell,h\) with weights \(2,1,2\), treating \(v,w\)
as coefficients. The leading parts of \(S_7,S_8,S_9\)
in that filtration have weights \(6,7,8\). They are obtained by
setting \(v=w=0\) and keeping the highest weights.
Over both \(\mathbb Q\) and \(\mathbb F_{331}\), their
homogeneous ideal has only the origin as a geometric zero.
The checker verifies this by a Gröbner computation in three
variables. In characteristic zero the leading forms are a
regular sequence and their homogeneous quotient has dimension 84.
The same standard monomials give a free basis of the
filtered quotient over \(\mathbb Q[v,w]\). In fact, a
cancellation of the leading terms of a combination is a
homogeneous syzygy of the leading forms. Their Koszul
syzygies lift to the Koszul syzygies of the three actual rows,
so subtracting them lowers the weight of the cancellation.
This identifies the associated graded quotient and proves
freeness by induction on weight.

At \(v=w=0\) over \(\mathbb F_{331}\), the prefix and row 11
generate the unit ideal. If \(v,w\) are integral, a root of
the first three rows must have \(r,\ell,h\) integral.
Otherwise scale by a parameter of minimum weighted valuation,
after a finite ramified extension if necessary. Its reduction
gives a nonzero triple annihilating the three homogeneous
leading parts with \(v=w=0\). This contradicts the homogeneous
test. Thus at \(v=w=0\) the characteristic-zero prefix has
no root in common with row 11, by reduction at \(331\).
The norm of row 11 in the rank-84 algebra is consequently
nonzero at the origin. It vanishes on the projection of the
full residual locus, which proves the dimension assertion.
If \(v,w\) are integral at any place above \(331\), integrality
of the other three coordinates and the existing \(331\) unit
identity for the full system give a contradiction. The
case \(b=0\) was already detected by shift 1.

## Exact affine slices

The following identities are computed by taking monomial
multiples up to the stated degree, selecting pivots modulo
\(331\), and verifying all coefficients of the constant-row
identity over \(\mathbb Q\). On a slice \(\ell=\alpha v\),
divide odd-parity rows by \(v\) and set \(u=v^2\); the
case \(v=0\) of that slice uses the existing \(c=0\) result.
For these slices the cap is ordinary degree in the transformed
coordinates. For \(r=0\) and \(h=17w\) the cap uses the
remaining original weights.

| Slice | Last shift | Degree cap | Pivot rank |
|---|---:|---:|---:|
| \(r=0\) | 16 | 18 | 792 |
| \(h=17w\) | 16 | 18 | 1080 |
| \(\ell=8v\) | 13 | 8 | 408 |
| \(\ell=(17/2)v\) | 16 | 9 | 680 |
| \(\ell=(56/5)v\) | 16 | 9 | 682 |
| \(\ell=8v,\ h=17w\) | 11 | 6 | 60 |

These exact affine exclusions do not assert the same for
nearby projective boundary directions with unequal pole orders.

## Remaining weighted residues above 331

The [two normal charts on \(L_3\)](CM_NEWTON_331_L3_NORMAL_CHARTS.md)
exclude every point of the third weighted boundary component
where \(v\) and \(w\) are nonzero in residue. Together with the
earlier boundary analysis, any algebraic blind marker must
reduce, at every place above \(331\), to one of four weighted
points. They are the two endpoints of \(L_3\) and the two
off-component points on \(L_2\) with \(v/\ell=2/17\) or
\(5/56\). The degree-17 branch \(E=0\) can use only the first
two. The earlier transverse restrictions remain in force on
the two off-component points.

## Exact coherent blindness through shift 12

**Proposition (shift 13 is necessary in any uniform initial window).**
There is an algebraic phase in the degree-\((5,1)\) family
with a nonzero coherent input invisible to shifts
\(0,\ldots,12\). Shift \(13\) detects that input. On the
slice \(\ell=v=0\) this window suffices everywhere.

On \(\ell=v=0\) the odd-parity residuals \(S_8,S_{10},S_{12}\)
vanish identically. Take primitive integral representatives of
the specialized rows \(S_7,S_9,S_{11},S_{13}\).
At \(11\) the point \((r,h,w)=(2,8,2)\) annuls the first
three. Their Jacobian determinant is \(7\), and the value of
the primitive row 13 is \(10\). Hensel's lemma gives an
algebraic common zero, with \(r\ne0\), detected by row 13.
A rational constant-span computation of degree 12,
with rank 441, checks that these four rows already
generate the unit ideal on the slice. For \(b=0\), shift 1
again applies.

There is also an example with all five normalized
coordinates nonzero in which shift 12 is the first
detecting row.

For the normalized residual equations take primitive integral
representatives of \(S_7,\ldots,S_{12}\). In their reductions
at \(23\), the point
\[
                   (r,\ell,h,v,w)=(17,5,12,6,17)
\]
annuls \(S_7,\ldots,S_{11}\). Their five-variable Jacobian
determinant is \(18\), and the primitive row-12 value is \(21\).
Multivariate Hensel lifting gives a \(23\)-adic phase with exactly
those five residual equations zero. Its coordinates are algebraic
over \(\mathbb Q\): over their coordinate field the invertible
Jacobian forces all coordinate differentials to vanish.
Since \(r\) is a unit, take \(b^2=r\) and undo the three
angular pivots, as well as the correction pivots. This gives
the stated allowed phase and a formal coherent direction with
coordinate one. The row-12 unit certifies detection at that
shift. The lifting argument concerns marker inputs; it imposes
no condition identifying them with actual Frobenius inputs.

A shorter prefix admits a uniform count on one slice.
In the cubic algebra at \(h=v=w=1\) write \(P,Q\) for the
remainders of rows 8 and 9. The norm of \(P\) is degree 21
in \(\ell\), has 21 distinct nonzero roots, and is coprime
to the constant term of the cubic. At every root the
multiplication matrix of \(P\) has rank two, since rank at
most one would make the derivative of the determinant zero.
Hence exactly one root of the cubic is common with \(P\).
The determinant with columns \(P,rP,Q\) is coprime to the
norm of \(P\) by an exact rational Bézout identity. Thus
these 21 cases are invisible through row 8 and detected
by row 9.

## Replay and remaining conditions

Run the exact computations with SageMath 10.9.

    sage -python research/foundations/extended-geometry/verify_cm_newton_seventeen_row_locus.py

The two projected cuts and the finite plane projection give
necessary loci for an algebraic blind marker through shift 16.
Deciding the residual equations on those loci remains a
characteristic-zero calculation. Coherent detection on the
[fixed-boundary stratum](CM_NEWTON_SPHERE_ATTACHMENT_STRATUM.md)
controls the full response there, with a separate finite
cutoff. It does not prescribe the seventeen-shift cutoff.
