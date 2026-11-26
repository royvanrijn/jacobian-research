# The conductor-18 seventeen-row blind locus

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
invisible to all seventeen rows. The four-point boundary
restriction first implies finiteness. The exact local
determinants below then show that the locus is empty over
\(\overline{\mathbb Q}\).

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

## Candidate weighted residues above 331

The [two normal charts on \(L_3\)](CM_NEWTON_331_L3_NORMAL_CHARTS.md)
exclude every point of the third weighted boundary component
where \(v\) and \(w\) are nonzero in residue. Together with the
earlier boundary analysis, any algebraic blind marker must
reduce, at every place above \(331\), to one of four weighted
points. They are the two endpoints of \(L_3\) and the two
off-component points on \(L_2\) with \(v/\ell=2/17\) or
\(5/56\). Both [the \(v=0\) endpoint](CM_NEWTON_331_L3_VZERO_ENDPOINT.md)
and [the \(w=0\) endpoint](CM_NEWTON_331_L3_WZERO_ENDPOINT.md)
of \(L_3\) are excluded. Only the two off-component points
remain. In particular the degree-17 branch \(E=0\) is empty
for the full system, since it could only meet \(L_3\).
The earlier transverse restrictions remain in force on the
two off-component points.

## Finiteness of the actual blind locus

**Corollary (finite seventeen-row exceptional algebra).**
For the degree-\((5,1)\) family, only finitely many geometric
phases can fail formal coherent detection in rows
\(0,\ldots,16\) with correction degree at most two.
The normalized residual algebra has length at most \(5670\).
Its horizontal model at \(331\) is finite flat and splits after
completion into at most two local algebras, each of even rank.
At any place above \(331\), an algebraic seventeen-row blind
point must have both \(c\) and \(e\) of negative valuation.

Let \(R=\mathbb Z_{(331)}\), and put the normalized blind scheme
\(Z\) in the chart \(T=1\) of
\(\mathbb P_R(2,1,2,1,2,1)\) with coordinates
\((r,\ell,h,v,w,T)\). The exact \(r=0\) certificate makes
\(r\) a unit on \(Z\). Take the schematic closure of \(Z\)
over \(\mathbb Q\), then its schematic closure
\(\mathscr Z\) over \(R\). In homogeneous coordinates this is
the saturation of the weighted homogenized residual ideal by
\(331rT\), followed by \(\operatorname{Proj}\). Saturation by
\(rT\) takes the closure from the original chart. Contraction
from \(\mathbb Q\), or saturation by \(331\), makes the affine
coordinate rings \(331\)-torsion free. Thus \(\mathscr Z\) is
flat over the DVR \(R\), and its underlying space is the closure
of its generic fibre.

Take a closed point of the special fibre and a generic point
of a horizontal component specializing to it. The closed
generic-fibre points specializing to the chosen closed point
are dense in that generic component. In particular one can
choose such a point in \(Z\), since \(Z\) is dense and open
in the projective closure over \(\mathbb Q\).
[Stacks, Lemma 37.52.4](https://stacks.math.columbia.edu/tag/053U).
Its residue field is a number field. The specialization is
realized by a DVR of that number field at a place above \(331\).
[Stacks, Lemma 28.5.10](https://stacks.math.columbia.edu/tag/054F).
The four-point theorem applies at this place. Consequently the
closed points of the special fibre, and hence its support, lie
in the four allowed points. This argument uses the proper
horizontal closure, including possible nonreduced fibres.

The special fibre has dimension at most zero. Proper upper
semicontinuity gives the same bound for the generic fibre.
The finite scheme \(Z\) is then already closed in the projective
generic fibre. Both fibres of \(\mathscr Z\) are finite, so
properness gives finiteness over \(R\).
[Stacks, Lemmas 37.30.5](https://stacks.math.columbia.edu/tag/0D4I)
[and 37.44.1](https://stacks.math.columbia.edu/tag/02LS).
Solving the correction and angular pivots and adjoining
\(b\) with \(b^2=r\) gives precisely the exceptional phases.
The sign rotation of the top transfers this to the plus convention.

### Length bound

The prefix \(S_7,S_8,S_9\) is the rank-84 finite free algebra
over \(\mathbb Q[v,w]\) above. The same homogeneous
integrality test shows that all its points over \((v,w)=(1,1)\)
are integral above \(331\). On that fibre the prefix and \(S_{10}\)
generate the unit ideal after reduction at \(331\).
Thus every characteristic-zero root of the prefix on that
fibre has \(S_{10}\ne0\).

Each irreducible component of the finite flat prefix dominates
the plane. Indeed nonzero base polynomials act as
nonzerodivisors, so the minimal primes contract to zero, and
a finite dominant map is surjective. In particular \(S_{10}\)
does not vanish identically on any component. The first four
rows cut to dimension at most one. Since the full residual
locus is finite, a rational linear combination of \(S_{11}\)
through \(S_{16}\) can avoid every one-dimensional component
of this intersection. The five chosen equations have weighted
degrees at most \(6,7,8,9,15\).

To count scheme lengths, pull back by \(r=x_r^2,\ h=x_h^2,\ w=x_w^2\)
and leave \(\ell,v\) fixed. This finite free map multiplies
length by eight and converts the degree bounds to ordinary
degrees. The five equations form a regular sequence wherever
the final scheme is nonempty. The first three form a complete
intersection by finiteness over the plane, and the last two
avoid the components at the respective stages. Flatness of
the squaring map preserves the regular sequence. In the
total-degree filtration of a quotient by an initial part of
the sequence, multiplication by the next regular element
of degree at most \(d\) embeds \(F_{n-d}\) into the part of
its ideal in \(F_n\). Taking the leading term of the Hilbert
polynomial bounds the new degree by \(d\) times the old one.
This proves the affine Bézout bound with scheme lengths even
when the equations have extra projective components at infinity.
Thus \(\operatorname{length}Z\le 6\cdot7\cdot8\cdot9\cdot15/8=5670\).
Adjoining \(b\) multiplies that length by two.

### An affine algebra and its two candidate residue factors

The two remaining points have \(v,w\ne0\) in weighted
residue. Thus \(D_+(v)\) contains the horizontal closure:
its proper closed complement has empty special fibre.
In the original affine chart \(T=1\) this says that at any
place above \(331\), \(\nu(e)=2\nu(c)<0\). Consequently
the seventeen shifts suffice whenever either \(c\) or \(e\)
is integral at one such place, with arbitrary \(B\).

Use the weighted chart \(v=1\). The normalized coordinate
\(\ell/v\) has distinct residues \(17/2,56/5\) at the two
points. The characteristic polynomial of multiplication
by this coordinate consequently splits by Hensel lifting
into the two residue clusters. Each factor degree is its
special-fibre length, since each point is rational over
\(\mathbb F_{331}\).
[Stacks, Lemma 10.153.4](https://stacks.math.columbia.edu/tag/04GH).

The weighted homogenizations \(\widehat S_k\) use only even
powers of \(T\). Put \(\tau=T^2\) and define the explicit
affine presentation over \(R\) by
\[
 A_0=R[r,\ell,h,w,\tau]\big/
 \left(\bigl(\widehat S_7|_{v=1},\ldots,\widehat S_{16}|_{v=1}\bigr)
                    :(331r\tau)^\infty\right),
\]
where \(T^{2j}\) means \(\tau^j\) in the displayed rows.
The algebra of the horizontal closure is
\(A_0[T]/(T^2-\tau)\). To see the equality, before saturation
this quadratic presentation is free on \(1,T\) over the
ring with \(\tau\). Inverting \(331rT\) is equivalent to
inverting \(331r\tau\), and contraction acts separately on
the two coefficients. It follows that \(A_0\) is finite
flat of rank at most \(2835\). The local summands of the
horizontal algebra have even rank, twice the corresponding
ranks of \(A_0\). This presentation preserves the horizontal
closure and halves the degree of the algebra. The local
determinants below show that both factors of \(A_0\) are zero.

### An exact univariate eliminant

**Proposition (an eliminant for the normalized slope).**
There is an explicit nonzero polynomial \(\Delta(s)\) of degree
at most \(835\) such that \(\Delta(\ell/v)\) vanishes in the
characteristic-zero residual algebra. It can be given by a
determinant with \(\Delta(1)=12473\bmod65521\).

Set \(s=\ell/v,\ x=v^2\). For \(S_{7+i}\), let
\(\epsilon_i=i\bmod 2\). Substitute \(\ell=sv\) and divide the
row by \(v^{\epsilon_i}\), which is allowed since \(v\) is a unit
on the residual locus. The resulting row \(A_i\) lies in
\(\mathbb Q[s,r,h,w,x]\); its degree in \(r,h,w,x\) is at most
\(d_i=\lfloor(6+i)/2\rfloor\). Let \(H_i\) be its ordinary
homogenization to degree \(d_i\) with new variable \(t\).
Let \(M(s)\) have as rows the coefficients of all monomial
multiples of the \(H_i\) of total degree nine. Its \(715\)
columns are the degree-nine monomials in \(r,h,w,x,t\).

Over \(\mathbb Q(s)\) consider two local quotients on the
homogeneous system. In the chart \(x=1\), the reference point
has \(r=h=t=0\) and \(w=s(8-s)/17\). Solve for \(r\) and the
displacement of \(w\) by two linear pivots, keeping \(h,t\)
with weights \(1,2\) through weight eight. In the chart
\(w=1\), use \(r=x=t=0,\ h=17\), solve for \(r\), and retain
terms through ordinary degree four in the three free variables.
The implicit jets are computed over \(\mathbb Q[s]\) by
scaling free coordinates by the square of each pivot determinant
and verifying each adjugate division exactly. The resulting
truncated quotients have dimensions \(25-10=15\) and
\(35-20=15\) over \(\mathbb Q(s)\).
These relation ranks are preserved at \(s=1\bmod65521\).
Stacking the specialization of the degree-nine restriction
map with the two relation matrices gives rank \(60\) in
\(60\) columns, so the restriction onto their direct sum has
rank \(30\) over \(\mathbb Q(s)\).

The rows of \(M(s)\) annihilate under the restriction, as does
\(t^9\) by both jet cutoffs. Thus their span lies in a
\(685\)-dimensional kernel. At \(s=1\bmod65521\), select
row and column pivots from the transpose and the original
matrix, respectively. Order the columns by increasing
\(x\)-exponent, and the rows by increasing
\(\epsilon_i-2b\), with \(b\) the \(x\)-exponent of their
monomial multiplier. Break ties in the original order.
This gives \(685\) independent rows and columns and a
minor with determinant \(12473\) in the residue field.
Take \(\Delta(s)\) to be the same determinant over
\(\mathbb Q[s]\).

The selected rows must span the whole kernel over
\(\mathbb Q(s)\). Solving for the coordinates of \(t^9\)
in these rows and multiplying by \(\Delta\) gives
\(\Delta(s)t^9\) in the homogeneous ideal, by the adjugate
formula. In any nonzero matrix entry, the exponent of \(s\)
is at most \(2(a-b)+\epsilon_i\), with \(a\) the column's
\(x\)-exponent. The sum in each determinant term is
therefore bounded by the sum of these row and column
contributions. It is \(835\) for the selected minor.
Dehomogenizing at \(t=1\) gives the asserted ideal
membership. Equivalently, after putting \(x=1,t=\tau\)
in the horizontal chart and inverting \(\tau\) over
\(\mathbb Q\), the same identity adds \(\Delta(\ell)=0\)
to the explicit presentation of \(A_0\).

The final exclusion uses two other minors chosen for their
Taylor orders at the two candidate residues.

    sage -python research/foundations/extended-geometry/verify_cm_newton_slope_eliminant.py

## Excluding the final two residue discs

**Theorem (universal seventeen-row coherent detection).**
For the fixed top with two opposite fourth-order roots and
\(\deg B\le5,\ D=-z^{-1}+c+ez\), the response in shifts
\(0,\ldots,16\) with two coherent inputs and correction
degree at most two has no kernel vector with a nonzero
coherent component over any characteristic-zero field.
The normalized exceptional algebra is the zero ring.

Use \(s,x,t,H_i\) from the slope eliminant above. The leading
parts of \(H_0,H_2\) by total degree in \(r,h\) are homogeneous
of degrees three and four in \(\mathbb Q[r,h]\). They form a
regular sequence whose quotient has dimension twelve.
Take their monic Gröbner basis in degree reverse lexicographic
order and lift each basis element using its homogeneous
multipliers applied to \(H_0,H_2\). The leading monomials of
the lifts are
\[
 h^6,\quad rh^4,\quad r^2h^2,\quad r^3.
\]
The monic lifts have coefficients integral at \(331\).
For the filtration by degree in \(r,h\), the associated
graded ideal of \((H_0,H_2)\) is generated by the two leading
parts. Indeed, the syzygy of the two coprime leading parts is
generated by their Koszul relation. Subtracting the corresponding
relation of \(H_0,H_2\) lowers the maximal \(r,h\)-degree of
a representation with leading cancellation. Descending proves
the associated graded assertion. Hence reduction by the
lifts leaves a free \(\mathbb Q[s,w,x,t]\)-module with the
twelve normal monomials in \(r,h\) as a basis. The reductions
preserve total degree in \(r,h,w,x,t\).

In degree nine the quotient has dimension \(394\), so
\(321\) dimensions of the original row span come from
\(H_0,H_2\). Reduce the monomial multiples of the other \(H_i\)
to this quotient. Multipliers using only normal monomials
suffice by the same reduction. The resulting matrix has
\(492\) rows and \(394\) columns. Its rank over
\(\mathbb Q(s)\) is \(685-321=364\), and its row span contains
the reduction of \(t^9\). For any nonzero \(364\)-minor
\(\Delta_\alpha(s)\), the adjugate therefore gives
\(\Delta_\alpha(s)t^9\in(H_0,\ldots,H_9)\). Dehomogenizing
gives \(\Delta_\alpha(s)\in(A_0,\ldots,A_9)\), where the
\(A_i\) here are the divided residual rows.

Select one such minor in each of the two remaining discs.
At \(s=\alpha+u\), over \(\mathbb F_{331}[[u]]\) select a
maximal rank pivot block of the constant matrix, take its
Schur complement, and divide that complement by \(u\).
Repeat until \(364\) pivots have been selected. The following
counts record the number of pivots at each order of \(u\).
The selected square submatrix has these same counts.

| \(\alpha\) | Pivot counts in successive orders | Order of the selected determinant |
|---|---|---:|
| \(17/2\) | \(359,2,0,0,0,2,0,1\) | \(19\) |
| \(56/5\) | \(361,1,0,1,0,0,1\) | \(10\) |

The replay verifies the same pivot counts over \(\mathbb Q\)
on the selected submatrices. For the first pivot block \(A\)
and the remaining blocks \(B,C,D\), solve \(AX=B\) and
\(YA=C\) through order three with exact rational arithmetic.
The identity
\[
 CA^{-1}B-CX-Y(B-AX)=(C-YA)A^{-1}(B-AX)
\]
gives the exact Schur complement \(D-CA^{-1}B\) through
order seven. Taking its successive constant ranks computes
the orders in the table. The original entries and the
row scalings used to clear denominators are integral at
\(331\), with the scalings units. Thus both the
characteristic-zero order and the residue order of
\(\Delta_\alpha\) at \(\alpha\) equal the displayed integer
\(m_\alpha\). The quotient
\[
 U_\alpha(s)=\Delta_\alpha(s)/(s-\alpha)^{m_\alpha}
\]
has coefficients in \(\mathbb Z_{(331)}\) and
\(U_\alpha(\alpha)\) is a unit there.

Every remaining algebraic blind point would have
\(s\equiv\alpha\) at a place above \(331\).
Then \(U_\alpha(s)\) is a unit at that place, whereas
\(\Delta_\alpha(s)=0\) by ideal membership. This forces
\(s=\alpha\). The exact affine slice identities above
exclude both equalities. Together with the two endpoint
exclusions this eliminates every possible specialization.
The localized rational residual ideal is therefore the
unit ideal by the Nullstellensatz. The \(b=0\) branch is
detected by shifts zero and one, using the three-pole
reduction. Base extension gives the stated conclusion over
any characteristic-zero field. The sign rotation gives the
same conclusion for the plus top.

For fixed algebraic coefficients, these finite exact
response tests persist at all sufficiently large good
primes. In the difficult class \(p\equiv7\bmod36\), the
possible correction polynomial has degree at most two,
and a nonzero coherent input must therefore give a nonzero
moment among \(m+kp\) with \(0\le k\le16\). In the class
\(p\equiv25\bmod36\), shift zero already suffices.
The [shorter-prefix calculation](CM_NEWTON_FOURTEEN_ROW_LOCUS.md)
gives the exact minimum \(N_{\min}=14\) for the uniform initial window,
using the example below for necessity.

    sage -python research/foundations/extended-geometry/verify_cm_newton_slope_residue_discs.py

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

## Replay

Run the exact computations with SageMath 10.9.

    sage -python research/foundations/extended-geometry/verify_cm_newton_seventeen_row_locus.py

The slice and lower-bound calculations use the checker above.
The weighted boundary and the determinant arguments have their
own linked checkers. The [fixed-boundary stratum](CM_NEWTON_SPHERE_ATTACHMENT_STRATUM.md)
also has a geometric detection proof for the full response there.
