# Small Bernstein degrees and constant-leading orders through forty-one

**Theorem PolynomialWeylBernsteinDegreeThroughFortyOne.** Let \(k\)
have characteristic zero and let \(P,Q\in A_1(k)\) have a nonzero scalar
commutator. If either member has Bernstein degree at most \(41\), the
pair generates \(A_1(k)\). The same conclusion holds if either member
has a nonzero constant differential leading coefficient and order at
most \(41\).

**Corollary PolynomialWeylIndividualLeadingDegreeFortyTwo.** In a
nongenerating polynomial Weyl pair with nonzero scalar commutator,
each member \(P\) satisfies
\[
                  \operatorname{ord}P+\deg\operatorname{lc}P\ge42.
\]

The proof uses the
[Bernstein bound through thirty-three](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_THIRTY_THREE.md),
the [constant-leading bound through thirty-three](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md)
and the
[nonconstant-leading degree calculation](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#nonconstant-leading-factors-and-bernstein-degree).
The crossings below are for actual subrectangular polynomial pairs.
Their Bernstein-degree gcd need not be the global minimum.

## From a sharp degree bound to a crossing edge

For the induction let \(n=34,36,38,40\) in succession, and suppose
both assertions hold through \(n-1\). In any nongenerating pair,
a member with nonconstant leading differential coefficient then has
order plus leading-coefficient degree at least \(n\). Indeed, keep
the member fixed and subtract powers of it from its mate whenever
the mate's order is divisible by its order. For a nondividing
resulting pair use the infinity shear, which gives exactly that
sum as Bernstein degree without changing the order or leading
coefficient. If instead the mate's new positive order \(j\)
divides the member's order \(i>j\), its leading coefficient is
nonconstant by the leading commutator identity. That same identity
gives \(i+\deg\operatorname{lc}P=(i/j)(j+\deg\operatorname{lc}Q)\).
Induction on \(i\) gives the leading-degree bound in this case too.

Suppose a nongenerating pair has one member \(P\) of
Bernstein degree \(n\). Extend constants to an algebraic
closure and scale the mate so that \([Q,P]=1\).
Proper generation persists under the constant extension.
[Guccione–Guccione–Valqui, Theorem 4.1](sources/guccione-guccione-valqui-1111.6100v3.pdf)
at weight \((1,1)\) gives a homogeneous quadratic \(F\)
whose Poisson bracket with the leading symbol of \(P\)
is that symbol, up to sign. If \(F\) is a square of a
linear form, its Hamiltonian action on fixed homogeneous
degree is nilpotent. Making the linear form one of the
canonical coordinates verifies this by reducing the exponent
of the other coordinate at each step. This cannot yield a
nonzero eigenvalue. Otherwise the quadratic is a product
of independent linear forms. A linear Weyl change makes it
proportional to \(XY\), whose Hamiltonian action has distinct
eigenvalues on the monomials of a fixed total degree. Thus
the leading symbol of \(P\) becomes a single monomial.
Theorem 4.1(3) excludes a diagonal monomial.

If that monomial is pure, orient it as \(x^n\).
The differential order is then below \(n\). A constant
leading coefficient contradicts the induction hypothesis.
A nonconstant leading coefficient contradicts the
individual leading-degree bound, because its term has
total degree below \(n\). In the mixed case, orient the
axes so that the differential exponent in the leading
monomial is smaller than its \(x\)-exponent. The actual order
is also less than \(n\), because the leading monomial is
mixed. Again its leading coefficient is nonconstant.
The individual leading-degree bound forces its term
to be the unique Bernstein leading monomial, and gives
\[
                i=\operatorname{ord}P,\qquad
                i+\deg\operatorname{lc}P=n .
\]
Reduce the mate by powers of this fixed \(P\) in
differential order. A mate whose reduced order \(j<i\)
divides \(i\) is impossible: the preceding identity
would give \(i+\deg\operatorname{lc}P\ge(i/j)n\ge2n\).
For the resulting nondividing orders \(ga,gb\), write
\(h\) for the common leading factor and \(d=\deg h\).
The nonconstant-leading theorem gives
\[
                   n=a(g+d),\qquad g+d\ge16 .
\]
Since \(n<48\) and \(a\ge2\), it follows that \(a=2\).
Thus \(g+d=n/2,\ g<d\), and the leading coefficients
are proportional to \(x^{2d},x^{bd}\), with \(b\) odd
and at least three.

Apply the infinity shear from the nonconstant-leading
degree calculation. Its bounds
\(\deg_x p_i\le2(d-g)+i,\ \deg_x q_j\le b(d-g)+j\)
make the Fourier presentations \(X=D,\ Y=-x\)
subrectangular, with corners \((2g,2d),(bg,bd)\).
The entire face of maximal \(X\)-exponent in either
operator is its single corner, since the leading
differential coefficients are pure powers. In particular
the two endpoints at \((1,0)\) are aligned in ratio
\(2:b\) and have negative \((1,-1)\) weight.

The following verifies the crossing conditions for this
pair. Increase \(s\) from zero to one in the directions
\((1,-s)\). Use primitive integer directions to describe
the finite list of edges. Just to the right of zero the
leading points are the corners above. At each subsequent
edge the highest \(Y\)-exponent endpoints are the leading
points from just before the edge. They are aligned in
ratio \(2:b\). By Guccione–Guccione–Valqui, Lemma 7.1,
both weights are positive and the edge bracket is zero.
The coprime power comparison in their Theorem 2.11(2)
then aligns the lowest \(Y\)-exponent endpoints in the
same ratio, allowing continuation. Lemma 7.1(3) also
ensures that the edges of the two operators occur
together. At \((1,-1)\), Lemma 7.1(1) makes the
weights positive. The difference of the \(X\)- and
\(Y\)-coordinates increases as the support moves from
the highest to the lowest \(Y\)-exponent endpoint
of an edge with \(0<s<1\). No such endpoint can have
difference zero by Theorem 4.1(3). Thus some edge has
one endpoint of negative difference and one of positive
difference, with bracket zero and ratio \(2:b\).

This supplies the six conditions of Proposition 7.2.
The numerical conditions in Proposition 7.3 follow
by its proof from those six conditions. In detail, the
integrality of the normalized endpoints follows from
Theorem 2.11(2). The homogeneous bracket mate and endpoint
alignments are from Theorem 4.1 and Corollary 4.4.
For the last numerical condition in Proposition 7.3, its
proof applies Propositions 5.3 and 5.6, whose hypotheses
are the positive nondividing ratio, zero bracket and
endpoint signs supplied here. None of these steps
requires the pair to be minimal.

## The normalized crossings up to sum twenty

Write \(C_0=(u,v)\) for the endpoint of largest
\(Y\)-exponent divided by two, and \(C_1=(r_0,s_0)\)
for the other endpoint divided by two. The rectangle
gives \(u+v\le g+d=n/2\le20\). Proposition 7.3 gives
\[
 u<v,\quad u\ge3,\quad \gcd(u,v)>1,\qquad
 0\le s_0<r_0<u.
\]
The endpoints have equal \((\rho,\sigma)\)-weight,
where \((\rho,\sigma)\) is the primitive direction of
\((f_2-1,1-f_1)\) for an integer
\((f_1,f_2)=\tau (u,v)\) with \(f_1\ge2\) and \(0<\tau<1\).
If \(\gcd(f_1-1,f_2-1)=1\), set
\[
       C_2=C_1+\left(\frac{v-s_0}{\rho}-s_0\right)
                                       (-\sigma/\rho,1).
\]
Proposition 7.3(6) excludes \(C_2=(j-1/\rho,j)\)
for \(j\ge2\). The numerical check in the proof
of Corollary 7.4 rules out sums at most fifteen.
Sums seventeen and nineteen cannot have
\(\gcd(u,v)>1\). The remaining rows are

| Sum | \(C_0\) | \((f_1,f_2)\) | \((\rho,\sigma)\) | \(C_1\) |
|---:|---|---|---|---|
| 16 | \((4,12)\) | \((3,9)\) | \((4,-1)\) | \((1,0)\) |
| 16 | \((6,10)\) | \((3,5)\) | \((2,-1)\) | \((1,0)\) |
| 18 | \((6,12)\) | \((2,4)\) | \((3,-1)\) | \((2,0)\) |
| 20 | \((5,15)\) | \((2,6)\) | \((5,-1)\) | \((2,0)\) |
| 20 | \((8,12)\) | \((2,3)\) | \((2,-1)\) | \((2,0)\) |

For each row take integers \(1\le k<\gcd(u,v)\) for
\(\tau=k/\gcd(u,v)\) and \(1\le r_0<u\) for the equation
of equal weights, then test the stated conditions.
The two sum-sixteen rows give the exact face profiles
of the [order-thirty-two first-edge calculation](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_FACES.md#the-polynomials-on-the-two-edges).
The [successor exclusion](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md)
uses those profiles and the ratio \(2:b\). It applies
regardless of the original Bernstein degrees.

## The new faces and their successors

For each remaining row put \(C_1=(h,j)\) and
\(z=XY^\rho\). The common factor and homogeneous
bracket mate of weight \(\rho-1\) have the forms
\[
         R=X^hY^j L(z),\qquad
         \deg L=u-h,\qquad L(0)\ne0,\qquad F=XYG(z).
\]
The equation \(\{R,F\}=cR\) gives, with the sign
absorbed into \(c\),
\[
   (\rho h-j)zG'
       +\left(h-j-(\rho-1)z\frac{L'}{L}\right)G=c,
                         \qquad c\ne0.                  \tag{1}
\]
Since \(\deg L>0\), pole cancellation forces \(G\)
to vanish at every root of \(L\), so \(G\) has positive
degree. Comparison at infinity gives
\[
                \deg G=\frac{(\rho-1)\deg L-h+j}{\rho h-j}.
\]
In all three new rows this degree is one.
Thus \(L\) has one root \(u_*\ne0\).
In the face coordinate \(X^{1/\rho}Y\), the maximal
root multiplicity of the normalized common factor is
\(\mu=\deg L\), realized at a nonzero root. A zero
root, when present, has multiplicity \(j<\deg L\).

Apply [Guccione–Guccione–Valqui II, Proposition 4.3](https://arxiv.org/pdf/1205.6827v1)
in the fractional Weyl algebra, as in the
[order-thirty-two successor proof](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md#successors-after-a-root-shift).
The crossing provides the hypotheses: both faces
are nonmonomial, their positive weights have ratio
\(2:b\) with neither quotient integral, their bracket
is zero, and their highest \(Y\)-exponent endpoints
have negative diagonal difference. The commutator is
one and the direction has negative second coordinate.
In the normalized weights \((1,-s)\), the successor
has larger \(s\), positive weights \(2w,bw\) and highest
endpoints \(2C,bC\), where
\[
 C=(h-j/\rho+\mu/\rho,\mu),\qquad
 w=C_x-\mu s .
\]
The field of coefficients is now \(k[X^{\pm1/\rho}]\),
and the successor slopes lie below \(C_x/\mu\).

For a subsequent endpoint \(C=(C_x,\mu)\) on the
fractional lattice of denominator \(l\), write
\(\beta=l(\mu-C_x)>0\). A scalar leading bracket
would require
\[
   (b+2)(C_x-\mu s)=1-s,\qquad
   ls=\frac{(b+2)l C_x-l}{(b+2)\mu-1}.
\]
The cancelled gcd divides \(\beta\). All cases below
have \(\beta\le\mu\). The reduced denominator is
then at least \(((b+2)\mu-1)/\beta>b+1\)
and coprime to \(\mu\). Therefore neither \(2\mu\)
nor \(b\mu\) is divisible by it. Both face polynomials
are divisible by \(\eta\) in \(\eta=YX^s\), since all
deficits in their \(Y\)-exponents are multiples of
that denominator. Their leading bracket vanishes at
\(\eta=0\). This excludes a scalar leading bracket.

At a commuting successor write the common factor
\(X^w f(YX^s)\), where \(\deg f=\mu\). If \(f\) has only
one distinct root, that root is nonzero because this is
an edge. Its adjacent coefficients force \(ls\) integral.
Shifting at the root keeps both \(C\) and \(l\), and
forces a strictly larger slope for the next successor.
If \(f\) has multiple distinct roots, the
[successor coefficient equation](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md#the-coefficient-equation-and-the-lattice-table)
gives
\[
    W'-\frac{M}{\mu}\frac{f'}f W=c_0,\qquad
    M=\deg W\ge2,\qquad ls=l-\frac{\beta M}{\mu(M-1)},
    \qquad c_0\ne0.                                   \tag{2}
\]
It follows from the homogeneous bracket mate in
Proposition 4.3(11), applicable by positivity,
nondividing ratio, zero bracket and the negative endpoint
differences. The polynomial \(W\) has a simple zero at
each distinct root of \(f\). A root of multiplicity
\(m_0\) is impossible if \(M m_0=\mu\).
If \(D\) is the reduced denominator of \(ls\), the
least exponent \(j_0\) of \(f\) satisfies
\(D\mid\mu-j_0\), and \(f=\eta^{j_0}J(\eta^D)\)
with \(J(0)\ne0\). The lattice equation gives
\(M-1\mid\beta(\mu-j_0)\), so \(M\le\beta\mu+1\).

The lattice rows below are written \((M,D)\).
They are obtained by looping over these finitely many
integers \(M\), retaining slopes strictly between the
previous slope and \(C_x/\mu\), and requiring \(D\le\mu\).

### Sum eighteen

At \(C_0=(6,12)\) the first shift gives \(C=(10/3,4)\),
\(l=3,\ \beta=2\), with previous slope \(1/3\).
The lattice rows are \((2,1),(3,4),(4,3)\).
In the last two, respectively, there are four distinct
roots which exceed \(M\), or all roots are simple
with multiplier \(M/\mu=1\).
The first row, of slope \(2/3\), permits multiple roots
only with multiplicities \(3,1\). A double root would
violate \(M m_0=\mu\), and three distinct roots exceed
\(M=2\). A single-root edge is possible only at the
same slope \(2/3\), and would preserve \(C,l\) after
shifting, with no listed successor of larger slope.
Thus the only continuation would shift at the triple
root at \(2/3\), giving \(C'=(8/3,3),l=3,\beta=1\).
The slope is now greater than \(2/3\) and smaller
than \(8/9\), so no single-root edge is possible.
The lattice rows \((2,3),(3,2)\) respectively give
three roots exceeding \(M\), or all simple roots with
multiplier one. This excludes the corner.

### Sum twenty

For \(C_0=(5,15)\), the first shift gives \(C=(13/5,3)\),
\(l=5,\beta=2\). The lattice rows are \((2,3),(3,1)\).
The first requires three roots exceeding \(M=2\).
The second has multiplier one, and every cubic with
more than one distinct root has a simple root.
A single-root successor keeps \(C,l\), and only finitely
many slopes \(s\) with \(ls\) integral lie strictly
between \(1/5\) and \(C_x/\mu\). Thus continued shifts
cannot avoid a contradiction.

For \(C_0=(8,12)\) with \(C_1=(2,0)\), the new
corner is \(C=(5,6),l=2,\beta=2\).
Again no single-root slope is possible. The rows
\((2,3),(6,5)\) fail by the number of nonzero roots
and by multiplier one. The row \((3,2)\), of slope
\(3/4\), has at most three distinct roots and
multiplier \(1/2\). The nonzero roots occur in
opposite pairs. Thus the only possibilities are
\[
    (\eta^2-u)^3,\qquad \eta^2(\eta^2-u)^2,\qquad
                 \eta^4(\eta^2-u),\qquad u\ne0 .
\]
The middle polynomial has a double root, forbidden
by multiplier \(1/2\). For the first, a triple-root
shift gives \(C'=(11/4,3),l=4,\beta=1\).
The new slope is strictly between \(3/4\) and
\(11/12\), with no integral \(ls\). Its two lattice
rows \((2,3),(3,2)\) fail as in the last sum-eighteen
step. For the third polynomial, the maximal root
is zero with multiplicity four. Proposition 4.3
then uses the identity shift. Although the proposition
may be applied in a larger fractional ring, the
actual operators stay in \(W(2)\) and retain its
support lattice. Their new endpoint is
\(C'=(7/2,4),l=2,\beta=1\), with slope strictly
between \(3/4\) and \(7/8\). There is no single-root
slope, and the sole lattice row \((4,3)\) has
all roots simple with multiplier one.

Every new corner used above has \(C_x-\mu<0\).
Together with the preserved positive nondividing
ratio, multiple-term faces, zero bracket and scalar
commutator, this verifies the remaining hypotheses of
Proposition 4.3 at every step. The finite increases of
the slope at a fixed endpoint, and the explicit
exclusions at the new endpoints, exhaust these rows.

## Closing the degree induction

The absence of all crossings proves that Bernstein degree \(n\)
is impossible in a proper pair. To handle \(n+1\), perform
Bernstein mate descent as in the degree-thirty-three proof.
When one degree divides the other, subtract the scalar
multiple of the power given by the top commuting symbols.
The sum of degrees strictly decreases, and the smaller
degree never increases. At the end the smaller degree of
a proper pair is at least twice its gcd, which is at least
sixteen by Corollary 7.4. It cannot be at most \(n\).
It cannot be exactly \(n+1\) either, since none of \(35,37,39,41\)
has a proper divisor at least sixteen. This closes the
Bernstein step and updates the individual nonconstant
leading-degree bound to \(n+2\).

For the constant-leading assertions at orders \(n,n+1\),
use the three branches of the
[monic weighted argument](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#the-monic-bound-and-the-first-residual).
A dividing mate left after fixed-member descent would
be a smaller constant-leading member of order at most
twenty, already excluded. In the nondividing case the
integer \(m=1/r\) of that proof satisfies \(m\ge1\).
The \(m=1\) branch would have Bernstein degree equal to
the prescribed order after depression, now impossible.
The binomial coordinate change gives a smaller
constant-leading order at most twenty. The intermediate
coordinate change gives order \(i\le20\). If its leading
coefficient is constant the same lower-order bound applies.
Otherwise the sum of order and leading degree in that
presentation is at most the original order by the
transformed weight bound, whereas the updated invariant
makes it at least \(n+2\). This proves both assertions
through \(n+1\) and closes the induction at \(41\).

## Exact finite checks

The following Python snippet reproduces the corner and lattice tables
using integer and rational arithmetic.

~~~python
from fractions import Fraction
from math import gcd

crossings = []
for total in range(16, 21):
    for u in range(3, (total + 1) // 2):
        v = total - u
        divisor = gcd(u, v)
        for k in range(1, divisor):
            f1, f2 = k * u // divisor, k * v // divisor
            if f1 < 2:
                continue
            cut_divisor = gcd(f1 - 1, f2 - 1)
            rho, sigma = (f2 - 1) // cut_divisor, (1 - f1) // cut_divisor
            for r0 in range(1, u):
                for s0 in range(r0):
                    if rho * r0 + sigma * s0 != rho * u + sigma * v:
                        continue
                    offset = Fraction(v - s0, rho) - s0
                    next_x = r0 - offset * Fraction(sigma, rho)
                    next_y = s0 + offset
                    if cut_divisor == 1 and next_y >= 2 and next_x == next_y - Fraction(1, rho):
                        continue
                    crossings.append((total, (u, v), (f1, f2), (rho, sigma), (r0, s0)))

assert crossings == [
    (16, (4, 12), (3, 9), (4, -1), (1, 0)),
    (16, (6, 10), (3, 5), (2, -1), (1, 0)),
    (18, (6, 12), (2, 4), (3, -1), (2, 0)),
    (20, (5, 15), (2, 6), (5, -1), (2, 0)),
    (20, (8, 12), (2, 3), (2, -1), (2, 0)),
]

def lattice_rows(ell, mu, beta, previous):
    upper = 1 - Fraction(beta, ell * mu)
    rows = []
    for degree in range(2, beta * mu + 2):
        slope = 1 - Fraction(beta * degree, ell * mu * (degree - 1))
        denominator = (ell * slope).denominator
        if previous < slope < upper and denominator <= mu:
            rows.append((degree, denominator))
    return rows

assert lattice_rows(3, 4, 2, Fraction(1, 3)) == [(2, 1), (3, 4), (4, 3)]
assert lattice_rows(3, 3, 1, Fraction(2, 3)) == [(2, 3), (3, 2)]
assert lattice_rows(5, 3, 2, Fraction(1, 5)) == [(2, 3), (3, 1)]
assert lattice_rows(2, 6, 2, Fraction(1, 2)) == [(2, 3), (3, 2), (6, 5)]
assert lattice_rows(4, 3, 1, Fraction(3, 4)) == [(2, 3), (3, 2)]
assert lattice_rows(2, 4, 1, Fraction(3, 4)) == [(4, 3)]
~~~
