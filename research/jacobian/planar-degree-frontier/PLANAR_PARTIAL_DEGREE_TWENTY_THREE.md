# Common planar partial degrees through twenty-three

**Theorem PlanarCommonPartialDegreeTwentyThree.**
Let \(P,Q\) be a Keller pair over a characteristic-zero field.
If its two partial degrees in the same source variable are at most
twenty-three, then it is invertible. With a common bound of
twenty-four, it also suffices that every component of partial
degree twenty-four has a leading coefficient of degree less
than \(84\) in the other variable.

The new reduced pair is \((14,21)\). We use a second polynomial
Jacobian identity to exclude one of its three possible primitive
corners and determine the right-hand Newton constraints in the
other two. A direct left-hand normalization puts both survivors
into the exact finite systems of Niedbala Giraudin, without needing
the minimal-pair normalization in those papers. The systems are the
ones in Sections 5 of the
[degree-126 preprint](sources/niedbala-giraudin-2611.00710v1.pdf)
and the [degree-147 preprint](sources/niedbala-giraudin-dz147-2611.00897v1.pdf).

## Reduction to three primitive corners

Work over the algebraic closure, with inverse descent by faithful
flatness. As in the [degree-twenty proof](PLANAR_PARTIAL_DEGREE_TWENTY.md),
subtract target powers until the positive partial degrees \(au,bu\)
do not divide one another, where \(2\le a<b\) are coprime. A zero
partial degree is triangular. For \(u\le5\), the
[small-gcd result](SMALL_PARTIAL_GCD_STANDARD_CORNERS.md)
excludes \(u<4\) and gives normalized larger total degree \(16b\)
or \(25b\), respectively, within Moh's inclusive degree-100 bound
under the cap twenty-three. For \(u=6\), only \((12,18)\) fits and
it is already excluded. Since \(b\ge3\), under the cap twenty-three
the sole remaining pair is \((14,21)\).

Let \(h\) be its common monic leading factor in the other variable
and \(d=\deg h\). The
[simultaneous triangular normalization](PARTIAL_DEGREE_FOURTEEN_AND_WEYL_TRANSFER.md#simultaneous-triangular-normalization)
makes the total degrees \(2(7+d),3(7+d)\) without altering the
two leading rows. The total-degree gcd floor of sixteen gives
\(d>7\). The [axial root lemma](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#one-translation-for-the-two-leading-coefficients)
provides a root of multiplicity greater than seven. Before
translating it to zero, we may make the two total highest forms
monomials. Indeed, their factors include the high power of the
other variable, since the total degrees equal the sums of the
two corner coordinates. The
[positive-weight slice classification](../../dixmier/dixmier-formal-lifting/WEIGHTED_PLANE_HAMILTONIAN_SLICES.md#positive-weights-and-compatibility-of-exposed-edges)
then gives either a monomial or, in the frame with partial degrees
in \(y\), the forms proportional to
\(x^{2d}(y-\lambda x)^{14},x^{3d}(y-\lambda x)^{21}\).
The common \(\lambda\) follows from the zero highest bracket
and homogeneity. The shear removing \(\lambda\) preserves both
partial degrees and both leading rows.

Translate the high root to zero and swap the axes, making the
pair standard in Guccione–Guccione–Valqui's Definition 4.3.
Use \(x,y\) for these swapped axes. Write \(A_0=(i,j)\) for
the primitive corner and \(p_0/q_0\) for the reduced proportionality
of its bracket mate. As in the small-gcd result, \(4\le i\le7\)
and \(2\le f_1<i\), so \(q_0>1\) divides \(i\).
If \(i<7\), there is a later corner at the end of the right
face, and GGV Theorem 7.6(7)–(8) gives \(q_0\mid7\), impossible.
Thus \(i=7\), which gives \(q_0=7\). The primitive corner is
the bottom of the right face. GGV Corollary 7.10 then collapses
that face to one endpoint. Hence \(d=j\) and \(h\) has a single
root of multiplicity \(j\). GGV Propositions 6.4 and 6.7 give
\(7\mid j\) and \(j\le42\). Moh's bound excludes \(j\le21\).
[GGV, primary text](https://arxiv.org/pdf/1401.1784v3)

For \(j=28,35,42\), apply Proposition 6.4 to every
\(2\le f_1\le6\), with \(f_2=jf_1/7\).
For \(e=\gcd(f_1-1,f_2-1)\) the direction is
\(((f_2-1)/e,(1-f_1)/e)\).
This gives the following possibilities with \(\rho\le7\).
A possible other endpoint \((i',j')\) must have
\(0\le j'<i'<7\) and the same weight as \(A_0\).

| \(j\) | \((f_1,f_2)\) | Direction | Other endpoint |
|---:|---:|---:|---:|
| 28 | \((2,8)\) | \((7,-1)\) | \((3,0)\) |
| 28 | \((4,16)\) | \((5,-1)\) | none |
| 35 | \((3,15)\) | \((7,-1)\) | \((2,0)\) |
| 35 | \((5,25)\) | \((6,-1)\) | none |
| 42 | \((6,36)\) | \((7,-1)\) | \((1,0)\) |

In each actual row put \(k=7-j/7\), so \(k=3,2,1\),
respectively. By homogeneity and the vanishing highest bracket,
the two initial forms at \((7,-1)\) are scalar multiples of
the square and cube of \(x^k H(xy^7)\), where \(H\) is monic
of degree \(7-k\) and \(H(0)\ne0\).

## Second Jacobian identity and the right-hand constraints

Temporarily use the other order of the axes \(x,y\), so that
the partial degrees in \(y\) are \(14,21\).
The two initial forms at \((-1,7)\) are proportional to
\(R^2,R^3\), with \(R=x^{-7k}r(z)\), \(z=x^7y\),
and \(r(z)=z^k H(z)\). Scale to make \(h^2,h^3\) the two
\(y\)-leading coefficients, with \(h=x^{7(7-k)}\).
Write \(J(P,Q)=\kappa\in K^\times\).

The same construction as in the degree-twenty proof now
uses \(T=P^{1/14}\), with leading \(y\)-term \(h^{1/7}y\).
In the Laurent expansion \(Q=\sum c_jT^j\), all
\(c_j\) with \(j>-13\) are constant, since
\((\partial Q/\partial x)_T=-\kappa/P_y\) has degree \(-13\).
Remove \(c_0+c_{14}P\) and put \(F=T^7(P+c_7)\).
If \(Q-F\) has nonnegative degree \(j\), then \(7\nmid j\).
The polynomial \(Q^2-F^2\) has \(y\)-degree \(21+j\)
and leading coefficient \(2c_jh^{(21+j)/7}\). Its
\((-1,7)\)-weight is at least \(k(21+j)\). Whenever that
weight is greater than \(7k+6\), the identity
\(J(P,Q^2-F^2)=2\kappa Q\) forces its top to commute with \(R\).
For a form \(x^{-w}U(z)\) this means
\(7k rU'=w r'U\). At zero, of multiplicity \(k\) in \(r\),
this requires \(w/7\) to be integral. The degree of \(U\)
is \(w/k\). Consequently \(w\ge k(21+j)\) and the degree
bound force equality, contradicting \(7\nmid j\).
In particular all nonnegative terms have been removed by \(F\).

Set \(F_1=F+c_{-7}T^{-7}\) and
\[
  W=Q^2-P(P+c_7)^2-2c_{-7}(P+c_7),\qquad
  J(P,W)=2\kappa Q .
\]
The degree \(j\) of \(Q-F_1\) lies in \([-13,-1]\)
and differs from \(-7\). The degree of \(W\) in \(y\)
is \(21+j\), since \(F_1^2\) differs from the subtracted
polynomial by \(c_{-7}^2 P^{-1}\).
If \(k(21+j)>7k+6\), then \(j>-13\), and the preceding
argument gives the same contradiction. Therefore
\(\deg_y W\le7+6/k\). A weight of \(W\) below \(7k+6\)
would fail to give the required weight in \(J(P,W)\).
A greater weight would force a commuting top whose
degree in \(y\) exceeds \(7+6/k\).
Its weight is exactly \(7k+6\).

Write the top of \(W\) as \(x^{-(7k+6)}w(z)\).
Up to a nonzero scalar the second Jacobian identity gives
\[
  7k r w'-(7k+6)r'w=C r^2,\qquad C\ne0 .
\]
At a root of multiplicity \(l\) in \(r\), an order \(v\le l\)
in \(w\) yields a nonzero leading multiplier and order below
\(2l\) on the left. An order \(v>l+1\) yields order above \(2l\).
Thus \(w\) has order \(l+1\). In particular \(w/(zr)\) is
a polynomial of degree at most \(6/k-1\) with a zero at
every root of \(H\). Also \(w/r\) has degree \(6/k\). Indeed,
it has at least two distinct zeros, and in the equation for
\(w/r\) the top term at infinity can cancel only at that degree.
Both facts will be used below.

Return to the standard axes where \(x\) has degrees \(14,21\).
Let \(m\) be the maximal multiplicity of a root of \(H\).
The preceding bound gives \(m=4\) for \(k=3\),
\(3\le m\le5\) for \(k=2\), and \(2\le m\le6\) for \(k=1\).
In \(K[x^{1/7},x^{-1/7}][y]\) let
\(\eta=x^{1/7}y\) at slope \(s_0=1/7\) for weight \((1,-s)\).
Shift \(y\) by \(\lambda x^{-1/7}\), where \(\lambda^7\)
is a root of multiplicity \(m\). The two forms and the form
of \(W\) now have smallest-\(y\) terms with exponents
\[
  x^{2k}\eta^{2m},\qquad x^{3k}\eta^{3m},
  \qquad x^{k+6/7}\eta^{m+1},
\]
respectively. For larger slopes, while the first two terms
remain highest, their weights are \(a=2(k-m\Delta), b=3a/2\),
where \(\Delta=s-1/7\). The third term bounds the weight
of \(W\) below by \(\omega_0=k+6/7-(m+1)\Delta\).
Above slope \(s_0\), terms of \(y\)-degree greater than \(m+1\)
cannot match the weight of that third monomial.

As in the degree-twenty proof, a first tie must occur before
the anchor weights vanish at \(1/7+k/m<1\), since the
scalar bracket would otherwise have negative weight there.
At a commuting tie the two top forms are powers of a
common polynomial \(f(x^s y)\) of degree \(m\), times their
respective powers of \(x\). The top of \(W\) has degree at
most \(m+1\) in \(x^s y\). The second Jacobian identity
shows its weight is \(\omega_0\): a larger weight would imply
commutation and require a weight of at most \(a(m+1)/(2m)\),
strictly less than \(\omega_0\) since
\[
 \omega_0-\frac{a(m+1)}{2m}=\frac67-\frac{k}{m}>0 .
\]
At that weight its coefficient polynomial \(U\) satisfies
\(U'-(2\omega_0/a)(f'/f)U=C_1 f\), \(C_1\ne0\).
The degree bound and the same root-order argument force
\(f\) to have one root and \(U\) to vanish there to order
\(m+1\). This root is nonzero, since a tie occurred.
Adjacent terms in the binomial power show \(7s\) integral.
The next translation preserves the three anchor points and
raises the next slope, preserving the degree bound for \(W\).
There are finitely many such integral \(7s\) in \((0,7)\).

At a scalar tie we must have \(a+b=1-s\), hence
\[
                        7s=\frac{35k+5m-7}{5m-1}.
\]
Its reduced denominator \(D\) divides every difference
of \(y\)-exponents at equal weight. Therefore \(D\le3m\)
at a tie. If \(D\) divides neither \(2m\) nor \(3m\),
both top coefficient polynomials in \(x^sy\) have zero
constant term and their bracket vanishes at \(y=0\). Here one can
compute the bracket in \((x,x^s y)\) after adjoining the
denominator of \(s\). The top scalar is a nonzero multiple
of \(a p q'-b p' q\), where \(p,q\) are the two coefficient
polynomials. The denominator tests leave the scalar ties in the table.

| \(k\) | Possible \(m\) | Corresponding \(D\) | Surviving scalar tie |
|---:|---|---|---|
| 3 | 4 | 19 | none |
| 2 | 3, 4, 5 | 7, 19, 3 | \(m=5,\ s=11/21\) |
| 1 | 2, 3, 4, 5, 6 | 9, 14, 19, 24, 1 | \(m=6,\ s=2/7\) |

There must be a scalar tie, since successive commuting ties
terminate. Thus \(k=3\) is excluded, and in the other cases
\(H(t)=(t-c)^{7-k}\) with \(c\ne0\). For \(k=2\), all
intermediate shifts after the first one have slope \(2/7\)
or \(3/7\). At the final slope \(11/21\), put \(X=x^{1/7}\)
and use the integer weight \(3i-11j\) on \(X^i y^j\).
The weights of the shifted \(P,Q\) are \(4,6\).
The highest \(y\)-exponents for \(P\) are congruent to
\(1\bmod3\) and those for \(Q\) to \(0\bmod3\).
Since the bracket is scalar and nonzero, their smallest
highest exponents are \(1,0\). The corresponding terms are \(X^5 y,X^2\).
For \(k=1\) there is no intermediate shift, and the
weights for \(i-2j\) are \(2,3\).

## Polynomial admission from the other Newton side

We give this step for both \(k=1,2\). In the axes with
\(y\)-degree \(14,21\) used in the second identity, the
form \(x^{-6}(w/r)(x^7y)\) is, up to scalar, a bracket
mate for the \((-1,7)\)-initial form of \(P\).
The difference of any two mates would commute with \(R\).
In \(K[x^{\pm1},y]\) that would require a weight divisible
by seven at the root zero of \(r\), so the mate is unique.
Thus it is the homogeneous mate in GGV Theorem 2.6
after scaling. Its start at the high-\(y\) endpoint
is \(((7-k)M,M)\), where \(M=6/k\), since \(w/r\)
has degree \(M\). The proportionality to the start of
\(R\) is \(M/7\), of reduced denominator seven.

In the standard axes, \(x\)-degree is bounded by \(14,21\),
and the high-\(x\) anchors \(x^{14}y^{14(7-k)}\) and
\(x^{21}y^{21(7-k)}\) are the single highest monomials
for ordinary total degree. Consider the other side in
directions \((-s,1)\) of these axes, starting at \(s=-1\).
The GGV Corollary 7.4 hypotheses in the flipped axes
are satisfied: the total-degree and second-axis degree
ratios are \(2:3\), at \((-1,7)\) the starts have normalized
coordinate \((7(7-k),7)\) with first coordinate larger
than second, and the mate denominator is seven as just
proved. At a first tie on this other side for \(s<7-k\),
all weights between the first face and the tie give
positive values because the anchors remain in the support.
Corollary 7.4 therefore makes the tied initial form of
\(P\) a fourteenth power in the ring allowing Laurent
powers of \(y\). The highest forms commute at the tie
because the sum of their weights exceeds \(1-s\).
They have weights in ratio \(2:3\), so a tie for \(Q\)
also gives one for \(P\).

Since the fourteenth root has degree one in \(x\), its
form is proportional to \(xy^{7-k}+\gamma y^{7-k-s}\).
Thus \(s\) is an integer and a shift of \(x\) by a
multiple of \(y^{-s}\) removes the tie simultaneously
for \(P,Q\). Each such step strictly increases the next
slope. These shifts have \(0\le s<7-k\). They preserve
the original \((7,-1)\)-forms, both total degrees, and
the two \(x\)-degrees with their leading anchors. In
particular Corollary 7.4 still applies after the shifts
in the same Laurent ring with the same bracket mate.
After finitely many steps at \(s=0,\ldots,6-k\) the anchors
are highest up to \(s=7-k\). At that weight they have
weight zero, so all supported monomials satisfy
\(j\le(7-k)i\) in \(x^i y^j\).

Consequently, putting \(v=y^{-1}\), the Laurent shifts
give
\[
 \widehat P(z,v)=P\left(zv^{7-k}+\sum_{s=0}^{6-k}g_s v^s,
                        v^{-1}\right)
\]
and the identical expression for \(\widehat Q\), both
polynomial. Their Jacobian is \(-\kappa v^{5-k}\).
In weights \((k,1)\) the top forms are proportional
to the square and cube of
\(z^k(z-c v^k)^{7-k}\), and the respective weights are
\(14k,21k\). The shifts at \(s=0,1\) are allowed here.
We check their effect explicitly when combining the two
sides.

For \(k=2\), substitute \(x=X^7\),
\(y=Y=\eta+\lambda X^{-1}+\lambda_2 X^{-2}+\lambda_3 X^{-3}\),
with \(\lambda^7=c\), from the right-hand constraints.
With weight \(3i-11j\) on \(X^i\eta^j\), \(v=Y^{-1}\)
starts with \(X/\lambda\). In these substitutions
\(z=X^7Y^5-\sum g_s Y^{5-s}\). The terms in that sum
have weight at most \(-3\).
Put
\[
 \zeta=z-cv^2-7\lambda^5\lambda_2 v
              -7\lambda^3(\lambda\lambda_3+2\lambda_2^2).
\]
Its leading term in the right chart is \(7\lambda^4\eta X^3\).
Indeed the terms of weights \(3,0\) in
\(X^7Y^5-cY^{-2}\) are
\(7\lambda^4\lambda_2 X\) and
\(7\lambda^4\lambda_3+7\lambda^3\lambda_2^2\).
The two subtractions cancel them, leaving the displayed
term of weight \(-2\). Leading monomials of
\(\zeta^a v^b\) are distinct, with weight \(3b-2a\).
The right-hand bounds therefore give
\[
 \begin{array}{ll}
 \widehat P:\ 2a+b\le28,\quad 3b-2a\le4,\\
 \widehat Q:\ 2a+b\le42,\quad 3b-2a\le6.
 \end{array}
\]
The coefficient of \(\zeta v^2\) in \(\widehat P\)
and that of \(v^2\) in \(\widehat Q\) are nonzero
by the lowest right-edge terms above. The top of
\(\widehat P\) for weights \((2,1)\) is proportional to
\(\zeta^{10}(\zeta+c v^2)^4\), and its bracket with
\(\widehat Q\) is still \(-\kappa v^3\).
These are exactly (S1)–(S3) of the degree-126
preprint, excluded there by Lemma 5.1 and Proposition 5.2.

For \(k=1\), the right chart uses
\(Y=\eta+\lambda X^{-1}\) with \(\lambda^7=c\).
For \(i-2j\) on \(X^i\eta^j\), \(v\) starts with
\(X/\lambda\), and \(\zeta=z-cv\) starts with
\(7\lambda^5\eta X^2\). This follows by expanding
\(X^7Y^6-cY^{-1}\) to first order in
\(\eta X/\lambda\), which has weight \(-1\).
Every \(g_sY^{6-s}\) has weight at most \(-1\).
The leading monomials of \(\zeta^a v^b\) are
distinct and have weight \(b\). Thus in \((\zeta,v)\)
the two polynomials have \(v\)-degree at most \(2,3\),
ordinary degrees \(14,21\), and highest forms
proportional to \([\zeta^6(\zeta+cv)]^2\) and
\([\zeta^6(\zeta+cv)]^3\). Their bracket is
\(-\kappa v^4\). Rescale \(v\), rescale the outputs and
exchange them to obtain exactly (S1)–(S2) of the
degree-147 preprint. Its Lemma 5.1 classifies the
normalized leading pairs, and Proposition 5.2
excludes every specialization. No completeness
claim about the total-degree census enters either
finite-system application.

## Partial degree twenty-four and the transfers

Under the cap twenty-four, the only new reduced pairs
outside the cap twenty-three that are not settled by
the small-gcd theorem and Moh are \((16,24)\) and
\((18,24)\). For the gcd-four possibility the larger
normalized total degree is at most \(16\cdot6\), within
Moh. For \((16,24)\), normalization gives \(3(8+d)>100\),
so \(d\ge26\). GGV Corollary 7.9 applied in the
[standardized axial frame](SMALL_PARTIAL_GCD_STANDARD_CORNERS.md#planar-pairs)
gives \(\gcd(8,d)>2\), since \(d>8\). Thus \(4\mid d\),
so \(d\ge28\), and the degree-twenty-four member
has leading coefficient of degree \(3d\ge84\).
For \((18,24)\), the larger total degree is \(4(6+d)>100\).
The gcd-six result gives either \(d=30\) or a multiple
of six in this range, including \(d=24\) as the smallest
possibility. Its leading coefficient has degree \(4d\ge96\).
Under target power subtraction, if a member of partial
degree twenty-four remains, one original member of that
degree has remained fixed. Once a degree drops below
twenty-four it cannot return. This proves the
coefficient-sensitive statement.

**Corollary SexticQuarticGenericFibreChart.**
A characteristic-zero sextic Keller triple \((P,Q,R)\) is
invertible if the generic fibre \(R=t\) over \(K(t)\) has
a polynomial parametrization \(Z(u,v)\) of geometric rank
two everywhere, with \(\deg_v Z_i\le4\) and the following
condition. For the coefficient vector \(A(u)\) of \(v^4\),
each nonzero value \(P_6(A),Q_6(A)\) has degree below \(84\)
in \(u\), where \(P_6,Q_6\) are the ordinary degree-six parts.
In particular \(\deg_u A\le13\) suffices.

The [polynomial generic-fibre argument](PARTIAL_DEGREE_SLICES_AND_HESSIAN_TRANSFER.md#a-three-variable-generic-fibre-criterion)
gives a planar Keller pair of common \(v\)-degree at
most twenty-four. Its only possible degree-twenty-four
leading coefficients are the stated evaluations. The
theorem applies and inversion descends. The existence
of such a parametrization is an additional hypothesis.

**Corollary HessianDirectionalDegreeTwentyThree.**
Assume Tang's complex structural alternative. For a
four-variable complex polynomial \(H\) of nonzero constant
Hessian determinant, its gradient is invertible if three
independent constant directions \(w\) satisfy
\(D_w^{24}H\in\mathbb C\). It also suffices that
\(D_w^{25}H\in\mathbb C\) and \(\deg D_w^{24}H\le84\)
in three independent directions.

As in the [degree-twenty proof](PLANAR_PARTIAL_DEGREE_TWENTY.md#consequences),
one direction projects nontrivially to Tang's base plane.
The two stated conditions there give common planar partial
degree at most twenty-three, or partial degree at most
twenty-four with each possible leading coefficient of
degree at most eighty-three. Apply the theorem.

**Corollary WeylCommonDifferentialOrderTwentyThree.**
A characteristic-zero first Weyl pair with nonzero scalar
commutator generates if its two differential orders are
at most twenty-three. In a proper pair with a member of
order fourteen, fixed-member reduction gives a mate whose
order is an odd multiple of seven at least thirty-five.

As in the [order-twenty proof](SMALL_PARTIAL_GCD_STANDARD_CORNERS.md#transfer-to-polynomial-weyl-pairs),
a positive nondividing Weyl order pair transfers to a
noninvertible planar pair with the same partial degrees.
Target power subtraction preserves generation, and order
zero forces generation. This gives the new cap. For a
fixed order-fourteen member, use the
[existing lower bound of twelve on both proper orders](../../dixmier/dixmier-formal-lifting/POLYNOMIAL_WEYL_PRIMITIVE_LEADING_FACTOR_OBSTRUCTION.md#consequences-for-centered-rectangles-and-low-orders).
Every proper divisor of fourteen is below that bound,
so the reduced mate and fourteen do not divide one another.
The planar gcd bound forces their gcd to be seven.
The mate then has odd quotient by seven, and the new cap
excludes quotient three.


## Exact finite-system inputs

For the degree-126 system, the [deposited source archive](sources/niedbala-giraudin-2611.00710v1.zip)
and the [internal replay receipt](SHAPE_EXCLUSION_REPLAY.json)
record the initial coefficient normalization and both ideal checks
in Proposition 5.2.

For the degree-147 system, the
[deposited source archive](sources/niedbala-giraudin-dz147-2611.00897v1.zip)
has SHA-256 `2501b0c2354d469ae4e1070a9a53a1173150a46370cdea9d0c385aea859ca97a`.
Its `anc/verify_147.py` and `anc/leading_pair_147.txt` have respective
hashes `a50c872e7b16f71df4a0cbe8fe88019416ae9fbe7f1c2570102c7b8885e859ee`
and `1ccae410930c502a4ba67ee9a031e191847a1a96967bb825154167dcbd99b272`.
On Ubuntu 24.04, Python 3.12.3, SymPy 1.14.0, python-flint 0.9.0,
cypari2 2.1.4, PARI 2.15.4 and Singular 4.3.2, the prefix through
check 7g passed. It verified the representative and its field of
degree six, the count of six rotation classes of triangulations,
the exact linear eliminations, and the unit ideal over the number
field. Lemma 5.1 supplies the mathematical upper bound on the
number of leading pairs via the weighted-plane-tree bijection.
The modular counts serve only as controls.

The full script under this cypari2 build stopped at the
supplementary check 7h while assigning a matrix entry. The
required prefix was replayed with the original code and a clean
exit before that extra check. For a byte-preserving replay,
run the following from the extracted `anc` directory with the
dependencies on the Python path and `Singular` on the executable path.

```python
import os
from pathlib import Path

script = Path("verify_147.py").resolve()
source = script.read_text()
stop = "mons = sorted(set(m for cc in res['c1'] for m in cc))"
assert source.count(stop) == 1
scope = {"__file__": str(script), "__name__": "__main__"}
exec(compile(source.split(stop)[0], str(script), "exec"), scope)
assert scope["NPASS"][0] == 23
print("REQUIRED CHECKS THROUGH 7g PASSED", flush=True)
os._exit(0)
```
