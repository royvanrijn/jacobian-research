# Bernstein degree through thirty-three for polynomial Weyl pairs

**Theorem PolynomialWeylBernsteinGcdSixteenQuotients.** Let \(P,Q\)
be polynomial Weyl operators over a characteristic-zero field,
with nonzero scalar commutator, generating a proper subalgebra.
If their Bernstein-degree gcd is sixteen, the unordered pair
of Bernstein degrees is
\[
                     \{16(2j+3),16(3j+4)\},\qquad j\ge0.
\]
In the standard minimal presentation the normalized corner is
\((4,12)\). The first face has the single-root profile
\(L(z)\) proportional to \((z-u)^3\) in the notation of
[the first-face calculation](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_FACES.md).
A necessary next scalar bracket, possibly after one intervening root
shift that preserves the corner, has normalized endpoint \((7/4,3)\).

**Corollary PolynomialWeylBernsteinDegreeThroughThirtyThree.**
If a member of a polynomial Weyl pair with nonzero scalar commutator
has Bernstein degree at most thirty-three, the pair generates.
In particular, in a nongenerating pair every member with nonconstant
leading differential coefficient satisfies
\[
                        \operatorname{ord}P+
                        \deg\operatorname{lc}P\ge34.
\]

## Standard minimal presentation and the commuting successors

Scale the commutator to one and extend the constants to an algebraic
closure. By [Guccione–Guccione–Valqui, Corollary 7.4](sources/guccione-guccione-valqui-1111.6100v3.pdf)
a proper pair of Bernstein-degree gcd sixteen is minimal. Their
Corollary 6.11 gives a standard minimal presentation preserving both
Bernstein degrees. Write them as \(16a,16b\), with \(a,b\) coprime.
If one quotient were one, repeatedly subtracting leading powers of
the degree-sixteen member by the top-symbol rule below would force
the mate's degree below sixteen or give a positive nonmultiple of
sixteen, either contrary to Corollary 7.4. Thus \(a,b\ge2\).

In the subrectangular presentation, divide the Bernstein vertex
coordinates by \(a\). By the commuting leading-symbol relation these
are integers and agree with the vertex for \(Q\) divided by \(b\).
Their sum is sixteen. The standard minimal first-face calculation
applies with \(a,b\) in place of \(2,b\): the two endpoints
in Propositions 7.2–7.3 are normalized by \(a\), and the support
rectangle bounds their coordinates by the normalized Bernstein
vertex. The exclusion of a normalized sum at most fifteen forces
the upper endpoint to equal that vertex. The only corners and slopes
are \((4,12)\) with slope \(1/4\), and \((6,10)\) with slope
\(1/2\). The common-power profiles are unchanged, since dividing the
homogeneous bracket identity for \(R^a\) by \(aR^{a-1}\) again
gives \(\{R,F\}=cR\).

In the [order-thirty-two successor proof](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md),
replace the positive weight ratio \(2:b\) by \(a:b\).
The root shifts keep that ratio, and the commuting-edge differential
equation (2), derived from the homogeneous bracket mate, depends
only on the common factor and its normalized weight. Thus the
lattice table, root obstructions, and possible single-root
intermediate edges remain the same. In particular, the results
for commuting successors are:

- At \((\ell,\mu)=(2,4),(2,5)\) no successor can commute.
- At \((4,3)\) a commuting successor can only be the single-root
  shift of slope \(1/2\) that keeps the same lattice and endpoint.
- At \((4,2)\) at most one such shift can precede the commuting
  binomial edge of slope \(5/8\). Its simple-root shift gives
  the endpoint \((7/8,1)\) on the lattice of denominator eight.
  There can then be at most one more commuting shift, at slope
  \(3/4\), preserving that endpoint and lattice.

Each alternative requires a scalar leading bracket at some successor,
since the successor proposition at a commuting edge supplies a strictly
larger slope with the prescribed positive weights and nonmonomial faces.

## The necessary scalar successor

Write \(N=a+b\ge5\). At the first shifted endpoint
\(C=(1+\mu/\ell,\mu)\), before or after a single-root shift
that preserves this endpoint and lattice, a scalar leading
bracket requires
\[
  N(C_x-\mu s)=1-s,\qquad
  \ell s=\frac{N(\ell+\mu)-\ell}{N\mu-1}.
\]
Let \(D\) be the reduced denominator of \(\ell s\).
The cancelled gcd divides \(\beta=\mu(\ell-1)-\ell\), and
\(\gcd(D,\mu)=1\). In the variable \(\eta=YX^s\), the two
face polynomials have \(Y\)-degrees \(a\mu,b\mu\), respectively,
and all their exponents are congruent to the associated degree
modulo \(D\). If both polynomials vanish at \(\eta=0\), their
Poisson bracket vanishes there: for each symbol the \(X\)
derivative at fixed \(Y=0\) is zero. A scalar nonzero bracket
therefore requires \(D\mid a\mu\) or \(D\mid b\mu\), hence
\(D\mid a\) or \(D\mid b\).

For \((\ell,\mu)=(2,4),(2,5),(4,2)\) we have \(\beta\le\mu\).
Consequently \(D\ge(N\mu-1)/\beta>N-1\), greater than both
\(a,b\). At the later endpoint \((7/8,1)\) the scalar equation
would be \(8s=7-1/(N-1)\); its reduced denominator \(N-1\)
is again greater than both face degrees at that endpoint,
violating edge membership.

Thus only \((\ell,\mu)=(4,3)\) can give a scalar successor.
There \(\beta=5\), so the cancelled gcd must be five, and
\(N\equiv2\bmod5\). Write \(N=5j+7\), with \(j\ge0\).
Then \(D=3j+4>N/2\), so divisibility by \(D\) forces one of
\(a,b\) to equal \(D\). The other is \(N-D=2j+3\).
This is the asserted necessary degree pair. The corner and
profile follow from \((\ell,\mu)=(4,3)\). If \(a=3j+4\), the
scalar slope is \(s=(7j+9)/(4(3j+4))\). The polynomials in
\(\eta\) then have the respective forms \(A(t),\eta B(t)\)
with \(t=\eta^a\), \(\deg A=3,\deg B=2\) and
\(A(0)B(0)\ne0\). The scalar leading bracket also requires
\[
                   b t A'B-A(B+a t B')=c,\qquad c\ne0 .
\]
This follows by taking the bracket of the symbols with prefactors
\(X^{aw},X^{bw}\) and dividing by its nonzero constant factor.

The scalar-face equation has a solution already for \(j=0\).
For a root \(v\) of \(54v^2-36v+5=0\), one can take
\[
 \begin{aligned}
 A(t)&=1+\frac52t+\frac{9v+5}{5}t^2+\frac{3v(9v+5)}{20}t^3,\\
 B(t)&=1+t+v t^2 .
 \end{aligned}
\]
Substitution modulo the quadratic equation gives \(c=-1\).
This realizes the necessary scalar bracket of the two symbols.
The corresponding Weyl operator equations remain open.

## Small Bernstein degree and a nonconstant leading coefficient

Suppose a proper pair has one Bernstein degree at most thirty-three.
When one degree divides the other, the highest Bernstein symbols
Poisson commute because their potential bracket weight is positive.
The homogeneous Euler identity then makes the higher symbol a
scalar multiple of the appropriate power of the lower. Subtract
that scalar times the operator power from the higher member.
The pair is still proper and the sum of the positive degrees
strictly decreases; a scalar member would contradict its
commutator. The smaller degree never increases. At the end the
degree quotient in either direction is nonintegral, so the
smaller degree is at least twice the degree gcd, which is at
least sixteen. The only possibility with smaller degree at most
thirty-three is smaller degree thirty-two and gcd sixteen.
Its smaller quotient would be two, which the theorem excludes.

For a differential nondividing pair with a nonconstant common
leading factor, the
[polynomial-shear calculation](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#nonconstant-leading-factors-and-bernstein-degree)
makes the Bernstein degree equal to order plus leading-coefficient
degree for each member. Hence each such sum is at least thirty-four.
For any member \(P\) with nonconstant leading coefficient, keep
it fixed and reduce its partner by its powers in differential
order. If the reduced pair is nondividing, apply the preceding
case. Otherwise the partner has order \(j<i=\operatorname{ord}P\)
dividing \(i\), and the leading commutator identity gives a
nonconstant leading coefficient for that partner and
\[
 i+\deg\operatorname{lc}P
             =(i/j)(j+\deg\operatorname{lc}Q).
\]
Induction on \(i\) gives the desired bound for every member.
