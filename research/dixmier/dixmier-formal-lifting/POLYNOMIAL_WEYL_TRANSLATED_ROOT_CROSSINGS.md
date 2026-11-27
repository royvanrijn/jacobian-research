# Newton crossings at roots of Weyl leading factors

## The finite list for a fixed differential-order gcd

**Theorem PolynomialWeylTranslatedRootCrossing.** Let \(P,Q\) be
polynomial differential operators over a characteristic-zero field with
nonzero scalar commutator and nondividing orders \(ga,gb\), where
\(a,b\ge2\) are coprime. Write their leading coefficients as
\(\alpha h^a,\beta h^b\). For each geometric root of \(h\) of
multiplicity \(m>g\), translating that root to zero and taking Fourier
coordinates gives a common Newton edge crossing the diagonal. Its
normalized endpoints \(C_0=(u,v)\), \(C_1=(r,t)\), ordered by decreasing
second coordinate, satisfy
\[
 3\le u\le g,\qquad u<v\le m,\qquad 0\le t<r<u .
\]
If \(u=g\), then \(v=m\). For an integer \(k\) with \(2\le k<u\),
the bracket-mate endpoint is \((k,kv/u)\) and
\[
        v=\frac{u\bigl(u-r-(k-1)t\bigr)}{u-kr},
        \qquad u-kr>0 .
\]
The primitive direction and the successor restriction are those of
Guccione–Guccione–Valqui, Proposition 7.3(5)–(6).

Extend constants to an algebraic closure and scale to \([Q,P]=1\).
For the chosen root \(x_0\), put \(X=D,\ Y=-(x-x_0)\), so
\([Y,X]=1\). In normal order \(X^iY^j\) the maximal \(X\)-exponents
are \(ga,gb\), with coefficient polynomials proportional to
\(h(x_0-Y)^a,h(x_0-Y)^b\). Reordering contractions lowers the
\(X\)-exponent. For small positive \(s\), the endpoints at weights
\((1,-s)\) are \(a(g,m),b(g,m)\), of negative coordinate difference.

As \(s\) increases toward one, the largest \(Y\)-exponent endpoint at
each edge is inherited from the preceding slope. The two incoming
endpoints start aligned in ratio \(a:b\).
[Guccione–Guccione–Valqui, Lemma 7.1](sources/guccione-guccione-valqui-1111.6100v3.pdf)
gives positive weights, zero leading bracket and simultaneous edges
from their alignment. Theorem 2.11(2) gives powers of a common
homogeneous polynomial, so the smallest \(Y\)-exponent endpoints
remain aligned and the induction continues. At \(s=1\) the entering
endpoints still align and Lemma 7.1 gives positive \((1,-1)\)-weights.
At every intervening edge, the difference \(i-j\) increases from the
largest to the smallest \(Y\)-exponent, since \(0<s<1\).
Theorem 4.1(3) excludes an endpoint on the diagonal. Hence a
common edge crosses it strictly.

The proof of Proposition 7.3 applies to this edge via its positive
weights, coprime nondividing ratio, zero bracket and opposite
diagonal signs, as in the
[arbitrary-pair crossing calculation](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md#from-a-sharp-degree-bound-to-a-crossing-edge).
In particular the normalized endpoints are integral and the
bracket-mate endpoint is a positive proper fraction of \((u,v)\).
The first coordinate \(u\) is at most \(g\) by the differential
order bound. The second coordinate of the leading point can only
decrease as the positive penalty \(s\) increases, so \(v\le m\).
At any positive \(s\), the only possible endpoint on the maximal
\(X\)-face is its point of smallest \(Y\)-exponent. This gives
\(v=m\) when \(u=g\). Conversely if \(v=m\) and \(u<g\), the
point \((g,m)\) on the normalized support has strictly larger
weight, so \(u=g\) if and only if \(v=m\).

With \((f_1,f_2)=(k,kv/u)\), the direction is proportional to
\((kv/u-1,1-k)\). Equal endpoint weights give
\[
 (u-kr)\frac{v}{u}=u-r-(k-1)t=(u-kr)+(k-1)(r-t).
\]
Since \(v/u>1\) and \(r>t\), the coefficient \(u-kr\) is positive.
This proves the formula and the bound \(v\le u(u-1)\).

For \(u\le5\), the formula, integrality and \(v>u\) give the seven
rows below before Proposition 7.3(6). They can be obtained by taking
\(2\le k<u\), \(1\le r<u\), and \(0\le t<r\).
For the rejected rows the direction denominator
\(c=\gcd(f_1-1,f_2-1)\) is one and the successor from Proposition
7.3(6) is \(C_2=(j-1/\rho,j)\) with \(j\ge2\).

| \(C_0\) | \(C_1\) | \((f_1,f_2)\) | Direction | Forbidden \(C_2\), if any |
|---|---|---|---|---|
| \((3,6)\) | \((1,0)\) | \((2,4)\) | \((3,-1)\) | \((5/3,2)\) |
| \((4,6)\) | \((1,0)\) | \((2,3)\) | \((2,-1)\) | \((5/2,3)\) |
| \((4,12)\) | \((1,0)\) | \((3,9)\) | \((4,-1)\) | |
| \((5,10)\) | \((2,1)\) | \((2,4)\) | \((3,-1)\) | \((8/3,3)\) |
| \((5,10)\) | \((1,0)\) | \((3,6)\) | \((5,-2)\) | \((9/5,2)\) |
| \((5,15)\) | \((2,0)\) | \((2,6)\) | \((5,-1)\) | |
| \((5,20)\) | \((1,0)\) | \((4,16)\) | \((5,-1)\) | |

The [exact checker](verify_weyl_translated_root_crossings.py) compares
the formula with a separate scan of the bounded equal-weight
equations and evaluates the successor restriction for these rows.
The bound in the theorem supplies the scan cutoff.

## Gcd four and differential order eight

**Corollary PolynomialWeylGcdFourNonconstantLeadingFactor.** In a proper
polynomial Weyl pair of nondividing orders \(4a,4b\), if \(h\) is
nonconstant, all its roots have multiplicities in \(\{1,2,3,12\}\)
and at least one has multiplicity twelve. The unordered quotients are
\[
                  \{a,b\}=\{2j+3,3j+4\},\qquad j\ge0.
\]
If there is more than one geometric root, at least two have
multiplicity less than four. For \(j=0\) one has \(\deg h\ge24\).

The [nonconstant-leading degree bound](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#nonconstant-leading-factors-and-bernstein-degree)
gives \(d=\deg h\ge12\). If \(h\) has one root, it has multiplicity
\(d>4\). If there are multiple roots, the
[residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
supplies a root of multiplicity greater than four. At any such
root the table permits only \(C_0=(4,12)\), so \(u=g\) and
the root has multiplicity twelve. The residue theorem also excludes
multiplicity four. In the multiple-root case, the
[perfect-power obstruction](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#perfect-power-leading-factors)
forces at least one multiplicity less than four, and these
multiplicities sum to a positive multiple of four since \(4\mid d\).
Thus there are at least two.

At the edge \((4,-1)\) the shape is \(R=X L(XY^4)\) with
\(\deg L=3,\ L(0)\ne0\). The
[bracket-mate calculation](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_FACES.md#the-polynomials-on-the-two-edges)
applies with \(a\) in place of two by dividing the homogeneous bracket
identity by the corresponding multiple of \(aR^{a-1}\).
Its equation gives a polynomial of degree two vanishing at every
root of \(L\), so the maximal multiplicity \(\mu\) is two or three.
The [successor calculation for arbitrary coprime quotients](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_THIRTY_THREE.md#the-necessary-scalar-successor)
uses only this face, its endpoint signs, positive weights in ratio
\(a:b\), and the scalar commutator. These are all supplied by the
translated-root crossing. In particular the calculation applies
even when \(4+d>16\).

In that calculation the \(\mu=2\) row can only reach a scalar
edge at \(C=(3/2,2)\) on the lattice of denominator four or at
\(C'=(7/8,1)\) on the lattice of denominator eight. At either
endpoint the scalar-edge denominator is greater than both
quotients, eliminating the row. In the \(\mu=3\) row the scalar
edge, possibly after a single-root shift preserving the endpoint,
has \(C=(7/4,3)\) and
\[
     4s=\frac{7N-4}{3N-1},\qquad N=a+b .
\]
The cancelled gcd divides five. The reduced denominator must
divide one of \(a,b\), since otherwise both face polynomials vanish
at the origin and their bracket cannot be scalar. Thus the gcd
is five and \(N=5j+7\). The denominator is \(3j+4>N/2\), so one
quotient equals \(3j+4\), and the other is \(2j+3\).
For \(j=0\) the sheared Bernstein degrees are \(3(d+4),4(d+4)\).
The possibilities \(d=12,16,20\) give \((48,64),(60,80),(72,96)\).
The [exact degree-pair transfer](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FORTY_NINE_FROM_F7.md#two-inputs-from-polynomial-degree-pairs)
and [Moh's planar degree-100 theorem](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds)
exclude all three. Here \(4\mid d\) also holds in the one-root case,
where \(d=12\).

**Corollary PolynomialWeylOrderEightGeneration.** A characteristic-zero
polynomial Weyl pair with a member of differential order eight
and nonzero scalar commutator generates.

Suppose it is proper. The
[degree-59 theorem](POLYNOMIAL_WEYL_BERNSTEIN_FIFTY_SIX.md)
excludes constant leading coefficient at order eight. Keep the
order-eight member fixed and subtract its powers from its mate
until the mate's order is no longer a multiple of eight. The
[prior order bound](POLYNOMIAL_GCD_THREE_LEADING_FACTOR.md#consequences-for-reduced-orders)
makes the mate's order at least eight, hence the two orders are
nondividing. Their gcd is four by the exclusions for gcd one and
two. The first quotient is then two, absent from the family above.
It follows that each differential order of a proper pair is
composite and at least ten.

## Gcd five

**Corollary PolynomialWeylGcdFiveNonconstantLeadingFactor.** In a proper
polynomial Weyl pair of nondividing orders \(5a,5b\) with nonconstant
common factor \(h\), there is a geometric root of multiplicity greater
than five. At each such root the translated crossing and unordered
quotient pair satisfy one of the following alternatives.

| Crossing | Root multiplicity of \(h\) | Maximal root multiplicity of the one-variable face factor | \(\{a,b\}\) |
|---|---:|---:|---|
| \((4,12)\) | at least 13 | 3 | \(\{2j+3,3j+4\}\), \(j\ge0\) |
| \((5,20)\) | 20 | 2 | \(\{2j+3,j+2\}\), \(j\ge0\) |
| \((5,20)\) | 20 | 3 | \(\{3j+2,4j+3\}\), \(j\ge0\) |
| \((5,20)\) | 20 | 3 | \(\{6j+4,j+1\}\), even \(j\ge2\) |
| \((5,20)\) | 20 | 4 | \(\{4j+5,7j+9\}\), \(j\ge0\) |
| \((5,20)\) | 20 | 4 | \(\{8j+10,3j+4\}\), odd \(j\ge1\) |

The infinity bound gives \(d\ge11\). There is a root of multiplicity
above five, either by this bound in the one-root case or by the
residue theorem in the multiple-root case. The crossing list is
\((4,12),(5,15),(5,20)\). For \((4,12)\) the calculation above
is local, so both the quotient family and its multiplicity-three
face apply. Since \(u<g=5\), the high root has multiplicity \(m>12\).

For \((5,15)\), the [sum-twenty calculation](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md#sum-twenty)
gives the shifted endpoint \(C=(13/5,3)\) and lattice denominator
five. The commuting-successor tests use the equation of the
common factor and are independent of the quotients. At a scalar
successor the reduced denominator of \(5s\) is at least
\((3N-1)/2>N-1\), where \(N=a+b\). Indeed the cancelled gcd
divides \(\beta=5(3-C_x)=2\). The denominator is coprime to three
and must divide one of \(a,b\) to allow a nonzero bracket at the
origin. This extends the scalar exclusion to every coprime
nondividing ratio. Hence the row \((5,15)\) is excluded.

For \((5,20)\), apply the
[sum-twenty-five calculation](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md#degree-fifty-and-the-sum-twenty-five-crossing)
at the local edge. Its common factor is \(X L(XY^5)\) with
\(\deg L=4,\ L(0)\ne0\), and the bracket mate shows that \(L\)
has at most three roots. The maximal root multiplicity \(\mu\)
is therefore two, three or four. The shifted corner is
\(C=(1+\mu/5,\mu)\) on the lattice of denominator five. The
commuting-edge lattice and the integer-multiplier obstruction
depend on \(C\) and the lattice, since the bracket-mate equation
is divided by the common power. They exclude every multiple-root
commuting successor. Single-root shifts preserve the endpoint and
lattice while strictly increasing the slope in a finite list.
Thus a scalar edge must occur at \(C\). There
\[
     5s=5-\frac{\beta N}{\mu N-1},\qquad \beta=4\mu-5\in\{3,7,11\}.
\]
The cancelled gcd divides \(\beta\). As the reduced denominator
\(D\) is coprime to \(\mu\), a nonzero scalar bracket requires
\(D\mid a\) or \(D\mid b\). If the gcd were one, \(D=\mu N-1\)
would exceed both quotients. Thus
\[
             D=\frac{\mu N-1}{\beta},\qquad
             \{a,b\}=\{kD,N-kD\},\qquad
             1\le k<\beta/\mu .
\]
Here \(kD\le N-2\) gives the strict bound on \(k\), and
\(\gcd(D,N)=1\) reduces coprimality to \(\gcd(k,N)=1\).
Solving \(\mu N\equiv1\bmod\beta\) gives
\(N=3j+5,\ 7j+5,\ 11j+14\) for \(\mu=2,3,4\), respectively.
The restriction on \(k\), coprimality and the lower bound two
on each quotient give the last five rows. Since \(u=g\) here,
the high root has multiplicity exactly twenty.

## The order-ten case

**Corollary PolynomialWeylOrderTenResidual.** In a proper
characteristic-zero polynomial Weyl pair with one member of
differential order ten, keeping that member fixed and subtracting
powers from its mate gives orders \(10,15\). Their nonconstant common
leading factor has degree \(d\ge30\) divisible by five. There is
at least one root of multiplicity twenty and at least two other
roots of multiplicities in \(\{1,2,3,4\}\). Every root of
multiplicity greater than five has multiplicity twenty. At each
such root the common translated Fourier face at weight \((5,-1)\)
is \(X L(XY^5)\), with \(L\) of degree four, nonzero constant
term and maximal root multiplicity two or three.

The cutoff through 59 excludes constant leading coefficient at
order ten. After reducing its mate, all proper divisors of ten
are below the lower bound on the mate's order. Thus the orders are
nondividing, their gcd is five and one quotient is two.
The gcd-five table gives only \((5,20)\) at each high root and
allows exactly the other quotient three, with face multiplicity
two or three.

The sheared degrees are \(2(d+5),3(d+5)\). The choice \(d=20\)
is excluded by the Weyl cutoff, and \(d=25\) by transfer of
\((60,90)\) and Moh's bound. Since \(d=20\) is the only one-root
possibility and the multiple-root case has \(5\mid d\), this gives
\(d\ge30\). The perfect-power obstruction then gives at least one
root of multiplicity below five, and their sum is a multiple of
five, so at least two such roots occur. At \(d=30\) the transferred
degrees are \((70,105)\), one of the cases retained by Nguyen's
degree-105 reduction ([Nguyen, Proposition 3.7](https://arxiv.org/abs/1902.05923)).
Under the planar bound through 107,
this case is excluded and \(d\ge35\).

One boundary pattern at \(d=30\) illustrates what the
[residue ramification condition](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#a-rational-map-from-the-residue-equation)
can determine about the roots. For multiplicities
\((20,2,2,2,2,2)\), the degree of \(B\) is six,
\(\Delta=3\), and \(S\) has degree five with a unique pole
at the multiplicity-twenty root \(x_0\). In the coordinate
\(z=1/(x-x_0)\), it is a degree-five polynomial, and its
fibre at \(S(0)\) is totally ramified at \(z=0\).
Hence \(S(z)=c_0+c_1z^5\) with \(c_0c_1\ne0\). In the
original variable \(t=x-x_0\), the residue equation therefore
forces
\[
       h=C\,t^{20}(t^5-\lambda)^2,\qquad
       B=\frac{t(t^5-\lambda)}{3\lambda},
       \qquad C\lambda\ne0 .
\]
Conversely these polynomials satisfy the residue equation
by direct substitution.

## The order-fifteen case

**Corollary PolynomialWeylOrderFifteenResidual.** A proper
characteristic-zero polynomial Weyl pair with a member of
differential order fifteen has a fixed-member reduction whose
other order is \(10,20,25\) or \(80\).

The leading coefficient is nonconstant by the cutoff through 59.
The reduced order is at least ten, so the two orders are
nondividing because no proper divisor of fifteen reaches ten.
The reduced gcd is five by the exclusions through gcd three.
In the gcd-five table, quotient three in the \((4,12)\) row
forces quotient four. In the \((5,20)\) rows it allows exactly
two, five or sixteen. Multiplying by five gives the list.
For mate order twenty, every high root must have the \((4,12)\)
crossing. For mate orders \(10,25,80\), every high root must
have the \((5,20)\) crossing and multiplicity exactly twenty.

## Primary comparisons

Guccione–Guccione–Valqui's proof of Proposition 7.3 is stated in
the standard minimal presentation. Its edge inputs apply after
translation at a high-multiplicity root as verified above, giving
a bound on the normalized first coordinate by the differential
order gcd. The fractional successors use their
[Part II, Proposition 4.3](sources/guccione-guccione-valqui-1205.6827v1.pdf).
The local corner conditions alone leave the displayed gcd-four
quotient family and the order-ten pair of orders \(10,15\). The
global polynomial Weyl equations for those residual cases remain
to be settled.
