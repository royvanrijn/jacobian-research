# Weyl generation through Bernstein degree and constant-leading order fifty-three

**Theorem PolynomialWeylBernsteinThroughFiftyThree.**
A characteristic-zero polynomial Weyl pair with nonzero scalar
commutator generates if one member has Bernstein degree at most
\(53\), or has nonzero constant differential leading coefficient
and differential order at most \(53\). In a proper pair each member
has differential order plus leading-coefficient degree at least \(54\).

Start from the
[bound through fifty-one](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md)
and its direct degree-forty-two exclusion.
Apply its axis orientation, fixed-member descent, infinity shear
and first-crossing argument to a hypothetical member of Bernstein
degree \(52\). This gives a shaped pair of nondividing ratio \(a:b\)
and normalized outer corner \((g,d)\), with
\(52=a(g+d)\), \(g+d\ge16\), \(a,b\ge2\) coprime.
Thus \(a=2\), \(g+d=26\).
In its Fourier presentation the support of the first member is
in \([0,2g]\times[0,2d]\), and the maximal \(X\)-exponent face is
the single corner \((2g,2d)\). The corresponding statements hold
for the mate with factor \(b\).

## The obstruction from the preceding edge

The earlier ratio-\(2:b\) calculations exclude all first-crossing
sums through twenty-five except a potentially embedded sum-twenty-one
crossing. In particular a sum-twenty-five crossing
forces ratio \(2:3\) by the local successor test, and the
[nondividing-degree transfer](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md)
followed by [Moh's degree bound](https://doi.org/10.1515/crll.1983.340.140)
excludes the actual degrees \((52,78)\).
In that case
\[
                  C_0=(6,15),\qquad(\rho,\sigma)=(3,-1)
\]
by the
[degree-forty-two crossing calculation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md).
The outer rectangle must have \(g\ge6\), and its corner's
\((3,-1)\)-weight satisfies \(4g-26\le3\), so \(g\le7\).
If \(g=6\), the term at \(C_0\) contradicts the pure maximal
\(X\)-exponent face. Thus \((g,d)=(7,19)\).

Use weights \((1,-s)\) for \(0<s<1/3\).
The actual term at \(C_0\) gives weights at least
\(2(6-15s)>2\) and \(b(6-15s)>b\) for the two members.
Their leading symbols therefore have zero Poisson bracket
throughout this interval, since their weight sum minus the
contraction loss \(1-s\) is strictly positive.
Except at the finitely many edge weights, both leading forms
are monomials. Their zero bracket makes their exponent vectors
proportional, so the ratio of the two weights is constant on
each such interval. Continuity and positivity give the same
ratio on the entire interval, and the pure maximal \(X\)-faces
at \(s=0\) give that ratio as \(2:b\).
Thus the normalized maximizing exponent vectors at generic
weights are integral, since \(\gcd(2,b)=1\).

For small positive \(s\) this common normalized vertex is
\((7,19)\), and for \(s<1/3\) sufficiently close to \(1/3\)
it is \((6,15)\), the highest \(Y\)-endpoint of the crossing.
Every change of normalized maximizing vertex as \(s\)
increases strictly lowers the \(X\)-coordinate. Indeed the
maximizers of \(i-sj\) at increasing \(s\) have nonincreasing
\(j\), and at any tie with \(s>0\) distinct such vertices
with different \(j\) also have different \(i\), ordered in
the same direction. Hence there are no intervening vertices
between these two. They are joined by an edge of slope
\(s=1/4\), of weight \((4,-1)\).

At that edge the leading forms are proportional to \(R^2,R^b\),
with
\[
                     R=X^6Y^{15}L(XY^4),
              \qquad \deg L=1,\quad L(0)\ne0.
\]
For example, the common power form follows from the zero
bracket and weighted Euler differentiation, followed by
polynomial factorization and \(\gcd(2,b)=1\).
The endpoints give the indicated monomial factor and the
unit lattice segment in the variable \(z=XY^4\).
By [Guccione–Guccione–Valqui, Theorem 4.1(4)](sources/guccione-guccione-valqui-1111.6100v3.pdf),
the leading form of \(P\) has a polynomial homogeneous bracket mate
of weight three. Since that leading form is proportional to \(R^2\),
the homogeneous mate \(F\) satisfies \(\{R,F\}=cR\) for
a nonzero scalar \(c\), up to the sign of the convention.
Every monomial of weight three has the form
\(XY(XY^4)^k\), \(k\ge0\), so \(F=XYG(z)\).
The bracket equation becomes
\[
                 9zG'+\left(-9-3zL'/L\right)G=c.          \tag{1}
\]
At the nonzero root of \(L\) the equation forces \(G\) to
vanish. Thus \(\deg G>0\).
For \(N=\deg G\), comparison of the leading term in (1)
at infinity gives \(9N-12=0\), impossible for an integer \(N\).
This rules out the embedded sum-twenty-one crossing.

## The sum-twenty-six corner test

For a first normalized corner \(C_0=(u,v)\) of sum twenty-six,
Guccione–Guccione–Valqui, Proposition 7.3(1)–(5), requires
\(3\le u<v\), \(\gcd(u,v)>1\), a point
\((f_1,f_2)=\mu(u,v)\) with \(0<\mu<1\) and integral entries,
and a direction \((\rho,\sigma)=\operatorname{dir}(f_1-1,f_2-1)\).
The other endpoint \(C_1=(r_0,s_0)\) must satisfy
\(0\le s_0<r_0<u\) at the same weight. The only \(u\)
are \(4,6,8,10,12\), each with \(\gcd(u,v)=2\), so
\((f_1,f_2)=(u/2,v/2)\). The direction and target weight
are as follows.

| \(u\) | \(v\) | \((\rho,\sigma)\) | Weight of \(C_0\) |
|---:|---:|:---:|---:|
| 4 | 22 | \((10,-1)\) | 18 |
| 6 | 20 | \((9,-2)\) | 14 |
| 8 | 18 | \((8,-3)\) | 10 |
| 10 | 16 | \((7,-4)\) | 6 |
| 12 | 14 | \((6,-5)\) | 2 |

In the first three rows \(r_0=1\) gives strictly smaller
weight, whereas \(r_0\ge2\) gives weight at least
\(2(\rho+\sigma)-\sigma\), strictly larger.
In the last two rows every \(r_0\ge1\) gives weight at
least \(\rho\), strictly larger.
Thus there is no \(C_1\), before the further endpoint
test of Proposition 7.3(6).
This exhausts the first-crossing possibilities and
excludes Bernstein degree \(52\).

For Bernstein degree \(53\), perform the same mate descent
as in the through-fifty-one proof. It would end at a proper
nondividing pair with smaller degree at most \(53\).
The smaller degree is not at most \(52\), by the preceding
bounds, and is not \(53\), since this is prime and its degree
gcd would have to be at least \(16\) and a proper divisor.
This proves the Bernstein assertion.

For starting constant-leading orders \(52,53\), the
same weighted descent gives a dividing-mate case or a
binomial coordinate change of smaller order at most
\(26\), already excluded, or an intermediate change to a
member whose order plus leading-coefficient degree is at most
the starting order. The intermediate member either has
constant leading coefficient and smaller order, or contradicts
the individual leading-degree bound from infinity shear
and the dividing-mate argument using the new Bernstein
cutoff. The weight-one branch has Bernstein degree equal
to the starting order after depression. This proves the
constant-leading and individual leading-degree assertions.
As in the preceding proof strengthening, finite coefficient
fields embed into \(\mathbb C\), and generation descends
along field extensions, giving the stated characteristic-zero
scope.
