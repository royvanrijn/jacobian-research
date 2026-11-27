# Weyl generation through Bernstein degree and constant-leading order fifty-one

**Theorem PolynomialWeylConditionalThroughFiftyOneFromF7.**
A characteristic-zero polynomial Weyl pair with nonzero scalar
commutator generates when one member has Bernstein degree at most
\(51\), or has nonzero constant differential leading coefficient
and order at most \(51\). Each member of a proper pair then has
order plus degree of its differential leading coefficient at least \(52\).

We start from the
[bound through forty-nine](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FORTY_NINE_FROM_F7.md).
The axis orientation, fixed-member mate descent and infinity shear
in that proof apply to the new Bernstein degrees \(n=50,51\).
They give a shaped pair with order ratio \(a:b\), coprime and
nondividing, and normalized outer corner \((g,d)\) with
\(n=a(g+d)\), \(g+d\ge16\).
The maximal \(X\)-face in its Fourier presentation is its single
corner. We use the same arbitrary-pair crossing argument based
on [Guccione–Guccione–Valqui, Lemma 7.1 and Proposition 7.3](sources/guccione-guccione-valqui-1111.6100v3.pdf)
as in the degree-forty-one proof.
Its first normalized crossing lies within that rectangle.
Its integer corner tests are independent of the ratio,
apart from the positive nondividing hypotheses.

## An integer-multiplier obstruction

Suppose \(f\) is a polynomial of degree \(d_f\) with \(r\ge2\)
distinct roots and \(j\ge1\) is an integer. Then
\(1/f^j\) has no rational primitive.
If it has a simple pole, that already prevents a primitive.
Otherwise a primitive would have, at every finite root of
multiplicity \(\nu\), a pole of order \(j\nu-1\), and it
would have no other pole. Indeed it is finite at infinity
because its derivative decays there as \(z^{-j d_f}\).
Its degree would therefore be \(j d_f-r\). But the primitive
has a first nonconstant term at infinity of order
\(z^{-(j d_f-1)}\), so its ramification index over its value
there is \(j d_f-1>j d_f-r\), impossible for that degree.

In particular a commuting successor with at least two distinct
roots and with integer multiplier \(M/\mu\) is impossible.
In the notation of the
[degree-forty-two successor test](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#testing-the-commuting-successors),
its equation is \(W'-(M/\mu)(f'/f)W=c\), \(c\ne0\).
It would make \(W/f^{M/\mu}\) a rational primitive of
\(c/f^{M/\mu}\), contrary to the preceding argument.

## Degree fifty and the sum-twenty-five crossing

For \(n=50\) the normalization forces \(a=2\) and \(g+d=25\).
Through twenty, the existing successor exclusions for
ratio \(2:b\) apply. The sum-twenty-one crossing
\(C_0=(6,15)\), of weight \((3,-1)\), cannot fit in the
sum-twenty-five rectangle. Its bound on the outer corner
is \(4g-25\le3\), while \(g\ge6\). If \(g=6\), the
crossing endpoint contradicts the pure maximal \(X\)-face.
If \(g=7\), the outer corner \((7,18)\) lies on the
same edge and has larger \(Y\)-exponent than \(C_0\),
which is defined as its highest \(Y\)-endpoint.
The sum-twenty-two and sum-twenty-four successor exclusions
used in the bound through forty-nine remain
applicable, and sum twenty-three is prime.

The full integer table at sum twenty-five gives
\[
 \begin{array}{c|c|c|c}
 C_0&C_1&(f_1,f_2)&(\rho,\sigma)\\ \hline
 (5,20)&(1,0)&(4,16)&(5,-1)\\
 (10,15)&(3,1)&(2,3)&(2,-1)
 \end{array}                                                  \tag{1}
\]
using Proposition 7.3(6) as in the earlier crossing tables.
An exact replay for both new sums is below.
For \(C_1=(h,j)\), use the common factor
\(X^hY^j L(XY^\rho)\) and the homogeneous bracket-mate
degree formula in the sources for the degree-forty-nine bound.
For the first row \(\deg L=4,\deg G=3\), so
\(L\) has at most three distinct roots and its largest
root multiplicity \(\mu\) is \(2,3,\) or \(4\).
For the second row \(\deg L=7,\deg G=1\), so
its unique root has multiplicity seven.
Here \(L(0)\ne0\). In the common face variable
\(X^{1/\rho}Y\) the zero root, when present, has
multiplicity \(j<\mu\).

Use [Guccione–Guccione–Valqui II, Proposition 4.3](https://arxiv.org/pdf/1205.6827v1).
The incoming faces are nonmonomial and commuting,
with positive weight ratio \(2:b\), a nonzero scalar
commutator, and negative endpoint diagonal difference.
The first fractional shifts have data
\[
 \begin{array}{c|c|c|c|c}
 C_0&\mu&C&l&\beta=l(\mu-C_x)\\ \hline
 (5,20)&2,3,4&(1+\mu/5,\mu)&5&4\mu-5\\
 (10,15)&7&(6,7)&2&2 .
 \end{array}                                                  \tag{2}
\]
Indeed the initial common factor is
\(t^{\rho h-j}E^jL(E^\rho)\), with
\(t=X^{1/\rho}\), \(E=tY\). The maximal root
shift gives normalized endpoint
\(((\rho h-j+\mu)/\rho,\mu)\).
All endpoints in (2) again have negative diagonal
difference.

For a scalar successor of slope \(s\), put \(N=b+2\).
The normalized endpoint equation gives
\[
             ls=l-\frac{\beta N}{N\mu-1}.
\]
Let \(D\) be the reduced denominator of \(ls\).
The lattice must permit a nonmonomial smaller face,
so \(D\le2\mu\), and at least one of \(2\mu,b\mu\)
must be divisible by \(D\). Indeed, otherwise both
face polynomials have zero constant term and their
leading Poisson bracket cannot be scalar.
Since \(D=(N\mu-1)/\gcd(N\mu-1,\beta)\) and is
coprime to \(\mu\), the row with \(\beta=2\le\mu\)
cannot have a scalar successor: \(D>N-1\) there.
In the first row the inequality \(D\le2\mu\)
bounds the integer \(N\) by \(2\beta\).
The finite scalar test leaves only
\[
              (\mu,b,s,D)=(2,3,2/3,3),\quad(3,3,1/2,2).
\]
In either case the actual Bernstein degrees are
\((50,75)\). The nondividing-degree center transfer
in the through-forty-nine proof gives a plane Keller
counterexample of those degrees, contrary to
[Moh's published degree-one-hundred bound](https://doi.org/10.1515/crll.1983.340.140).

For commuting successors with multiple roots,
the same lattice tests as in the earlier proofs give
\[
 \begin{array}{c|c|c|c}
 \mu & M & D & j_0\\ \hline
 2&2,4&1&0,1\\
 3&3,15&2&1\\
 4&4&3&1\\
 4&12&1&0,1,2,3\\
 7&7&3&1,4
 \end{array}
\]
where \(j_0\) is the least exponent in the one-variable
face polynomial. Each multiplier \(M/\mu\) is an
integer, excluded by the rational primitive argument.
A commuting single-root face must have a nonzero
root, since zero alone would give a monomial face.
Its binomial-power form forces \(ls\) integral.
Shifting at that root preserves \(C,l\), and strictly
increases \(s\); the positive-weight and ratio
hypotheses of Proposition 4.3 continue to hold.
There are finitely many such slopes between
\(1/l\) and \(C_x/\mu\). They therefore cannot avoid
the exclusions of all remaining successor types.
This excludes degree fifty.

## Degree fifty-one and the order bounds

For \(n=51\) the normalization instead forces
\(a=3\), \(g+d=17\).
There is no crossing of sum at most fifteen by the
integer tests behind Guccione–Guccione–Valqui,
Corollary 7.4. Sum seventeen is prime. At sum sixteen
the full corner table gives \(C_0=(4,12)\) of
weight \((4,-1)\) and \(C_0=(6,10)\) of
weight \((2,-1)\). Neither fits a sum-seventeen
outer vertex with pure maximal \(X\)-face.
Indeed for \(C_0=(u,v)\) one has \(g\ge u\) and
\((\rho+1)g-17\le \rho u-v\), forcing \(g=u\)
in both rows. Then \(d=v+1\), so that \(C_0\)
contradicts the single maximal \(X\)-face.
This excludes degree fifty-one without a
quotient-three successor classification.

The constant-leading weighted descent cited in the
through-forty-nine proof now applies at orders
\(50,51\). A residual dividing mate or a binomial
change gives smaller order at most \(25\), already
excluded. The intermediate change either gives a
constant-leading member of that smaller order, or a
member whose order plus nonconstant leading degree
is bounded by the original order. The latter
contradicts the updated Bernstein bound by infinity
shear and dividing-mate descent. The weight-one
branch has Bernstein degree equal to the starting
order after depression. This proves the claimed
constant-leading cutoff. With the Bernstein cutoff proved,
the individual leading-degree argument in the
[degree-forty-one proof](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md#from-a-sharp-degree-bound-to-a-crossing-edge)
applies with threshold \(52\). The nondividing case is put into
that Bernstein bound by infinity shear. In the dividing case
the leading-degree sum is a positive integer multiple of the
smaller member's sum, so induction on order gives the same
lower bound. For constant leading coefficient the lower bound
is the constant-leading cutoff itself. The field extension and
descent are as in the through-forty-nine theorem.

## Exact arithmetic replay

The [exact checker](verify_weyl_f7_through_fifty_one.py) enumerates the integer
corners and the scalar and commuting successor lattices.

