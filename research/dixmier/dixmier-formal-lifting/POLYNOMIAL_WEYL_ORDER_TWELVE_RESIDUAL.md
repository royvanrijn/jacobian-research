# The differential-order-twelve Weyl residual

**Theorem PolynomialWeylOrderTwelveResidual.** Let \(P,Q\) be a
proper polynomial Weyl pair over a characteristic-zero field with one
member of differential order twelve. Keep that member fixed and subtract
powers of it from its mate until neither order divides the other.
The resulting mate and nonconstant common leading factor \(h\) satisfy
the following necessary conditions. The degree \(d=\deg h\) is divisible
by the reduced differential-order gcd.

| Reduced mate order | Reduced gcd | Multiplicities of roots above the gcd | Lower bound for \(d\) |
|---:|---:|---|---:|
| 16 | 4 | 12 | 24 |
| 18 | 6 | 26 or 30 | 30 |
| 30 | 6 | 30 | 30 |
| 42 | 6 | 15 | 24 |
| 78 | 6 | 30 | 30 |

There is at least one root above the gcd in every row. In the
\((12,18)\) row with \(d=30\), either \(h\) has a single root of
multiplicity thirty or it has a root of multiplicity twenty-six and
all other multiplicities sum to four. The normalized Bernstein
degrees in this boundary case are \((72,108)\).

The [primitive leading-factor obstruction](POLYNOMIAL_WEYL_PRIMITIVE_LEADING_FACTOR_OBSTRUCTION.md#consequences-for-centered-rectangles-and-low-orders)
requires a common multiplicity divisor greater than one in every
row. In the \((12,18)\) row, the
[next-face exclusion](POLYNOMIAL_WEYL_ORDER_TWELVE_HIGH_ROOT_FRACTIONAL_EXCLUSION.md#exclusion-by-the-next-newton-face)
rules out every multiplicity-26 root. The
[multiplicity-thirty argument](POLYNOMIAL_WEYL_ORDER_TWELVE_HIGH_ROOT_FRACTIONAL_EXCLUSION.md#the-multiplicity-thirty-root-forces-a-sixth-power)
then forces a single root of multiplicity thirty over the
algebraic closure. The coefficient-field test in the same
source excludes this entire row over \(\mathbb Q\).

## The crossings at gcd six

The [Weyl bound through 59](POLYNOMIAL_WEYL_BERNSTEIN_FIFTY_SIX.md)
excludes constant leading coefficient at order twelve. Fixed-member
subtraction does not change that coefficient. The mate retains
differential order at least ten by the
[order-eight exclusion and prior order bounds](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#gcd-four-and-differential-order-eight).
After subtraction it is not a multiple of twelve, and no proper
divisor of twelve reaches ten. Thus the orders are nondividing. The
exclusions through gcd three leave gcd four or six. At gcd four,
quotient three in the
[gcd-four family](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#gcd-four-and-differential-order-eight)
forces the pair of quotients \(\{3,4\}\), giving the first row.

For gcd six the quotients have ratio \(2:b\), with \(b\ge3\)
odd. The
[nonconstant-leading normalization](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#nonconstant-leading-factors-and-bernstein-degree)
gives \(6+d\ge16\). There is a root of multiplicity above six,
either by this bound in the one-root case or by the
[residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
in the multiple-root case. Apply the
[translated-root crossing theorem](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#the-finite-list-for-a-fixed-differential-order-gcd)
at any such root. Besides the three rows through first coordinate five,
the same integer tests at first coordinate six give four rows.
The sole additional equal-weight candidate \((6,9)\) has
\(C_2=(7/2,4)\), excluded by Proposition 7.3(6).

| \(C_0\) | \(C_1\) | \((f_1,f_2)\) | Direction |
|---|---|---|---|
| \((6,10)\) | \((1,0)\) | \((3,5)\) | \((2,-1)\) |
| \((6,12)\) | \((2,0)\) | \((2,4)\) | \((3,-1)\) |
| \((6,15)\) | \((1,0)\) | \((4,10)\) | \((3,-1)\) |
| \((6,30)\) | \((1,0)\) | \((5,25)\) | \((6,-1)\) |

The prior local successor exclusions remove
[\((4,12)\) at quotient two](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#gcd-four-and-differential-order-eight),
[\((5,15)\)](POLYNOMIAL_WEYL_TRANSLATED_ROOT_CROSSINGS.md#gcd-five),
[\((6,10)\)](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_THIRTY_THREE.md#the-necessary-scalar-successor)
and [\((6,12)\)](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md#sum-eighteen).
The \((5,20)\) row requires \(b=3\) by the gcd-five note's
local calculation. The
[\((6,15)\) successor calculation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#testing-the-scalar-successors)
requires \(b=7\). Its commuting-edge tests use only the shifted
corner, the lattice and the ratio, so this conclusion does not use
the actual Bernstein degrees \((42,147)\). Here its root
multiplicity is exactly fifteen, since \(u=g\).
The \((6,30)\) calculation below leaves \(b=3,5,13\) and
forces root multiplicity thirty, again by \(u=g\).

## The crossing at \((5,20)\) inside the gcd-six support

For \(C_0=(5,20)\) the initial point at small positive slope is
\((6,m)\) with \(m>20\). On an edge with slope \(s>0\), both
coordinates strictly decrease when passing from the upper to the
lower endpoint. Every intervening vertex is integral after dividing
by the quotient, since the commuting faces are powers with coprime
exponents. Thus there is exactly one preceding edge, from
\((6,m)\) to \((5,20)\), at direction \((\rho,-1)\) with
\(\rho=m-20\ge6\), because its slope is below \(1/5\).

The common factor on this edge is \(X^5Y^{20}L(z)\) with
\(z=XY^\rho\), \(\deg L=1\) and \(L(0)\ne0\). Its homogeneous
bracket mate has symbol \(XYG(z)\), since the polynomial monomials
of weight \(\rho-1\) have exponents \((1+i,1+\rho i)\) for
\(i\ge0\). Dividing the bracket identity by the common power, as in the
[sum-twenty calculation](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md),
gives
\[
      (5\rho-20)zG'+
      \left(-15-(\rho-1)z\frac{L'}{L}\right)G=c,\qquad c\ne0.
\]
The pole at the root of \(L\) forces \(G\) to vanish there.
The top-degree coefficient then gives
\[
                          \deg G=\frac{\rho+14}{5\rho-20}.
\]
For \(\rho\ge9\) this is below one. Of \(\rho=6,7,8\), only
\(\rho=6\) gives an integer. Hence \(m=26\).

## The local successors at \((6,30)\)

**Proposition PolynomialWeylSixThirtyScalarSuccessors.** A
strict polynomial Weyl crossing with \(C_0=(6,30)\),
\(C_1=(1,0)\), direction \((6,-1)\), and positive weights in
ratio \(2:b\) for odd \(b\ge3\), with zero leading bracket and
nonzero scalar commutator, requires \(b\in\{3,5,13\}\).

The common power factor is \(X L(XY^6)\), with \(\deg L=5\)
and \(L(0)\ne0\). Its bracket mate is \(XYG(XY^6)\).
The equation \(6zG'+(1-5zL'/L)G=c\ne0\) gives
\(\deg G=4\), and \(G\) vanishes at each root of \(L\).
The maximal root multiplicity \(\mu\) is thus in \(\{2,3,4,5\}\).

Use [Guccione–Guccione–Valqui II, Proposition 4.3](sources/guccione-guccione-valqui-1205.6827v1.pdf).
The crossing has negative second direction coordinate, both faces
have multiple terms, their weights are positive in the nondividing
ratio \(2:b\), and both upper endpoints are below the diagonal in
the sense that their first coordinate is smaller than the second.
The scalar commutator and zero leading bracket give the remaining
hypotheses. A shift at a maximal root produces the endpoint
\(C=(1+\mu/6,\mu)\) on lattice denominator \(l=6\), and its
successor has slope \(s>1/6\). At a general shifted endpoint
\(C=(\mu-\beta/l,\mu)\) in the calculation, the weight
\(w=C_x-\mu s\) is positive and \(\beta>0\).

At a scalar edge, \(Nw=1-s\), where \(N=2+b\). Thus
\[
            ls=l-\frac{\beta N}{\mu N-1}.
\]
The denominator \(D\) divides every exponent deficit in both
one-variable face polynomials and is coprime to \(\mu\).
The first face needs two terms, so \(D\le2\mu\).
The bracket is zero at the origin if both polynomials vanish there,
so \(D\mid 2\) or \(D\mid b\).
The cancelled gcd divides \(\beta\), giving
\(\mu N-1\le2\mu\beta\). We include \(N=2\beta+1\)
in the search when \(\mu=1\); this equality contributes to the table.

At a commuting edge the [bracket-mate equation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#testing-the-commuting-successors)
for the common one-variable factor \(f\) of degree \(\mu\) gives,
if \(f\) has multiple roots,
\[
 W'-\frac{M}{\mu}\frac{f'}{f}W=c,\qquad
 M=\deg W,\qquad ls=l-\frac{\beta M}{\mu(M-1)},\qquad c\ne0 .
\]
Here \(M\ge2\). With \(D\) the reduced denominator of \(ls\),
write \(f=\eta^j J(\eta^D)\), \(J(0)\ne0\). Then
\(D\mid\mu-j\), \(M-1\mid\beta(\mu-j)\), and
\(M\equiv1\bmod D\) by the character projection of \(W\).
The possible nonzero root multiplicities are the parts of
a partition of \((\mu-j)/D\), each repeated \(D\) times.
Include a root of multiplicity \(j\) if \(j>0\). There must be
between two and \(M\) distinct roots. A multiplicity \(m_0\)
with \(Mm_0=\mu\) is impossible by the simple zero condition
on \(W\). If \(\mu\mid Mm_0\) for all roots, then there is a
polynomial \(A\) with \(A'/A=(M/\mu)f'/f\) and the equation
becomes \((W/A)'=c/A\), excluded by the
[multi-root rational-primitive obstruction](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md#an-integer-multiplier-obstruction).

The next maximal root has multiplicity \(\nu\). The successor
proposition replaces \(C\) by
\((C_x-s(\mu-\nu),\nu)\); the lattice becomes \(l'=\mathrm{lcm}(l,\operatorname{den}(s))\)
for a nonzero shift. When zero has maximal multiplicity we choose
that root, and the operator remains on its original lattice \(l'=l\).
If the edge has a single root, it must be nonzero and \(ls\) is
integral by adjacency of the exponents. Its shift preserves
the endpoint and lattice and strictly increases the slope. Thus
there are only finitely many such steps at a corner, and the tests
at its initial lower slope also cover all larger slopes.

The [exact finite test](verify_weyl_order_twelve.py) bounds
\(M\) by \(\beta\mu+1\) and \(N\) by
\(\lfloor(2\mu\beta+1)/\mu\rfloor\).
It gives no scalar edge at any of the four initial shifted
endpoints. The table gives every permitted multiple-root shift
from them. A successor state is listed as \((l',\nu,\beta')\),
where \(\beta'=l'(\nu-C'_x)\). Each has \(\beta'>0\), so preserves
the required strict diagonal sign. There is no further
multiple-root successor from any state in the table.

| Initial \(\mu\) | \(M\) | \(D\) | Shift slope | Next state | Subsequent scalar \((b,s)\) |
|---:|---:|---:|---|---|---|
| 2 | 3 | 1 | \(1/2\) | \((6,1,1)\) | |
| 2 | 5 | 2 | \(7/12\) | \((12,1,3)\) | \((5,17/24)\) |
| 3 | 4 | 1 | \(1/3\) | \((6,1,1)\) | |
| 3 | 4 | 1 | \(1/3\) | \((6,2,5)\) | |
| 3 | 7 | 2 | \(5/12\) | \((6,1,2)\) | \((3,7/12)\) |
| 3 | 10 | 3 | \(4/9\) | \((18,1,7)\) | \((13,7/12)\) |
| 4 | 22 | 3 | \(7/18\) | \((6,1,3)\) | \((5,5/12)\) |

An actual scalar commutator therefore requires one of the listed
scalar edges. The computation gives necessary conditions on the
Weyl pair, including when the edge calculations have not tested all
coefficients of its lower terms.

## The degree bounds and a real refinement

In the gcd-six cases with all high multiplicities thirty, one has
\(d\ge30\), and \(6\mid d\) either from the residue theorem or from
the only possible single-root degree. For \(b=3\), the single-root
case \(d=26\) would have sheared degrees \((64,96)\), excluded by
the [planar degree-100 bound and exact degree-pair transfer](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FORTY_NINE_FROM_F7.md#two-inputs-from-polynomial-degree-pairs).
All other cases have \(d\ge30\) with \(6\mid d\).
When \(b=7\), the single-root case \(d=15\) is excluded by the
Bernstein cutoff through 59. The same cutoff requires
\(2(d+6)\ge60\), giving \(d\ge24\). In the multiple-root case
the residue theorem makes \(d\) a multiple of six. This proves
the order table.

In the \((12,18)\) boundary profile \((26,1,1,1,1)\), the
[ramification bound for the residue equation](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#a-rational-map-from-the-residue-equation)
has \(M=5,\Delta=5,N=4\), with a unique pole at the high root.
In the reciprocal coordinate \(z=1/(x-x_0)\), the map \(S\)
is a polynomial of degree four and is totally ramified at \(z=0\).
Therefore \(S(z)=c_0+c_1z^4\) with \(c_0c_1\ne0\).
In the coordinate \(t=x-x_0\), the residue equation forces
\[
         h=A\,t^{26}(t^4-\lambda),\qquad
         B=\frac{3t(t^4-\lambda)}{10\lambda},
         \qquad A\lambda\ne0 .
\]
These polynomials satisfy the residue equation by direct substitution.

**Corollary PolynomialWeylOrderTwelveFortyTwoRootField.** Over a
characteristic-zero field \(k\), each multiplicity-fifteen root \(x_0\)
in the reduced \((12,42)\) branch satisfies
\(\sqrt{-15}\in k(x_0)\).
If \(-15\) is nonsquare in \(k\), the common factor degree is
at least thirty. At degree thirty the common factor is a scalar
multiple of the fifteenth power of an irreducible quadratic with
splitting field \(k(\sqrt{-15})\).

At such a root the \((6,15)\) crossing is defined over \(k(x_0)\).
Make both powers monic at the upper endpoint by multiplying by
scalars in that field. The quotient \(R^7/(R^2)^3\) gives the
common factor \(R=XL(XY^3)\) over \(k(x_0)\).
In the [local two-root calculation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#the-two-roots-and-real-coefficients),
\(L\) has roots of multiplicities four and one, whose ratio \(c\)
satisfies \(2c^2-5c+5=0\). Both roots belong to \(k(x_0)\) because
their multiplicities are distinct, so \((4c-5)^2=-15\) in that
field. If \(-15\) is nonsquare in \(k\), the root \(x_0\) has
even degree over \(k\). Its conjugates have the same multiplicity,
giving the degree bound and the quadratic assertion at equality.

**Corollary PolynomialWeylOrderTwelveRealFortyTwoResidual.** Over
\(\mathbb R\), the reduced \((12,42)\) branch requires
\(d\equiv6\pmod{12}\), with \(d\ge30\), and at most one
real root of \(h\), whose multiplicity is even if it occurs.

All real roots of \(h\) have multiplicity below six by the
root-field condition. The
[real residue obstruction](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#real-roots-and-parity)
then gives the root count and makes \(d/6\) odd.
