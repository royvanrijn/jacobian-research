# Weyl generation through Bernstein degree and constant-leading order fifty-five

**Theorem PolynomialWeylBernsteinThroughFiftyFive.**
A characteristic-zero polynomial Weyl pair with nonzero scalar
commutator generates if it has a member of Bernstein degree at
most \(55\), or a member of nonzero constant differential leading
coefficient and order at most \(55\).
In a proper pair each member has differential order plus
leading-coefficient degree at least \(56\).

Start from the
[bound through fifty-three](POLYNOMIAL_WEYL_BERNSTEIN_THROUGH_FIFTY_THREE.md).
For a hypothetical proper member of Bernstein degree \(54\), axis
orientation, differential mate descent and the infinity shear
give a shaped pair of coprime nondividing ratio \(a:b\), with both
\(a,b\ge2\), outer corner \((g,d)\), and \(54=a(g+d)\),
\(g+d\ge16\). Its maximal \(X\)-exponent face is the single
outer corner in Fourier orientation.
Thus \((a,g+d)\) is either \((2,27)\) or \((3,18)\).

## First crossings

For outer sum twenty-seven, the earlier ratio \(2:b\) arguments
exclude all lower crossing sums except the potentially embedded
sum-twenty-one row. In particular a sum-twenty-five crossing
forces ratio \(2:3\) by the local successor test, and the
[nondividing-degree transfer](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md)
followed by [Moh's degree bound](https://doi.org/10.1515/crll.1983.340.140)
excludes the actual degrees \((54,81)\).
For its corner \(C_0=(6,15)\),
\(g\ge6\) and \(4g-27\le3\), so \(g\le7\).
As in the degree-fifty-two proof, \(g=6\) contradicts the pure
maximal \(X\)-face. With \(g=7,d=20\), the same positive-weight
argument from that proof gives the unique preceding edge
from \((7,20)\) to \((6,15)\), now of weight \((5,-1)\).
Its common factor is \(X^6Y^{15}L(XY^5)\) with
\(\deg L=1,\ L(0)\ne0\). The polynomial bracket mate of
[Guccione–Guccione–Valqui, Theorem 4.1(4)](sources/guccione-guccione-valqui-1111.6100v3.pdf)
has the form \(XYG(z)\) at \(z=XY^5\), and its equation is
\[
                  15zG'+(-9-4zL'/L)G=c\ne0 .
\]
The root of \(L\) forces \(\deg G>0\), while the coefficient
at infinity requires \(15\deg G=13\). This excludes the
embedded sum-twenty-one crossing.

For outer sum eighteen, the general crossing tests exclude
sums at most fifteen, and sum seventeen is prime.
The two sum-sixteen rows \((4,12),(6,10)\) have weights
\((4,-1),(2,-1)\), respectively.
If one of these is \(C_0=(u,v)\), the outer corner satisfies
\(g\ge u\) and \((\rho+1)g-18\le \rho u-v\).
Since \(u+v=16\) and \(\rho+1>2\), it follows that \(g=u\).
Then \(d=v+2\), contradicting the pure maximal \(X\)-face.

It remains to test normalized crossing sum \(18\) at ratio
\(3:b\), and sum \(27\) at ratio \(2:b\), each saturating the
outer sum. The integer corner tests of
Guccione–Guccione–Valqui, Proposition 7.3, give
\[
 \begin{array}{c|c|c|c|c}
 \text{sum}&C_0&C_1&(f_1,f_2)&(\rho,\sigma)\\ \hline
 18&(6,12)&(2,0)&(2,4)&(3,-1)\\
 27&(9,18)&(3,0)&(2,4)&(3,-1).
 \end{array}                                                  \tag{1}
\]
The exact replay is below. The common factors have the form
\(X^h L(XY^3)\), with \((h,\deg L)=(2,4),(3,6)\) respectively.
Their polynomial bracket mates have
\(\deg G=(2\deg L-h)/(3h)=1\) in the earlier
[bracket equation](POLYNOMIAL_WEYL_BERNSTEIN_DEGREE_THROUGH_FORTY_ONE.md#the-new-faces-and-their-successors).
As \(G\) vanishes at all roots of \(L\), each \(L\) has one
nonzero root, with multiplicity \(\mu=4,6\), respectively.

## Successor tests

Apply [Guccione–Guccione–Valqui II, Proposition 4.3](https://arxiv.org/pdf/1205.6827v1)
at the two crossings. Its hypotheses hold by the original first
crossing: the incoming edge has nonmonomial, positive-weight
forms, with nondividing ratio and zero leading bracket, a scalar
commutator, and negative diagonal differences of their highest
\(Y\)-endpoints. The first fractional shifts give the endpoints
\[
 \begin{array}{c|c|c|c|c}
 C_0&C&\mu&l&\beta=l(\mu-C_x)\\ \hline
 (6,12)&(10/3,4)&4&3&2\\
 (9,18)&(5,6)&6&3&3 .
 \end{array}
\]
Here \(l\) is the actual support denominator after the shift.
Both have previous slope \(1/3\).

We use the lattice and coefficient tests of the
[earlier successor analysis](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#testing-the-commuting-successors).
For clarity, the scalar exclusion in the present argument
does not require ratio \(2:b\). At a scalar successor of
ratio \(a:b\), write \(N=a+b\). Its slope satisfies
\[
                 ls=l-\frac{\beta N}{\mu N-1}.
\]
The reduced denominator \(D\) divides every difference of
\(Y\)-exponents within a face, and is coprime to \(\mu\).
When \(\beta\le\mu\), it satisfies
\(D\ge(\mu N-1)/\beta>N-1\).
It therefore divides neither \(a\mu\) nor \(b\mu\),
so neither face polynomial has a constant term.
Their bracket then cannot be scalar. Every endpoint
used below has \(\beta\le\mu\).

For a commuting multiple-root successor, write \(j\) for the
least exponent of the normalized polynomial \(f(\eta)\) of
degree \(\mu\). The tests use
\[
 ls=l-\frac{\beta M}{\mu(M-1)},\quad
              D\mid\mu-j,\quad M\equiv1\pmod D,\quad
              2\le M\le\beta\mu+1 ,
\]
where \(D\) is the reduced denominator of \(ls\) and \(M\) is
the degree in the coefficient equation
\(W'-(M/\mu)(f'/f)W=c\), \(c\ne0\).
As in the bound through fifty-one, an integral \(M/\mu\)
is impossible when there are multiple distinct roots:
\(W/f^{M/\mu}\) would be a rational primitive of
\(c/f^{M/\mu}\).
Also every distinct root of \(f\) is a simple zero of \(W\).
If the multiplicity of such a root is \(m_0\),
evaluation there gives \(c=(1-(M/\mu)m_0)W'\).
In particular the product \((M/\mu)m_0\) cannot be one.
A commuting single-root shift preserves \(C,l\), and
requires an integral \(ls\) and a strict increase in slope.

### Sum eighteen

For \((l,\mu,\beta)=(3,4,2)\), the multiple-root rows are
\(M=2,D=1,0\le j<4\) at slope \(2/3\), and
\(M=4,D=3,j=1\) at slope \(7/9\).
The latter has integral multiplier.
At \(M=2\) exactly two distinct roots are allowed.
Their multiplicities are \(3,1\), since \(2,2\)
would make the multiplier times a root multiplicity one.
Shift at the triple root to \(C'=(8/3,3)\),
\(l=3,\beta=1\), at previous slope \(2/3\).
The denominator \(3\) suffices for a nonzero root
shift at slope \(2/3\). If the triple root is zero,
the shift is the identity and the actual denominator
is still \(3\).
There is no single-root continuation and the only multiple-root
row has \(M=3,D=2,j=1\), an integral multiplier.
A single-root shift at the original data could only have
slope \(2/3\), which leaves the corner unchanged and only
the rejected \(M=4\) row at a larger slope.

### Sum twenty-seven

For \((l,\mu,\beta)=(3,6,3)\), the multiple-root rows are
\[
 \begin{array}{c|c|c|c}
 M&D&j&s\\ \hline
 2&1&0,\ldots,5&2/3\\
 4&3&0,3&7/9\\
 6&5&1&4/5 .
 \end{array}
\]
The last has integral multiplier.
For \(M=2\), the two distinct roots have multiplicities
\(5,1\) or \(4,2\), since \(3,3\) is forbidden
by the root-multiplicity test.
For \(M=4\), the two values of \(j\) give
\(f=(\eta^3-u)^2\) or \(f=\eta^3(\eta^3-u)\)
up to scale, with \(u\ne0\).
These four cases give the following successor data
after shifting at a maximal root. In the first two
rows the denominator \(3\) suffices whether the root
is zero or nonzero at slope \(2/3\).
In the third row the nonzero root shift at \(7/9\)
enlarges the denominator to \(9\).
In the last row the root is zero and the identity
shift preserves the actual denominator \(3\).

| Previous slope | \(C'\) | \(l'\) | \(\beta'\) | Multiple-root rows \((M,D,j)\) | Single-root slopes |
|---:|:---:|---:|---:|:---|:---|
| \(2/3\) | \((13/3,5)\) | 3 | 2 | \((5,2,1),(5,2,3)\) | none |
| \(2/3\) | \((11/3,4)\) | 3 | 1 | \((4,3,1)\) | none |
| \(7/9\) | \((17/9,2)\) | 9 | 1 | \((2,1,0),(2,1,1)\) | \(8/9\) |
| \(7/9\) | \((8/3,3)\) | 3 | 1 | \((3,2,1)\) | none |

All multiple-root rows have integral multiplier.
In the third row the only single-root shift has
the same slope as its rejected multiple-root rows,
so it leaves no larger continuation.
At the original \((3,6,3)\) data a single-root
shift could only have slope \(2/3\), preserving
the corner and leaving the already considered
\(M=4,6\) rows at larger slopes.
This excludes the sum-twenty-seven crossing.
Each application of Proposition 4.3 is either
the original crossing or a commuting successor
with the preserved edge, positive-weight,
nondividing-ratio and negative-diagonal conditions.

## The final degrees and orders

The two branches exclude Bernstein degree \(54\).
For \(55\), Bernstein mate descent would give a
nondividing pair of smaller degree \(55\),
since the lower degrees are excluded.
But no proper divisor of \(55\) is at least
sixteen, contradicting the degree-gcd bound
of Guccione–Guccione–Valqui, Corollary 7.4.
The constant-leading weighted descent then
works at orders \(54,55\) as in the proof
through fifty-three. The smaller-order cases
have order at most \(27\), the nonconstant
intermediate case contradicts the new individual
leading-degree bound, and the weight-one case
has Bernstein degree equal to its starting order.
The extension and descent to arbitrary
characteristic-zero fields are as in the
preceding bounds.

## Exact lattice replay

The [exact checker](verify_weyl_bernstein_through_fifty_five.py) enumerates the integer corners
and the seven successor lattices used above.
