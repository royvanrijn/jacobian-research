# A weighted angular collision satisfying quadratic tanh packing

Positive row weights admit an angular collision at orders eight and nine
even if the homogeneous matched tanh inequalities are supplemented by the sharp
capped-list bound for the second moment at zero skew. The example uses
one fractionally weighted mirrored row pair and a tail of ordinary
mirrored unit rows. The same weights are used at both orders.

## The two linked polynomials

Put \(H=250\), \(\ell=199\), and \(m=95\). Give each of
\((1,1+H),(1+H,1)\) weight \(\kappa/H>0\), and add \(m\) copies
of each of \((\ell,\ell+1),(\ell+1,\ell)\) with weight one.
The coefficient of \(z^k\) in the endpoint-normalized moment of degree
\(j\) is
\[
 c_{j,k}(\kappa)=\binom jk
 \left(\frac{(j+Hk)!(j(1+H)-Hk)!}{j!(j(1+H))!}\right)^{\kappa/H}
 \left(\frac{(j\ell+k)!(j(\ell+1)-k)!}{(j\ell)!(j(\ell+1))!}\right)^m.
\]
These polynomials are reciprocal. The matched Gauss factors and
their cross-order map are exact. The homogeneous rowwise moment
and kernel estimates apply by summing with the positive weights.
[Matched tanh transport](MULTIVARIATE_MATCHED_TANH_MOMENTS.md)
and [full kernel transport](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md)

For one orientation let \(C_j\) and \(B_j\) be its weighted first
and second tanh moments and \(A_j\) its maximal tanh factor. The
other orientation has the same values, so the signed skew is zero.
For every \(\kappa\) in
\[
                 I=[9734800,9734801]/10^8
\]
and \(j=8,9\), one has \(C_j<A_j\) and \(C_j^2-B_j>1/100\).
Thus the example passes the [sharp quadratic packing
test](MULTIVARIATE_SHARP_SKEW_CONDITIONED_TANH_MOMENTS.md)
with the actual cap \(A_j\) and combined first-moment budget \(2C_j\).
In this range the packed list has one entry \(C_j\) in each orientation.
The weighted cubic moment exceeds \(C_j^3\) by more than \(1/25\).
The example therefore leaves open what follows from the higher-power
packing constraints together with full row integrality.

## Separation fails inside the parameter interval

Set \(x=z+z^{-1}\). As in the
[weighted two-row example](MULTIVARIATE_WEIGHTED_TANH_ANGLE_OBSTRUCTION.md#an-exact-resultant-sign-change),
division by \(z^4\) in degree eight and by \(z^4(1+z)\) in degree nine
gives monic quartics \(E(x),O(x)\). The resultant
\(\operatorname{Res}_x(E,O)\) is strictly negative at the left
endpoint of \(I\) and strictly positive at the right endpoint.
The sign brackets below put every root of each quartic in \((-2,2)\)
for the whole interval. The intermediate value theorem provides a
common root of the two quartics at an interior parameter. It gives a
simple common conjugate pair of unit-circle roots of the original
reciprocal polynomials. In particular, the example fails even when the
two arguments are compared at the same radius and common endpoint scale.

The checks below use rational interval endpoints, integer factorials,
and rigorous Arb balls for logarithms, exponentials, and the determinant.
The packing checks use only rational arithmetic.

```python
from fractions import Fraction
from math import comb, factorial

from flint import arb, arb_mat, ctx, fmpq

ctx.prec = 256
gap, base, count = 250, 199, 95
left, right = Fraction(9734800, 10**8), Fraction(9734801, 10**8)


def coefficient(j, k, parameter):
    hard = fmpq(
        factorial(j + gap*k) * factorial(j*(1 + gap) - gap*k),
        factorial(j) * factorial(j*(1 + gap)),
    )
    tail = fmpq(
        factorial(j*base + k) * factorial(j*(base + 1) - k),
        factorial(j*base) * factorial(j*(base + 1)),
    )
    correction = parameter * arb(hard).log() / gap + count * arb(tail).log()
    return comb(j, k) * correction.exp()


def quartics(parameter):
    c4, c5, c6, c7 = [coefficient(8, k, parameter) for k in range(4, 8)]
    d5, d6, d7, d8 = [coefficient(9, k, parameter) for k in range(5, 9)]
    return (
        [1, c7, c6 - 4, c5 - 3*c7, 2 + c4 - 2*c6],
        [1, d8 - 1, d7 - d8 - 3, d6 - d7 - 2*d8 + 2,
         d5 - d6 - d7 + d8 + 1],
    )


def resultant(parameter):
    matrix = arb_mat(8, 8)
    for block, polynomial in enumerate(quartics(parameter)):
        for row in range(4):
            for col, coefficient_value in enumerate(polynomial):
                matrix[4*block + row, row + col] = coefficient_value
    return matrix.det()


def evaluate(polynomial, x):
    value = arb(0)
    for coefficient_value in polynomial:
        value = value*x + coefficient_value
    return value


a = fmpq(left.numerator, left.denominator)
b = fmpq(right.numerator, right.denominator)
assert -2 < 10**13 * resultant(arb(a)) < -1
assert 9 < 10**12 * resultant(arb(b)) < 10
parameter = arb((a + b)/2, (b - a)/2)
isolations = (
    ((-193, -190), (-133, -129), (-76, -73), (-3, 0)),
    ((-172, -169), (-102, -99), (-76, -73), (3, 7)),
)
for polynomial, brackets in zip(quartics(parameter), isolations):
    for index, (u, v) in enumerate(brackets):
        sign = (-1)**index
        assert sign*evaluate(polynomial, arb(fmpq(u, 100))) > 0
        assert sign*evaluate(polynomial, arb(fmpq(v, 100))) < 0


def moments(j, p, parameter):
    numerator = (j - 1)*gap
    hard = [Fraction(numerator, numerator + 2*j + 2*h)
            for h in range(1, gap + 1)]
    tail = Fraction(j - 1, j*(2*base + 1) + 1)
    return parameter/gap * sum(t**p for t in hard) + count*tail**p


for j in (8, 9):
    cap = Fraction((j - 1)*gap, (j - 1)*gap + 2*j + 2)
    assert Fraction(j - 1, j*(2*base + 1) + 1) < cap
    assert moments(j, 1, right) < cap
    assert moments(j, 1, left)**2 - moments(j, 2, right) > Fraction(1, 100)
    assert moments(j, 3, left) - moments(j, 1, right)**3 > Fraction(1, 25)
```
