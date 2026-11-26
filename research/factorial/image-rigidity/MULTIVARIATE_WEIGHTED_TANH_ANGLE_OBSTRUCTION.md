# A shared angular root for positive weights on matched tanh blocks

The homogeneous matched tanh moment and kernel transport inequalities
remain valid with their actual row caps when each complete row block
has a positive real weight.
They do not imply angular separation in that weighted class,
even with exact sampling and exponentiation of the linked
factorial symbols. The following example has simple unit-circle
roots in both degrees and satisfies the midpoint constraint
\(\delta\le S_*\). The unweighted sharp packing constraint
distinguishes it from the ordinary profile class.

## The weighted example

Take the mirrored rows \((1,31),(31,1)\), each with weight
\(\mu=\kappa/30>0\). Use the normalized polynomials
\[
 F_j(z;\kappa)=\sum_{k=0}^j c_{j,k}(\kappa)z^k,\qquad
 c_{j,k}(\kappa)=\binom jk
       \left(\frac{(j+30k)!(31j-30k)!}{j!(31j)!}\right)^{\kappa/30}.
\]
These are the exact row-weighted coefficients. Both endpoint
coefficients equal one. The exponent pair will remain fixed
as \(\kappa\) varies over the interval
\[
                         I=[1021/2000,\ 2553/5000].
\]
There is a \(\kappa\) in the interior of \(I\) at which \(F_8\)
and \(F_9\) have a common conjugate pair on the unit circle.
All the roots of each polynomial are simple on the entire
interval \(I\).

The weights multiply each lower and upper row contribution by
the same positive scalar. The individual Gauss grid points
and the exact adjacent map for their tanh values are unchanged.
The [matched moment comparison](MULTIVARIATE_MATCHED_TANH_MOMENTS.md#exact-tanh-transport-and-moment-errors-at-every-order),
the [integer-grid excess-loss inequality](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md),
and the [full kernel comparison](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md)
follow by summing their rowwise bounds with weight \(\mu\).
The Jensen and pointwise cap bounds for the weighted moments
also hold by taking probability weights proportional to the
weighted tanh values. The signed quantities cancel identically
by reflection. The positive Laplace curvature kernel and
the [midpoint curvature hierarchy](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md)
are compatible with positive row weights. In particular the
weighted midpoint curvature sum is
\[
          \delta=30/16,\qquad S_*=2\mu(30^2/16)=2\kappa\delta,
\]
so \(\delta<S_*\) throughout \(I\).

These polynomials also arise from genuine positive random
variables. In [Berg, equation (10) and Lemma 2.1](https://web.math.ku.dk/~berg/manus/potens1.pdf)
the measure \(e_\mu\) satisfies
\(\mathbb E X^a=(a!)^\mu\) for nonnegative integers \(a\).
Take independent \(X,Y\) of that law and the monomials
\(U=XY^{31}\), \(V=X^{31}Y\).
The polynomial \(\mathbb E[(U+zV)^j]\) is \(F_j\) multiplied
by its positive constant coefficient.

## An exact resultant sign change

Put \(x=z+z^{-1}\), abbreviate \(c_k=c_{8,k}\), \(d_k=c_{9,k}\),
and use reciprocity to write
\[
 F_8(z;\kappa)=z^4 E(x),\qquad
 F_9(z;\kappa)=z^4(1+z)O(x).
\]
The two monic quartics are
\[
\begin{aligned}
 E(x)&=x^4+c_7 x^3+(c_6-4)x^2+(c_5-3c_7)x+2+c_4-2c_6,\\
 O(x)&=x^4+(d_8-1)x^3+(d_7-d_8-3)x^2
         +(d_6-d_7-2d_8+2)x+d_5-d_6-d_7+d_8+1.
\end{aligned}
\]
For instance, \(z^j+z^{-j}\) generates the even expression.
For the odd quotient the analogous polynomials are
\((z^{j+1}+z^{-j})/(1+z)\), with the same Chebyshev recurrence
and initial terms \(1,x-1\).

The Sylvester determinant \(\operatorname{Res}_x(E,O)\) is
in \((-15,-14)\cdot10^{-6}\) at the left endpoint of \(I\)
and \((4,5)\cdot10^{-6}\) at the right endpoint.
The interval sign checks below isolate one root of each
quartic in each of its four intervals throughout \(I\).
All eight intervals are in \((-2,2)\), and the intervals
for the two polynomials overlap only at the fourth pair.
Continuity gives a zero of the resultant in the interior of
\(I\). The shared \(x\) lies in the fourth interval and gives
the two simple unit roots solving \(z^2-xz+1=0\).
The odd polynomial has the additional simple root \(z=-1\),
since \(O(-2)\ne0\).

The checks use rigorous Arb ball arithmetic with exact
rational input and integer factorials. The endpoint tests
bound the determinant as stated, while the root isolation
tests bound the signs over the whole parameter interval.

```python
from math import comb, factorial

from flint import arb, arb_mat, ctx, fmpq

ctx.prec = 192
lower = fmpq(1021, 2000)
upper = fmpq(2553, 5000)


def coefficient(order, index, kappa):
    moment = fmpq(
        factorial(order + 30 * index) * factorial(31 * order - 30 * index),
        factorial(order) * factorial(31 * order),
    )
    return comb(order, index) * (arb(moment).log() * kappa / 30).exp()


def quartics(kappa):
    c4, c5, c6, c7 = [coefficient(8, k, kappa) for k in range(4, 8)]
    d5, d6, d7, d8 = [coefficient(9, k, kappa) for k in range(5, 9)]
    return (
        [1, c7, c6 - 4, c5 - 3 * c7, 2 + c4 - 2 * c6],
        [1, d8 - 1, d7 - d8 - 3, d6 - d7 - 2 * d8 + 2,
         d5 - d6 - d7 + d8 + 1],
    )


def resultant(kappa):
    even, odd = quartics(kappa)
    matrix = arb_mat(8, 8)
    for row in range(4):
        for index in range(5):
            matrix[row, row + index] = even[index]
            matrix[row + 4, row + index] = odd[index]
    return matrix.det()


def evaluate(coefficients, x):
    value = arb(0)
    for coefficient in coefficients:
        value = value * x + coefficient
    return value


assert -15 < 10**6 * resultant(arb(lower)) < -14
assert 4 < 10**6 * resultant(arb(upper)) < 5
parameter = arb((lower + upper) / 2, (upper - lower) / 2)
isolations = (
    (("-39/20", "-9/5"), ("-6/5", "-1"), ("-1/20", "1/10"), ("9/20", "3/4")),
    (("-7/4", "-31/20"), ("-9/10", "-7/10"), ("1/5", "7/20"), ("9/20", "3/4")),
)
for polynomial, intervals in zip(quartics(parameter), isolations):
    for index, (left, right) in enumerate(intervals):
        sign = (-1) ** index
        assert sign * evaluate(polynomial, arb(fmpq(left))) > 0
        assert sign * evaluate(polynomial, arb(fmpq(right))) < 0
```

## A rational replay of the interval checks

Each table entry \(b\) brackets the corresponding coefficient between
\(b/10^{18}\) and \((b+1)/10^{18}\). The columns use the two endpoints
of \(I\).

| \(j\) | \(k\) | Left \(b\) | Right \(b\) |
|---:|---:|---:|---:|
| 8 | 4 | 6376052590192719150 | 6373060799072591520 |
| 8 | 5 | 5752847079127953382 | 5750283202891099896 |
| 8 | 6 | 4187169757267832575 | 4185611501934222172 |
| 8 | 7 | 2381340998134509137 | 2380775804476987382 |
| 9 | 5 | 8716165605402117626 | 8711606216210748070 |
| 9 | 6 | 7211306921756310822 | 7207839594995384261 |
| 9 | 7 | 4863636463616575692 | 4861729746503390535 |
| 9 | 8 | 2577146824346924326 | 2576515593241504952 |

For \(\kappa/30=p/q\), write the factorial ratio as \(A/B\) and
\(M=10^{18}\binom jk\). The two checks for each entry are integer
comparisons
\[
                    b^q B^p < M^q A^p < (b+1)^q B^p .
\]
Since \(A/B<1\), the right lower endpoint and left upper endpoint
bound the coefficient throughout \(I\).
The table verifies both determinant bounds by rational interval
Gaussian elimination on the same Sylvester matrix, pivoting on the first
available interval that excludes zero and tracking row swaps. Horner
evaluation of the quartics with the full-interval bounds verifies all
the root endpoint signs. This replay uses only integer powers and
rational interval operations.

## Where the ordinary row count enters

The obstruction concerns the homogeneous matched inequalities.
It even retains the exact same row weighting at both degrees.
A proof can use further integer-row constraints.

One such constraint is the
[sharp unweighted packing bound](MULTIVARIATE_SHARP_SKEW_CONDITIONED_TANH_MOMENTS.md).
In the example the \(30\) tanh values in either orientation
carry total mass \(\kappa<1\). If \(V_1^{(+)}\) and \(V_2^{(+)}\)
are their weighted first and second moments, Cauchy–Schwarz gives
\[
                    V_2^{(+)}\ge (V_1^{(+)})^2/\kappa
                                      >(V_1^{(+)})^2.
\]
The sharp packed-list bound with zero skew and the same total
first-moment budget gives \(V_2^{(+)}\le (V_1^{(+)})^2\) for
an unweighted list. Indeed the first-moment budget of each
orientation in its extremal list is \(V_1^{(+)}\), which is
smaller than the maximal tanh value since \(\kappa<1\).
Thus the weighted collision respects the homogeneous transport
but fails a discrete packing restriction. It also distinguishes the
ordinary curvature-only cap \(T_*\le(S+1)/(S+3)\) from the actual
row caps: the weighted \(S=60\kappa/31<1\), while the factor
at \(j=8,h=1\) is \(T=35/38\). That curvature-only cap uses the
ordinary row count as well. The ordinary integer-row angular
conjecture remains open.

