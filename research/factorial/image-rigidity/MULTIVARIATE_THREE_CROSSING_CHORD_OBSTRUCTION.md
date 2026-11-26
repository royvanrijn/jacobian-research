# An ordinary obstruction to the unrestricted three-crossing chord test

The [two-crossing comparison](MULTIVARIATE_TWO_SECTOR_COEFFICIENT_CHORD.md)
subtracts a linear multiple of the lower moment polynomial and
uses a signed chord of its actual coefficients. A natural attempt
at a third crossing is to subtract a quadratic multiple and
interpolate at three sine zeros. One ordinary fully shared profile
shows that the corresponding mixed minors have both signs.
For the same profile, they also have both signs when the
interpolation points are three active sine cuts at a prescribed
angle. These are failures of the proposed sign criteria.

## The mixed determinants

Take sixteen copies of the mirrored row pair
\((20,21),(21,20)\). Positive constant normalization gives
the reciprocal coefficients
\[
 C_{m,k}=\binom mk
       \left(\frac{(20m+k)!(21m-k)!}{(20m)!(21m)!}\right)^{16}.
\]
The profile has signed gap sum zero, \(V=8/5\),
\(\eta=1/20\) and \(J=4/105\), so it meets all three
hypotheses of the two-crossing theorem.

For adjacent degrees \(n,N=n+1\), put \(A_k=C_{n,k}\),
\(B_k=C_{N,k}\), and extend each array by zero outside its
support. Writing
\[
              R(k)=(A_k,A_{k-1},A_{k-2},B_k)
\]
amounts to choosing a quadratic multiplier of \(A(z)\) and
comparing its coefficients with \(B(z)\).
For \(n=8\), every three-row determinant of the first three
columns, among rows \(0,\ldots,8\) in increasing order, is
positive. But the mixed determinants of the four columns
have the following signs:
\[
 \begin{array}{c|cc}
  \text{row indices} &(0,1,2,3)&(0,1,5,8)\\ \hline
  \operatorname{sgn}\det &+&-
 \end{array}
\]
Dividing the rows by their positive \(A_k\) leaves all the
signs unchanged. Hence the good three-column interpolation
basis alone does not give a uniform sign for the next mixed
divided differences.

## Three active cuts with a padded tail

Use the same ordinary profile at \(n=4,N=5\), now keeping rows
through \(6\) so as to contain the quadratic multiple of \(A\).
Linearly interpolate \(R(k)\) between consecutive integer indices.
Take \(\theta=2\pi/3\). The three cuts within the support occur at
\[
                    x_1=3/2,\qquad x_2=3,\qquad x_3=9/2 .
\]
The ending index \(6\) has zero sine.
The matrix of the first three columns at the cuts is nonsingular
with positive determinant. Thus there is a unique \(c_0,c_1,c_2\)
such that the interpolated residual
\[
                B_k-c_0A_k-c_1A_{k-1}-c_2A_{k-2}
\]
vanishes at these cuts. For each of the nonzero sine indices,
insert \(R(k)\) into the three cut rows in increasing index order.
The determinant signs are
\[
 \begin{array}{c|rrrr}
    k&1&2&4&5\\ \hline
    \operatorname{sgn}\det &+&+&-&-
 \end{array}
\]
Indeed the determination of the interpolant expresses
each sorted determinant as the residual times the positive
anchor determinant times \((-1)^j\), where \(j\) is the
number of cuts after \(k\). Multiplying the residuals by
\(\sin(k\theta)\) gives two negative and two positive
terms. Consequently this interpolation does not yield
a termwise signed sine comparison even at all the active
cuts of one angle.

Here \(\theta\) is a comparison angle. It is not a root
argument of \(A\) at \(n=4\), as is seen from
\(C_{4,2}-C_{4,1}-1\ne0\) after the reciprocal projection
\(u=z+z^{-1}=-1\). The example leaves open an interpolation
criterion restricted to a lower root, or an estimate using
the aggregate size of the mixed signed terms.

## Exact arithmetic replay

~~~python
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, factorial, prod

def coefficients(m):
    denominator = factorial(20*m)*factorial(21*m)
    return [comb(m,k)*
            F(factorial(20*m+k)*factorial(21*m-k),denominator)**16
            for k in range(m+1)]

def matrix_row(A, B, k):
    def at(C, j):
        return C[j] if 0 <= j < len(C) else F(0)
    return [at(A,k-j) for j in range(3)] + [at(B,k)]

def det(M):
    n = len(M)
    return sum(
        (-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        * prod(M[i][p[i]] for i in range(n))
        for p in permutations(range(n))
    )

A,B = coefficients(8), coefficients(9)
assert all(det([matrix_row(A,B,k)[:3] for k in s]) > 0
           for s in combinations(range(9),3))
assert det([matrix_row(A,B,k) for k in (0,1,2,3)]) > 0
assert det([matrix_row(A,B,k) for k in (0,1,5,8)]) < 0

A,B = coefficients(4), coefficients(5)
assert A[2]-A[1]-1 != 0

def interpolate(x):
    k = x.numerator // x.denominator
    return [(k+1-x)*a+(x-k)*b for a,b in
            zip(matrix_row(A,B,k),matrix_row(A,B,k+1))]

cuts = [F(3,2),F(3),F(9,2)]
assert det([interpolate(x)[:3] for x in cuts]) > 0
for k,sign in [(1,1),(2,1),(4,-1),(5,-1)]:
    assert sign*det([interpolate(x) for x in sorted(cuts+[F(k)])]) > 0
~~~
