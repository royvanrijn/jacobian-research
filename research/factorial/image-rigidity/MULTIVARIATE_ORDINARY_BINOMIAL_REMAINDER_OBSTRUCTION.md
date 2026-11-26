# Failure of binomially normalized remainder transfer for ordinary shared rows

The endpoint descent theorem applies to positive coefficient arrays
with linked interlaced ratios. For fully shared factorial moments,
the gamma arrays obtained by removing the binomial coefficients
still have that ratio order after signed degree compensation.
One possible extension would show that the ordinary endpoint
remainder, when divided by the binomial coefficients of its own
degree, inherits the gamma ratio order. This fails for an actual
balanced ordinary profile with just two shared rows.
[Endpoint descent theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)

Write the lower and upper ordinary arrays as
\(A_k=\binom nk a_k\) and \(B_k=\binom{n+1}k b_k\),
with \(a_0=b_0=1\). Put \(d=b_{n+1}/a_n\).
The endpoint remainder \(D=((1+dz)A-B)/z\) has degree at most
\(n-1\). Its binomially normalized coefficients are exactly
\[
 e_k=\frac{[z^k]D}{\binom{n-1}{k}}
       = n\left(\frac{a_{k+1}-b_{k+1}}{k+1}
                   +\frac{d\,a_k-b_{k+1}}{n-k}\right),
                   \qquad 0\le k\le n-1.                     \tag{1}
\]
The identity uses Pascal's rule and the endpoint normalization.
If \(q_k=b_{k+1}/b_k\le r_k=a_{k+1}/a_k\le q_{k+1}\),
both differences in (1) are nonnegative by telescoping from
the corresponding endpoints.

Take the two exponent rows \((11,10),(2000,2001)\), so the
signed total gap is zero. Take \(n=11\) and \(N=12\).
Scale the arguments of both polynomials by \(1/2000\) and
divide by their constant terms. This gives rational gamma ratios
\[
 r_k=\frac{22001+k}{2000(121-k)},\qquad
 q_k=\frac{24001+k}{2000(132-k)} .
\]
All the linked inequalities are strict at \(0\le k<11\),
and (1) has strictly positive coefficients through index ten.
The proposed next comparison would require
\(e_{10}/e_9\le r_{10}=7337/74000\).
Instead the exact excess is
\[
 \frac{e_{10}}{e_9}-r_{10}
 =\frac{57502869484256356600677210547845668393}
        {202205654427548900888929919713910570127974000}>0 .
\]
The following replay computes just the exhibited ratio chain
and the obstruction. It uses coefficient ratios and no
factorials or polynomial root computations.

~~~python
from fractions import Fraction as F

n = 11
r = [F(22001 + k, 2000*(121 - k)) for k in range(n)]
q = [F(24001 + k, 2000*(132 - k)) for k in range(n + 1)]
assert all(q[k] < r[k] < q[k+1] for k in range(n))

a, b = [F(1)], [F(1)]
for t in r:
    a.append(a[-1] * t)
for t in q:
    b.append(b[-1] * t)

d = b[-1] / a[-1]
e = [n*((a[k+1]-b[k+1])/(k+1) + (d*a[k]-b[k+1])/(n-k))
     for k in range(n)]
assert all(value > 0 for value in e)
assert r[10] == F(7337, 74000)
assert e[10]/e[9] - r[10] == F(
    57502869484256356600677210547845668393,
    202205654427548900888929919713910570127974000,
)
~~~

This calculation isolates failure of the proposed inheritance
of binomially normalized ratios in an ordinary profile. The
[relaxed binomial obstruction](MULTIVARIATE_BINOMIAL_RATIO_DESCENT_OBSTRUCTION.md)
gives a different failure, where positive abstract moment arrays
with both linked quotients can acquire an actual shared factor
after binomial weighting.
