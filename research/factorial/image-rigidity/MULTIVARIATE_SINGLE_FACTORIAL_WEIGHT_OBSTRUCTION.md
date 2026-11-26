# A single reciprocal factorial can destroy ratio-descent coprimality

Let \(T\) be the coefficient map
\[
                         T\left(\sum_{k=0}^m a_kz^k\right)
                                  =\sum_{k=0}^m\frac{a_k}{k!}z^k.
\]
There are strictly positive rational arrays in degrees six and
seven whose successive ratios strictly interlace, but whose
images under \(T\) share an irreducible quadratic factor.

This is relevant to one-private-side profiles. After reserving
one private unit factorial, the remaining combinatorial multiplier
is \(1/(m-k)!\) up to a scalar. Reversing the coefficient orders
gives the map \(T\) above. Reversal also preserves interlacing of
the unweighted ratios. Thus the bare ratio invariant of the
[endpoint descent](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
cannot justify insertion of even one reciprocal factorial.
The arrays here are abstract rational arrays. The test leaves
open adding the full linked Mellin and Gamma restrictions of
ordinary factorial arrays.

## Exact rational construction

Take
\[
 C(z)=z^2-\frac{19}{12}z+\frac{27}{2},
 \qquad
 U(z)=1+\frac34\sum_{k=1}^6
                     \left(1+\frac{k(k-1)}{1000}\right)
                     \frac{z^k}{k!}.
\]
Let \(J\) be the polynomial quotient of \(U\) by \(C\). Use
\(a_k=k![z^k](JC)\) as the degree-six unweighted array, and
write \(r_k=a_{k+1}/a_k\). All \(a_k>0\) and all the \(r_k\)
are strictly increasing, as verified below.

Start with \(b_0=1\) and positive successive ratios
\[
 \begin{aligned}
 q_0&=r_0, & q_1&=(99r_0+r_1)/100,\\
 q_k&=(r_{k-1}+r_k)/2 &&(2\le k\le5),\qquad q_6=r_5 .
 \end{aligned}
\]
This determines the degree-seven array \(b\).
Let \(B(z)=\sum_{k=0}^7 b_k z^k/k!\). In the quotient by \(C\),
write \(B(z)=s_0+s_1z\) and \(z^7/7!=t_0+t_1z\).
Now increase \(b_7\) by \(\delta_7=-s_1/t_1\) and \(b_0\)
by \(\delta_0=-s_0-\delta_7t_0\). The exact values satisfy
\(6<\delta_7<7\) and \(12<\delta_0<13\). These positive
changes decrease the first ratio, increase the last ratio,
and preserve all interior ratios. The final unweighted
arrays therefore satisfy \(q_k<r_k<q_{k+1}\) at every index.
By construction the final \(B\) and \(JC\) are divisible by
\(C\), which is irreducible over the rationals by its negative
discriminant. Everything is in rational arithmetic.

The following replay uses just polynomial division and
rational coefficient ratios.

~~~python
from flint import fmpq, fmpq_poly
from math import factorial

C = fmpq_poly([fmpq(27, 2), fmpq(-19, 12), 1])
U = fmpq_poly([1] + [
    fmpq(3, 4) * (1 + fmpq(k*(k-1), 1000)) / factorial(k)
    for k in range(1, 7)
])
J, _ = divmod(U, C)
A = J * C
a = [A[k]*factorial(k) for k in range(7)]
assert all(v > 0 for v in a)
r = [a[k+1]/a[k] for k in range(6)]
assert all(r[k]<r[k+1] for k in range(5))

q = [r[0], (99*r[0]+r[1])/100]
q += [(r[k-1]+r[k])/2 for k in range(2, 6)]
q += [r[5]]
b = [fmpq(1)]
for rate in q:
    b.append(b[-1]*rate)
B = fmpq_poly([v/factorial(k) for k,v in enumerate(b)])
s = B % C
t = fmpq_poly([0]*7+[fmpq(1, factorial(7))]) % C
delta7 = -s[1]/t[1]
delta0 = -s[0]-delta7*t[0]
assert 6<delta7<7 and 12<delta0<13
b[-1] += delta7
b[0] += delta0
rates = [b[k+1]/b[k] for k in range(7)]
assert all(rates[k]<r[k]<rates[k+1] for k in range(6))
B = fmpq_poly([v/factorial(k) for k,v in enumerate(b)])
assert A % C == B % C == 0
assert C[1]*C[1]-4*C[0]<0
~~~
