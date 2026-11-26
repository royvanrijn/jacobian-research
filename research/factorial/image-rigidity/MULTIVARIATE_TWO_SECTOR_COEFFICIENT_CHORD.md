# Two-sector separation from a linked factorial coefficient chord

For distinct exponent columns \(a,b\) with only shared unequal
rows, write \(d_i=b_i-a_i\), \(\ell_i=\min(a_i,b_i)\), and set
\[
 V=\sum_i d_i^2/\ell_i,\qquad \eta=\max_i |d_i|/\ell_i,\qquad
 J=\left|\sum_i d_i/a_i\right|,\qquad \Delta=\sum_i d_i.
\]
Let \(P_m(z)=\sum_{k=0}^m\binom mk F_{m,k}z^k\), where
\(F_{m,k}=\prod_i((m-k)a_i+kb_i)!\).

**Theorem.** Suppose
\[
                 V\le 8/5,\qquad \eta\le1/8,\qquad J\le V/8.
\]
For every \(n\ge2\), \(N=n+1\), the polynomials \(P_n(z)\) and
\(P_N((n/N)^\Delta z)\) have no common root at
\(|\arg z|\le3\pi/N\). For balanced columns \(\Delta=0\) this
is a physical comparison. No low-order scan or large-order
cutoff is needed. The exact coefficient test below applies
without caps on \(V,\eta,J\).

## The two-crossing and one-crossing tests

Let \(A,B\) be arbitrary positive coefficient arrays of degrees
\(n,N=n+1\). For \(0\le k\le n\) set
\[
 D_k=A_{k-1}/A_k,\quad C_k=B_k/A_k,\qquad A_{-1}=0,\qquad
 L=B_N/A_n .
\]
Suppose \(D_k\) strictly increases, the consecutive chord
slopes of \(C\) as a function of \(D\) strictly increase,
and the last chord slope is strictly smaller than \(L\).
Suppose also the \(C_k\) are strictly increasing. Then
the two polynomials have no common zero at
\(0<|\arg z|<\pi,\ N|\arg z|\le3\pi\).

The
[one-crossing test](MULTIVARIATE_DIRECT_RATIO_FIRST_SECTOR.md#one-sine-crossing-for-the-direct-coefficient-ratio)
covers \(0<N\theta\le2\pi\). For \(2\pi<N\theta\le3\pi\)
and \(0<\theta<\pi\), set \(x_1=\pi/\theta\) and
\(x_2=2\pi/\theta\). Interpolate \(D_k\) linearly in the index
\(k\), after adding \(D_N=D_n+1\). Form a piecewise linear convex
function of \(D\) through the actual points \((D_k,C_k)\),
then extend it to \(D_N\) with slope \(L\).
The values of this function at \(D(x_1),D(x_2)\) determine
a chord line \(c+dD\) with \(d\le L\).
The differences \(C_k-c-dD_k\) are nonnegative outside
\((x_1,x_2)\) and nonpositive inside. These signs agree
with \(\sin(k\theta)\) at all \(0\le k\le n\).
The inequality is strict at \(k=1\): \(x_1>1\),
\(x_2-x_1>1\), and there is a strict slope increase between
the chord points. Thus for \(z=r e^{i\theta}\),
\[
 \Im\{B(z)-(c+dz)A(z)\}
    =\sum_{k=0}^n A_k r^k(C_k-c-dD_k)\sin(k\theta)
              +(B_N-dA_n)r^N\sin(N\theta)>0.
\]
The last term is nonnegative and the \(k=1\) term is
positive. Conjugation treats negative \(\theta\).
For ordinary factorials a positive real root is impossible,
and a negative real root cannot be common since the even
moment is strictly positive there.

## Arithmetic and harmonic interpolation of the stripped ratios

For \(A=P_n,\ B=P_N(\kappa z)\), take
\(a_k=F_{n,k}\) and \(b_k=\kappa^k F_{N,k}\). Define
\[
 r_k=a_{k+1}/a_k,\quad t_k=b_{k+1}/b_k,\quad
 G_k=b_k/a_k,\quad H_k=b_{k+1}/a_k .
\]
The actual binomial factors give \(C_k=G_k+D_kH_{k-1}\)
at \(1\le k\le n\), \(C_0=G_0\), and
\(D_k=k/[(N-k)r_{k-1}]\) for \(k\ge1\).
In particular \(H_n=L\).

Suppose \(D_k\) strictly increases, \(t_0\le r_0\) and
\(r_{n-1}\le t_n\). For \(1\le k<n\) assume
\[
 \frac{(n-k)r_k+(k+1)r_{k-1}}{N}
       < t_k <
       \frac{N}{(N-k)/r_k+k/r_{k-1}}.                         \tag{1}
\]
Then the chord slope from index \(k\) to \(k+1\) is
strictly between \(H_{k-1}\) and \(H_k\).
Indeed, subtracting respectively \(H_{k-1}\) and \(H_k\)
times \(D_{k+1}-D_k\) from \(C_{k+1}-C_k\) gives
\[
  (G_{k+1}-G_k)+D_{k+1}(H_k-H_{k-1}),\qquad
  (G_{k+1}-G_k)+D_k(H_k-H_{k-1}).
\]
Divide by \(G_k\) and substitute
\(G_{k+1}/G_k=t_k/r_k,\ H_k/H_{k-1}=t_k/r_{k-1}\).
The resulting signs are exactly (1).
The first chord is at most \(H_0\) by \(G_1\le G_0\),
and the last is below \(H_{n-1}\le H_n=L\).
It remains only to check the first chord is positive to
ensure the monotonicity of the \(C_k\).

## Exact factorial verification in the stated range

Take \(\kappa=(n/N)^\Delta\). The first-coefficient
real-order estimate on each shared row gives \(t_0\le r_0\)
and, in the reverse direction, \(r_{n-1}\le t_n\).
In detail, a positive row gap \(d_i\) contributes
\(\sum_{h=1}^{d_i} a_i/(m a_i+h)\le d_i/m\)
to the real-order derivative of the log first ratio.
A negative gap contributes the negatives of the
fractions \(a_i/(m a_i-h)\) for \(0\le h<|d_i|\),
giving the same upper bound \(d_i/m\).
The compensating factor removes \(\Delta\log(N/n)\).
Likewise the bound \(t_0/r_0\ge e^{-\mathcal E}\) holds with
\(\mathcal E\le V/n^2\). To see this latter bound, compare
each signed term in the log first ratio directly with \(d_i/m\).
At real \(m\in[n,N]\), its deficit is at most
\[
 \sum_{h=1}^{|d_i|}\frac{h}{m^2\ell_i}\le d_i^2/(m^2\ell_i).
\]
For a decreasing row the sum over \(h\) starts at zero, and
the bound uses \(m a_i-h\ge m\ell_i\).
Thus \(C_1/C_0=(N/n)t_0/r_0>1\), using
\(V/n^2\le (8/5)/n^2<1/N<\log(N/n)\) for \(n\ge3\).
At \(n=2\), the stronger deficit
\(\mathcal E\le V\int_2^3 m^{-2}\,dm=V/6\) gives the
same result. This integral bound works at every \(n\ge2\).

Let \(\delta=\log r_k-\log r_{k-1}\) for \(1\le k<n\).
Grouping factors in a row by integer offsets, its contribution
to a log stripped ratio at order \(m\) is
\[
 d_i\log m+\sum_{\sigma\in\mathcal S_i}
        \operatorname{sgn}(d_i)\log(a_i+d_i(k+\sigma)/m),
\quad
 \mathcal S_i=
 \begin{cases}
 \{1/d_i,\ldots,1\}&d_i>0,\\
 \{0,1/|d_i|,\ldots,(|d_i|-1)/|d_i|\}&d_i<0 .
 \end{cases}
\]
Write \(f_i(s)=\operatorname{sgn}(d_i)\log(a_i+d_i s)\).
For each slot, the contribution to \(\delta\) is
\(\delta_{i,\sigma}=f_i((k+\sigma)/n)-f_i((k-1+\sigma)/n)>0\).
Since \(|d_i|/\max(a_i,b_i)\le f'_i\le|d_i|/\ell_i\),
\[
                     V/[n(1+\eta)]\le\delta\le V/n .
\]
Consequently \(D_k\) strictly increases. For \(k\ge1\) its
successive log increment is
\(\log((k+1)/k)+\log((N-k)/(n-k))-\delta\);
the two logs are each larger than \(1/(k+1)\) and \(1/(N-k)\)
respectively, whose sum is at least \(4/(n+2)\).
This exceeds \(V/n\) under the stated caps.

Put \(\omega=(k+1/2)/N\). The upper sampling point
\((k+\sigma)/N\) is the convex combination of the two lower
points with weight \((k+\sigma)/N\) on the earlier point.
For each slot the interpolation error is at most
\(\eta(1+\eta)\delta_{i,\sigma}/(8n)\), since
\(|f_i''|\le\eta f_i'\) and the maximum \(f_i'\) is at most
\(1+\eta\) times its minimum. To replace its weight by
\(\omega\), sum the signed offsets. The mean of \(\sigma-1/2\)
on a row is \(\operatorname{sgn}(d_i)/(2|d_i|)\).
The contribution from that mean differs from
\(d_i/(2n a_i)\) by at most \(\eta d_i^2/(2n\ell_i)\),
by comparing the reciprocal of \(a_i+d_i s\) with \(1/a_i\).
The centered positive deviations of the equally spaced offsets
have mean at most \(1/8\). The range of the slot differences
on each row is at most \(\eta\) times their mean. The
centered part therefore costs at most \(\eta\delta/8\)
after summing over the rows. It follows that
\[
 \left|\log t_k-(1-\omega)\log r_k-\omega\log r_{k-1}\right|
 \le \frac{\delta}{N}
 \left(\frac{(1+\eta)(J/V+\eta)}2+\frac{\eta}{8}
                  +\frac{N\eta(1+\eta)}{8n}\right).
 \tag{2}
\]

For two numbers whose log difference is \(\delta\),
the log of a weighted arithmetic mean deviates from the
weighted log mean by at most \(\delta^2/8\).
In fact the second derivative of the log exponential mean
is a variance of a random variable supported on two points
at distance \(\delta\), at most \(\delta^2/4\).
Apply the same bound to reciprocals for the harmonic mean.
The two weights on the earlier ratio in (1) are
\(\omega+1/(2N)\) and \(\omega-1/(2N)\). Therefore (2)
implies both strict inequalities in (1) if the parenthesis
in (2) plus \(N\delta/8\) is strictly below \(1/2\).
Under the stated caps the sum is at most
\[
         \frac9{64}+\frac1{64}+\frac{27}{1024}+\frac3{10}
                           =\frac{2473}{5120}<\frac12.
\]
The coefficient chord and the two-crossing test now give
the result at every order.

## Exact control

The following check uses a nonmirrored balanced profile made of
seven copies of \((9,10),(11,10)\). All inequalities are rational.
The checks illustrate the boundary and interior comparisons.

~~~python
from fractions import Fraction as Q

V = 7*(Q(1,9)+Q(1,10))
eta = Q(1,9)
J = 7*(Q(1,9)-Q(1,11))
assert V <= Q(8,5) and eta <= Q(1,8) and 0 < J <= V/8
bound = Q(9,64)+Q(1,64)+Q(27,1024)+Q(3,10)
assert bound == Q(2471,5120) < Q(1,2)

def ratios(m):
    return [Q(9*m+k+1,11*m-k)**7 for k in range(m)]

for n in (2, 13, 24):
    N = n+1
    r, t = ratios(n), ratios(N)
    assert t[0] <= r[0] and r[-1] <= t[-1] and N*t[0] > n*r[0]
    D = [Q(0)]+[Q(k,N-k)/r[k-1] for k in range(1,n+1)]
    assert all(D[k] < D[k+1] for k in range(n))
    for k in range(1,n):
        a = ((n-k)*r[k]+(k+1)*r[k-1])/N
        h = N/((N-k)/r[k]+k/r[k-1])
        assert a<t[k]<h
~~~
