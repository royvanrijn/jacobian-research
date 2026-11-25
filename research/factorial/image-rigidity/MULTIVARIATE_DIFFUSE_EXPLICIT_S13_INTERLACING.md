# All-degree Curie–Weiss interlacing for \(s\ge13\)

Let \(B_j(z;s)=\sum_{k=0}^j\binom jk e^{-s k(j-k)/(2j)}z^k\).
For every \(s\ge13\) and \(j\ge1\), the roots of \(B_j\) and
\(B_{j+1}\) strictly interlace in argument. This is the
slowest physical upper time in the
[Gaussian reference comparison](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md).
Outward monotonicity and the reference inner inequalities give
strict interlacing throughout the physical time interval
\(\tau'/\tau\in[j/(j+1),1]\) whenever \(2j\tau\ge13\).

For \(j\ge22\), we subdivide the far horizontal contour in the
[analytic phase comparison for \(s\ge24\)](MULTIVARIATE_DIFFUSE_EXPLICIT_S24_INTERLACING.md).
The estimates use sixteen rational endpoint bounds and hold
uniformly at every higher degree. For \(2\le j\le21\), the
[exact diffuse resultants](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#exact-exclusions-through-terminal-22-in-the-limiting-model)
exclude shared roots at every \(s>0\). The roots lie on the
unit circle and are simple. Their order thus follows by
continuity from \(s\to\infty\), where \(B_j\to1+z^j\) and
the alternating order is strict. At \(j=1\) the lower root is
\(-1\) and the two roots of \(B_2\) are a nonreal conjugate
pair on the circle.

## Bounds for \(j\ge22\)

Put \(a=s/2\), \(a_0=13/2\) and \(\varepsilon=9/200\).
Use the same half-Gaussian integral \(I_j\), phase \(F\),
stationary value \(F_*\), and normalized Gaussian \(g_0\)
as in the \(s\ge24\) proof. On \(|z|\le2\) the stationary
equation has one simple solution \(|x_*|<\varepsilon\).
Indeed \(q_0=2e^{\varepsilon-a_0}\) gives
\[
 \frac{2a_0q_0}{1-q_0}<\varepsilon,\qquad
 \frac{\varepsilon^2}{4a_0}+\frac{q_0}{1-q_0}<\frac1{250}.
\]
Monotonicity of \(ae^{-a}\) extends the stationary-disk bound
to all \(a\ge a_0\). Rouché and the analytic implicit function
theorem give the same holomorphic saddle construction, and the
second inequality shows \(|F_*|<1/250\).

On the horizontal through \(x_*\) with \(\Re x\le a/5\),
the curvature perturbation from \(-1/(2a)\) is bounded by
\(q/(1-q)^2\), where \(q=2e^{-4a/5}\).
Thus the relative quadratic error constant is
\(\delta=2aq/(1-q)^2\le3/20\), by checking at \(a=a_0\)
and using the decrease of \(ae^{-4a/5}\).
The normalized near-horizontal error after matching the Gaussian
is at most \(\delta/[2(1-\delta)^{3/2}]<12/125\).

For a far subinterval \(\alpha\le\Re x/a\le\gamma\), put
\(t=x-x_*\in\mathbb R\) and \(g(t)=e^{-F_*}G(x_*+t,z)\).
In the exact ratio \(g/g_0=L(t)(1+ze^{x_*+t-a})\) of the
\(s\ge24\) proof, use
\[
\begin{split}
 \ell &=\frac{\exp(\varepsilon(\gamma+\varepsilon/a_0)/2)-1+q_0}{1-q_0},
 &|L-1|&\le\ell,\\
 G_0 &=\exp(-a_0\alpha^2/4+\varepsilon\alpha/2),\\
 D&=\ell G_0+2(1+\ell)\exp[-a_0(1-\gamma/2)^2+\varepsilon\gamma/2],
 &R&=G_0+D .
\end{split}
\]
In fact \(|t|\le\gamma a+\varepsilon\).
The normalized \(g_0\) is at most
\(G_0e^{-(a-a_0)\alpha^2/4}\).
The exponent \(\Re x-a-t^2/(4a)\) increases up to
\(\Re x=\gamma a\), where it is at most
\(-a(1-\gamma/2)^2+\varepsilon\gamma/2\).
Since \((1-\gamma/2)^2\ge\alpha^2/4\), the same factor
\(e^{-(a-a_0)\alpha^2/4}\) multiplies the upper bounds \(D,R\)
on \(|g-g_0|\) and \(\max(|g|,|g_0|)\).

Take \(\alpha=i/20,\ \gamma=(i+1)/20,\ 4\le i\le19\).
The difference of powers bounds the cost of each subinterval by
\[
 \frac{\sqrt a}{40\sqrt\pi}\,
        D j^{3/2}R^{j-1}e^{-j(a-a_0)\alpha^2/4}
 \le \frac{37}{1000} D j^{3/2}R^{j-1}.
\]
For the last inequality the \(a\)-dependent part decreases for
\(a\ge a_0\), since \(j\alpha^2/4\ge22/100>1/(2a_0)\).
We also used \(\sqrt{a_0}/(2\sqrt\pi)<37/50\).
The table gives outward rational \(D_i,R_i\) and \(M_i\), where
\(M_i\) bounds \(j^{3/2}R_i^{j-1}\) for \(j\ge22\).

| \(i\) | \(10^6D_i\) | \(10^6R_i\) | \(10^3M_i\) |
| ---: | ---: | ---: | ---: |
| 4 | 22443 | 963736 | 59910 |
| 5 | 27751 | 936272 | 25909 |
| 6 | 34212 | 904005 | 12395 |
| 7 | 42115 | 868093 | 5291 |
| 8 | 51801 | 829824 | 2053 |
| 9 | 63668 | 790589 | 743 |
| 10 | 78166 | 751846 | 259 |
| 11 | 95799 | 715085 | 91 |
| 12 | 117117 | 681795 | 34 |
| 13 | 142708 | 653426 | 14 |
| 14 | 173182 | 631359 | 7 |
| 15 | 209159 | 616874 | 5 |
| 16 | 251245 | 611120 | 4 |
| 17 | 300009 | 615086 | 4 |
| 18 | 355951 | 629575 | 7 |
| 19 | 419477 | 655178 | 15 |

In each row \(R_i^2(61/60)^3<1\), so it suffices to compare
\(j^3 R_i^{2(j-1)}\) to \(M_i^2\) for \(22\le j\le60\).
Summing the far costs gives
\((37/1000)\sum_i D_i M_i=107754776657/10^{12}<11/100\).

The Gaussian completion begins at
\(b=(a-\Re x_*)\sqrt{j/(4a)}>5\), since the lower bound
\((a-\varepsilon)\sqrt{j/(4a)}\) increases in \(a,j\) and already
exceeds \(5\) at \((a_0,22)\). Thus its cost is at most
\(e^{-25}/(10\sqrt\pi)<10^{-6}\).
On the connector \(\Re x=a\) and
\(\left|e^{-F_*}G(x,z)\right|\le
3\exp(1/250+\varepsilon^2/(4a_0)-a_0/4)<3/5\).
Its cost is at most \(\varepsilon\sqrt j(3/5)^j<10^{-5}\),
checked by squaring at \(j=22\) and using the decreasing ratio.

It follows that \(I_j=e^{jF_*}(1+\eta_j)\) with
\(|\eta_j|<21/100\) on \(|z|\le2\).
For \(r=21/100\) put \(C=r+r^2/(2(1-r))\).
The logarithmic series bounds \(|\log(1+\eta_j)|\) by \(C\).
Cauchy's estimate now gives
\[
 |H_j'|\le j/250+C<j/4,\qquad
 |H_{j+1}'-H_j'|\le1/250+2C<1/2,\quad |z|\le1.
\]
These are the same two strict inequalities required by the
[alternating phase comparison](MULTIVARIATE_DIFFUSE_EXPLICIT_S24_INTERLACING.md#strict-alternating-phases).
They prove the claim for every \(j\ge22\).

## Exact replay of the endpoint and subinterval bounds

For any rational \(x\), use the Taylor sum \(S=\sum_{k=0}^{36}|x|^k/k!\).
When \(x\le0\), the reciprocal \(1/S\) bounds \(e^x\) above.
Otherwise \(S+|x|^{37}/(37!(1-|x|/38))\) bounds it above.
All required \(x\) satisfy \(|x|<38\). The following standard Python
uses those rational bounds to recover the table.

~~~python
from fractions import Fraction as Q
from math import factorial, isqrt

def exp_up(x):
    y = abs(x)
    terms = sum(y**k / factorial(k) for k in range(37))
    if x <= 0:
        return 1 / terms
    return terms + y**37 / (factorial(37) * (1 - y/38))

def ceil_scaled(x, m):
    y = m * x
    return (y.numerator + y.denominator - 1) // y.denominator

a = Q(13, 2)
e = Q(9, 200)
q = 2 * exp_up(e-a)
assert 2*a*q/(1-q) < e and e*e/(4*a)+q/(1-q) < Q(1,250)
q_near = 2 * exp_up(-4*a/5)
assert 2*a*q_near/(1-q_near)**2 < Q(3,20)
assert (Q(3,20)/(2*Q(12,125)))**2 < Q(17,20)**3
assert ((a-e)**2)*22/(4*a) > 25
assert 3*exp_up(Q(1,250)+e*e/(4*a)-a/4) < Q(3,5)

cost = Q(0)
for i in range(4, 20):
    lo, hi = Q(i,20), Q(i+1,20)
    ell = (exp_up(e*(hi+e/a)/2)-1+q)/(1-q)
    g0 = exp_up(-a*lo*lo/4+e*lo/2)
    d = ell*g0 + 2*(1+ell)*exp_up(-a*(1-hi/2)**2+e*hi/2)
    di, ri = ceil_scaled(d, 10**6), ceil_scaled(g0+d, 10**6)
    D, R = Q(di, 10**6), Q(ri, 10**6)
    assert R**2 * Q(61,60)**3 < 1
    m2 = 10**6 * max(j**3 * R**(2*(j-1)) for j in range(22,61))
    mi = isqrt(m2.numerator // m2.denominator) + 1
    cost += Q(37,1000) * D * Q(mi,1000)
    print(i, di, ri, mi)
assert cost == Q(107754776657,10**12) < Q(11,100)
~~~

