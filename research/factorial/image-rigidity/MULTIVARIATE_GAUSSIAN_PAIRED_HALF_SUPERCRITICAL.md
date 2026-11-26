# All-degree sharp Gaussian interlacing for \(t\ge26/5\)

For every \(t\ge26/5\) and every \(n\ge1\), the consecutive Gaussian
polynomials
\[
 Q_j(z;t)=\sum_{k=0}^j\binom jk e^{-t k(1-k/j)}z^k,\qquad j=n,n+1,
\]
have strictly interlacing unit-circle roots. The same statement holds
when the lower raw time is \(t/n\) and the upper raw time is anywhere
from \(t/(n+1)\) to \(t/n\), by outward root motion and the
[inner-gap comparison](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md).
In the convention \(s=2t\) of the
[preceding uniform contour theorem](MULTIVARIATE_DIFFUSE_EXPLICIT_S13_INTERLACING.md),
the new bound is \(s\ge52/5\).

At \(n\ge22\) the proof controls the difference of successive
normalized half-integral errors. At \(n\le21\) it uses the existing
[all-coupling diffuse resultants through terminal 22](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#exact-exclusions-through-terminal-22-in-the-limiting-model).
There are no new tests in individual degrees. In these lower degrees
the existing coprimality, circle location, and the limit
\(Q_j(z;t)\to1+z^j\) as \(t\to\infty\) give the required ordering.
For \(n=1\) strict interlacing also follows from the two conjugate
roots of \(Q_2\) and the root \(-1\) of \(Q_1\).

## The analytic half integral

Write \(a=t\), \(a_0=26/5\), and take
\[
 \rho=51/40,\qquad \varepsilon=1/12,\qquad
 G(x,z)=e^{-x^2/(4a)}(1+ze^{x-a}),\qquad
 I_j(z)=\sqrt{j/(4\pi a)}\int_{-\infty}^{a}G(x,z)^j\,dx .
\]
The exact half decomposition is
\(Q_j(z;t)=I_j(z)+z^j I_j(1/z)\).
In the disk \(|z|\le\rho\), the saddle of
\(F(x,z)=-x^2/(4a)+\log(1+ze^{x-a})\)
lies in \(|x_*|<\varepsilon\).
Indeed the stationary equation
\[
 x=\frac{2az e^{x-a}}{1+ze^{x-a}}
\]
has its right side of modulus less than \(\varepsilon\) there.
For \(q_0=\rho e^{\varepsilon-a_0}\), the necessary bound is
\(2a_0q_0/(1-q_0)<\varepsilon\).
Rouché's theorem and monotonicity of \(ae^{-a}\) give a unique
simple saddle, analytic in \(z\). Its phase satisfies
\[
             |F_*|=|F(x_*,z)|<9/1000
\]
by \(\varepsilon^2/(4a_0)+q_0/(1-q_0)<9/1000\).

Move the half integral to the horizontal through \(x_*\) and
then close vertically at \(\Re x=a\).
Set \(u=x-x_*\in\mathbb R\) on the horizontal and
\[
 g(u,z)=e^{-F_*}G(x_*+u,z),\qquad g_0(u)=e^{-u^2/(4a)}.
\]
Let \(\eta_j=e^{-jF_*}I_j-1\).
The error is the difference from the full Gaussian integral,
which is one after the same normalization. We prove
\[
                      |\eta_j|<12/25,\qquad
                 |\eta_{j+1}-\eta_j|<11/100                       \tag{1}
\]
for \(j\ge22\), throughout \(|z|\le51/40\).

## The paired near-horizontal error

On \(\Re x\le 11a/40\), the exact logarithm gives
\[
 b=u^2/(4a),\qquad g=e^{-b+E},\qquad
 |E|\le\delta b,\qquad \delta=13/40 .
\]
Indeed
\(2a\rho e^{-29a/40}/(1-\rho e^{-29a/40})^2\le\delta\)
by monotonicity and the endpoint check at \(a_0\).
Integrate the second-derivative bound twice from the saddle.
The usual normalized Gaussian difference bound is
\[
        \frac{\delta}{2(1-\delta)^{3/2}}<295/1000 .
\]

The near segment also permits real order interpolation for
each fixed \(z\). For \(m\in[j,j+1]\), differentiate its
normalized integrand
\[
                  \sqrt{m/(4\pi a)}
                            \left(e^{-m(b-E)}-e^{-mb}\right).
\]
The differentiated absolute value is at most
\[
 \sqrt{m/(4\pi a)}\delta
           \left(\frac32 b+m b^2\right)e^{-m(1-\delta)b}.
\]
To see this, differentiate the prefactor and separate the
exponent derivative into \(-b\) times the original difference
and \(E e^{-m(b-E)}\).
Integrating the majorant over the entire horizontal line bounds
the paired error by
\[
     \frac{3\delta}{4j}
       \big((1-\delta)^{-3/2}+(1-\delta)^{-5/2}\big)
         <1/20 .
\]

## The paired far-horizontal error

Use fifty-eight intervals \(\alpha a\le\Re x\le\gamma a\) with
\(\alpha=k/80,\ \gamma=(k+1)/80,\ 22\le k<80\).
For \(e=a-a_0\ge0\), write
\(L=e^{-e\alpha^2/4}\).
The exact horizontal ratio is
\[
 g/g_0=\mathcal L(u)(1+ze^{x_*+u-a}),\qquad
 \mathcal L(u)=\frac{e^{-x_*u/(2a)}}{1+ze^{x_*-a}}.
\]
Define the endpoint bounds
\[
 \ell=\frac{\exp(\varepsilon(\gamma+\varepsilon/a_0)/2)-1+q_0}
                   {1-q_0},\quad
 G_0=e^{-a_0\alpha^2/4+\varepsilon\alpha/2},\quad
 D=\ell G_0+\rho(1+\ell)e^{-a_0(1-\gamma/2)^2+\varepsilon\gamma/2}.
\]
Then \(|\mathcal L-1|\le\ell\), \(g_0\le G_0 L\),
\(|g-g_0|\le D L\), and
\(\max\{|g|,g_0\}\le (G_0+D)L\).
For the term containing \(ze^{x-a}\), the real exponent
\(\Re x-a-(\Re x-\Re x_*)^2/(4a)\) increases up to
\(\Re x=\gamma a\), where it gives the displayed bound.
The factor \(L\) is valid since
\((1-\gamma/2)^2\ge\alpha^2/4\).

Round \(D\) and \(G_0+D\) outwards to \(D_i,R_i\), respectively.
The single-degree cost of this interval is at most
\[
 F_i=\frac{81}{10000} D_i M_i,\qquad
 M_i\ge\sup_{j\ge22}j^{3/2}R_i^{j-1}.
\]
This uses the difference of powers and the length \(a/80\).
The \(\sqrt a\,L^j\) ratio relative to \(a_0\) is bounded by
\[
 \sqrt{a/a_0}\,L^j\le e^{-A_i e},\qquad
 A_i=22\alpha^2/4-1/(2a_0)>0.
\]
Also \(\sqrt{a_0}/(160\sqrt\pi)<81/10000\).

To bound the adjacent difference, let
\[
 B_i\le e^{-a_0\gamma^2/4-\varepsilon\gamma/2-\varepsilon^2/(4a_0)},
 \qquad b_i=\gamma^2/4
\]
with \(B_i\) a positive rational lower bound.
The lower Gaussian estimate is \(g_0\ge B_i e^{-b_i e}\),
so \(|g-1|\le1+D_i-B_i e^{-b_i e}\).
The replay verifies
\((A_i+b_i)B_i\le A_i(1+D_i)\).
It follows by differentiating that
\(e^{-A_i e}(1+D_i-B_i e^{-b_i e})\) decreases on \(e\ge0\).
The paired far cost of the interval is therefore at most
\[
                (1+D_i-B_i+3/44)F_i .
\]
Indeed
\[
 (g^{j+1}-g_0^{j+1})-(g^j-g_0^j)
       =(g-1)(g^j-g_0^j)+(g-g_0)g_0^j .
\]
The first term uses the decreasing product bound just proved.
The second costs at most \(F_i/j\), since \(R_iL\le1\).
The change in the square-root prefactor costs at most \(F_i/(2j)\).
The rational sums verify
\[
         \sum_i F_i<182/1000,\qquad
         \sum_i (1+D_i-B_i+3/44)F_i<58/1000 .
\]

The Gaussian completion begins at
\((a-\Re x_*)\sqrt{j/(4a)}>5\), with cost \(<10^{-6}\).
On the vertical connector, the normalized modulus of \(G\) is
less than
\((1+\rho)e^{9/1000+\varepsilon^2/(4a_0)-a_0/4}<13/20\).
Its length is at most \(\varepsilon\), so its normalized cost
is less than \(1/25000\).
Twice these two individual bounds suffice for their paired costs.
Combining all pieces proves (1).

## Phase comparison

Put \(H_j=jF_*+\log(1+\eta_j)\) with the logarithm from its
power series. Equation (1) and integration along the straight
segment from \(\eta_j\) to \(\eta_{j+1}\) give
\[
 \|H_{j+1}-H_j\|_{|z|\le51/40}
    < 9/1000+(11/100)/(1-12/25).
\]
Schwarz--Pick on the disk of radius \(\rho\) bounds the derivative
of a holomorphic function with modulus at most \(M\) by
\(M\rho/(\rho^2-1)\) on the unit disk.
The factor here is \(2040/1001\), which gives
\(|H_{j+1}'-H_j'|<1/2\).
The individual bound obtained the same way is
\(|H_j'|<(2040/1001)(j\cdot9/1000+12/13)<j/4\).
Both \(H_j\) have conjugation symmetry.
The [strict alternating phase argument](MULTIVARIATE_DIFFUSE_EXPLICIT_S24_INTERLACING.md#strict-alternating-phases)
applied to the exact half decomposition now proves interlacing
at every \(j\ge22\).

## Exact replay of the analytic bounds

The script checks the fixed contour bounds on the entire
half-line \(a\ge26/5\) and all orders \(j\ge22\).
It does not check resultants or individual terminal cases.

~~~python
from fractions import Fraction as Q
from math import factorial, isqrt


def exp_up(x):
    y = abs(x)
    terms = sum(y**k / factorial(k) for k in range(37))
    if x <= 0:
        return 1 / terms
    return terms + y**37 / (factorial(37) * (1 - y / 38))


def round_up(x, scale):
    y = x * scale
    return Q(-(-y.numerator // y.denominator), scale)


a = Q(26, 5)
z_radius = Q(51, 40)
epsilon = Q(1, 12)
q0 = z_radius * exp_up(epsilon-a)
assert 2*a*q0/(1-q0) < epsilon
assert epsilon**2/(4*a)+q0/(1-q0) < Q(9,1000)

q_near = z_radius * exp_up(-Q(29,40)*a)
delta = Q(13,40)
assert 2*a*q_near/(1-q_near)**2 < delta
assert (delta/(2*Q(295,1000)))**2 < (1-delta)**3
assert (3*delta*(2-delta)/(4*22*Q(1,20)))**2 < (1-delta)**5

assert 22*(a-epsilon)**2/(4*a)>25
assert (1+z_radius)*exp_up(Q(9,1000)+epsilon**2/(4*a)-a/4)<Q(13,20)
assert epsilon**2*22*Q(13,20)**44<Q(1, 25000**2)
assert a/(Q(25,8)*160**2)<Q(81,10000)**2

far = paired_far = Q(0)
for i in range(22,80):
    lo, hi = Q(i,80), Q(i+1,80)
    rate, b = Q(22,4)*lo**2-1/(2*a), hi**2/4
    assert rate > 0

    ell=(exp_up(epsilon*(hi+epsilon/a)/2)-1+q0)/(1-q0)
    g0=exp_up(-a*lo**2/4+epsilon*lo/2)
    d=ell*g0+z_radius*(1+ell)*exp_up(-a*(1-hi/2)**2+epsilon*hi/2)
    D,R=round_up(d,10**6),round_up(g0+d,10**6)
    assert R**2*Q(81,80)**3<1
    m2=10**6*max(j**3*R**(2*(j-1)) for j in range(22,81))
    M=Q(isqrt(m2.numerator//m2.denominator)+1,1000)
    F=Q(81,10000)*D*M
    far += F

    B=1/exp_up(a*hi**2/4 + epsilon*hi/2 + epsilon**2/(4*a))
    assert (rate+b)*B<=rate*(1+D)
    paired_far += (1+D-B+Q(3,44))*F

assert far < Q(182,1000)
assert paired_far < Q(58,1000)
assert Q(295,1000)+Q(182,1000)+Q(41,10**6) < Q(12,25)
assert Q(1,20)+Q(58,1000)+Q(82,10**6) < Q(11,100)
scale = z_radius/(z_radius*z_radius-1)
assert scale * (Q(9,1000)+Q(11,100)/(1-Q(12,25))) < Q(1,2)
assert scale * (Q(22*9,1000)+Q(12,13)) < Q(22,4)
~~~
