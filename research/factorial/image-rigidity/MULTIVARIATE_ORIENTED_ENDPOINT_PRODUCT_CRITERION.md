# An oriented endpoint criterion without mirrored pairs

Consider any ordinary two-monomial profile and a common Gaussian
strength \(\tau\ge0\). In degree \(m\) the Gaussian factor at index
\(k\) is \(e^{-\tau k(m-k)/(2m)}\) in the endpoint normalization.
Delete equal rows. Write \(S\) for the total absolute curvature
\(\sum (b_i-a_i)^2/\max(a_i,b_i)\) plus \(\tau\), and \(V\) for
\(\sum (b_i-a_i)^2/\min(a_i,b_i)\) on shared rows plus the sum
of private exponents and \(\tau\). Put \(W=1+V\).

For either orientation, retain any subset of the shared unequal
rows and put \(s_0\) equal to the curvature of the other rows
plus \(\tau\). One may take every shared row. As shown below, doing so gives
the smallest endpoint bounds among all choices of retained rows
in each orientation. On a retained row starting at \(a>0\) and ending at \(b>0\),
use \(d=b-a\) with its sign and \(x(w)=a+dw\). Define the
endpoint chord differences
\[
 \begin{split}
 A(w)&=x(w)\log x(w)-(1-w)a\log a-wb\log b,\\
 G(w)&=\log x(w)-(1-w)\log a-w\log b .
 \end{split}
\]
For \(0<w\le1/2\), define \(Q\) separately in the two orientations by
\[
             \log Q(w)=-s_0(1-w)/2
                   +\sum_{\mathrm{retained}}
                              \left(A(w)/w+G(w)/2\right).
\]
Use the continuous limit at zero, and let \(q\) be the larger
of the two values \(Q(0)\).

**Theorem.** The consecutive moment polynomials have disjoint
root arguments in every degree if both orientations satisfy
\[
           Q(1/30001)\le10^{-6},\qquad Q(1/2)\le1/55,
                    \qquad Wq\le10^{-4}.                       \tag{1}
\]
The ordinary pair at \(\tau=0\) is then coprime over every
characteristic-zero field. A profile certified by (1) remains
certified after appending any ordinary rows or increasing
\(\tau\). There is no symmetry or degree restriction.
In particular, one shared row with gap \(d\) and smaller exponent
\(\ell\) certifies the whole profile if \(d\ge18\ell\), or if
\(d\ge20\) and \(d\ge7\ell\). Several smaller rows can instead
pass the two tests collectively.

## A Gamma chord bound in both orientations

A retained row contributes \(mL_m(w)\) to the normalized
coefficient logarithm, with
\(L_m(w)=m^{-1}[\log\Gamma(1+mx(w))]_{\rm chord}\), where
the subscript subtracts the endpoint chord. The exact trigamma
integral in the
[matched curvature hierarchy](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md#proof-of-the-curvature-hierarchy)
gives
\[
 L_m''(w)=d^2\int_0^\infty e^{-x(w)u}h(u/m)\,du,
 \qquad h(t)=t/(e^t-1).
\]
Since \(h(t)\ge1-t/2\), by \((t/2)\coth(t/2)\ge1\),
this curvature is at least
\(d^2/x(w)-d^2/(2m x(w)^2)\).
That is the second derivative of \(A(w)+G(w)/(2m)\).
Both functions have zero endpoint values, hence
\(L_m(w)\le A(w)+G(w)/(2m)\).
For an integer \(k=mw\ge1\), \(G(w)\ge0\) by concavity,
so \(L_m(w)/w\le A(w)/w+G(w)/2\).
The unretained rows and Gaussian supply \(e^{-s_0 k(1-w)/2}\)
by the [discrete curvature estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate).

The function \(w(\log Q)'(w)\) is nonnegative and increasing
on \(0\le w\le1/2\). For one row put \(\delta=d/a\) and
\(I(x)=\int_0^1 u(1+xu)^{-2}du\). Differentiation under this
finite integral gives
\[
 2\frac{d}{dw}\!\left[w\frac{d}{dw}
                  \left(A(w)/w+G(w)/2\right)\right]
   =2a\delta^2 I(\delta w)
         +\frac{\delta}{(1+\delta w)^2}-\log(1+\delta).
\]
For \(\delta>0\) each nonconstant term decreases in \(w\).
For \(\delta=-\zeta<0\), replace \(a\) by \(1\) for a lower
bound, and put \(y=\zeta w\). The derivative in \(y\) of
\(2\zeta^2 I(-y)-\zeta/(1-y)^2\) is at most
\((4\zeta^2/3-2\zeta)/(1-y)^3<0\) since \(0<\zeta<1\).
Thus in either case it suffices to bound the expression at
\(w=1/2,\ a=1\), where it equals
\[
 E(\delta)=8\log(1+\delta/2)-\frac{8\delta}{2+\delta}
                  +\frac{4\delta}{(2+\delta)^2}-\log(1+\delta).
\]
Its derivative is
\(\delta(7(\delta+1)^2+1)/[(\delta+2)^3(\delta+1)]\).
Since \(E(0)=0\), \(E(\delta)\ge0\) on \(\delta>-1\).
Adding the residual term proves the asserted monotonicity.
In particular \(Q\) increases, since \(w(\log Q)'(w)\) vanishes at zero.

Also
\[
 0\le\log Q(w)-\log Q(0)\le (3V/4)w .
\]
In fact \((A(w)/w)'=d^2\int_0^1u/x(wu)\,du
\le d^2/(2\min(a,b))\). The concave correction \(G(w)/2\)
is at most \(\xi w\) with
\[
 \xi=\tfrac12(d/a-\log(b/a))
                      \le d^2/(4a\min(a,b)).
\]
The latter estimate follows by integrating \(t/(1+t)\) between
\(0\) and \(d/a\), with its orientation, and bounding the denominator
below by \(\min(1,b/a)\). Finally \(s_0\) is bounded by the
corresponding portion of \(V\).

These estimates give for the normalized inclusive coefficients
in either orientation
\[
 c_{m,k}^{1/k}\le (e/w) Q(w),\qquad
          c_{m,k}\le \frac{(mq)^k}{k!}\exp(3V k^2/(4m)).
                                                               \tag{2}
\]
Multiplying the logarithmic derivative of \(eQ(w)/w\) by \(w\)
gives an increasing function. Its maximum on a closed
\(w\)-interval is therefore attained by \(eQ(w)/w\) at an endpoint.
Also \(mQ(1/m)\) has no interior maximum on an interval in
\(m\ge2\), because \(m\) times its logarithmic derivative is
\(1-w(\log Q)'(w)\) at \(w=1/m\), increasing as \(m\) increases.

## Retaining every shared row improves the bounds

For \(0<w\le1/2\), let \(U\) have density
\[
       \frac{2\min(u,w)(1-\max(u,w))}{w(1-w)}
                          \quad (0\le u\le1).
\]
The two elementary integrals over \(u\le w\) and \(u\ge w\)
give \(\mathbb E U=(1+w)/3\). Put \(X=a+dU\),
\(M=\max(a,b)\), and \(e_0=|d|\). The Green integral for
the two chord differences gives exactly
\[
  A(w)/w+G(w)/2
   =-\frac{(1-w)d^2}{2}
                  \mathbb E\left[\frac1X-\frac{w}{2X^2}\right].
\]
In both orientations \(\mathbb E X\le M-e_0/3\).
The function \(f(x)=1/x-w/(2x^2)\) is convex and decreasing
for \(x\ge1\). Jensen's inequality bounds its expectation
below by \(f(M-e_0/3)\ge1/M\).
For the last inequality use
\[
 (e_0/3)(M-e_0/3)\ge (M-1/3)/3>M/4\ge wM/2.
\]
The first bound uses \(1\le e_0<M\) and the increase of
the quadratic over this range. Thus
\[
            A(w)/w+G(w)/2
                    \le-\frac{(1-w)d^2}{2M}.
\]
The inequality at zero follows by continuity. This shows
that keeping a shared row explicitly improves on using
its curvature alone. In particular the best version of (1)
uses every shared row, and the unused curvature is then
just the private curvature plus \(\tau\).

## Uniform comparison across degrees

Use the two oriented halves and the
[short phase comparison](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#short-halves-and-large-absolute-curvature),
which needs only the coefficient bounds given below.
Write \(N=n+1\) and assume \(n\ge2\).
For \(N\le30001\), (1) and the two endpoint maxima
give \(c_{m,k}\le10^{-k}\) for \(k\le m/2\).
At \(w=1/2\) use \(2e/55<1/10\), and at \(w=1/m\)
use the endpoint bound for \(mQ(1/m)\), namely
\(e\max(2/55,30001\cdot10^{-6})<1/10\).
For \(N>30001,\ S\ge2\log N+5\), the
[large-curvature bound](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
gives \(c_{m,k}\le8^{-k}\); its calculation uses only the discrete
coefficient bound with \(S\), which includes \(\tau\).

Otherwise \(N>30001,\ S<2\log N+5\).
Put \(D=Nq\) and \(L=3V/4\). Then
\[
 L/n\le \frac{3(2\log N+5)(2\log N+7)}{4(N-1)}<1/32.
\]
The upper expression decreases for \(N\ge30001\), by
differentiating its numerator as a function of \(\log N\).
At \(N=30001\), use \(\log30001<11\).
If \(D\le1/10\), (2) gives
\(c_{m,k}\le(1/8)^k/k!\) for \(1\le k\le7\), since
\((1/10)e^{7/32}<1/8\).
For \(8\le k\le m/2\) its first bound gives \(10^{-k}\),
testing \(k=8\) with \((e/80)e^{1/4}<1/10\) and the
midpoint. These are the three short ranges of the cited
phase comparison. In each, every root is in the annulus
\(|\log|z||<1/(2m)\) and its phase correction has modulus
less than \(\theta/2\) at angle \(0<\theta\le\pi\).
Subtracting phases at a common argument is impossible.

In the long range \(D>1/10\) we have \(n\ge100000\), since
\(q\le10^{-6}\). Put \(K=\lceil32D\rceil+14\). As
\(Wq\le10^{-4}\),
\[
 D/n<1/720000,\quad WD/n<1/9000,\quad
 K\le182D,\quad WK/n<182/9000,\quad K+1<n/2 .
\]
Testing (2) at \(k=K\) using
\((e/32)e^{(3/4) WK/n}<1/8\), and then at the midpoint,
gives \(c_{m,k}\le8^{-k}\) for \(K\le k\le m/2\).
At a positive starting endpoint the interpolated first factor
\(B_m\) has \(|(\log B_m)'|\le2/n\) by the same
[asymmetric Gamma derivative estimate](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#long-halves-and-the-change-of-order).
Indeed \(\log N>11\), \(S<2\log N+5\), and the numerator
\(5V+S(3+\log N)\) is at most \(40(\log N)^2\).
Equation (2) then gives \(B_m\le De^{L/n+2/n}<101D/100\),
using \(D>1/10\) and \(WD/n<1/9000\).
The same connected logarithm estimate applies. Its inputs are
\(\alpha=2\max(1,\max_i|d_i|/a_i)/n\le2W/n\),
\(\beta=8eW/n\), the first-factor bounds, and the coefficient
tails above. In particular \(\alpha K<1/20,\ \beta K<1/2\) and
\(\beta B_m<1/300\).

If a half has a private starting row, cancel one factorial against
the binomial coefficient as in the
[private-start endpoint estimate](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation).
That row was included in \(s_0\). Removing one curvature unit
gives \(c_{m,k}^{1/k}\le e^{3/2} Q(w)\).
For \(1\le k\le K\), \(Q(w)\le q e^{LK/n}\);
on radius \(6/5\) the resulting geometric ratio is less than
\(6\cdot10^{-6}\). Beyond \(K\) the ratio is at most \(3/20\).
Hence \(\|H_m-1\|_{6/5}<1/200\) and
\(\|\log H_m\|_{6/5}<1/190\).
In both the positive and private starting cases we thus obtain
\[
       \|\log H_m\|_{6/5}\le3(D+1),\qquad
                   \|\log H_N-\log H_n\|_{6/5}<1/55 .
\]
These are exactly the bounds used in the
[paired radius and phase argument](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#the-two-radii-do-not-affect-the-phase-contradiction).
It uses \(D>.1,\ D/n<1/720000\) and the common angle.
It again excludes a common root argument.
Order one versus two follows from strict positivity of the
quadratic factorial moment on the real line, preserved by the
positive Gaussian averaging.

## Preservation under arbitrary added rows

The definition at zero and the endpoint Green integral give
\[
  -\log Q(0)=s_0/2+
       \sum_{\rm retained}\int_0^1(1-u)\frac{d^2}{x(u)}\,du
                      \le V/2.
\]
Consequently a qualifying profile has \(W\ge5\).
For added curvature \(\sigma\ge0\), keep the retained rows
and add \(\sigma\) to the unused curvature. Each \(Q(w)\)
is multiplied by \(e^{-\sigma(1-w)/2}\), and
\(W_{\rm new}\le W+\sigma(\sigma+2)\) by
\(V\le S(S+2)\) for the added profile.
The same comparison holds for increasing \(\tau\) by \(\sigma\),
which adds \(\sigma\) to \(V\).
The function \((W+\sigma(\sigma+2))e^{-\sigma/2}\) is
nonincreasing for \(W\ge5\), since the sign of its derivative
is that of \(\sigma+2-\sigma^2/2-W/2\le0\).
Thus each of the tests is preserved.

## A single shared row with relative gap at least eighteen

Let the retained row be \((\ell,\ell+d)\) with \(d\ge18\ell\).
It suffices by preservation to check the profile consisting
of that row, at zero Gaussian strength. Write \(r=d/\ell\)
and \(P=\log(1+r)\). At zero the two endpoint quantities are
\(e^{d-(\ell+d)P}\) and \(e^{-d+\ell P}\).
The decreasing orientation has the larger one because
\(\log(1+r)\ge2r/(2+r)\), by differentiation.
As \(P/r\) decreases in \(r>0\), it gives
\[
 \bigl(1+d^2/\ell\bigr)q
 \le (1+d^2)\exp[-d(1-\log(19)/18)]
 \le 325\cdot19e^{-18}<10^{-4}.
\]
The middle expression decreases for \(d\ge18\), using
\(1-\log(19)/18>5/6\).

For the other two tests use \(G(w)/2\le \xi w\).
In the increasing orientation, the logarithm is at most
\[
  -\ell P+d\left(1+\frac w2+
                \int_0^1\log\frac{1+rwu}{1+r}\,du\right).
\]
The integral decreases in \(r\), so evaluate at \(r=18\).
At \(w=1/2\) the bracket is
\(1/4+(10/9)\log10-\log19<-1/9\), using
\(\log10<231/100,\ \log19>294/100\).
It follows that \(Q(1/2)\le e^{-2}/19<1/55\).
At \(w\le1/100\), Jensen's inequality bounds the integral by
\(\log(109/1900)<-14/5\), which makes \(Q(w)<e^{-18}<10^{-6}\).

For the decreasing orientation set
\(\lambda=d/(\ell+d)\). The upper logarithm is
\[
 -d+\ell P -d\int_0^1\log(1-\lambda wu)\,du
                                      +w(P-\lambda)/2.
\]
Use \(-\log(1-\lambda wu)\le wu/(1-w)\) and
\(\ell P\le d/6\). At \(w=1/2\) this gives
\(-d/2+(5/4)\ell P\le-7d/24\le-21/4\).
At \(w\le1/1000\) it gives
\(-d+d/6+d/1998+d/12000<-4d/5\le-72/5\).
These prove the tests since \(e^{21/4}>55\) and
\(e^{72/5}>10^6\).

## A second one-row condition and a collective example

Suppose \(d\ge20\) and \(r=d/\ell\ge7\). If \(d\ge31\),
then \(S\ge7d/8>27\), so the
[Gaussian curvature theorem](MULTIVARIATE_MIRRORED_CORE_GAUSSIAN_ENDPOINT_TRANSFER.md#absolute-curvature-for-a-general-ordinary-core-plus-gaussian)
handles the full profile directly. For \(20\le d\le30\)
it suffices to check the one row by preservation.
At fixed \(d\), the logarithmic derivative of
\((1+d^2/\ell)q\) with respect to \(\ell\) is
\[
    \log(1+r)-\frac{r}{1+r}-\frac{d^2}{\ell(\ell+d^2)}>0,
\]
using \(\log8-7/8>1\) and \(\ell\ge1\).
For \(d=20\) take \(\ell=2\) by integrality, giving
\(Wq\le201\cdot121e^{-20}<10^{-4}\).
For \(d\ge21\) use \(\ell\le d/7\), giving
\(Wq\le(1+7d)8^{d/7}e^{-d}\le148\cdot512e^{-21}<10^{-4}\).
The upper bound decreases for \(d\ge21\).
The same monotonicity of \((1+d/\ell)^\ell\) bounds
\(q\) by the larger of \(121e^{-20}\) and \(512e^{-21}\).
As \(V\le900\), the increment bound gives \(Q(1/30001)<10^{-6}\).

At the midpoint, \(2A(1/2)/d\) decreases with \(r\), by its
Green integral with integrand proportional to \(r/(1+ru)\).
At \(r=7\) it is \((9\log(9/2)-8\log8)/7<-11/25\);
equivalently \((9/2)^9 e^{77/25}<8^8\).
In either orientation \(G(1/2)/2\le(r-\log(1+r))/4\):
for the decreasing orientation use
\(r+r/(1+r)\ge2\log(1+r)\) by differentiation.
It follows that
\(\log Q(1/2)<-19d/100-\log8/4\le-19/5-\log8/4\),
so \(Q(1/2)<1/55\).

For another example take four copies of \((1,7)\), all with
the same orientation. They have \(S=144/7<27\) and \(W=145\).
The two endpoint tests and the product test in (1) follow
from the rigorous evaluations below. There are no mirrored
pairs and the individual rows do not pass (1) on their own.
A lower-curvature example is the two rows \((1,9),(8,1)\),
with \(S=953/72\) and \(W=114\). The same three tests
pass in both orientations. Arbitrary additional rows are
allowed in either example.

## Scalar replay

The numerical margins use rational Taylor bounds.

~~~python
from fractions import Fraction as Q
from math import factorial

def lo(x):
    return sum(Q(x)**j / factorial(j) for j in range(81))

def hi(x):
    x = Q(x)
    return lo(x) + x**81 / (factorial(81)*(1-x/82))

assert hi(Q(294,100)) < 19 and lo(3) > 19 and lo(Q(231,100)) > 10
assert Q(1,4)+Q(10,9)*Q(231,100)-Q(294,100)<-Q(1,9)
assert hi(Q(14,5))<Q(1900,109)
assert lo(18) > 325*19*10**4 and lo(18)>10**6
assert 19*lo(2)>55 and lo(Q(21,4))>55 and lo(Q(72,5))>10**6
assert -1+Q(1,6)+Q(1,1998)+Q(1,12000)<-Q(4,5)
assert lo(11)>30001 and hi(11)<10**5
assert Q(3,4)*27*29/30000<Q(1,32)
assert hi(Q(7,32))/10<Q(1,8) and hi(Q(5,4))/80<Q(1,10)
assert hi(1)*Q(2,55)<Q(1,10)
assert hi(1)*Q(30001,10**6)<Q(1,10)
assert hi(1+Q(3,4)*182/9000)/32<Q(1,8)
assert hi(Q(3,4)/900 + Q(2,100000))<Q(101,100)
assert 8*hi(1)*Q(101,100)/9000<Q(1,300)
assert 8*hi(1)*182/9000<Q(1,2)
assert Q(2*182,9000)<Q(1,20)
assert Q(6,5)*hi(Q(3,2)+Q(3,4)*182/9000)*Q(1,10**6)<Q(6,10**6)
assert hi(Q(15,8))<8 and lo(3)>8
assert lo(20)>201*121*10**4 and lo(21)>148*512*10**4
assert hi(Q(675,30001))*max(121/lo(20),512/lo(21))<Q(1,10**6)
assert Q(9,2)**9 * hi(Q(77,25)) < 8**8
assert lo(Q(19,5))*Q(3,2)>55 and Q(3,2)**4<8
~~~

The collective example uses rigorous Arb intervals.

~~~python
from flint import arb, ctx, fmpq
ctx.prec = 160

def log_q(a, b, w):
    a, b = arb(a), arb(b)
    if w == 0:
        return b-a-b*(b/a).log()
    w = arb(w)
    x = a+(b-a)*w
    A = x*x.log()-(1-w)*a*a.log()-w*b*b.log()
    G = x.log()-(1-w)*a.log()-w*b.log()
    return A/w+G/2

for W, rows in [(145, [(1,7)]*4), (114, [(1,9),(8,1)])]:
    for oriented in [rows, [(b,a) for a,b in rows]]:
        def Q(w):
            return sum(log_q(a,b,w) for a,b in oriented).exp()
        assert Q(fmpq(1,30001)) < arb(fmpq(1,10**6))
        assert Q(fmpq(1,2)) < arb(fmpq(1,55))
        assert W*Q(0) < arb(fmpq(1,10**4))
~~~
