# All-degree interlacing with one private unit pair and a mirrored core

Consider an ordinary two-monomial support with a distinguished
private unit pair \((1,0),(0,1)\). Its other nonconstant
exponent rows occur in mirrored pairs
\[
        (\ell_i+d_i,\ell_i),\qquad
        (\ell_i,\ell_i+d_i),\qquad \ell_i\in\mathbb Z_{\ge0},\quad d_i\in\mathbb Z_{>0}.
\]
Rows with equal exponents contribute only a scalar.

**Theorem.** For every such support, the consecutive factorial
moment polynomials strictly interlace on the unit circle in
every degree, with angles listed from the positive real axis.
Consequently SFC(2) holds for this family over every
characteristic-zero field.

The same angular conclusion holds for arbitrary nonnegative
real \(\ell_i\), positive real \(d_i\), and equal positive
Gamma shape \(\gamma_i\) on the two coordinates of each
shared pair. The private unit pair retains Gamma shape
one. One may also multiply the coefficient at \(k\) in
degree \(j\) by \(\exp[-\tau k(j-k)/(2j)]\), for a common
\(\tau\ge0\) at all degrees.

The private pair cancels the binomial coefficient. The
matched law of the shared log ratio then gives a uniform
bound on the change of the exact root phase. The argument
uses no size restriction on the core.

## The exact paired angular equation

Use the
[resolvent formula and angle boxes](MULTIVARIATE_ONE_PRIVATE_MIRRORED_ANGLE_BOXES.md#the-exact-resolvent-and-angle-boxes)
and the
[folded likelihood-ratio comparison](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md#the-folded-likelihood-ratio-comparison).
Let \(\rho_j^+\) be the folded midpoint law. In the ordinary
case this is the law of the absolute core log ratio under
the order-\(j\) geometric mean tilt. Write \(c=\cos\theta\) and
\[
 R_j(\theta)=
 \frac{\displaystyle\int
       \frac{\cosh((j/2+1)s)}{\cosh s-c}\,d\rho_j^+(s)}
      {\displaystyle\int
       \frac{\cosh(js/2)}{\cosh s-c}\,d\rho_j^+(s)} .
\]
The ratio exceeds one and increases with \(\theta\) on
\((0,\pi)\). The positive angles are simple and the \(k\)-th
angle is in the open box
\[
             \left(\frac{(2k-1)\pi}{j},\ \frac{2k\pi}{j+1}\right).
                                                               \tag{1}
\]
There it solves \(R_j(\theta)=\cos\theta-\tan(j\theta/2)\sin\theta\),
whose right side is strictly decreasing within the box.
Odd orders also have the root \(-1\). If the core is empty,
the polynomials are geometric sums and already strictly
interlace when \(\tau=0\).

Set \(N=n+1\), hold \(\theta\in(0,\pi/2]\) fixed, and put
\[
 Q(t)=
 \frac{\displaystyle\int
        \frac{\cosh((n/2+t)s)}{\cosh(ts)-c}\,d\rho_n^+(s)}
      {\displaystyle\int
        \frac{\cosh(ns/2)}{\cosh(ts)-c}\,d\rho_n^+(s)},
                      \qquad n/N\le t\le1 .
\]
The matched likelihood-ratio order says that dilation of
\(\rho_n^+\) by \(n/N\) is below \(\rho_N^+\). A common
positive reweighting preserves this order, and the ratio of
the two resolvent integrands increases in \(s\). Therefore
\[
                    R_N(\theta)\ge Q(n/N),\qquad Q(1)=R_n(\theta).
                                                               \tag{2}
\]

## A bound independent of the core

Write \(x=ts\), \(D=\cosh x-c\), and
\[
 f=\cosh((n/2+t)s)-\cosh(ns/2)
     =2\sinh(ts/2)\sinh((n+t)s/2).
\]
Let \(B\) and \(\nu\) be the two probability laws obtained
from \(\rho_n^+\) by reweighting with \(\cosh(ns/2)/D\)
and \(f/D\), respectively. The latter dominates the former
in likelihood-ratio order, since
\[
           \frac{f}{\cosh(ns/2)}
               =\cosh x-1+\tanh(ns/2)\sinh x
\]
increases in \(s\). These integrals may be differentiated in
\(t\) on the indicated compact interval. The numerator weight
has growth \(e^{ns/2}\) at infinity, and the corresponding
endpoint tilt has a positive exponential moment. In a modified Gamma pair, the
midpoint density before folding has exponential decay
\(\gamma_i/d_i+n\ell_i/d_i+n/2>n/2\).
Finite convolution and an optional Gaussian factor preserve
a strict exponential moment margin. The
[Gaussian augmentation of the folded comparison](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md#adding-the-gaussian-symbol)
gives (2) under that modification.

Set \(\delta=1-c\), \(h=\sin(\theta/2)\), and
\[
       p=\frac{\cosh x-1}{\cosh x-c},\qquad d_0=s\coth(x/2).
\]
The definitions extend at \(s=0\) by continuity.
We have \(\partial_t\log D=d_0p\), and the factorization
of \(f\) gives \(\partial_t\log f\le d_0\). Thus
\[
  \partial_t\log(Q-1)
            \le \mathbb E_\nu[d_0(1-p)]+\mathbb E_B[d_0p].
\]
Using \(d_0\le2/t+s\), the terms with \(2/t\) sum to at
most \(2/t\), since \(p\) increases in \(s\) and \(\nu\)
dominates \(B\). The other two terms admit the bounds
\[
       s(1-p)\le\frac{s\delta}{\delta+x^2/2}\le\frac{h}{t},
                 \qquad
       \mathbb E_B[sp]\le\mathbb E_Bs\le\frac{\operatorname{arcosh}Q}{t}.
\]
The last inequality is Jensen's inequality for \(\cosh\), since
\(Q\ge\mathbb E_B\cosh(ts)\). It follows that
\[
              \partial_t\log(Q-1)
                     \le\frac{2+h+\operatorname{arcosh}Q(t)}{t}.
                                                               \tag{3}
\]
If \(Q(n/N)\ge Q(1)\), the next estimate follows immediately.
Otherwise integrate (3) from \(n/N\) to the first point at
which \(Q(t)=Q(1)\). In that interval \(Q(t)\le Q(1)\).
Writing \(r=R_n(\theta)-1\), (2) then gives
\[
           R_N(\theta)-1 \ \ge\
           r\,(n/N)^{2+h+\operatorname{arcosh}(1+r)} .
                                                               \tag{4}
\]

## A scalar phase comparison

For \(n\ge3\), \(r>0\), and
\(\pi/n\le\theta\le\pi/2\), we prove
\[
 r\,(n/N)^{2+h+\operatorname{arcosh}(1+r)}
       >\frac{2(1-h^2)r}{r+2}-2h^2,\qquad h=\sin(\theta/2).
                                                               \tag{5}
\]
After division by \(r\), the right side is
\(b=2((r+1)c-1)/(r(r+2))\). If \(b\le0\) the result is
immediate. Put \(L=\log(1+1/n)\). For \(b>0\),
\[
       -\partial_\theta\log b
           =\frac{(r+1)\sin\theta}{(r+1)c-1}>\tan\theta .
\]
The derivative of \(L(2+h+\operatorname{arcosh}(1+r))\)
is at most \(1/(2n)\). Thus the logarithmic margin in (5)
is increasing in \(\theta\), and it is enough to use
\(\theta=\pi/n\).

First take \(n\ge4\). At \(\theta=\pi/n\), we have
\(n^2\delta\ge16(1-2^{-1/2})>23/5\) and \(h<2/5\).
The first bound follows from the decrease of \(\sin u/u\)
for \(u\in(0,\pi/8]\), and the second from the value
of \(\sin(\pi/8)\). Since
\(b=(1+r/2)^{-1}[1-\delta(1+1/r)]>0\),
\[
 \begin{split}
 n(-\log b)&\ge n\log(1+r/2)+\frac{23}{5n}(1+1/r)\\
            &\ge 2\sqrt{\frac{46(r+1)}{5(r+4)}}.               \tag{6}
 \end{split}
\]
The last step combines \(\log(1+r/2)\ge2r/(r+4)\) with
the arithmetic–geometric mean inequality. It suffices
to bound \(n(-\log b)\) above \(12/5+\operatorname{arcosh}(1+r)\)
because \(nL<1\).

For \(0<r\le5\), set \(y=\operatorname{arcosh}(1+r)<5/2\).
Squaring the desired comparison of (6) with \(12/5+y\)
reduces it to
\[
        \bigl(184/5-(12/5+y)^2\bigr)\cosh y
                       -3(12/5+y)^2>0.
\]
The coefficient of \(\cosh y\) is positive in this interval.
Replacing \(\cosh y\) by \(1+y^2/2+y^4/24+y^6/720\)
still gives a positive expression. In order from zero to eight,
the Bernstein coefficients of the resulting polynomial on
\([0,5/2]\) are
\[
 \frac{344}{25},\quad\frac{194}{25},\quad\frac{758}{175},\quad
 \frac{7853}{2800},\quad\frac{14173}{4800},\quad
 \frac{623969}{134400},\quad\frac{24189281}{3225600},\quad
 \frac{2320729}{230400},\quad\frac{5421491}{921600}.
\]
They follow by substituting \(y=5z/2\) and replacing each
power \(z^j\) by its Bernstein coefficients
\(\binom{k}{j}/\binom{8}{j}\) for \(k\ge j\), with zero
in the preceding positions. The basis weights are
nonnegative and sum to one on \([0,1]\).

For \(r\ge5\), use the first line of (6) and
\[
              4\log(1+r/2)>12/5+\operatorname{arcosh}(1+r).
\]
Indeed the derivative of the difference is
\(4/(r+2)-1/\sqrt{r(r+2)}>0\) for \(r\ge5\).
At five, \(\operatorname{arcosh}6<5/2\) follows by the
same positive Taylor lower bound for \(\cosh(5/2)\),
and \(\log(7/2)>5/4\) follows from the first four terms
of \(2\operatorname{artanh}(5/9)\). This gives (5) for
\(n\ge4\).

For \(n=3\), we may use \(c\le1/2\) and \(h\le1\).
Then \(b\le(r-1)/(r(r+2))<1/6\), the latter inequality
being equivalent to \((r-2)^2+2>0\).
For \(0<r\le8\), the divided left side of (5) exceeds
\((3/4)^6>1/6\), using
\(\operatorname{arcosh}(1+r)\le\log(2+2r)\) and
\(\log18<3\). For \(r\ge8\), the undivided left side
is at least \(r\exp[-L(3+\log(2+2r))]\), an increasing
function of \(r\) since \(L<1\). At eight it exceeds
\(8(3/4)^6>1\). The right side of (5) is at most
\(2-4h^2=2c\le1\). This completes the scalar comparison.

## Strict alternation

Consider the \(k\)-th positive lower root, of angle
\(\theta\). If \(\theta>\pi/2\), (1) gives
\(4k>n+1\), hence \(4k\ge n+2\). Since the lower angle
exceeds \((2k-1)\pi/n\), the upper \(k\)-th angle is
already earlier by its upper box endpoint \(2k\pi/(n+2)\).
For \(\theta\le\pi/2\), the lower box gives
\(\theta>\pi/n\), and adding \(\theta/2\) to the lower
phase equation yields
\[
        \cos\theta-\tan(N\theta/2)\sin\theta
                      =c+\frac{(1+c)r}{r+2},
                   \qquad r=R_n(\theta)-1 .
\]
The lower angle exceeds the left endpoint of the upper
\(k\)-th box. If it is at or beyond the upper endpoint,
that upper root already comes first. Otherwise (4) and
(5) show \(R_N(\theta)\) strictly above the displayed
right side. The two functions cross in the upper box
with \(R_N\) increasing and the right side of its root
equation decreasing. This again puts the upper root first.

The reverse inequality between the lower \(k\)-th root
and the upper \((k+1)\)-st follows from (1).
For \(n=1,2\), both sides follow directly from (1) and
the odd structural root. Reflection and that structural
root give the stated strict interlacing in every degree.
The ordinary moment polynomials are rational, so their
Bézout identity transfers coprimality to every
characteristic-zero field.

Only the lower matched likelihood-ratio comparison was
used. More generally, the proof applies to even midpoint
laws \(\rho_j\) with positive folded densities,
\(\int (1+s)e^{js/2}\,d\rho_j^+(s)<\infty\), and
\((j/(j+1))\rho_j^+\preceq\rho_{j+1}^+\) at each step.
The corresponding order-\(j\) coefficients are
\(\int \cosh((k-j/2)s)\,d\rho_j^+(s)\) for
\(0\le k\le j\). The root equation and boxes (1) follow
from summing that finite geometric series. The stated
moment condition justifies differentiation of \(Q\)
on its order-comparison interval.
