# A monotone curvature for the direct factorial coefficient ratio

Fix \(N=n+1\), \(n\ge2\), and ordinary exponent columns \(a,b\).
Write \(d_i=b_i-a_i\) and let \(P_j^\tau\) be the order-\(j\)
factorial moment polynomial with each coefficient at \(z^k\)
multiplied by \(e^{-\tau k(j-k)/(2j)}\), where \(\tau\ge0\).
Compare \(P_n^\tau(z)\) and \(P_N^\tau(\kappa z)\) for a
positive real \(\kappa\). For \(0\le k\le n\), the ratio
of the upper coefficient to the lower coefficient is
\(e^{f(k)}\), where
\[
 \begin{aligned}
 g(x)&=\prod_{i:a_i>0}\ \prod_{h=1}^{a_i}(na_i+d_i x+h),\\
 f(x)&=\log \frac{N}{N-x}+\log g(x)+x\log\kappa
                                      -\frac{\tau x^2}{2nN}.
 \end{aligned}                                                     \tag{1}
\]
The factorial and binomial ratios are those of the
[direct first-sector test](MULTIVARIATE_DIRECT_RATIO_FIRST_SECTOR.md#the-exact-factorial-ratio).
The difference of the two Gaussian exponents at the same \(k\)
is exactly \(-\tau k^2/(2nN)\).

Assume \(a_i\le N b_i\) for each \(i\) with \(a_i>0\).
Then the factors in \(g\) remain positive through \(x=N\).
The following exact identity holds on \(0\le x<N\):
\[
 \begin{aligned}
 \Phi(x):=(N-x)^2 f''(x)
   &=1-\sum_{i:a_i>0}\sum_{h=1}^{a_i}
        \left(\frac{d_i(N-x)}{na_i+d_i x+h}\right)^2
                        -\frac{\tau}{nN}(N-x)^2.                 \tag{2}
 \end{aligned}
\]
The function \(\Phi\) is nondecreasing, and it is strictly
increasing if \(\tau>0\) or a contributing row is unequal.
Indeed
\[
   \frac{d}{dx}
       \frac{N-x}{na_i+d_i x+h}
       =-\frac{N b_i-a_i+h}{(na_i+d_i x+h)^2}<0 .
\]
Thus each unequal row decreases the corresponding squared
term. In particular \(f''\) can change sign only once,
from negative to positive, and \(f''>0\) near \(N\).

The hypothesis of (2) holds in at least one orientation of
every profile with \(N\ge14\), curvature
\(S=\sum_i d_i^2/\max(a_i,b_i)\le12\), and private rows on
at most one side. Orient any private rows so \(a_i=0\).
For shared rows the ratio of larger to smaller exponent is
at most \(13\). Indeed if the smaller exponent is one,
a larger exponent at least \(14\) gives curvature at least
\(169/14>12\). If the smaller exponent is at least two and
the ratio is greater than \(13\), the row contributes more
than \(2\cdot144/13>12\). This proves the hypothesis in
both orientations for the shared rows. The Gaussian
strength and bias do not affect this curvature-shape conclusion.

## The discrete shape on the exact matched coefficients

Let \(w(u)=\min(u,2-u)\) for \(0\le u\le2\), and zero
elsewhere. For integer \(0\le k\le n-2\), put
\[
 I_k=\int_k^{k+2}\frac{w(x-k)}{(N-x)^2}\,dx>0 .
\]
Then \(\Delta^2 f(k)/I_k\) is nondecreasing in \(k\).
It is strictly increasing in the strict case of (2).
Indeed the integral representation of the second
difference makes this ratio the mean of \(\Phi\) against
density proportional to \(w(x-k)/(N-x)^2\).
On the overlap of two adjacent windows the ratio of the
later and earlier weights is
\((x-k-1)/(k+2-x)\), which increases in \(x\).
Outside their overlap only the earlier density is supported
on the left and only the later one on the right.
Consequently the means of an increasing function are ordered.
Explicitly, the difference of the two means is the integral
over \(x<y\) of \(\Phi(y)-\Phi(x)\) times the nonnegative
cross product of the later weight at \(y\) and earlier at
\(x\) minus the reversed product, divided by \(I_k I_{k+1}\).
This also gives strictness. In particular the discrete log
ratios can switch from negative to positive second differences
only once. This assertion includes the actual raw integer
coefficients and makes no claim about the signs of angular
velocities at arbitrary roots.

## A second-sector test from the valley of the ratio

Suppose \(f'(0)<0\) under the hypothesis of (2).
As \(f'(x)\to+\infty\) at \(N\), the increasing curvature
sign in (2) makes \(f\) strictly decrease and then strictly
increase, with a unique minimum.
Write \(t=\pi/\theta\) for a candidate upper-half-plane angle.
The two parents cannot have a common root of argument
\(\theta\) under either of these tests:
\[
 \begin{array}{ll}
 N/3<t<N/2,& f(N-t)< f(N-2t)\quad\text{and}\quad f(2t)>f(t),\\[2pt]
 N/2\le t<N,& f(t)<f(0).
 \end{array}                                                       \tag{3}
\]
For the first line continuity gives an \(x\in(N-2t,t)\)
such that \(f(x)=f(x+t)\). For the second line use
\(0<x<N-t\), since the difference \(f(x+t)-f(x)\) starts
negative and tends to \(+\infty\) by the pole at \(N\).
In both cases \(x<t,\ x+t<N<x+2t\). The valley property
makes \(e^{f(k)}-e^{f(x)}\) nonnegative outside the two
cuts \(x,x+t\) and nonpositive between them.

Suppose \(z=R e^{i\theta}\), \(R>0\), were a common root.
Write \(A_k,B_k\) for the lower and scaled upper coefficients.
Subtract \(e^{f(x)}P_n^\tau(z)\) from \(P_N^\tau(\kappa z)\)
and take imaginary parts after multiplication by
\(e^{-i(\theta x-\pi)}\).
For \(0\le k\le n\) the terms are
\[
       (e^{f(k)}-e^{f(x)}) A_k R^k
                                   \sin(k\theta-\theta x+\pi).
\]
The sine is nonnegative outside the two cuts and
nonpositive between them, by \(x<t\) and \(N<x+2t\).
Their sum is nonnegative. The last term has strictly
positive sine because \(x+t<N<x+2t\), and its coefficient
is positive. This contradicts the two root equations.
The proof is a two-cut analogue of the first-sector sine
test and includes the boundary \(\theta=2\pi/N\)
when the second line applies.

## A uniform collar under the surviving curvature bound

Assume \(N\ge14\), \(S+\tau\le12\), and orient private rows
to \(a_i=0\) when present. Take
\(\kappa=(N/n)^h e^{-\lambda/(2nN)}\) with \(|\lambda|\le\tau\),
as in the physical Gaussian reference. For the ordinary pair
use \(\tau=\lambda=h=0\).
Write
\[
 L=\log(N/n),\qquad
 D=h+\sum_{i:a_i>0}d_i,\qquad b_n=2-6/n .
\]
Suppose \(D\le-2\). Then the two parents have no common root
when
\[
           \frac{\pi}{N s_+}<|\arg z|<
           \frac{\pi}{N s_-},                                    \tag{4}
\]
where \(s_+\in(1/2,1)\) is the unique solution of
\(-\log(1-s_+)=b_n s_+\), and one may take
\[
 \begin{aligned}
 B&=\max(0,-\log\kappa-(\log g)'(N)+\tau/n),\\
 E&=e^{NB/2},\\
 s_-&=\max\left(\frac13,\frac{\log2}{b_n},
                                \frac{E-1}{2E-1}\right)<\frac12.
 \end{aligned}                                                     \tag{5}
\]
This is an explicit nonempty collar around \(2\pi/N\).
For a tolerance uniform in the core and tail at fixed \(n,D\)
one can replace \(B\) in (5) by
\[
                    -D L+\frac{78}{nN}
                               +\frac{156N}{n(n-12)} .            \tag{6}
\]
In particular the ordinary collar applies to every fully
shared profile of curvature at most \(12\) and absolute
total degree difference at least two, after reversing
columns to make the degree difference negative if needed.
The root inversion from reversal preserves absolute angle.

Here are the interval estimates. For each \(a_i>0\), the
right sum of the decreasing convex function
\(1/(na_i+t)\) for \(0\le t\le a_i\) differs from its
integral \(L\) by at most \(1/(2nN a_i)\).
The upper bound at \(x=0\) uses the error only on
decreasing rows, where \(|d_i|/a_i\le d_i^2/a_i\) and
\(a_i=\max(a_i,b_i)\). With \(G=\log g\), this gives
\[
  \log\kappa+G'(0)
        \le D L+\frac{S+|\lambda|}{2nN}
        \le-\frac{b_n}{N},
\]
using \(L\ge1/N\).
In particular \(f'(0)<0\).
The function \(x\mapsto G(x)+x\log\kappa-\tau x^2/(2nN)\)
is concave, so for \(s=t/N\) the same bound gives
\[
 \begin{aligned}
   f(N-t)-f(N-2t)&\le\log 2-b_n s \qquad(s<1/2),\\
   f(t)-f(0)&\le-\log(1-s)-b_n s .
 \end{aligned}
\]
Its derivative is everywhere at least \(-B\),
since \(G'\) is decreasing and \(x\le N\).
For \(s<1/2\) this gives
\(f(2t)-f(t)\ge\log((1-s)/(1-2s))-NB/2\).
These estimates imply (3) over (4). Since
\(b_n/2\ge10/13>\log 2\),
the quantities in (5) satisfy the required strict bounds,
and \(-\log(1-s)/s\) increases from \(1\) to infinity.

For (6), the lower bound at the origin uses the
right-sum error only on increasing rows and
\(d_i/a_i\le d_i^2/a_i\le13 d_i^2/\max(a_i,b_i)\).
Also the exact difference from \(G'(0)\) to \(G'(N)\) is
\[
 -N\sum_{i:a_i>0}d_i^2
           \sum_{h=1}^{a_i}
                \frac1{(na_i+h)(N b_i-a_i+h)} .
\]
Its magnitude is at most \(13 N S/[n(n-12)]\).
For increasing rows each denominator contains two factors
at least \(na_i\). For decreasing rows the second factor
is at least \((n-12)b_i\), and \(a_i/b_i\le13\).
Adding \(\log\kappa-\tau/n\), and using
\(13S+|\lambda|\le156\) and \(S+\tau\le12\), gives (6).

The exact shape (2) also applies when the signed degree
condition for (4) fails. The collar does not exclude other
angles or all signed degree differences, and it does not
settle the consecutive coprimality problem.

As an exact instance of the second-sector test, use the mixed
shared rows \((30,20),(62,70)\), \(n=50,N=51\) and the
ordinary physical scale. Their curvature is \(446/105\) and
their total degree difference is \(-2\). At \(t=25\) the
ratios \(e^{f(N-t)-f(N-2t)}\) and \(e^{f(2t)-f(t)}\) lie in
\((1/3,2/5)\) and \((13/10,4/3)\), respectively, by direct
products of the factors in (1). Thus (3) excludes the angle
\(\pi/25>2\pi/51\), with a neighborhood by strictness.
