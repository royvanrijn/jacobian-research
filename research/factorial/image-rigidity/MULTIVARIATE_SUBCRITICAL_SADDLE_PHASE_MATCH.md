# Subcritical equal-degree saddle phase matches

For every integer \(m\ge3\) there is a mirrored two-row support with
\(a=(m,m+1)\), \(b=(m+1,m)\) and two distinct regular upper saddles at
the same target \(X\in\mathcal T\) with equal \(e^{F_a}\). Here \(X,F_a,\mathcal T\)
are as in the [regular saddle comparison](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md).
Writing \(M=2m+1\), the total degree difference is zero, the
[total curvature](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md) is \(S=2/(m+1)\),
and the stationary curvature satisfies
\[
 \rho(q)=\frac{M q(1-q)}{m(m+1)+q(1-q)},\qquad
 \max_{0\le q\le1}\rho(q)=1/M.
\]
These matches therefore persist at arbitrarily small total curvature.

Set \(\omega=2q-1\). The exponential of the stationary target and
the exponential of its phase are
\[
 R(\omega)=e^{X(q)}
   =\frac{(M+\omega)(1-\omega)}{(M-\omega)(1+\omega)},\qquad
 \Psi(q)=e^{F_a(q)}
   =2^{-2m}\frac{(M+\omega)^m(M-\omega)^{m+1}}{1-\omega}.
\]
Since
\[
 \frac{1-R}{1+R}=\frac{(M-1)\omega}{M-\omega^2},\qquad
 \Im\frac{\omega}{M-\omega^2}
   =\frac{\Im\omega\,(M+|\omega|^2)}{|M-\omega^2|^2},
\]
the image of the upper half-plane under \(R\) lies in the lower
half-plane. The principal-log expression for \(X\) is continuous
and has small negative imaginary part near \(q=1/2+i\varepsilon\)
for small positive \(\varepsilon\), since \(X'(1/2)<0\). Consequently
\(-\pi<\Im X<0\) throughout the upper half-plane. The involution
\(\omega'=-M/\omega\) preserves the upper half-plane and leaves
\(R\), hence \(X\), unchanged. The only upper ramification point
is its fixed point \(\omega^2=-M\), as follows by differentiating
the rational Cayley transform above.

Put \(Y=\omega^2\). Substituting \(\omega'\) into \(\Psi\) shows that
the two phases match precisely when
\[
 P(Y)=Y^m(M^2-Y)^{m+1}-(-1)^m M^{2m+1}(1-Y)^{m+1}=0.
\]
For \(Y=-z<0\), put \(p=m/(m+1)\). The equation reduces to
\(h(z)=M^{p+1}\) with \(h(z)=z^p(M^2+z)/(1+z)\). Its logarithmic
derivative is
\[
 \frac{h'(z)}{h(z)}=\frac{p(z-M)^2}{z(1+z)(M^2+z)}.
\]
Thus the only negative root is \(-M\) with multiplicity three,
since \(h(M)=M^{p+1}\) and \(h'\) vanishes to order two there.
For odd \(m\) there are no positive roots by the signs of the
two terms. For even \(m\) there are no roots in \((1,M^2)\) by
the same sign comparison. There is one simple root in each
of \((0,1)\) and \((M^2,\infty)\), because the respective functions
\[
 Y^p\frac{M^2-Y}{1-Y},\qquad Y^p\frac{Y-M^2}{Y-1}
\]
have positive derivative and increase from zero to infinity
on their stated intervals. The points \(0,1,M^2\) are never roots.
The degree is \(2m+1\), so a nonreal root exists for every \(m\ge3\).
Its upper square root \(\omega\) and the point \(-M/\omega\) give
distinct regular saddles with the required equalities.

## A match inside the conditional integration strip

In this family the
[conditional coordinate](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md)
satisfies \(e^{s(q)}=(M+\omega)/(M-\omega)\). It maps the upper
\(\omega\)-half-plane into the upper half-plane and the
logarithmic lift has \(0<\Im s(q)<\pi\). The real part of the
ratio has the sign of \(M^2-|\omega|^2\). Both saddles of a matched
fiber therefore lie in \(0<\Im s<\pi/2\) exactly when
\(1<|Y|<M^2\).

For \(m=9\), divide \(P\) by \(19^{19}\) to obtain
\(Q(Y)=19Y^9(1-Y/361)^{10}+(1-Y)^{10}\). Let
\(c=(6555+105165i)/100000\), \(r=1/1000\) and
\(Q(c+rZ)=\sum_{k=0}^{19}T_kZ^k\). The following strict rational
bounds use \(\|z\|_1=|\Re z|+|\Im z|\).
\[
 \|T_0\|_1<1/500,\qquad
 \sum_{k=2}^{19}\|T_k\|_1<1/500,\qquad
 \Re T_1>1/10.
\]
For direct verification use
\(T_k=r^k\sum_{j=k}^{19}[Y^j]Q(Y)\binom jk c^{j-k}\),
with the coefficients of \(Q\) given by the two binomial
expansions above. Rouché's theorem on \(|Z|=1\) gives one
root in \(|Y-c|<r\). This disk has \(\Im Y>1\) and
\(|Y|<|\Re Y|+|\Im Y|<2\), so it lies within the required
annulus.

These equalities concern candidate saddle phases. The
moment polynomials of the same support strictly interlace
at every consecutive order by
[the single mirrored unit-gap theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#consecutive-interlacing-for-a-single-mirrored-unit-gap).
A saddle argument for the subcritical class needs information
beyond target, exponential phase, curvature and strip membership.

## The squared Gaussian amplitudes are separated

The
[regular contour comparison](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#consequences-away-from-the-real-cusp)
requires the ratio of the unweighted squared amplitudes at two
active saddles to be positive real if the contour-weighted
squared amplitudes agree. For the unit-gap family above the
unweighted invariant and its ratio on a regular upper fiber are
\[
 H(q)=\frac{x_1(q)x_2(q)}{1-\rho(q)}
     =\frac{(M^2-Y)^2}{4(M-1)(M+Y)},\qquad
 \frac{H(q')}{H(q)}
     =\frac{M^3(Y-1)^2}{Y(M^2-Y)^2},
 \qquad Y=(2q-1)^2.
\]
The ratio follows by replacing \(Y\) with \(M^2/Y\) at \(q'\).

Suppose the ratio were \(\kappa>0\) and set
\(T=M(Y-1)/(M^2-Y)\). Then
\[
 \frac YM=\frac{T^2}{\kappa}=\frac{MT+1}{M+T},
 \qquad T^3+MT^2-\kappa MT-\kappa=0 .
\]
The last cubic has three distinct real roots, since its signs
at \(-\infty,-M,0,+\infty\) are respectively negative,
positive, negative and positive. Thus \(T\) and \(Y\) are real.
Since \(2q-1\) lies in the upper half-plane, this forces
\(Y<0\), where the displayed amplitude ratio is negative.
This contradiction holds on every regular upper fiber in the
family, including all the exact phase matches constructed above.
