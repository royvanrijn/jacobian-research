# Regular fields in the subcritical thermal approach to a shared cusp

Use the positive bounded \(W\)-band, \(\delta=\max\rho-1\),
and the centered quantities \(F,q_*,w_*,a_2,g\) from the
[uniform moving quartic proof](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
Here the regime is \(\delta\to0^-\) with
\(\sqrt n\,|\delta|\to\infty\). Let
\[
 r^2=a_2/(2g)>0,\qquad c=-4g r^4>0,qquad
 \beta_j=w_*-w_j,\qquad b_j=r\beta_j/c,qquad M_j=jc .
\]
These scales use the same profile at both adjacent orders.

**Proposition.** Along any consecutive common-root sequence
in this regime for which one of the centered parameters
tends to \(w_*\) in difference, both \(b_n,b_{n+1}\) must tend,
up to subsequences, to the same one of
\[
                      b=\pm i\frac{2}{3\sqrt3}.
\]
This localizes the possible common roots to the two simple
fold windows of the scaled phase. The conclusion places no
bound on the total-degree difference \(\Delta\).

## Reduction to bounded thermal fields

For each approaching profile
\(q_*-q_m=O(|\delta|)\) at the curvature maximum \(q_m\),
so \(\rho(q_*)-1=\delta+O(\delta^2)\). Hence \(-a_2\)
is positive and comparable to \(|\delta|\). In particular
\(r\to0\), \(M_j\to\infty\) and \(M_j=o(j)\).
Extract a retained and diffuse limit. The other centered
parameter also tends to the limiting real cusp: the
[uniform root annulus](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md)
bounds it, its imaginary part is shared, and a noncusp
real limit is excluded by the uniform finite-interval
one-saddle estimate in the moving quartic proof.
Its finite integral representation with a fixed physical
gap on the exterior thus applies at both orders.

Suppose \(|b_j|\to\infty\) at one order. With
\(\epsilon_j=|\beta_j|^{1/3}\) we have
\(\epsilon_j/r\to\infty\), \(a_2/\epsilon_j^2\to0\)
and \(j\epsilon_j^4\to\infty\). Therefore the
[escaping-field argument](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md#escaping-fields-at-bounded-scaled-temperature)
applies in its field scaling. Its one-saddle estimate
forces the argument of such a root's field to tend to
\(\pi/2\) or \(-\pi/2\). A bounded \(b_j\) at the other
order would then contradict the equality of their physical
imaginary parts. If both escape, both have the same
imaginary direction, and the adjacent comparison in
that argument excludes their simultaneous cancellation.
That comparison uses only the vanishing scaled quadratic
coefficient, the diverging saddle scale, and the real
shift \(\beta_{n+1}-\beta_n=\Delta\log(1+1/n)\).
It therefore also covers the present thermal regime.

For bounded \(b_j\), set \(q=q_*+rz\). After centering
and division by \(c\), the phase becomes
\[
        \Psi_j(z)=-z^4/4-z^2/2+b_jz+O(r|z|^5).
 \tag{1}
\]
On fixed scaled disks this approximation holds with
derivatives. For bounded scaled heights, large real
parts and \(|rz|\) in a sufficiently small fixed physical
disk, the real part is at most \(-c_1|\Re z|^4\).
Indeed the remainder in (1) is then bounded by a small
multiple of \(|z|^4\).

## The regular thermal contours

Consider \(G_b(z)=-z^4/4-z^2/2+bz\) with \(b=B+iY\).
For \(B>0\) and any real \(v\), the equation
\[
 u^3+(1-3v^2)u=B,qquad u>\sqrt{\max(0,3v^2-1)}
\]
has exactly one solution. The derivative in \(u\) of
its left side is \(A=3u^2+1-3v^2>0\).
Implicit differentiation gives \(u'(v)=6uv/A\).
Thus the function \(v\mapsto v(3u(v)^2+1-v^2)\)
is strictly increasing, with derivative
\(A+6uv u'(v)>0\), and tends to opposite infinities.
There is exactly one choice of \(v\) giving the
imaginary equation \(Y=v(3u^2+1-v^2)\). This makes
\(z=u+iv\) a stationary point. On its horizontal the
real phase is a constant plus
\[
              -x^4/4+(3v^2-1)x^2/2+B x.
\]
For positive \(x\) this is strictly larger than at
the reflected point. Its unique positive stationary
point is \(u\), with negative curvature \(-A\).
It is consequently the unique strict maximum on that
horizontal. Reflection gives the corresponding result
for \(B<0\).

Write \(Y_e=2/(3\sqrt3)\). For \(B=0\) and
\(0\le Y<Y_e\), solve \(Y=v-v^3\) with
\(0\le v<1/\sqrt3\). On the horizontal through
the simple saddle \(iv\) the real phase has its
unique maximum there, with curvature \(3v^2-1<0\).
Conjugation covers \(-Y_e<Y<0\).
At all these one-saddle parameters, the strict real
loss and nondegenerate saddle persist on the shifted
actual contour for (1), as do the quartic tails.
The connectors at the fixed physical endpoints keep
the fixed phase gap of the finite local integral.
A uniform Gaussian comparison in a neighborhood of
each such limiting parameter makes the individual
polynomial nonzero at all sufficiently large orders.

For \(B=0\), \(Y>Y_e\), put
\(u=\sqrt{3v^2-1}\) and solve
\(Y=8v^3-2v\) with \(v>1/\sqrt3\).
The two saddles are \(z_R=u+iv\), \(z_L=-u+iv\).
On their horizontal each half-line has a unique strict
maximum with negative curvature there. The midpoint
is lower than both maxima by \(u^4/4\).
Fix a bounded limiting \(iY\) with \(Y>Y_e\).
In a small neighborhood of it, shift and split the
finite contour at the midline and then move each
half to its actual saddle height. The short midpoint
connectors keep a strict scaled gap and the fixed
physical connectors retain their original gap. This
gives two actual Gaussian packets with separate
nonzero amplitudes, uniformly bounded above and away
from zero after extracting their stationary
exponentials and the real length factor \(rM_j^{-1/2}\).
The relative errors in this approximation are
\(O(M_j^{-1/2})\) by the local quadratic comparison.
The polynomial remainder from the fixed physical
exterior is exponentially negligible relative to
either centered packet, since \(M_j=o(j)\) and both
stationary scaled heights are bounded.
The case \(Y<-Y_e\) is its conjugate.

These estimates show that any bounded pair of
scaled common-root fields must limit to the same
pure imaginary point with \(|Y|\ge Y_e\), since their
imaginary parts agree. It remains in this proposition
to exclude the strictly larger arguments \(|Y|>Y_e\).

## The regular two-packet comparison

By conjugation take the limiting point \(iY\) with
\(Y>Y_e\), and put \(M=M_n\). Let \(D(\beta)\) be the
right-minus-left physical stationary phase difference
for the same profile. Modulus cancellation at both
orders gives \(\Re D(\beta_j)=O(n^{-1})\).
On the short real segment between the two fields in
their scaled coordinate, stationarity gives
\[
 D_\beta=q_R-q_L,qquad
            \Re(q_R-q_L)\asymp r .
\]
Integrating the real difference on the
segment of length \(|\Delta|\log(1+1/n)\) therefore gives
\(|\Delta|r=O(1)\).

Keep the lower-order contours at the upper order.
As in the [adjacent packet comparison](MULTIVARIATE_SHARED_SUPERCRITICAL_THERMAL_APPROACH.md#adjacent-comparison-and-small-imaginary-fields),
the centered multiplier is the amplitude ratio times the
exponential of
\[
 \Phi_{w_n}(q)-\Phi_{w_n}(q_*)
           +\Delta(n+1)\log(1+1/n)(q-q_*).
\]
Locally the exponent differs from its value at the
saddle by \(O(c|z-z_j|^2+|\Delta|r|z-z_j|)\), using
exact stationarity. The packet ratio thus equals its
stationary multiplier times one plus
\[
                 O(c/M+|\Delta|r/\sqrt M+n^{-2}).
]
The two Gaussian moments give the first two terms,
and the analytic order derivative of the amplitude
gives the third. The strict nonlocal gaps and the
quartic majorant justify the same bound on the
entire contour: \(|\Delta|r=O(1)\) makes its
additional scaled growth at most exponential in
a linear term and in \(c|z|^4\). The fixed physical
connectors have a gap that absorbs \(e^{O(|\Delta|+1)}\).

For \(\Delta\ne0\), the real difference of the two
stationary multiplier exponents has magnitude
comparable to \(|\Delta|r\), since
\(\Re D(\beta_n)=O(n^{-1})\). This exceeds the ratio
error by a factor tending to infinity. For \(\Delta=0\),
the imaginary difference of their exponents has
magnitude comparable to \(c\).
In fact the limiting stationary difference of
\(G_b\) vanishes at \(iY_e\), and the derivative of its
imaginary part along \(iY\), \(Y>Y_e\), is \(2u>0\).
Since \(Y\) is separated from \(Y_e\), it supplies the
positive margin after rescaling by \(c\).
Again the ratio error is smaller.
At a common root both packets would cancel at each
order, forcing their ratios to agree up to the
exponentially small physical remainders. This is
impossible in either case and finishes the localization.
\(\square\)
