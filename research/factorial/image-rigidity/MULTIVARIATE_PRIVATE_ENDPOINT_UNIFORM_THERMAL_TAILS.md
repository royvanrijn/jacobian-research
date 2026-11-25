# Uniform thermal tails at the one-private linear endpoint

Orient the ordinary exponent columns so the only private row
is \((0,1)\). Use the shared rows \(a_i,b_i>0\), gaps
\(d_i=b_i-a_i\), curvature \(S\), endpoint curvature \(\sigma\),
and fields \(\xi_j\) from the
[moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
Set \(N=n+1\) and \(\mu=1-\sigma\).

**Theorem.** Fix an upper bound on \(S\). Suppose the profiles vary
along a sequence \(n\to\infty\), \(\xi_n,\xi_N\to0\), and
\(\sigma\) tends to a limit in \([0,1]\). If these are consecutive
common roots, then
\[
 \sigma\to1,\qquad n^{1/3}|\mu|=O(1),\qquad
 n^{2/3}(|\xi_n|+|\xi_N|)=O(1).
\]
The result allows every integer total-degree difference.
In particular it gives a uniform exclusion when both fields lie
in a sufficiently small fixed disk on a strictly subcritical margin.

## Linked moments and rescaled phases

Write the normalized moment at integer \(j\) as
\[
 I_j(\xi)=\int_0^L A_j(q)e^{j\Phi(q,\xi)}\,dq+B_j(\xi)+O(e^{-\eta j}),
 \qquad
 \Phi(q,\xi)=\xi q-\mu q^2/2-c_2q^3/3+q^4E(q),
 \tag{1}
\]
with the variables and normalization in equation (5) of the
[moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md#from-the-lattice-to-a-linked-analytic-integral).
On fixed disks the analytic \(A_j\) has \(A_j(0)=1\),
bounded derivatives, and \(A_N-A_n=O(n^{-2})\).
The coefficient \(c_2\) is real and stays above a positive
bound, and \(E\) is uniformly bounded and analytic.
At limits strictly below one the same representation holds.
Indeed \(2c_2\ge1-R\sigma\) with a fixed \(R<1\) as in that
proof, and \(R(q)=(1-q)\sum d_i^2/(a_i+d_iq)\) decreases.
Thus \(\Phi(q,0)<0\) throughout \(q>0\) and there is uniform
loss beyond \(L/2\) on each strictly subcritical margin.
The cotangent argument uses this loss and the same fixed
analytic neighborhoods. Its sectorial estimate also gives
\[
 |B_j(\xi)|=O(n^{-1}),\qquad
 |\partial_\xi B_j|=O(n^{-1}),\qquad
 |\partial_j B_j|=O(n^{-2})
 \tag{2}
\]
on smaller disks, by Cauchy on order disks of radius
proportional to \(n\). The even-order moment is positive
at real \(t\), so by conjugation it suffices to take
\(\Im\xi_n=\Im\xi_N>0\).

Consider a scale \(q=r z\) with \(r\to0\).
Put \(b=r\xi/c\), \(M=nc\), and
\(\Psi_b(z)=\Phi(rz,\xi)/c\), where \(c>0\).
Three cases with \(M\to\infty\) are used below.

| Scale | \(r\) and \(c\) | Limiting phase |
| --- | --- | --- |
| Cubic field | \(r=\sqrt{|\xi|/c_2}\), \(c=c_2r^3\), \(|\mu|/r\to0\) | \(bz-z^3/3\), \(|b|\to1\) |
| Thermal field | \(r=|\mu|/c_2\), \(c=c_2r^3\) | \(bz-a z^2/2-z^3/3\), \(a=\operatorname{sgn}\mu\), \(b\) bounded |
| Quadratic field | \(r=|\xi|/\mu\), \(c=\mu r^2\), \(\mu>0\), \(r/\mu\to0\) | \(bz-z^2/2\), \(|b|\to1\) |

The convergence is analytic on compact \(z,b\) sets.
A convenient feature is that at a fixed profile the phase
depends on \(b\) only through the term \(bz\).
In each scale \(c=O(r^2)\). Bounded horizontal shifts in
the \(z\) coordinate are allowed in (1).
To check uniform outer tail bounds, the real parts of
\(z^2,z^3\) on such lines are at least half the corresponding
powers of \(x=\Re z\) once \(x\) is large. The term
\(r^4z^4E(rz)\) is at most a small multiple of \(c_2r^3x^3\)
for \(r x\le L\), by decreasing \(L\).
In the quadratic field both the quadratic and cubic have
negative signs on these tails. In the two cubic scales
the quadratic coefficient is bounded and the cubic absorbs
it for large \(x\). Thus \(\Re\Psi_b\le-C x^p\) for large positive
\(x\), with \(p=2\) in the quadratic field and \(p=3\)
in the two cubic scales. The changes at \(q=L\) cost
\(O(e^{-\eta' n})\) by the fixed real loss there.

## The contours and the limiting equal-modulus fields

First take a limiting cubic phase
\(P(z)=bz-a z^2/2-z^3/3\), with \(a\in\{-1,0,1\}\)
and \(\Im b>0\). Write its right saddle
\(z_s=(-a+\sqrt{a^2+4b})/2=x+iy\) with \(y>0\)
and \(2x+a>0\). Then
\[
 \Re b=x^2+ax-y^2,\qquad \Im b=(2x+a)y,\qquad
 D=P(z_s)=az_s^2/2+2z_s^3/3 .
 \tag{3}
\]
When \(x>0\), deform the real interval into the vertical leg
from zero to \(iy\) and the rightward horizontal line through
\(z_s\). The real phase on the leg at \(iv\) is
\(-(\Im b)v+a v^2/2\). It is strictly negative for \(0<v\le y\):
for \(a\le0\) both terms are nonpositive, and for \(a=1\)
use \(\Im b=(2x+1)y>y\).
The horizontal line starts at strict negative height. Its
real derivative at \(s+iy\) is
\[
                       (x-s)(x+s+a),\qquad s\ge0 .
\]
Consequently any horizontal maximum is at the starting
point or the saddle \(s=x\). Its Hessian at the saddle
has real part \(-2x-a<0\). If \(\Re D\) is near zero, the
saddle is the unique maximum with a strict loss elsewhere
outside its neighborhood. These strict bounds persist
on the nearby actual contours through the continued
saddle. If \(\Re D<0\), the horizontal part is exponentially
small on the scale \(M\). If \(\Re D>0\), the continued
saddle dominates the endpoint.

For \(x\le0\), the condition \(2x+a>0\) forces \(a=1\)
and \(-1/2<x\le0\). Choose any fixed
\(0<v<\min(y,\Im b)\) and shift to height \(v\).
On the vertical leg there is strict linear loss.
The horizontal phase starts strictly below zero and
has real derivative
\(x^2+x-y^2+v^2-s-s^2<0\).
Thus only the nonzero endpoint term contributes there.

On the real \(b\)-axis, the undeformed real line gives the
following alternatives. For \(a\in\{0,1\}\), the endpoint is
strict and dominates at \(b<0\), and a unique nondegenerate
interior maximum dominates at \(b>0\). The point \(b=0\)
is excluded by \(|b|\to1\) in the cubic field and is treated
separately below for \(a=1\). For \(a=-1\) the endpoint
dominates below \(b=-3/16\), the unique interior maximum
dominates above it, and their heights agree there.
Indeed the only possible strictly positive local maximum
is at \(x>1/2\), with height \(x^2(2x/3-1/2)\).
At coexistence \(x=3/4\), and the endpoint derivative
\(b=-3/16\) is strictly negative. Use a small fixed real
segment from zero for the endpoint. The right contour
starts at that segment's endpoint and goes through the
continued saddle, with both its start and far tail lying
at strict loss. These contour statements on the real
axis persist for nearby complex fields by the strict
losses and local deformations at the maximum.

For the quadratic limit, the real-axis alternatives are
the same as for \(a=1\), with the origin excluded by
\(|b|\to1\). For a fixed nonreal \(b=x+iy\) with \(x>0\)
shift to height \(y\). The vertical leg has negative
real phase \(-yv+v^2/2\) for \(0<v\le y\), and the
horizontal real phase has its unique maximum at \(x\).
If \(x\le0\), the horizontal real part instead decreases
from its starting point, which is strictly below zero.
The quadratic stationary height when \(x>0\) is \(D=b^2/2\).
Thus a quadratic-field zero limit has \(\arg b=\pi/4\).
In the cubic-field limit \(a=0\), (3) similarly gives
\(\arg b=\pi/3\) at a zero.

At every regular equality just described the limiting
saddle has positive real coordinate and the endpoint
has nonzero phase derivative. At a fixed positive
imaginary height there is at most one limiting equality.
In fact \(D'(b)=z_s\) by stationarity, so its real part
increases strictly with \(\Re b\) on the candidate region
\(x>0\). That region is an interval at fixed height, since
\(\Re b=x^2+ax-(\Im b)^2/(2x+a)^2\) increases with \(x\)
when \(2x+a>0\). At height zero the only thermal equality
for \(a=-1\) is the coexistence just described, and for
\(a=1\) the only possible root limit is \(b=0\).

The contour assertions give these localizations for
actual zeros. Specifically, a fixed first leg with strict
negative linear slope at zero has a nonzero term of size
\(r/M\), by the substitution \(z=v/M\) and domination by
\(e^{-c|v|}\). Adding \(B_j=O(1/n)\) changes this term by a
relative \(O(c/r)=o(1)\). An active saddle with a strict
quadratic maximum gives a term of size \(r/\sqrt M\)
times its stationary exponential. This follows by the
local Gaussian substitution and the uniform outer losses
proved above. The upper and lower bounds on each leading
amplitude are uniform on a fixed neighborhood of each
limiting point under consideration. Thus off an equality
the exponentially smaller contribution cannot cancel
the larger one.

## The adjacent comparison at a regular equality

Suppose both orders have fields tending to the same regular
equality on one common scale \(r,c\). Write \(b=b_n\),
\(b'=b_N\), \(M=nc\), and \(M'=Nc\).
The actual saddle \(z_s(b)\) and height \(D(b)=\Psi_b(z_s(b))\)
are analytic near the limiting point and depend only on
the profile and the chosen scale. In particular \(D'(b)=z_s(b)\).
For each of the two orders the local contour decomposition,
including \(B_j\) in the endpoint term, has the form
\[
 E_j(b)=\frac{r}{jc}\,e_j(b)(1+o(1)),\qquad
 R_j(b)=\frac{r}{\sqrt{jc}}\,s_j(b)e^{jcD(b)}(1+o(1)).
 \tag{4}
\]
Here \(E_j\) and \(R_j\) are the exact local integrals with
amplitude \(A_j\). The coefficients \(e_j,s_j\) in (4) may
be taken independent of order for a given profile.
Use the endpoint slope and the leading saddle amplitude
\(A_0(rz_s)\) divided by the square root of the negative
Hessian along the horizontal contour. These functions are
bounded, bounded away from zero, and analytic with bounded
derivatives on a smaller field disk. The discarded
physical tails are \(O(e^{-\eta'n})\).
Both scales \(r\) and \(r/M\) are bounded below by a
negative power of \(n\), since \(M\to\infty\),
\(c=O(r^2)\), and \(M=nc\). Thus those tails are negligible
even after division by the endpoint term.

At an actual zero, (4) gives
\[
 M\Re D(b)=-\tfrac12\log M+\log|e_n(b)/s_n(b)|+o(1),
 \tag{5}
\]
with its analogue at the next order. Subtracting the two
relations gives \(M'\Re D(b')-M\Re D(b)=o(1)\), since
the leading coefficients are linked, \(b'-b\to0\), and
\(M'/M\to1\). The difference \(b'-b\) is real, and
\(\Re D'(b)\) stays above a fixed positive bound.
With \(c\Re D(b)=o(1)\) it follows that
\[
 d=M'(b'-b)=r N\Delta\log(N/n)=o(1).
 \tag{6}
\]
This preliminary comparison handles arbitrarily large
possible total-degree differences before any expansion
using the increment \(d\).

Retain the lower-order contours for the upper-order
integrals. The integrands, apart from the amplitude ratio
\(1+O(n^{-2})\), differ by the exact factor
\(\exp(c\Psi_b(z)+d z)\).
The endpoint's linear decay makes the error from this
factor \(O((c+|d|)/M)\) relative to its integral.
For the additive correction \(B_j\), (2) gives a
difference of size \(O(n^{-2}+n^{-1}|\xi_N-\xi_n|)\).
After division by the endpoint integral this is
\(O(c/(rn)+c|d|/(r^2n))\), because
\(|\xi_N-\xi_n|=|d|/(rN)\). Thus
\[
 \frac{E_N(b')}{E_n(b)}
 =1+O\left(\frac{c+|d|}{M}+n^{-2}
               +\frac{c}{rn}+\frac{c|d|}{r^2n}\right).
 \tag{7}
\]
The Gaussian width at the saddle is \(O(M^{-1/2})\).
At \(z=z_s(b)\), the part \(c\Psi_b(z)\) is stationary,
so
\[
 \frac{R_N(b')}{R_n(b)}
 =e^{cD(b)+d z_s(b)}
     \left(1+O\left(\frac{c}{M}
                    +\frac{|d|}{\sqrt M}+n^{-2}\right)\right).
 \tag{8}
\]
To justify both error bounds, subtract the respective
multiplier values and use
\(|e^u-1|\le |u|e^{|u|}\).
On the endpoint ray the exponent increment is
\(O((c+|d|)|z|)\). Near the saddle its difference from the
central value is \(O(c|z-z_s|^2+|d||z-z_s|)\).
The indicated linear or Gaussian absolute-integral
bounds have the same scales as the leading moduli in
(4). On bounded complementary pieces, the strict
scaled gap absorbs the increment since \(c,|d|\to0\).
On unbounded horizontal pieces, the original tail
has loss at least a positive multiple of \(M x^3\) in the
cubic scales. In the quadratic field the loss is at least
a positive multiple of \(M(x^2+(c_2r/\mu)x^3)\). This extra loss absorbs the same cubic
part in \(c\Psi_b\), since \(c/M=1/n\).
Polynomial factors from the difference are then
integrable with exponentially small relative error.

If both orders vanished, the ratios in (7) and (8) would
agree up to errors exponentially small in \(n\).
For \(\Delta\ne0\), (6) gives \(|d|\ge r\).
The logarithm of the modulus in (8) has leading term
\(d\,\Re z_s(b)\), since \(c=o(r)\).
Every error in (7) or (8) is \(o(|d|)\), using
\(c/(r^2n)=O(1/n)\). This contradicts \(\Re z_s(b)>c_0>0\).
For \(\Delta=0\), (7) has error \(O(1/n)\).
The logarithm of the modulus in (8), by (5), is
\[
                    -\frac{\log M}{2n}+O(1/n).
\]
Again the two ratios cannot agree.

## The bounded quadratic window

It remains to treat \(\mu>0\), \(n\mu^3\to\infty\), and both
\(H_j=\sqrt{j/\mu}\,\xi_j\) bounded. Set
\(r_n=(n\mu)^{-1/2}\). Use \(q=r_n y\) in the analytic
interpolation (1) for \(|j-n|\le\varepsilon n\).
For \(H=\sqrt{j/\mu}\,\xi\) the first two terms of
the rescaled exponent are
\(\sqrt{j/n}\,H y-(j/n)y^2/2\).
The cubic coefficient \(j c_2 r_n^3\) tends uniformly
to zero. For real \(y\ge0\) the cubic has nonpositive
real contribution, and the term \(q^4 E(q)\) is absorbed
by half its loss for \(q\le L\).
Thus \(\exp(-c_0 y^2+C|y|)\) dominates the integrals
uniformly on the complex order disk and bounded
\(H\) sets. The amplitude tends to one on compact
\(y\)-sets, and the normalized \(B_j\) tends to zero by
\(\sqrt{j\mu}\,B_j=O(\sqrt{\mu/n})\).
The normalized interpolation converges uniformly to
\[
                     G(H)=\int_0^\infty e^{-y^2/2+Hy}\,dy .
\]
For the complex square root used in the normalization,
the Gaussian model integral has this value after
multiplication by \(\sqrt{j/n}\); this follows for real
\(j\) by rescaling and on the disk by holomorphy.

The zeros of \(G\) are simple, since \(G'(H)=1+HG(H)\)
by integration by parts. Moreover they have distinct
imaginary coordinates. For \(\Re H>0\), Fourier inversion
of the Gaussian and integration of \(e^{-Hy}\) give
\[
 G(-H)=\frac1{\sqrt{2\pi}}\int_{\mathbb R}
                  \frac{e^{-u^2/2}}{H-iu}\,du .
\]
It has strictly positive real part there. At \(\Re H=0\)
the real part of \(G(H)\) is a positive Gaussian
by the cosine transform. These facts put every zero
in \(\Re H>0\). The identity \(G(H)+G(-H)=
\sqrt{2\pi}e^{H^2/2}\) follows by the whole-line Gaussian
integral. At a zero it forces the ratio of the two
moduli \(|G(-H)|/|e^{H^2/2}|\) to equal \(\sqrt{2\pi}\).
At fixed \(\Im H\) its logarithmic \(x=\Re H>0\)
derivative is \(-\Re(1/G(-H))<0\), by the differential
identity above. Hence each height contains at most one zero.

At a proposed common root the two bounded rescaled fields
have the same limiting height, since their physical
imaginary coordinates agree and \(N/n\to1\).
They therefore approach the same simple zero \(H_0\)
of \(G\). It is nonreal, since \(G\) is positive
on the real line. Rouché gives a single analytic
zero branch \(H(j)=H_0+o(1)\) uniformly on the complex
order disk. Cauchy gives \(H'(j)=o(1/n)\) on a smaller
disk. In physical coordinates this yields
\[
 \xi'(j)=\sqrt{\mu}\,j^{-3/2}(-H_0/2+o(1)).
\]
The adjacent imaginary displacement has strict size
\(\sqrt{\mu}\,n^{-3/2}\). It exceeds the exponentially
small integer-order errors in (1), by \(n\mu^3\to\infty\).
This excludes a common root because the exact drift is real.

## Exhaustion of the endpoint scales

If at either order both
\(|\xi|/\mu^2\to\infty\) and \(n^{2/3}|\xi|\to\infty\),
use the cubic field scale, with the first condition read
as automatic when \(\mu=0\). Its zero angle tends to
\(\pi/3\), so its imaginary part is comparable to its
modulus. The next field has the same imaginary part
and hence belongs to the same escaping regime. By its own
angle limit its modulus is asymptotic to the first one.
On the common scale the two fields reach the same
regular limiting equality, already excluded.

If \(n^{1/3}|\mu|\) is unbounded and \(\sigma\to1\),
pass to \(n^{1/3}|\mu|\to\infty\) with fixed sign of \(\mu\).
The previous case rules out unbounded \(|\xi|/\mu^2\)
at either order. Both thermal fields are now bounded
and \(M\to\infty\). Their imaginary parts are the same.
At positive limiting height the regular equality is
unique, and for \(\mu<0\) this also holds on the real
boundary. Thus the only remaining thermal limit has
\(\mu>0\) and both \(b,b'\to0\).

This last limit gives \(|\xi_j|/\mu^2\to0\) and
\(n\mu^3\to\infty\), as it also does at a strictly
subcritical limit of \(\sigma\). If at either order
\(\sqrt{n/\mu}|\xi|\to\infty\), use the quadratic field
scale. Its zero angle tends to \(\pi/4\). The common
imaginary coordinate ensures the other order belongs
to the same escaping regime and their moduli are
asymptotically equal. The comparison on one common
regular scale applies. Otherwise both quadratic fields
stay bounded, and the bounded quadratic window gives
the contradiction. These cases leave only the bounded
moving cubic fields in the theorem.
