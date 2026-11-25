# Effective bounded supercritical thermal comparison at the shared cusp

Delete equal rows. Use the ordinary fully shared rows, chord-centered phase \(F\),
curvature maximum \(\rho(q_m)=1+\delta\), and real center
\(F'''(q_*)=0\) of the
[moving shared quartic proof](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
Write \(w_*=F'(q_*)\), \(a_2=F''(q_*)/2\), \(g=F''''(q_*)/24\), and
for a common physical coefficient ratio put \(\beta_j=w_*-w_j\).
The fields for \(n,N=n+1\) have the exact relation
\(\beta_N-\beta_n=\Delta\log(N/n)\).
The logarithms are taken on the branches near the positive real
scales, as is automatic for the bounded fields below when the
curvature cap is sufficiently small.

**Theorem.** Given positive rational \(W_-\le W_+\) and rational
\(B\ge1\), one can compute \(\delta_0>0\), \(M_0\ge1\), and an
integer \(n_0\) with the following property. For \(W_-\le W\le W_+\)
and \(0<\delta\le\delta_0\), put
\[
 \ell=(-4g)^{-1/3},\qquad r=(2a_2\ell)^{1/2},\qquad
 A=\ell r^4,\qquad b_j=\beta_j/r^3,\qquad M_j=jA.
\]
Then \(r>0\), and the orders \(n,N\) have no common physical root
with \(|b_n|,|b_N|\le B\) if \(n\ge n_0\) and \(M_n\ge M_0\).
The constants are uniform in the exponent rows and the integer
\(\Delta\). One can equivalently use a computable lower bound
on \(\sqrt n\,\delta\) in place of the lower bound on \(M_n\).
The result allows the imaginary parts of the bounded fields to
approach zero.

## Effective physical reduction

The [effective center bounds](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md#effective-center-and-localization-on-the-real-lattice)
give a fixed physical disk, \(q_*-q_m=O(\delta)\), uniform
bounds on the needed derivatives, and positive upper and lower
bounds for \(-g\). Taylor's theorem at \(q_m\) gives
\(|\rho(q_*)-(1+\delta)|\le C\delta^2\). Decrease the cap on
\(\delta\) so that \(\rho(q_*)-1\ge\delta/2\).
Using \(2a_2=(\rho(q_*)-1)/(q_*(1-q_*))\) then bounds \(a_2\)
above and below by positive computable multiples of \(\delta\).
In particular \(\ell\) is bounded above and away from zero and
\(A\) is comparable to \(\delta^2\) with computable constants.

Write the exactly centered local phase for \(h=q-q_*\) as
\[
 f_\beta(h)=\beta h+a_2h^2+gh^4+h^5E(h).
 \tag{1}
\]
Choose a fixed positive rational physical radius \(h_0\) well
inside the available disk. It may be reduced by any prescribed
positive amount in the following construction. Shrink it so
\(h^5E(h)\) and its real derivative are small relative to the
quartic term and its derivative. Decrease the caps on \(\delta\)
and \(|\beta|\) so at \(h=\pm h_0/2\) the real phase has a strict
negative gap and its real derivative points inward. Outside the
interval between these points, \(\rho-1<0\) follows from the
strict real concavity of \(\rho\) and \(q_*-q_m=O(\delta)\).
Hence the phase gap extends to the two physical endpoints.

Use the uniform
[cotangent integral reduction](MULTIVARIATE_EFFECTIVE_SHARED_ESCAPING_QUARTIC_FIELDS.md#effective-finite-contour-near-the-center).
After dividing the exact moment by its nonzero scalar center
and common order factor, its value is
\[
 I_j(\beta)+{\cal E}_j(\beta),\quad
 I_j(\beta)=\int_{-h_0}^{h_0}A_j(q_*+h)e^{j f_\beta(h)}\,dh,\quad
 |{\cal E}_j(\beta)|\le C_{\rm ph}e^{-\gamma j}.
 \tag{2}
\]
Here \(A_j\) is analytic on a fixed physical rectangle,
positive on its real interval, uniformly bounded there
together with its needed derivatives, and bounded away
from zero. The adjacent ratio satisfies
\[
                  |A_{n+1}/A_n-1|\le C_{\rm amp}/n^2 .
 \tag{3}
\]
The constants in these statements are computable. For (2), the
strict outer phase gap bounds the discarded summands using the
binomial entropy bound. Choose the cotangent height \(\nu>0\)
rational with \(C\nu^2\) below one quarter of the endpoint loss
and below \(\pi\nu/4\). Decrease the field and curvature caps so
these inequalities also absorb the imaginary-field term
\(|\Im\beta|\nu\) and the real-axis increase near the center.
On the half-grid walls the cotangent is bounded by one. Its
difference from the constant on the two horizontal sides is
bounded by \(2e^{-2\pi j\nu}/(1-e^{-2\pi j\nu})\). Moving
the half-grid endpoints back to the fixed endpoints has the
same strict loss. These give (2) after absorbing the
polynomial factors in the order into the exponential bound.
Spatial amplitude bounds and (3) follow from the relative
Stirling remainder on the fixed rectangle, using
\(\sum_i1/m_i\le W_+\) and \(\sum_i |d_i|/m_i\le W_+\).
For example, Binet's remainder in \(\Re z>0\) has kernel
\((\tfrac12\coth(t/2)-1/t)/t\) in its Laplace integral,
between zero and \(1/12\). Thus its magnitude and the magnitude
of its derivative in \(z\) are bounded by
\(1/(12\Re z)\) and \(1/(12(\Re z)^2)\), respectively.
Applying this in the row, binomial, and endpoint terms gives
(3) directly at positive adjacent orders. The undiscarded
half-log row factors in the limiting amplitude are bounded
in logarithm by a constant times \(\sum_i|d_i|/m_i\).

Make the substitution \(h=\ell r z\), put \(L_0=h_0/(\ell r)\),
and write \(a_j(z)=A_j(q_*+\ell r z)\). Dividing (1) by \(A\)
gives the scaled phase
\[
 P_b(z)=-z^4/4+z^2/2+bz+r\ell^4z^5 E(\ell r z).
 \tag{4}
\]
The integral (2) is thereby \(\ell r\) times its scaled
integral with exponent \(M_j P_{b_j}\). On any fixed scaled
disk, the perturbation and its needed derivatives have
computable \(O(r)\) bounds. For bounded \(b\) and bounded
imaginary height of \(z\), on the scaled rectangle with
\(|\Re z|\le L_0\), there is a computable \(R\) beyond which
\[
    \Re P_b(z)\le-c|\Re z|^4,\qquad
                  |P_b(z)|\le C(1+|\Re z|^4).
 \tag{5}
\]
Indeed, choose \(h_0\) small to absorb the perturbation, using
\(r\ell |z|\le2h_0\), then take \(R\) large enough to absorb
the bounded-height lower-degree terms of the quartic. After
shrinking the thermal \(r\)-cap, the entire required
bounded-height rectangle maps into the physical disk.
This order of choices keeps \(h_0\) fixed as \(r\to0\).

## Effective horizontal contours and the analytic seam

The following tests are on the pure model
\(G_b(z)=-z^4/4+z^2/2+bz\) and its small perturbation (4).
All tests use a fixed bound larger than \(B+2\) for
\(|b|\). Tail constants from (5) are used outside a fixed
scaled compact set.

For \(b=iY\), solve \(Y=2v+8v^3\) for the unique real \(v\),
and put \(u=(1+3v^2)^{1/2}\). The right and left saddles
are \(u+iv\) and \(-u+iv\), respectively. Their real
curvatures on their horizontals are \(-2u^2\). At the same
height the real phase loss from either saddle to \(x+iv\)
is exactly \((x^2-u^2)^2/4\). In particular on its own
half-line there is strict loss away from the saddle,
and the midpoint at \(iv\) is below both saddles by at
least \(1/4\). The third stationary point is \(-2iv\)
and is separated from both by at least one.

Here is a finite choice of uniform constants for these
tests. Enclose \(v\) on \(|Y|\le B+2\) by the monotone
cubic, or use \(|v|\le(B+2)/2\) for upper bounds. Bound
\(u\), all three saddle positions and all model
derivatives on a larger compact disk. Choose a small
rational saddle radius \(\varepsilon>0\) so that on a
disk of twice this radius the Taylor error in \(G_{iY}'\)
is less than one quarter of its linear part on the
circle of radius \(\varepsilon/2\). Also take it small
enough to keep negative real horizontal curvature
on that disk. On a fixed compact block outside the
short horizontal saddle segments, the displayed
quartic loss gives a lower bound by a positive
multiple of \(\varepsilon^2\). Near the midpoint
use its loss of at least \(1/4\) and the derivative
bounds. Choose \(0<\eta<1/4\) and then the thermal
\(r\)-cap so that \(|\Re b|\le4\eta\) and the perturbation
of (4) preserve all these losses, with half the margins
if necessary. These choices use just bounds for
polynomials, their derivatives, and the \(O(r)\) bound
on fixed disks. Equivalently, interval subdivisions on
the compact block can certify the strict margins after
each decrease of the perturbation caps. The displayed
strict model inequalities ensure termination.

Rouché's theorem on the \(\varepsilon/2\) circles,
using the linear parts of \(G'_{iY}\), locates one
actual saddle near each of \(u+iv,-u+iv\). Take the
analytic branches \(z_R(b),z_L(b)\) so obtained on
\[
                 |\Re b|<4\eta,\qquad |\Im b|<B+2.
\]
The branches on overlapping disks agree by uniqueness.
For example, use also the test on the larger
\(\varepsilon\) circles. The center locations are
continuous in \(Y\), and the smaller disks keep the
root inside the larger disk under a sufficiently
small change of \(Y\). These local identifications
give single-valued branches on the strip. With the
above choices they satisfy
\[
  \Re(z_R-z_L)\ge c_d>0,\qquad
  |z_s|+|P_b(z_s)|\le C_d,
 \tag{6}
\]
and are real at real \(b\), by conjugation and the
distinct real saddle disks there.

Use the analytic split point
\[
                       z_c(b)=\tfrac12(z_R(b)+z_L(b)).
\]
The integral from \(-L_0\) to \(L_0\) can be split
exactly at \(z_c(b)\), since the integrand has a
primitive on the scaled physical rectangle. Denote
the two integrals, including the factor \(\ell r\), by
\(J_{L,j}(b)\) and \(J_{R,j}(b)\). These are analytic
in \(b\) on the strip. To estimate each integral,
deform it to the horizontal through its own saddle,
with a short connector at the real part of \(z_c\)
and a vertical connector at its physical end.
The midpoint lies close to the model midpoint
\(iv\), so the short connectors keep the strict
midpoint loss. Near the midpoint the horizontal
may extend a little past the original imaginary
axis, also within this loss neighborhood.
The distant horizontals and the outer end edges
are covered by (5). This construction avoids
assigning a connector from zero to \(iv\) to either
packet when \(v\) is not small. At real \(b\), the
exact split integrals have positive integrands on
the real halves with real split point \(z_c(b)\).

For \(|\Re b|\ge\eta\) and \(|b|\le B\), a single
horizontal contour has a strict global maximum
at its actual saddle. Here are uniform model
bounds for this claim. With \(X=\Re b>0\),
\(Y=\Im b\), solve for \(u>\sqrt{1+3v^2}\) and real
\(v\) from
\[
 X=u(u^2-3v^2-1),\qquad Y=v(3u^2-v^2-1).
\]
For fixed \(v\) the first equation is strictly
increasing in \(u\) on the stated range. Its
solution is even in \(v\) and increases with
\(|v|\). The right side of the second equation
increases strictly in \(v\) and has opposite
infinite limits; its derivative is
\(3u^2-3v^2-1+6uv u'(v)>0\).
The bound \(|Y|\ge2|v|\) bounds \(v\), and the
first equation bounds \(u\). On that horizontal
the real phase is a constant plus
\(V(x)=-x^4/4+(1+3v^2)x^2/2+Xx\).
For \(x\ge0\),
\[
 V(u)-V(x)=\frac{(x-u)^2}{4}
                  \big((x+u)^2+2(u^2-1-3v^2)\big).
\]
For \(x=-t\le0\), the loss is \(V(u)-V(t)+2Xt\).
If \(t\ge u/2\) the latter is at least \(\eta\);
otherwise the first term is at least \(1/16\).
Thus there is a uniform strict gap away from
the positive saddle, with strictly negative local
curvature. Rouché and the same compact and tail
bounds preserve these tests for the perturbed
saddle after decreasing \(r\). Reflection handles
\(X<0\). The contour is the full horizontal with
its two outer vertical connectors.

## Uniform packet bounds

For each of the horizontal packets there are computable
constants \(c_1,C_1>0\) such that, for all \(M_j\) above
a computable \(M_1\),
\[
 c_1r M_j^{-1/2}e^{M_j\Re P_b(z_s)}
       \le |J_{s,j}(b)|
       \le C_1r M_j^{-1/2}e^{M_j\Re P_b(z_s)} .
 \tag{7}
\]
The same estimate applies to the full single-saddle
integral. To see this, take a short real segment
\(h=z-z_s\) on the shifted horizontal and use
\(P_b(z)-P_b(z_s)=P_b''(z_s)h^2/2+O(|h|^3)\).
The real parts of the quadratic exponent and the
actual exponent are both at most \(-c h^2\).
The integrated difference of their exponentials is
bounded by \(C M_j\int |h|^3e^{-cM_jh^2}\,dh=O(M_j^{-1})\).
Variation of \(a_j\) on this segment costs
\(C\int |h|e^{-cM_jh^2}\,dh=O(M_j^{-1})\).
The remaining bounded pieces have a fixed loss
below the saddle; the long pieces obey (5).
Increase \(R\) there so their loss below the saddle
is at least \(\tau+c'|\Re z|^4\) for some positive
\(\tau,c'\). The Gaussian coefficient
\(a_j(z_s)\sqrt{2\pi/(-P_b''(z_s))}\) is uniformly
nonzero, taking the square root in the right
half-plane for \(-P_b''\). Its relative error is
bounded by \(C M_j^{-1/2}\), proving (7). Each
Gaussian moment and each tail moment in these
bounds has an explicit computable upper bound.
For example one can bound the tails by integrals
of polynomials times \(e^{-c'x^4}\). For a computable tail
bound beyond \(R\), use \(x^4\ge R^3x\) and integrate the
resulting polynomial times a decaying exponential. In particular
the estimates do not multiply an \(e^{-\tau M_j}\)
error by the full scaled length \(O(1/r)\).

Shrink the \(r\)-cap so \(C_d A\le\gamma/4\), with
\(\gamma\) from (2), and also \(A\le1\).
By \(M_n=n\ell r^4\ge1\),
\(r^{-1}\le(\ell_+ n)^{1/4}\), where \(\ell_+\) is
an upper bound for \(\ell\). Consequently the physical
error in (2) relative to either packet in (7) is
bounded by
\[
                    C n^{3/4}e^{-\gamma n/2}.
 \tag{8}
\]
The same bound holds at order \(N\), using \(N\le2n\).
The one-saddle comparison and (8) show that a zero
with \(|b|\le B\) must have \(|\Re b|<\eta\), above
computable thresholds in \(M_n,n\).

For these nearby fields put
\[
             d(b)=P_b(z_R(b))-P_b(z_L(b)).
\]
Stationarity gives \(d'(b)=z_R(b)-z_L(b)\). At a zero,
(2), (7), and (8) imply that the ratio of the
two packet magnitudes is between \(1/2\) and \(2\),
after increasing \(n_0\) if needed. Thus at the two
orders of a putative common root,
\[
             |\Re d(b_j)|\le C_2/M_j .                  \tag{9}
\]
The two fields have the same imaginary part. Integrate
\(\Re d'\ge c_d\) along their joining real segment in the
strip, using (9). This gives
\[
  |b_N-b_n|\le 2C_2/(c_d M_n),\qquad
  |T|\le 4C_2/c_d,\qquad T=M_N(b_N-b_n).
 \tag{10}
\]

## Adjacent multiplier and the two kinds of drift

Write \(b=b_n\) and \(M=M_n\). Compare the upper-order
integral on the exact lower-order contours, with the
same split point \(z_c(b)\). There is no replacement
error from moving a seam. Its multiplier relative to
the lower integrand is exactly
\[
                     \frac{a_N(z)}{a_n(z)}
                                      e^{A P_b(z)+Tz}.
 \tag{11}
\]
The upper packets used in obtaining (9) need not have
these same contours. Only the exact sum of the two
upper integrals is used in the comparison.

Let \(m_s\) be the ratio of the upper integral on the
chosen half to the exact lower packet on that half.
For each half, (3), (7), and (10) give
\[
 m_s=e^{A P_b(z_s)+Tz_s}(1+\varepsilon_s),\qquad
 |\varepsilon_s|\le C_3\bigl(A/M+|T|/\sqrt M+n^{-2}\bigr).
 \tag{12}
\]
For clarity, the estimate includes the whole contour.
After factoring out the saddle multiplier, the exponential
in (11) differs from one by the exponential of
\[
       A(P_b(z)-P_b(z_s))+T(z-z_s)
\]
minus one. Bound this difference by the absolute
value of its exponent times the larger of one
and the exponential of its real part. Near the
saddle stationarity bounds the factor by
\(C A|z-z_s|^2+|T||z-z_s|\). Integration with the Gaussian
majorant gives \(A/M+|T|/\sqrt M\) after division by
the packet size. The amplitude ratio error is bounded
by \(C n^{-2}\). On bounded complementary pieces
the same difference contains a factor \(O(A+|T|)\),
and the fixed loss gives \(O((A+|T|)\sqrt M e^{-\tau M})\)
relative to the packet. On the long pieces (5) gives
\(|P_b(z)-P_b(z_s)|\le C(1+|\Re z|^4)\) and their loss
absorbs the possible linear growth from \(Tz\), since
\(|T|\) is bounded by (10). Polynomial quartic moments
give the same bound as on the bounded complement,
with new computable constants. The vertical end
edges have bounded scaled length and the same loss.
Finally \(\sqrt M e^{-\tau M}\) is bounded by
computable multiples of both \(M^{-1}\) and \(M^{-1/2}\)
for \(M\ge1\). In particular the complementary pieces
of the multiplier comparison retain the small factors
\(A,|T|,n^{-2}\).

At a common root, subtract \(m_L\) times the lower
equation \(J_{L,n}+J_{R,n}=-{\cal E}_n\) from the upper
equation, on these common contours. The multipliers
in (12) are bounded above, so (8) gives
\[
                 |m_R-m_L|\le C_4 n^{3/4}e^{-\gamma n/2}.
 \tag{13}
\]

If \(\Delta\ne0\), the integer difference of degrees gives
\(|T|=N\ell r|\Delta|\log(N/n)\ge \ell_-r\), with
\(\ell_-\) a positive lower bound for \(\ell\).
For the difference of the two leading exponents in (12),
\[
 \Re(A d(b)+T d'(b))
     =T\Re d'(b)+O(1/n)
\]
by (9), and \(\Re d'(b)\ge c_d\) by (6).
The term \(1/n\), divided by \(|T|\), is at most
\(C r^3/M\). Increase \(M_0\) to make it less than
\(c_d/2\). The leading exponentials are bounded above
and away from zero, and their difference is at least
a fixed multiple of \(|T|\) by the real difference of
their exponents. The errors in (12) are smaller than
half this separation after increasing \(M_0,n_0\):
use \(A/M=1/n\), \(|T|\ge\ell_- r\) and
\(r^{-1}\le(\ell_+ n)^{1/4}\).
This contradicts (13) for large enough \(n\).

If \(\Delta=0\), the two fields agree and \(T=0\).
For all \(b\) in a fixed smaller strip around
\(|\Re b|\le\eta,\ |\Im b|\le B\), the exact packet
ratios on the common analytic seam give holomorphic
functions \(\varepsilon_s(b)\) in (12), and
\(|\varepsilon_s(b)|\le C_5/n\). The logarithms
\(\log(1+\varepsilon_s(b))\) use the branch near one,
and are real for real \(b\) because the exact real-half
integrals and the stationary multipliers are positive there.
Cauchy's bound on a bandwidth with fixed positive slack,
followed by integration from \(\Re b\) to \(b\), gives
\[
       |\Im\log(1+\varepsilon_s(b))|
                        \le C_6|\Im b|/n .               \tag{14}
\]
The function \(d\) is real on the real axis. Its derivative
has real part at least \(c_d\) throughout the vertical
segment from \(\Re b\) to \(b\). Consequently
\(|\Im d(b)|\ge c_d|\Im b|\), with the corresponding sign.
By (14) the imaginary part of the difference of the exact
log multipliers has magnitude at least \(c_d A|\Im b|/2\)
when \(M\ge4C_6/c_d\). Its magnitude is less than one
after shrinking the \(r\)-cap and increasing \(n_0\),
using the uniform bound on \(d\). Thus exponentiation
preserves this lower bound up to a positive uniform
factor and
\[
                         |m_R-m_L|\ge c_5 A|\Im b| .
 \tag{15}
\]

For this last comparison, a zero cannot have
\(|\Im b|\le c_0/M\) for sufficiently small computable
\(c_0>0\). Indeed at \(b=x\) real in the strip the
full retained real integral is positive. Its scaled
integrand, divided by that integral, is a probability
density on the retained interval with uniformly bounded
first moment of \(|z|\). To verify this, use (7) at
the larger real stationary value as a lower bound on the
total mass. The part with \(|z|\le R\) contributes at
most \(R\) to the first moment. The remainder is bounded
by (5) with a strict quartic loss relative to that
stationary value, giving a computable bound independent
of \(r,M\). Changing \(x\) to \(x+iY\) multiplies the
integrand by \(e^{iMYz}\), and changes the integral
relatively by at most \(M|Y|\) times the first moment.
Choose \(c_0\) so this is at most \(1/4\). The lower
bound on the real integral and (8) now exclude a zero
in this range. At a root, therefore,
\(A|\Im b|>c_0/n\). Inequality (15) contradicts (13)
also for \(\Delta=0\).

## Choice of the cutoff

The construction uses only the following decreasing caps
and increasing thresholds. First choose the physical
radius and the physical remainder constants from (1)–(3).
Choose the bounded-height quartic tail and the compact
horizontal contour margins for the specified \(B\), then
shrink the thermal \(r\)-cap to verify the two-saddle and
one-saddle tests. The Gaussian bounds supply the constants
for (7), (9), and (10). These give a fixed upper bound on
\(|T|\) before choosing the constants in (12).
Finally shrink the \(r\)-cap further for (8) and the phase
bound in (15), increase \(M_0\) to absorb the Gaussian and
multiplier errors, and increase \(n_0\) so the bound (8)
is smaller than the needed multiples of \(1/n\) and
\(n^{-1/4}\). The latter comparisons terminate by
exponential decay, and can use any verified monotone
polynomial-exponential upper bound beyond its maximum.
All interval tests on the bounded model parts have strict
margins already established above. The relation
\(r^2=2a_2\ell\) turns the thermal \(r\)-cap into a decrease
of \(\delta_0\), and \(A\asymp\delta^2\) gives the
equivalent threshold in \(\sqrt n\,\delta\).
No restriction on \(\Delta\) enters the initial caps.
