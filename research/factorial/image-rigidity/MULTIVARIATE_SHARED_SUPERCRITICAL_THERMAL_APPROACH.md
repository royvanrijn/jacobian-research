# Supercritical thermal approach to a varying shared cusp

Use \(W,\delta,C_0,q_*,w_*,a_2,g\) and the centered scales
\(w_j=-\operatorname{Log}(j^\Delta t)-C_0\) from the
[uniform moving quartic proof](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
The profile is fixed within each pair of adjacent orders.

**Theorem.** Fix a positive bounded \(W\)-band. There exist
\(\delta_0,r_*>0\) and a terminal cutoff such that for every
larger terminal and every fully shared ordinary profile with
\(0\le\delta\le\delta_0\), there is no consecutive common root
satisfying \(|w_j-w_*|<r_*\) at either order. The total-degree
difference \(\Delta\) is unrestricted. The cutoff and radii
are unevaluated.

The [uniform supercritical theorem](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md)
covers sequences with \(\delta\) bounded away from zero,
and the [moving quartic proof](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md)
covers bounded \(\sqrt n\delta\). Thus it suffices to exclude
sequences with \(\delta\to0^+\), \(\sqrt n\delta\to\infty\) and
one centered scale tending to \(w_*\) in difference, since
otherwise one could take the physical exclusion radius to
zero along a counterexample sequence.

The limiting centered scales are both the same real cusp:
the other scale stays compact by the uniform root annulus,
has the same imaginary part, and a noncusp real scale
has the one-saddle exclusion of the moving quartic proof.
For both actual orders, therefore, the finite local integral
and exponentially smaller physical remainder in that
proof apply. Write \(\beta_j=w_*-w_j\). The latter two fields
have the same imaginary part and differ by
\(\Delta\log(1+1/n)\) in their real parts.

## Two scales for the almost critical phase

The maximizer \(q_m\) of \(\rho\) satisfies \(q_*-q_m=O(\delta)\).
Taylor expansion at \(q_m\) gives \(\rho(q_*)-1=\delta+O(\delta^2)\),
so \(a_2=F''(q_*)/2>0\) is comparable to \(\delta\).
Put
\[
 \ell_0=(-1/(4g))^{1/3},\qquad r=\sqrt{2a_2\ell_0},\qquad
 b_j=\beta_j/r^3,\qquad M_j=j\ell_0 r^4 .
\]
Then \(r\to0\) and \(M_j\to\infty\). When \(b_j\) is bounded,
the centered phase under \(q-q_*=\ell_0 r z\), divided by
\(\ell_0 r^4\), is
\[
                   -z^4/4+z^2/2+b_j z+O(r|z|^5) .           \tag{1}
\]
The remainder is analytic on the fixed physical disk
and the local scaled convergence includes all derivatives.

If \(|b_j|\to\infty\), use instead
\(\epsilon_j=|\beta_j|^{1/3}\) and divide the phase by
\(\ell_0\epsilon_j^4\). This is the escaping-field phase (6)
of the [moving quartic proof](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md#escaping-fields-at-bounded-scaled-temperature).
Indeed \(a_2\ell_0/\epsilon_j^2=(r/\epsilon_j)^2/2\to0\)
and \(j\epsilon_j^4\to\infty\). Its one-saddle argument
forces the field of any such individual zero to have
argument tending to \(\pi/2\) or \(-\pi/2\).
Since the physical imaginary part is shared, a bounded
\(b_j\) at one order rules out an escaping \(b_j\) at the
other. If both escape, the two-saddle adjacent comparison
in that proof applies verbatim. Its hypotheses there
are precisely \(\epsilon_j\to0\),
\(a_2/\epsilon_j^2\to0\), \(j\epsilon_j^4\to\infty\) and
the common real parameter shift. We treat the remaining
case of bounded \(b_n,b_{n+1}\).

## The two saddle packets for bounded \(b\)

For \(b=B+iY\) with \(B>0\), a stationary point of
\(G_b(z)=-z^4/4+z^2/2+bz\) has the form \(u+iv\)
with \(u>\sqrt{1+3v^2}\). For each real \(v\), solve
\[
                  B=u(u^2-3v^2-1)
\]
uniquely for that positive \(u\). It is even and increases
with \(|v|\). The equation \(Y=v(3u^2-v^2-1)\) has
a unique solution \(v\), since its right-hand side has
positive derivative and tends to opposite infinities.
Indeed its derivative is
\(3u^2-3v^2-1+6uv u'(v)>0\).
On the horizontal through this saddle, the real phase
is a constant plus
\[
                  -x^4/4+(1+3v^2)x^2/2+B x.
\]
Positive \(x\) wins against its reflection, and the
derivative has exactly one positive zero, with negative
second derivative there. Thus the saddle is the
unique strict real-phase maximum of its horizontal.
Reflection handles \(B<0\). On bounded regions away
from \(B=0\) the same horizontal maximum and Gaussian
comparison hold for (1), so an individual zero can
only limit to \(B=0\).

On the imaginary line, the two saddles are
\(z_R=u+iv\), \(z_L=-u+iv\), where
\(u=\sqrt{1+3v^2}\) and \(Y=2v+8v^3\).
Each has negative real curvature on its horizontal
half-line. The midpoint is lower than both by
\(u^4/4\ge1/4\). In a small complex disk about any
fixed \(iY\), use the two continued saddles and the
same split horizontal contours as in the escaping
quartic proof. The midpoint connectors keep the
strict scaled loss. The distant horizontal tails
have quartic decay at bounded heights, because
\(r|z|^5\) is a small multiple of \(|z|^4\) on the fixed
physical disk. The physical end connectors and the
exterior retain their fixed phase gaps. These facts
give two actual Gaussian packets with amplitudes
uniformly bounded above and away from zero after
removing their stationary exponentials and the common
length factor \(rM_j^{-1/2}\). In particular all
these estimates hold for complex \(b\) in a fixed small
disk about zero, using the real halves with common split
point \(q_*\) as the initial contours. Their exact
integrals there are analytic in \(b\), conjugate-symmetric,
and positive for real \(b\) after removal of the common
real exponential. This last assertion holds on
the initial real halves, whose integrands are positive.

Let \(D(\beta)\) be the right-minus-left actual stationary
phase difference. At a common root, both scaled parameters
tend to the same \(iY\), because each has real part tending
to zero and their imaginary parts are the same. Cancellation
and the separated saddle estimates give
\(\Re D(\beta_j)=O(n^{-1})\) at both scales. As before,
\(D_\beta=q_R-q_L\) and its real part is comparable to \(r\)
on the segment joining the two fields. Integration across
that real segment gives \(|\Delta|r=O(1)\).

## Adjacent comparison and small imaginary fields

Use the same lower-order contours to compare the
integrals at the two orders. After removing the common
centered factors, their multiplier is
\[
 \frac{A_{n+1}(q)}{A_n(q)}
     \exp\bigl(\Phi_{w_n}(q)-\Phi_{w_n}(q_*)
          +\Delta(n+1)\log(1+1/n)(q-q_*)\bigr).
 \tag{2}
\]
The amplitude ratio is \(1+O(n^{-2})\), by the
complex-order Stirling derivative bound of the moving
quartic proof on this fixed physical disk. With
\(M=M_n\), the relative errors when evaluating the
exponential at each lower saddle are
\[
                 O\left(r^4/M+|\Delta|r/\sqrt M+n^{-2}\right).
 \tag{3}
\]
Indeed the phase in (2) changes from the saddle by
\(O(r^4|z-z_j|^2+|\Delta|r|z-z_j|)\) locally, by exact
stationarity, and these factors have the indicated
Gaussian moments. The other pieces of the contour
have the same gaps as above. On the long scaled tails
the growth of the multiplier is bounded by
\(\exp(Cr^4(1+|z|^4)+C|\Delta|r(1+|z|))\), so the
quartic majorant absorbs it. The fixed physical
connector loss absorbs the factor \(e^{O(|\Delta|+1)}\).
This gives (3) uniformly on the bounded scaled disks.
The exterior remainder is exponentially small
relative to either centered packet, whose modulus
is bounded below by \(c rM^{-1/2}e^{-C M}\).

When \(\Delta\ne0\), the real difference of the exponents
in (2) at the two saddles has modulus bounded below by
a constant times \(|\Delta|r\), by the derivative and
modulus-cancellation bounds above. The exponentials
at the saddles have bounded moduli with positive lower
bounds. The separation therefore dominates (3). Double
cancellation of the packets is impossible.

When \(\Delta=0\) and the common limit \(iY\) has \(Y\ne0\),
the imaginary stationary difference is bounded below
in modulus by a positive multiple of \(r^4\). Indeed the
limiting difference vanishes at \(b=0\), and along the
imaginary axis its derivative in \(Y\) has imaginary part
\(2u>0\) by the stationary-value derivative. Hence the
same comparison excludes this case.

Suppose finally that \(\Delta=0\) and \(b=x+i y\to0\).
On a fixed small disk about zero the ratios of the exact
centered packets from (2), divided by their stationary
multipliers, are analytic and nonzero. They tend to one
with error \(O(r^4/M+n^{-2})=O(n^{-1})\) there by (3).
Their logarithms are real on the real diameter, by the
positivity of the two real-half integrals. Cauchy's
estimate on a smaller disk and integration from \(x\) to
\(x+i y\) show that the imaginary error in each logarithm
is \(O(n^{-1}|y|)\).
On the other hand \(D\) is real at \(b=x\), and
\[
 \partial_b D=r^3(q_R-q_L)
]
has real part at least \(c r^4\) throughout a fixed
smaller disk. Thus the leading imaginary difference
has modulus at least \(c r^4|y|\) with a fixed sign
relative to \(y\). Since \(nr^4\to\infty\), it dominates
the imaginary errors.

The exponentially small exterior remainders are also
negligible at this precision. To check this even for
very small \(y\), use the exact centered real-half
integrals at \(b=x\). Each is positive and its positive
integrand, normalized to total mass one, has bounded
first moment of \(|z|\), uniformly at real \(x\) near zero.
This follows from its local Gaussian lower bound
and the strict outer loss and quartic bound. Changing
\(b=x\) to \(x+i y\) multiplies that integrand by
\(e^{i M y z}\). The relative change in each half is
therefore at most \(C M|y|\). For sufficiently small
fixed \(c_0>0\), if \(|y|\le c_0/M\) both perturbed
halves have positive real parts at least half their
original values. Their sum dominates the exponentially
small exterior, so the actual value cannot vanish.
Thus a root in this disk requires \(|y|>c_0/M\), and
\(r^4|y|\gtrsim 1/n\). This bound absorbs the exterior
errors in the adjacent logarithmic phase comparison.
The two exact packets cannot cancel at both orders,
which finishes the last thermal regime and proves
the theorem. \(\square\)
