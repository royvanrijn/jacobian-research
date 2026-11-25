# Actual boundary and saddle contributions for one private side

Let \(a,b\) be two distinct ordinary factorial exponent
profiles with private variables on exactly one side.
Orient them so
\[
 P=\sum_{a_i=0}b_i>0,\qquad b_i>0\ \text{for every varying row}.
\]
Remove equal rows, put \(d_i=b_i-a_i\), \(\Delta=\sum_i d_i\),
and set
\[
 K_a=\prod_{a_i>0}a_i^{d_i},\qquad
 \sigma_a=\sum_{a_i>0}d_i^2/a_i.
\]
All statements concern each fixed profile, with scaled
parameter \(\zeta=n^\Delta t\).

**Theorem.** If \(P\ge2\), or if \(P=1\) and \(\sigma_a>1\),
the ordinary consecutive moment polynomials are coprime
for all sufficiently large orders. The degree-one disjoint
case \(P=1,\ \sum_i|d_i|=1\) is coprime at every order.
For \(P=1,\ \sigma_a=1\), the argument below also excludes
limiting common roots at every scaled parameter except the
positive boundary transition; at that transition the
[bounded scaled windows](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#continuous-endpoint-windows)
are excluded separately.

The proof identifies the actual terms on the two vertical
sides of the [Rice contour](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md).
The right side is the single upper saddle. The pinching
left side gives the boundary contribution as an exact
Mellin–Barnes integral. For \(P=1,\ \sigma_a<1\), its direct
bound also covers a complex sector described below.

## The actual right contour

Put \(x_i(q)=a_i+d_iq\), \(D_\pm=\sum_{\pm d_i>0}|d_i|\),
\(D=D_++D_-\), and
\(p_R=\min_{d_i<0}(-a_i/d_i)>1\), or \(+\infty\)
if no gaps are negative. Assume \(D\ge2\); otherwise,
after removing equal rows, the two terms are a constant
and one variable and the normalized reciprocal is a
truncated exponential.

For \(\zeta\) in the upper half-plane use the principal
upper phase \(\phi\) and its continuation \(\phi_-\) from
the Rice formula. Thus
\[
 \phi'(q)=X(q)+\operatorname{Log}\zeta,\qquad
 \phi_-'(q)=X(q)+\operatorname{Log}\zeta-2\pi i .
\]
The cut domain is formed by the two open half-planes
joined through \((1,p_R)\). The other real points are
removed. It contains the right vertical contour through
any fixed \(c_R\in(1,p_R)\); there are no poles in the
domain. The real phase and its gradient agree across this
gap, with a possible imaginary constant \(2\pi i\) that
does not affect integer orders.
The [private sheet count](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#the-degree-on-every-open-horizontal-strip)
gives exactly one upper saddle \(q_*(\zeta)\). It is simple
and analytic in the upper parameter half-plane. The lower
phase has no critical point by the Rice inequality
\(\Im\phi_-'<\arg\zeta-\pi<0\).

The two descending arms of this saddle escape to infinity.
The only possible finite real exit of an upper arm
is through the right gap. Indeed \(\Im\phi'>0\) on all
other real intervals, by the boundary heights of \(X\).
This prevents both crossing and forward accumulation,
using the local backward flow from a hypothetical accumulation
point. At the real branch points other than zero the
[real-boundary argument](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#where-the-upper-descending-contours-can-meet-the-real-axis)
applies with the same signs: the endpoints of the right
gap have opposite imaginary derivative signs, and at
all other poles both signs are positive. At zero, if \(P=1\),
the upper phase is regular with positive vertical
velocity \(\arg\zeta\). If \(P>1\), it has logarithmic
coefficient \(P-1>0\) in the phase derivative and positive
boundary imaginary derivatives on both sides. The
farther-pole argument therefore applies at zero too.
In detail, a level curve approaching zero must have
the value of \(\Im\phi(0)\). Its imaginary phase difference
has leading form
\((P-1)y\log|q|+(\arg\zeta+(P-1)\arg q)x+O(y+|q|^2)\).
Vanishing forces \(x>0\), \(y/|x|\to0\); but the real
descent velocity then points strictly away from zero.
The same observation excludes subsequential accumulation,
by starting at the first entry on a fixed sufficiently
small semicircle, where \(|x|\) is at least half its radius.
No arm can accumulate in a compact interior set: the
strict phase loss leaves no stationary limit point.
Once an arm enters the lower half-plane, it keeps
descending vertically and cannot return.

For completeness the escape estimates use the
[general Rice tail argument](MULTIVARIATE_SHARED_SUPERCRITICAL_EVENTUAL_COPRIMALITY.md#descending-paths-and-stirling-bounds-at-infinity).
At infinity the two phase derivatives have form
\(\Delta\log q+C_\pm+O(|q|^{-1})\), with
\(\Im C_+=\arg\zeta+\pi(D_--1)\) and
\(\Im C_-=\arg\zeta-\pi(D_-+1)\).
If \(\Delta\ne0\) the horizontal descent velocity has
size comparable to \(\log|q|\), directed left for
\(\Delta>0\) and right for \(\Delta<0\), and the
vertical velocity is bounded. In that preferred
horizontal half-plane the upper velocity is bounded
below by a positive constant at infinity, and the lower
one above by a negative constant, since their leading
bounds are respectively
\(\arg\zeta+\pi(D/2-1)>0\) and
\(\arg\zeta-\pi(D/2+1)<0\).
Consequently an unbounded path eventually reaches one
of the two outward regions. In fact, if its height
becomes large before its real part is in the preferred
region, bounded vertical speed and the horizontal speed
directed to that region carry it there before its
height can return to a fixed bounded interval. After
entry, it stays there and escapes. Integrating the two
speed bounds gives
\(|\Re q|=O(|\Im q|\log(2+|\Im q|))\).
For \(\Delta=0\), both vertical velocities at infinity
are already outward in their whole half-planes and the
horizontal speed is bounded, giving the linear escape
bound. These are stable end regions, uniformly on
compact parameter paths.

As \(\arg\zeta\) approaches \(\pi\), the saddle lies on
the real right gap at the limiting parameter. In fact
\(X'\) is positive there and its real part ranges over
the line. Its two descending directions are vertical,
one in each half-plane, and at this negative-real
parameter the upper vertical velocity is positive
everywhere. The corresponding contour is homotopic
to the right vertical line by interpolation at the same
imaginary height. Its end types stay fixed throughout
the upper parameter half-plane, by finite-time continuous
dependence until reaching the stable end regions and
the absence of finite accumulation. The same
continuation provides the contour homotopy, truncating
at a common large height on compact parameter paths.

The exact Gamma expansions justify closing these
homotopies with horizontal connectors. For \(\Re z\le0\)
and large \(|\Im z|\), reflect \(\Gamma(1+z)\) using
\(-\pi/(\Gamma(-z)\sin\pi z)\) and apply right-half-plane
Stirling to \(-z\). The dominant sine exponential
recovers the original chosen phase and amplitude with
relative error \(O(e^{-2\pi|\Im z|}+1/|z|)\).
Use this also on the two arguments of the Gamma ratio,
applying \(\Gamma(z)=\Gamma(1+z)/z\) as needed.
This works at all large heights even when the
horizontal term is of order \(|\Im q|\log|\Im q|\).
At each such height the real phase on the fixed right
vertical line has negative linear rate in the height.
For \(\Delta\ne0\), its horizontal derivative has sign
\(\Delta\) at large height, so the maximum on the
connector is at that fixed line, since the descending end
leans in the opposite direction. For \(\Delta=0\)
use the linear phase \(C_\pm q+O(\log|q|)\): both ends
of the connector have linear loss, as does the entire
connector. The exact contour connector integrals therefore
vanish as the height goes to infinity.

It follows that the actual right period, relative to
the pure \(a\) moment, has the expansion
\[
 {\cal R}_n(\zeta)
   = n^{p/2}\,J(\zeta)e^{n S(\zeta)}
                    (1+O(n^{-1})),                       \tag{3}
\]
uniformly on compact parameter disks in the upper
half-plane. Here \(p=\#\{i:a_i=0\}\), \(J\) is analytic
and nonzero locally, and
\[
 S(\zeta)=F_a(q_*)-\Delta q_* - A_a-\sum_i a_i,\qquad
 F_a(q)=\sum_i a_i\log x_i(q)-\log(1-q),\qquad
 A_a=\sum_{a_i>0}(a_i\log a_i-a_i).
\]
Equivalently \(S\) is the stationary Rice phase minus
its finite value at \(q=0\); the displayed expression
uses that \(F_a(0)=A_a+\sum_i a_i\).
The power \(p/2\) is the ratio of the interior and pure
boundary Gaussian powers. It also follows directly from
the Rice prefactor: \(n!\) after cancellation of its
exponential and \(n^n\) leaves \(n^{1/2}\), canceled
by the saddle Gaussian \(n^{-1/2}\), and the ratio of the
row Stirling factors at \(m_i\) and the positive \(a_i\)
leaves \(n^{p/2}\). To see (3), expand the exact
phase in its analytic quadratic coordinate at the saddle.
The amplitude is nonzero and the oriented Gaussian
integral is nonzero. Local Stirling and Gaussian
expansions give the indicated error. Beyond a fixed
saddle neighborhood the strict phase loss is uniform on
a smaller parameter disk, and the large-height loss
just established controls the infinite tails. The
stationary derivative obeys
\(\partial_{\log\zeta} S=q_*\).

## The pinching left contour

Choose its real part to be \(-c/n\) with
\(0<c<\min_{a_i=0}1/d_i\), and set \(z=nq\).
Reflection of the rational Gamma ratio in the Rice
formula gives the exact decomposition, relative to the
pure \(a\) moment,
\[
 \frac{M_n(n^{-\Delta}\zeta)}{M_n(0)}
      ={\cal R}_n(\zeta)+{\cal L}_n(\zeta),
 \qquad
 {\cal L}_n=\frac1{2\pi i}\int_{\Re z=-c}
  \Gamma(-z)\prod_{a_i=0}\Gamma(1+d_i z)
        (-\zeta K_a n^{1-P})^z T_n(z)\,dz,                 \tag{4}
\]
where
\[
 T_n(z)=
   n^{-z}\frac{\Gamma(n+1)}{\Gamma(n+1-z)}
   \prod_{a_i>0}\left[
       (na_i)^{-d_i z}\frac{\Gamma(1+na_i+d_i z)}
                              {\Gamma(1+na_i)}\right].
 \tag{5}
\]
The sign in (4) comes from subtracting the upward left
vertical integral and the reflection factor
\((-1)^{n+1}\). The identity of the powers uses
\(\sum_{a_i>0}d_i=\Delta-P\).
Take the branch of \(\log(-\zeta)\) with argument
\(\arg\zeta-\pi\).

### Private degree at least two

For fixed \(\Re z=s_1\) and \(n\) sufficiently large,
all numerator moduli of the active Gamma ratios are
bounded by their values at the real parts, by the Gamma
integral. After normalization these are bounded
uniformly in \(n\). For \(\lambda=n+1-s_1\ge1\) the
Gamma Euler product gives
\[
 \frac{\Gamma(\lambda)}{|\Gamma(\lambda+i v)|}
   =\prod_{k=0}^\infty
       \left(1+\frac{v^2}{(\lambda+k)^2}\right)^{1/2}
   \le e^{\pi|v|/2}.
\]
Compare the product to its value at \(\lambda=1\),
namely \((\sinh(\pi v)/(\pi v))^{1/2}\).
Thus the Gamma factors outside \(T_n\) in (4) leave
exponential decay at least
\(e^{-\pi P|v|/2+(\pi-\arg\zeta)v}\), times a
fixed polynomial. This is integrable uniformly on
upper parameter compacts when \(P\ge2\).

Move the left line right across the fixed integer
poles at \(z=0,1,2\) to \(\Re z=5/2\), for \(n\)
large enough that no active Gamma pole is encountered.
The horizontal integrals vanish by the same estimate,
uniformly on the bounded real range of ordinates.
The residue at zero gives \(1\) in \({\cal L}_n\).
The residues at one and two are \(O(n^{1-P})\)
and \(O(n^{2(1-P)})\), with the first having an
ordinary expansion in inverse powers of \(n\)
from the factorial ratios. The remaining integral is
\(O(n^{-5(P-1)/2})\). In particular
\[
 {\cal L}_n(\zeta)=1+\frac{L_1(\zeta)}n+O(n^{-2})
 \qquad(P\ge2) .                                         \tag{6}
\]
This is the actual left contour, not just a local
real endpoint estimate.

### One private linear variable

When \(P=1\) the factors outside \(T_n\) in (4) are
\(\Gamma(-z)\Gamma(1+z)(-\zeta K_a)^z\).
The Gamma Euler product gives a sharper collective
bound for \(T_n\). Define
\[
 J(a,c)=\frac12\int_a^\infty\log(1+c^2/u^2)\,du,\qquad
 B(y)=\sum_{a_i>0}|d_i|\arctan(|d_i|y/a_i)-\arctan y .
\]
Comparing the logarithmic sums in the Euler products
with their integrals, and shifting the lower bounds of
integration by a bounded amount, gives
\[
 \log|T_n(s_1+iv)|
 \le O(\log(2+|v|/n))
        +n\left[J(1,|v|/n)
                -\sum_{a_i>0}J(a_i,|d_i v|/n)\right].
 \tag{7}
\]
Indeed \(u\mapsto\log(1+c^2/u^2)\) is decreasing;
the difference between its sum and integral is at
most one initial term. The real parts of the Gamma
arguments differ from \(n\) and \(na_i\) by bounded
amounts and remain comparable to them.
Also \(dJ(a,c)/dc=\arctan(c/a)\), so the derivative of
the bracket with respect to \(|v|/n\) equals \(-B\).

If \(\sigma_a\ge1\), then \(B\ge0\).
Use the subadditivity of \(\arctan\) on nonnegative
numbers, with \(|d_i|\) copies of the argument
\(|d_i|y/a_i\), to obtain
\(\sum |d_i|\arctan(|d_i|y/a_i)\ge\arctan(\sigma_a y)\).
More generally put \(\theta_*=\sup_{y\ge0}(-B(y))\).
When \(D\ge2\) this is less than \(\pi/2\):
\(B(0)=0\), the limit of \(B\) at infinity is
nonnegative, and at any finite \(y\) one has
\(-B(y)<\pi/2\).
For \(\arg\zeta>\theta_*\), the bound (7) and the
decay from \(\Gamma(-z)\Gamma(1+z)\) give a
uniformly integrable exponential majorant on every
compact parameter set. On \(|v|\le C\log n\) the
fixed-shift Gamma ratios in (5) give
\(T_n=1+t_1(z)/n+O(n^{-2}(1+|z|)^K)\), with
\(t_1\) polynomial. This follows by Taylor expansion
of the ordinary Stirling formula in \(z/n\).
The exponentially small remainder beyond \(C\log n\)
is \(O(n^{-2})\) for sufficiently large fixed \(C\).
Thus the Mellin–Barnes identity gives
\[
 {\cal L}_n(\zeta)=\frac1{1-\zeta K_a}
                     +\frac{L_1(\zeta)}n+O(n^{-2})
 \qquad(\arg\zeta>\theta_*).                             \tag{8}
\]
For the leading identity, first take \(|\zeta K_a|<1\)
and close right to obtain the geometric series from the
poles of \(\Gamma(-z)\); both sides continue throughout
the upper parameter half-plane. Equivalently this is
the Mellin inversion of \((1+u)^{-1}\) with
\(u=-\zeta K_a\).

## Eventual exclusion

On the imaginary \(q\) axis the height of \(X\) is
\[
 (P-1)\pi/2+
 \sum_{a_i>0}|d_i|\arctan(|d_i|y/a_i)-\arctan y .
\]
For \(P\ge2\) it is positive. For \(P=1\) it is at least
\(-\theta_*\). Hence on the parameter sets covered by
(6) or (8) the unique upper saddle has \(\Re q_*>0\).
Indeed it cannot cross the axis there, and starts on
the right side at the negative-real parameter boundary.
The relevant parameter sector is connected.

Suppose along a subsequence that a common scaled root
tends to an interior point of one of these sectors. At
a root (3) and (6) or (8) imply \( \Re S\to0\), since
the leading boundary factor never vanishes there.
At the next order the scaled parameter is
\(\zeta'=(1+1/n)^\Delta\zeta\). Dividing the two zero
equations yields the necessary condition
\[
 \exp(S+\Delta q_*)\longrightarrow1 .
\]
If \(\Delta\ne0\), its limiting modulus is
\(e^{\Delta\Re q_*}\ne1\). If \(\Delta=0\), division
instead forces \(n\Re S=O(1)\), by the equal parameter
and the \(1+O(n^{-1})\) errors. The first zero equation
requires \(n\Re S=-(p/2)\log n+O(1)\), also impossible.

For the global assertion, all scaled roots lie in a
fixed compact annulus once \(n\) is large. At small
\(|\zeta|\) the pure \(a\) boundary dominates, with
nonzero stable amplitude (4) of the
[boundary analysis](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#actual-boundary-contributions).
To justify localization, compare the modulus outside
a fixed boundary neighborhood with the positive real
integral at a small fixed radius, whose phase has a
strict gap there. At large \(|\zeta|\), reverse the
monomials and use \(\eta=1/\zeta\). At \(\eta=0\) the
pure \(b\) maximum is fully interior and nondegenerate
because every \(b_i>0\). This yields a single
nondegenerate saddle and no zeros on a fixed small
complex disk in \(\eta\), by the local contour and
compact upper-level-set estimates.

The negative real scales are regular here as well.
For \(\zeta=-r_0<0\), split the original real integral
at \(1-r_0e^s=0\). On the positive side its scalar
phase decreases strictly to the split, because the
logarithmic derivative of that factor is negative and
\(q(s)>0\). The pure \(a\) boundary is the unique
maximum on this half. On the negative side the phase
for the modulus has its unique maximum at the right
gap saddle \(q>1\), by the strict monotonicity and full
range of \(\Re X\) in that gap. Other points have a
strict gap in the phase. The gap persists in a complex
neighborhood. For \(|\zeta+r_0|\le\eta\), use the upper
envelope \(|y^a-r_0 y^b|+\eta y^b\), whose upper
phase level sets are uniformly bounded. If its maximum
outside the two retained neighborhoods failed to
converge to at most the negative-real maximum on that
complement as \(\eta\downarrow0\), a converging sequence
of maximizers would contradict the continuity of that
envelope before taking its logarithm. Its one-order
integral is finite, so the phase gap bounds the entire
discarded integral. On the retained neighborhoods
the stable boundary and regular interior expansions
continue analytically, with an extra sign \((-1)^n\)
for the negative half. At an equal-modulus point the
consecutive multiplier ratio is a negative real number,
so cannot equal one. Off that point one term dominates.

For positive real \(\zeta\), the
[real geometry and coexistence analysis](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md#real-geometry-and-consecutive-comparison-at-coexistence)
applies. If \(P\ge2\), or \(P=1\) with \(\sigma_a>1\),
there are only stable boundary or nondegenerate interior
global maxima and their regular coexistence. At the
spinodal the degenerate saddle is strictly below the
stable global maximum. If the \(P=1\) boundary stability
limit is reached, the interior already dominates there.
Thus every positive real scale has a neighborhood
excluding all sufficiently large common roots in these
cases. Compactness of the root annulus and conjugation
now complete the theorem. Rational Bézout identities
give the same consecutive coprimality over any
characteristic-zero field.

If \(P=1,\ \sigma_a=1\), equations (3) and (8) still
rule out all compact upper or lower sequences away
from the real axis. The only remaining accumulation
point is the positive endpoint transition. The known
bounded windows there leave the sequences for which
\(n^{2/3}|\log(\zeta K_a)|\to\infty\).
For \(P=1,\ \sigma_a<1\), the direct left-contour
estimate leaves the shallow sectors
\(0<|\arg\zeta|\le\theta_*\) as well as escaping
fields at the continuous endpoint.
