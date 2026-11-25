# Consecutive diffuse coprimality in a moving critical window

Write
\[
 B_n(z;s)=\sum_{j=0}^n\binom nj e^{-s j(n-j)/(2n)}z^j .
\]

**Theorem.** For every fixed \(0\le K<\infty\) there is \(n_0(K)\)
such that, for \(n\ge n_0(K)\) and \(s>0\) satisfying
\[
 |s-4|\le \frac{K}{\sqrt n},
\]
the polynomials \(B_n(z;s)\) and \(B_{n+1}(z;s)\) are coprime.
The two degrees use the same value of \(s\). This concerns the
[diffuse limiting polynomials](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#the-limiting-pair-when-each-relative-gap-tends-to-zero).
The threshold is obtained by compactness and is not evaluated.

The distinction between consecutive degrees in the bounded field
window is a root-speed inequality for a quartic Fourier integral.
Its proof uses a Rayleigh equation and a normalized Airy limit
on the subcritical side. The large-field part then follows from
two separated quartic saddles. All roots are on the unit circle
by the [unitary Hermite circle theorem](https://doi.org/10.5802/ahl.227).

## A quartic root-speed inequality

Consider the real Fourier integral
\[
 V(t,x)=\int_{\mathbb R}\exp\!\left(-u^4/4-tu^2/2+ixu\right)\,du,
 \qquad t\in\mathbb R .
\]
For each real \(t\), all its zeros in \(x\) are real by Theorem 2 of
Newman's *Fourier Transforms with Only Real Zeros*. In his notation
(1.8), use the empty product, quartic coefficient \(1/4\), and
quadratic coefficient \(t/2-1\). Taking the damping parameter in
(1.6) to be \(1\) gives this Fourier integral by \(z=ix\).
The measure is finite for every \(t\) because of the quartic term.
[Newman, Theorem 2](https://www.math.northwestern.edu/~auffing/papers/Newman.pdf)

The zeros are simple. The integral is analytic in both variables,
and integration by parts gives
\[
 V_t=\tfrac12 V_{xx},\qquad V_{xxx}=xV+tV_x .
\]
A root of multiplicity at least three contradicts uniqueness for
the third-order equation and \(V(t,0)>0\). At a double zero
\((t_0,x_0)\), the functions
\(\epsilon^{-2}V(t_0+\epsilon^2,x_0+\epsilon w)\) converge
locally uniformly in \(w\) to
\(\tfrac12 V_{xx}(t_0,x_0)(w^2+1)\). Rouché's theorem would
give nonreal zeros at real \(t_0+\epsilon^2\), contradicting
Newman's theorem. Also \(V(t,\cdot)\) has order at most \(4/3\)
by the quartic integral bound
\(|V(t,x)|\le C_t\exp(C_t|x|^{4/3})\). If there were only
finitely many zeros, Hadamard factorization would make it an
exponential of a linear polynomial times a polynomial. Such
a function cannot decay to zero on both ends of the real line,
as \(V\) does by integration by parts. Thus the positive zeros
form a strictly increasing infinite sequence. Each indexed
zero \(y_k(t)\) is real analytic in \(t\) by simplicity.
Locally the index is preserved by continuity on a bounded
interval containing the first \(k\) zeros and whose endpoints
are free of zeros.

For \(t>0\) we prove
\[
                    3y_k(t)-2t y_k'(t)>0.                    \tag{1}
\]
Fix the index and write \(y=y_k\), \(p=y'\).
Differentiating \(V(t,y(t))=0\) and using the two differential
identities yields
\[
 y''=2p^3-\tfrac12tp-\tfrac14y.                                \tag{2}
\]
Indeed \(V_{xx}/V_x=-2p\) at a root. Differentiating once
more gives
\(V_{tt}+2pV_{tx}+p^2V_{xx}+y''V_x=0\), with
\(V_{tx}=tV_x/2\) and
\(V_{tt}=(yV_x+tV_{xx})/4\) at that root. This is the
Rayleigh equation of Hernández-del-Valle, derived here
for the indexed simple roots.
[Hernández-del-Valle, Theorem 4.2](https://arxiv.org/pdf/1506.02502)

We first record the large-\(t\) bound
\[
 \limsup_{t\to\infty}\frac{y_k(t)}{t^{3/2}}
                \le \frac{2}{3\sqrt3}.                        \tag{3}
\]
Put \(h=\sqrt{t/3}\), \(x_e=th-h^3\), and \(c=(3h)^{1/3}\).
Shift the full integration line to \(u=ih+w\). The integrand
is entire and quartic decay kills the closing sides. At \(x=x_e+c v\),
after division by a positive scalar and substitution \(w=\xi/c\),
the integral becomes
\[
 \int_{\mathbb R}
    \exp\!\left(-\frac{\xi^4}{4c^4}+iv\xi-\frac{i\xi^3}{3}\right)
    d\xi .
\]
Rotate the right half to argument \(-\pi/12\) and the left
to argument \(-11\pi/12\). In the intervening wedges the
quartic and cubic terms have nonpositive real parts and the
quartic term gives decay on the closing arcs. On the two rays
the cubic term alone gives a uniform integrable majorant for
bounded \(v\) as \(c\to\infty\). The limit is
\(2\pi\,\mathrm{Ai}(-v)\), by the Airy contour after the change
of variable \(w=i\xi\) and rotation within its cubic decay
sectors. The simple negative zeros of \(\mathrm{Ai}\) give
at least \(k\) sign changes of the shifted integral for \(v\)
in some fixed bounded positive interval, for all large \(t\).
Thus \(y_k(t)\le x_e+O_k(c)\), proving (3).
[DLMF Airy integral and zeros](https://dlmf.nist.gov/9.5.E4)
[DLMF zero distribution](https://dlmf.nist.gov/9.9.i)

Put \(\ell=\log t\), \(z=y/t^{3/2}\), \(q=p/\sqrt t\), and
use a dot for \(d/d\ell\). Equation (2) gives
\[
 \dot z=q-\tfrac32z,\qquad
 \dot q=t^2(2q^3-q/2-z/4)-q/2 .
\]
Initially \(\dot z<0\) as \(t\downarrow0\), since \(y(0)>0\)
and \(p\) is bounded there. At a first zero of \(\dot z\),
one has \(\ddot z\ge0\) and
\[
 \ddot z=z\,[t^2(27z^2/4-1)-3/4] .
\]
It follows that
\[
 z\ge b(t):=\sqrt{\frac4{27}\left(1+\frac{3}{4t^2}\right)} .
\]
If equality held, both \(\dot z\) and \(\dot q\) would vanish.
Differentiating the equation for \(\dot q\) then gives
\(\dddot z=3z/2>0\), forcing \(\dot z>0\) immediately to the
left by Taylor's formula. This contradicts the choice of the
first zero. Thus \(z>b(t)\) there and
\(\ddot z>0\). While \(\dot z>0\) afterward, \(z\) increases
and \(b(t)\) strictly decreases. At any later first zero of
\(\dot z\) the same formula would give \(\ddot z>0\), which
is incompatible with approaching zero from positive values
of \(\dot z\). So \(z\) would increase thereafter, staying
strictly above its limiting upper bound \(2/(3\sqrt3)\)
from (3). This contradiction proves \(\dot z<0\), giving (1).

## The bounded moving critical field

Set \(a=s/2\), \(\lambda=n^{-1/4}\), and use
\[
 a=2+\tau\lambda^2,\qquad \theta=v\lambda^3 .
\]
In the normalized Gaussian integral from the
[critical-coupling proof](MULTIVARIATE_DIFFUSE_CRITICAL_COUPLING.md#the-bounded-quartic-window),
the phase after removal of \(e^{in\theta/2}\) is
\[
 F_a(x,\theta)=-a/4-x^2/(4a)+
                  \log\!\left(2\cosh\frac{x+i\theta}{2}\right).
\]
Use the same branch and divide the integral by
\(\lambda\exp(\lambda^{-4}F_a(0,\theta))\), omitting its
nonzero Gaussian prefactor. Denote it by \(J(\tau,v,\lambda)\),
where \(\lambda>0\) may be real. After the change \(x=\lambda w\),
its limit is
\[
 J(\tau,v,0)=\int_{\mathbb R}
   \exp(-w^4/192+\tau w^2/16+ivw/4)\,dw .                     \tag{4}
\]
In particular, for \(\kappa=48^{1/4}\), this is
\(\kappa V(-\sqrt3\,\tau/2,\kappa v/4)\).

The convergence of \(J,\partial_vJ,\partial_\tau J\) is uniform
on compact real \((\tau,v)\)-rectangles, with \(\partial_\lambda J\)
uniformly bounded. To verify this, the quadratic contribution
from the change in \(a\) in the scaled exponent is exactly
\(\tau w^2/(8a)\). The other terms have the Taylor bounds of the
fixed-critical proof. In particular, on a fixed local interval in
\(x\), the real phase difference from \(F_a(0,\theta)\) is bounded
above by
\(-c x^4+C|a-2|x^2+C|\theta||x|\).
In the scaled integral this gives the envelope
\(\exp(-c w^4+C|\tau|w^2+C|v||w|)\). Each differentiated exponent
is bounded by a polynomial in \(w\), uniformly there. For the
quadratic contribution this is immediate from \(\tau w^2/(8a)\)
with \(a=2+\tau\lambda^2\). For the other contributions expand the
logarithm as in the fixed-critical estimate. In particular the
scaled fifth-order remainder has the form
\(\lambda w^5 R(\lambda w,\lambda^3 v)\) with bounded first
derivatives of \(R\). Derivatives of the endpoints of the local
interval after scaling are exponentially small. The tails have
strict loss, persisting from the exactly critical phase for small
\(\lambda\), and Gaussian decay. Their differentiated contributions
are exponentially small as well. This proves the claimed
differentiated convergence and bounds.

Write \(v_k(\tau)\) for an indexed positive root of (4).
We have the strict inequality
\[
                    3v_k(\tau)-2\tau v_k'(\tau)>0
                      \qquad(\tau\in\mathbb R).               \tag{5}
\]
For \(\tau<0\) this is (1), since the corresponding parameters
are \(t=-\sqrt3\,\tau/2\) and \(v_k=4y_k(t)/\kappa\).
For \(\tau\ge0\), \(v_k'(\tau)\le0\) follows from fixed-degree
Curie–Weiss monotonicity. In fact at fixed \(n\) the coupling
in the spin convention \(\exp(t(\sum\sigma_i)^2)\) is
\(t=s/(8n)\). The principal positive zero arguments decrease
with \(s\). In the critical scaling they therefore decrease
with \(\tau\) as well. For fixed \(\tau\) the first \(k\) of
these scaled arguments converge to the first \(k\) positive
zeros of (4). The uniform finite-window convergence, the simple
real limit zeros, and the nonzero value at zero give this
identification by counting up to any fixed point between two
successive limit zeros. Taking the limit of the monotonicity
inequality between two values of \(\tau\) proves that \(v_k\)
is nonincreasing. It is differentiable by simplicity, giving
the claim.
[Hou–Jiang–Newman, Theorem 2](https://doi.org/10.1007/s10955-023-03066-x)

Fix bounded \(\tau,v\) ranges. Only finitely many limit zero branches
can enter the field range. To see this, the indexed model roots are
continuous, and for each fixed \(\tau\) there is some index beyond
the range. A finite cover of the compact \(\tau\)-interval supplies
a uniform bound on that index. Choose a larger bounded field
interval containing all the finitely many relevant branches.
Uniform convergence and simplicity then give nearby real roots
\(v_k(\tau,\lambda)\), accounting for every root in the original
bounded range. Branches have a uniform positive separation on the
compact \(\tau\)-interval. The implicit-function formula gives
\(\partial_\tau v_k(\tau,\lambda)\to v_k'(\tau)\) uniformly and
bounded \(\partial_\lambda v_k\) on these branches.

For consecutive orders with the same physical coupling \(a\),
the scaled temperature follows
\(\tau(\lambda)=(a-2)/\lambda^2\). Between
\(\lambda=n^{-1/4}\) and \(\lambda=(n+1)^{-1/4}\) it remains in a
bounded interval and changes by \(O(1/n)\). The physical angle on
each branch satisfies
\[
 \frac{d}{d\lambda}\bigl(\lambda^3 v_k(\tau(\lambda),\lambda)\bigr)
   =\lambda^2\left(3v_k-2\tau\,\partial_\tau v_k
                                    +\lambda\partial_\lambda v_k\right)
   >0
\]
for all sufficiently small \(\lambda\). The strict sign follows
from the positive compact minimum of (5) for the finitely many
branches. Thus two consecutive angles on the same branch are
distinct. A common physical angle would make the two scaled
field coordinates differ by a factor \((1+1/n)^{3/4}\). Their
boundedness, the \(O(1/n)\) temperature change and continuity
of the separated model branches exclude distinct branch indices.
This proves the exclusion for bounded \(n^{3/4}\theta\).

## Arguments outside bounded field windows

If a sequence in the stated temperature range had common roots
with principal absolute argument bounded away from zero, the
[regular-arc estimate](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#localization-of-possible-common-roots-to-the-edge)
would exclude them. It remains to consider sequences with
\(\theta\to0\), \(n^{3/4}\theta\to\infty\). For \(a<2\) the edge
argument is \(O((2-a)^{3/2})=O(n^{-3/4})\), so these arguments
are above the edge for large \(n\). In all cases use the
left saddle \(-r+iy\) and \(u=\theta+y\) from the regular-arc
construction. Here \(r\) is bounded because
\(a=P_{\cos u}(r)\ge r\tanh(r/2)\). A subsequential positive
limit of \(r\) at \(a=2,\theta=0\) would imply \(u\to0\) by
\(\theta=u-(r/\sinh r)\sin u\), and then \(P_1(r)>2\) at the
positive limit, an impossibility. Thus \(r\to0\).
The saddle identities give
\[
 \begin{gathered}
 \cos u=a\sinh r/r-\cosh r,\qquad
 1-\cos u=2-a+(3-a)r^2/6+O(r^4),\\
 u^2=O(r^2+|a-2|),\qquad
 \theta=u-\frac r{\sinh r}\sin u
             =O\bigl((r^2+|a-2|)^{3/2}\bigr).
 \end{gathered}
\]
Thus \(\sqrt n\,r^2\to\infty\) and \((a-2)/r^2\to0\).
It follows that
\[
 u/r\longrightarrow1/\sqrt3,\qquad
 \theta/r^3\longrightarrow2/(9\sqrt3).
\]

Center at \(iy\), and set
\[
 F_r(\zeta)=-a/4-\frac{(iy+\zeta)^2}{4a}
               +\log\!\left(2\cosh\frac{\zeta+iu}{2}\right),
 \qquad F_s=F_r(-r).
\]
The scaled Taylor polynomial has the same limit as in the
fixed-critical proof:
\[
 \frac{192}{r^4}\bigl(F_r(rw)-F_r(0)\bigr)
   \longrightarrow -w^4+2w^2+\frac{4i}{\sqrt3}(3w-w^3)
\]
locally uniformly in \(w\). In the linear coefficient, the
change from critical \(a=2\) contributes \(O(|a-2|\,y)=o(r^3)\).
In the quadratic coefficient it contributes \(O(|a-2|)=o(r^2)\).
The higher coefficients use \(u/r\to1/\sqrt3\).
In particular \(\Im F_s\sim-r^4/(24\sqrt3)\).

To control the whole horizontal contour, let \(Y=(\Re\zeta)^2\).
The real log modulus is a constant plus
\[
 -Y/(4a)+\tfrac12\log(\cosh\sqrt Y+\cos u).
\]
Its derivative vanishes at \(Y=r^2\), and its second derivative
at the limiting parameters \(a=2,u=Y=0\) is \(-1/96\).
On a fixed small neighborhood the real loss is therefore
at least \(c(\zeta^2-r^2)^2\). At fixed distance outside this
neighborhood the horizontal has uniform strict loss by
continuity from the critical maximum. The far tail has Gaussian
decay. These are the estimates used in the fixed-critical
packet proof. In detail, with \(m=nr^4\),
\[
 \Psi_r(w)=\frac{F_r(rw)-F_s}{r^4},\qquad
 K_r(m)=\int_{-\infty}^0 e^{m\Psi_r(w)}\,dw ,
\]
the saddle at \(w=-1\) has uniform positive real negative
curvature and bounded third derivative. For \(|rw|\) small,
the Taylor bounds give \(|\Psi_r(w)|\le C(1+|w|^4)\), while
the real loss is at least \(c(w^2-1)^2\). Gaussian comparison
near \(w=-1\) and the strict loss elsewhere yield
\[
 \sqrt m K_r(m)=C_r(1+\eta_r(m)),\qquad
 \eta_r=O(m^{-1/2}),\qquad \eta_r'=O(m^{-1}),
\]
with nonzero \(C_r\) bounded away from zero. The derivative
estimate follows by inserting \(\Psi_r\) in the integral.
It is \(O((w+1)^2)\) near the saddle, and the same loss controls
the other pieces. On the original far tail a fixed loss
multiplied by \(m/r^4\) absorbs the rescaled lengths and
powers of \(r^{-1}\).

The left packet after removing \(e^{in\theta/2}\) is proportional
to \(r\sqrt n\, e^{nF_s}K_r(m)\). Reflection about the midpoint
\(iy\) makes the full normalized polynomial twice its real
part. At fixed physical \(a,\theta\) the packet's change in
argument from \(n\) to \(n+1\) is therefore
\[
 \Im F_s+O(r^4/m).
\]
Here \(m=(\sqrt n\,r^2)^2\to\infty\). The change is strictly
between \(-\pi\) and zero for large \(n\), so both packet values
cannot be purely imaginary. This contradicts a common root
in the last field regime, completing the theorem.
\(\square\)
