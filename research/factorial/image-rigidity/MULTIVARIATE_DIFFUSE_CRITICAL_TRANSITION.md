# Diffuse coprimality through the critical transition

Write
\[
 B_n(z;s)=\sum_{j=0}^n\binom nj e^{-s j(n-j)/(2n)}z^j .
\]

**Theorem.** For each \(s_->0\), there is an order \(n_0(s_-)\)
such that \(B_n(z;s)\) and \(B_{n+1}(z;s)\) are coprime whenever
\(n\ge n_0(s_-)\) and \(s\ge s_-\). The same \(s\) is used in
both polynomials. The proof leaves the threshold unevaluated.

The [compact subcritical and supercritical theorems](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md)
and the [bounded moving critical window](MULTIVARIATE_DIFFUSE_MOVING_CRITICAL_WINDOW.md)
reduce the theorem to sequences
\[
 a=s/2\to2,\qquad n(a-2)^2\to\infty.
\]
Indeed, coupling \(s\ge64\) is already covered at every order, and
any remaining sequence of common roots has a subsequence in one of
these regimes or the earlier ones. By the circle theorem and
[regular-arc localization](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#localization-of-possible-common-roots-to-the-edge),
the principal absolute argument of any such sequence tends to zero.
It is positive because all coefficients are positive.

For small positive \(\theta\), use the real centered integral
for \(e^{-in\theta/2}B_n(e^{i\theta};s)\) with phase
\[
 f(t)=-a/4-t^2/(4a)+\log\!left(2\cosh\frac{t+i\theta}{2}\right).
\]
The prefactor is \(\sqrt{n/(4\pi a)}\). In the strip
\(|\Im t+\theta|<\pi\) the cosh factor has positive real part.
We use its principal logarithm and deform inside the strip,
with Gaussian decay on the closing sides. This also defines
the integrals at positive real orders. On a horizontal line
through \(iy\), set \(F(\zeta)=f(iy+\zeta)\) and \(u=y+\theta\).
The symmetry \(F(-\bar\zeta)=\overline{F(\zeta)}\) makes the
full integral twice the real part of its left half.

In any region with two separated saddles use the stationary
parameters from the localization proof,
\[
 t_L=-r+iy,quad a=\frac{r(\cosh r+\cos u)}{\sinh r},\quad
 y=\frac{r}{\sinh r}\sin u,\quad \theta=u-y ,
\]
with \(r>0\). The left saddle has value \(F_s=F(-r)\).
These parameters tend to zero along the sequences under
consideration. In fact \(r\) stays bounded by
\(a\ge r\tanh(r/2)\) and \(u\in[0,\pi]\). A subsequential positive
limit of \(r\) would force \(u\to0\) by the formula for \(\theta\)
and then \(a\to r\coth(r/2)>2\). Once \(r\to0\), the same formula
implies \(u\to0\).

Two estimates hold uniformly for small \(r,u\). First, with
\(Y=x^2\) on the horizontal line, the real part of \(F(x)\)
is a constant plus
\[
 -\frac{Y}{4a}+\frac12\log(\cosh\sqrt Y+\cos u).
\]
Its \(Y\)-derivative vanishes at \(Y=r^2\), and its second
derivative at \(a=2,u=Y=0\) is \(-1/96\). For some fixed small
\(\delta>0\), this gives
\[
 \Re(F(x)-F_s)\le-c(x^2-r^2)^2\qquad (|x|\le\delta).
 \tag{1}
\]
For \(|x|\ge\delta\) it is at most \(-c(1+x^2)\) with a possibly
smaller \(c>0\), by continuity from the strict critical maximum
and the uniform Gaussian bound at infinity. Second, stationarity
gives \(y/a=\sin u/(\cosh r+\cos u)\), so
\[
 \begin{split}
 \Im F_s&=-\frac{\sin u}{2}
   \int_0^r\left(\frac{1}{\cosh x+\cos u}
                -\frac{1}{\cosh r+\cos u}\right)\,dx,\\
 -\Im F_s&\asymp u r^3\qquad(u>0).
 \end{split}                                                    \tag{2}
\]
The comparison follows by bounding the two denominators above
and below and using \(\cosh r-\cosh x\asymp r^2-x^2\).
All comparisons below are uniform in their indicated limits.

## Separated saddles with bounded \(u/r\)

Suppose \(u/r\) stays bounded and \(m=nr^4\to\infty\).
Write \(v=u/r\) and temporarily vary \(v\) in a fixed symmetric
bounded interval at fixed \(r\), defining \(a,y,\theta\) by the
stationary formulas. This preserves the saddle at \(\zeta=-r\)
and gives \(F_{-v}(x)=\overline{F_v(x)}\) for real \(x\). Scale
the left horizontal integral by \(w=\zeta/r\) and put
\[
 \Psi(w)=\frac{F(rw)-F_s}{r^4},\qquad
 J(m,v)=\sqrt m\int_{-\infty}^0 e^{m\Psi(w)}\,dw .
\]
For \(|rw|\le\delta\), the functions \(\Psi,\partial_v\Psi\)
are \(O(1+|w|^4)\) and are both \(O((w+1)^2)\) near the saddle.
Indeed, \(F'(0)=O(r^3)\), \(F''(0)=O(r^2)\),
\(F'''(0)=O(r)\), and \(F^{(4)}\) is bounded on a fixed small
disk. These estimates hold after one \(v\)-derivative as well:
the first three derivatives are
\[
 \begin{split}
 F'(0)&=\frac{i\sin u}{2}\left(\frac1{1+\cos u}-\frac1{\cosh r+\cos u}\right),\\
 F''(0)&=\frac1{2(1+\cos u)}-\frac1{2a},\\
 F'''(0)&=-\frac{i\sin u}{2(1+\cos u)^2}.
 \end{split}
\]
Here \(a-(1+\cos u)=O(r^2)\) with the same derivative bound.
The local quadratic vanishing follows from stationarity at
fixed \(w=-1\) for the entire \(v\)-family.

By (1), the real negative curvature of \(\Psi\) at \(-1\)
is bounded below. Its third derivative is locally bounded,
so the Gaussian approximation gives
\[
 J(m,v)=\sqrt{2\pi/(-\Psi''(-1))}+O(m^{-1/2})
\]
with the square root of positive real part. The amplitude
is bounded away from zero. On the local segment this error
follows by replacing the exponential by its quadratic part,
with error bounded by a constant times
\(\sqrt m\int m|w+1|^3e^{-cm(w+1)^2}\,dw\).
Away from that segment, (1) supplies strict loss and then
quartic decay in \(w\) until \(|rw|=\delta\). Farther out,
the loss is \(-c(1+r^2w^2)/r^4\). The latter exponential
absorbs all rescaling factors and the at-most-quadratic
growth in \(rw\) of \(F,F_v\).

The same bounds after differentiation give
\[
 J_m=O(m^{-1}),\qquad J_v=O(1),\qquad J_{mv}=O(m^{-1}).
\]
In the local integral the new factors are \(\Psi\),
\(m\Psi_v\) and \((1+m\Psi)\Psi_v\), respectively, together
with the derivative of \(\sqrt m\). The quadratic vanishing
and the Gaussian bound give the stated orders. The nonlocal
bounds justify the differentiations and absorb those factors.
It follows that \(\partial_v(J_m/J)=O(m^{-1})\). At \(v=0\),
\(J_m/J\) is real by reflection, so
\(\Im(J_m/J)=O(v/m)\).

At the two physical orders the left integral is proportional
to \(e^{nF_s} J(nr^4,v)\) with a positive factor. Hold the
physical \(a,\theta\), hence \(r,v\), fixed in the comparison.
Its change in argument is
\[
 \Im F_s+O(r^4v/m).
\]
By (2) this lies strictly between \(-\pi\) and zero for
large \(m\) and small positive \(u,r\). The two left integrals
therefore cannot both be purely imaginary, excluding common
roots in this regime.

## Separated saddles with \(r/u\to0\)

Suppose now that \(r/u\to0\) and \(m=nu r^3\to\infty\).
The horizontal real curvature alone does not give the required
cubic scale. Use
\[
 A=\frac{\sin u}{12(1+\cos u)^2}\asymp u,\qquad
 M(\zeta)=iA(3r^2\zeta-\zeta^3),\qquad E=F-M .
\]
On the fixed small disk, \(E''(\zeta)=O(r^2+|\zeta|^2)\) because
the cubic term agrees exactly and \(F''(0)=O(r^2)\).
Since \(E'(-r)=0\), integration gives
\[
 |E(\zeta)-E(-r)|\le C|\zeta+r|^2(r^2+|\zeta|^2),\qquad
 |E(\zeta)-E(0)|\le Cr^3|\zeta|\quad (|\zeta|\le r).
 \tag{3}
\]

Take a small fixed \(\epsilon>0\) and set \(R=\epsilon u\).
In \(w=\zeta/r\) coordinates, approach \(-1\) on
\(w=-1-(1+i)\tau\) from its intersection with \(|\zeta|=R\).
Leave on \(w=-1+(1+i)\tau\), \(0\le\tau\le1/2\), then take
\(w=(-1+i)\tau\) with \(\tau\) decreasing from \(1/2\) to zero.
Join the first segment to the horizontal by the short circle
arc from \(-R\). This deforms the left half without changing
the full value at integer orders. On the approach and departure
segments the real model differences from the saddle are
\(-A r^3(6\tau^2+2\tau^3)\) and
\(-A r^3(6\tau^2-2\tau^3)\), respectively. The ratio of the
error in (3) to these losses is bounded by
\(C(r/u)(1+\tau)\) on the approach and \(C r/u\) on the
departure. These ratios are small by the choice of \(\epsilon\)
and \(r/u\). On the endpoint segment the real model difference
from zero is \(-A r^3\tau(3-2\tau^2)\), and the perturbation
is \(O(r^4\tau)\). Since \(\Re F(0)\le\Re F_s\), the last
segment has loss at least \(c u r^3\tau\) from the saddle.

For the arc, reflect to \(\zeta=Re^{i\phi}\) with
\(-\pi/4\le\phi\le0\). The linear coefficient of \(F\) is
imaginary and \(O(u r^2)\), its quadratic coefficient is
real and \(O(r^2)\), its cubic coefficient is \(-iA\),
and its quartic coefficient tends to \(-1/192\).
For \(-\pi/12\le\phi\le0\) the negative cubic term in the
real part controls the linear term, using
\(|\sin\phi|\le C|\sin3\phi|\) and \(r/R\to0\).
The quartic term then gives loss at least \(cR^4\),
absorbing the quadratic and fifth-order remainder.
On the rest of the arc \(\sin3\phi\le-\sqrt2/2\), so the
cubic loss controls the other terms when \(\epsilon\) is small.
The loss applies relative to \(F(0)\), hence to \(F_s\).
The normalized exponent \(\Psi=(F(rw)-F_s)/(u r^3)\) on the
arc has magnitude \(O((u/r)^3)\) and real loss at least
\(c(u/r)^3\).

On the retained horizontal \(R\le|x|\le\delta\), (1) gives
\(-\Re\Psi\ge c(r/u)|w|^4\ge c\epsilon|w|^3\).
The Taylor bound gives \(|\Psi|\le C(1+|w|^4)\) there.
For \(|x|\ge\delta\), the loss in (1) outside the small
interval, divided by \(u r^3\), absorbs the rescaled length
and the quadratic growth. On the diagonals (3) bounds
\(|\Psi|\) by \(C(\tau^2+\tau^3)\), and on the endpoint
segment it stays bounded.

Write \(J(m)=\sqrt m\int e^{m\Psi}\,dw\) on this contour.
Locally \(w=-1+(1+i)\tau\), the quadratic negative curvature
in \(\tau\) is \(12A/u+O(r/u)\), bounded above and with
positive real part bounded below. The scaled third derivative
is bounded. The Gaussian comparison thus gives
\(J(m)=C_{r,u}+O(m^{-1/2})\), where \(C_{r,u}\) is the
nonzero Gaussian amplitude bounded away from zero. The final
segment contributes \(O(m^{-1/2})\) after normalization by
its linear loss. The nonlocal diagonals, arc and horizontal
pieces are exponentially small by the displayed estimates.
Inserting \(\Psi\) in the same integral gives \(J'(m)=O(m^{-1/2})\).
Here the differentiated local contribution is \(O(m^{-1})\)
by quadratic vanishing, and the differentiated endpoint
piece is \(O(m^{-1/2})\) by its linear loss. All constants
are uniform as \(r/u\to0\) and \(u\to0\).

At fixed physical parameters the change in left argument from
\(n\) to \(n+1\) is now \(\Im F_s+O(u r^3m^{-1/2})\).
Equation (2) again rules out a common root.

## Subcritical angles in the shrinking Airy window

For \(a<2\), set \(b=arccos(a-1)>0\),
\(\theta_* = b-\sin b\) and \(\lambda=(n b^4)^{-1/3}\).
Then \(n b^4\to\infty\) in the remaining subcritical regime.
We first exclude common roots when
\[
 \eta=\theta-\theta_* = b^3\lambda^2 v
\]
with \(v\) bounded. Center the phase at \(i\sin b\) and set
\[
 \begin{split}
 f_b(\zeta,\eta)
 &=-a/4-(i\sin b+\zeta)^2/(4a)
      +\log\!left(2\cosh\frac{\zeta+i(b+\eta)}{2}\right),\\
 G_b(w,t)&=b^{-4}\bigl(f_b(bw,b^3t)-f_b(0,b^3t)\bigr).
 \end{split}
\]
On a fixed neighborhood of \((w,t)=(0,0)\) this extends
analytically to \(b=0\). The numerator is even in \(b\), and
its quadratic term vanishes by cancellation of the Gaussian
and log cosh terms at \(a=2\). The constant term is zero.
At \(t=0\) the first two \(w\)-derivatives of \(G_b\) vanish,
the mixed derivative \(\partial_w\partial_t G_b(0,0)=i/(2a)\),
and the cubic coefficient is \(-i A_b\), where
\[
 A_b=\frac{\sin b}{12b(1+\cos b)^2}\longrightarrow\frac1{48}.
\]
The quartic coefficient tends to \(-1/192\).

In the scaled \(w\)-plane, use the two rays of arguments
\(-5\pi/6,-\pi/6\) up to a small fixed radius \(\epsilon\),
connected to the retained horizontal by arcs of that radius.
At \(t=0\) the rays have real loss at least \(c|w|^3\).
The arcs have uniform strict loss. On the right arc the cubic
real part is \(A_b\epsilon^3\sin3\phi\), negative uniformly
away from \(\phi=0\). On \(-\pi/12\le\phi\le0\) it is
nonpositive and the quartic real part is negative with
\(\cos4\phi\ge1/2\). Uniform Taylor remainders and reflection
give the two arc bounds for small fixed \(\epsilon\).

On the horizontal part \(\epsilon\le|w|\le\delta/b\), the
real loss at \(t=0\) is at least \(c w^4\). This is (1)
with its critical point at \(x=0\) for \(u=b\), \(a=1+\cos b\).
Also there \(|G_b(w,0)|\le C(|w|^3+|w|^4)\) by Taylor
expansion in the physical displacement. The perturbation
by \(t\) is bounded by \(C|t||w|\), because
\[
 \partial_\eta(f_b(\zeta,\eta)-f_b(0,\eta))
  =\frac i2\left(\tanh\frac{\zeta+i(b+\eta)}2
                      -\tanh\frac{i(b+\eta)}2\right)
\]
is \(O(|\zeta|)\) in the strip neighborhood in use.
Farther out, for \(|\zeta|\ge\delta\), the loss is at least
\(c(1+|\zeta|^2)/b^4\) in \(G_b\), uniformly for small \(t\),
by the same strict critical bound and Gaussian tail as in (1).
There \(G_b\) grows at most quadratically in \(\zeta\) with
coefficient \(O(b^{-4})\).

Divide the shifted integral by
\(b\lambda e^{n f_b(0,\eta)}\), omitting its real Gaussian
prefactor. This defines a real function \(J_b(v,\lambda)\)
for positive real \(\lambda\), with \(n=(b^4\lambda^3)^{-1}\).
On the rays put \(w=\lambda\xi\).
The scaled exponent \(\lambda^{-3}G_b(\lambda\xi,\lambda^2v)\)
has limit \(iv\xi/(2a)-iA_b\xi^3\) as \(\lambda\to0\).
It is dominated in real part by \(-c|\xi|^3+C|v||\xi|\).
For clarity, its difference from the limiting exponent
is a sum of terms of the form
\[
 \lambda\bigl(v\xi^2 R_2+\xi^4 R_4\bigr)
    +\lambda^2\bigl(v^2\xi R_1+v\xi^3 R_3\bigr).
\]
The analytic remainders and their first partial derivatives
are bounded uniformly, with arguments \((\lambda\xi,\lambda^2v)\)
for \(R_4\) and \(\lambda^2v\) for the other terms. Hence the
\(\lambda\)- and \(v\)-derivatives of the exponent are bounded
by \(C(1+|\xi|^4)\) on \(|\lambda\xi|\le\epsilon\), for bounded \(v\).
The moving endpoint terms vanish exponentially. The arc and
horizontal estimates likewise give exponentially small terms
for the normalized integral and its two derivatives.
On the retained horizontal up to \(\delta/b\) this follows
using \(e^{-c|w|^4/\lambda^3}\), absorbing powers of \(w\)
and \(\lambda^{-1}\). Beyond \(\delta/b\), the loss
\(e^{-c(1+|\zeta|^2)/(b^4\lambda^3)}\) absorbs the powers
of \(b^{-1}\) as well.

It follows that \(J_b\) and \(\partial_v J_b\) converge
uniformly to the Airy integral and its derivative, and
\(\partial_\lambda J_b\) stays bounded. The limit for each
\(b\) is
\[
 J_b(v,0)=\frac{2\pi}{c_b}\mathrm{Ai}\!
                  \left(-\frac{v}{2a c_b}\right),\qquad
 c_b=(3A_b)^{1/3},
\]
by the [same ray integral as at fixed subcritical coupling](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#the-bounded-airy-window-in-the-subcritical-range).
The positive scaled zeros in any bounded interval are
simple and uniformly separated, since \(2a c_b\) tends to
a positive constant. Near each of these zeros the finite-order
branch \(v_k(b,\lambda)\) has bounded \(\lambda\)-derivative,
by the implicit function theorem and the derivative bounds.
The physical angle is strictly increasing in \(\lambda\)
on that branch, since its derivative is
\[
 b^3\lambda(2v_k+\lambda\partial_\lambda v_k)>0.
\]
Distinct branches cannot give the same physical angle at
consecutive orders, because the two scaled \(v\) coordinates
of a fixed angle differ by a factor \((1+1/n)^{2/3}\).
Taking the convergence interval larger than the window in
question accounts for both orders. This excludes common roots
in every bounded scaled window.

## Angular coverage

For \(a>2\), all positive arguments use the two separated
saddles. The stationary identities give
\[
 1-\cos u=2-a+(3-a)r^2/6+O(r^4).
 \tag{4}
\]
Consequently \(u=O(r)\) and \(r^2\ge c(a-2)\). Thus \(nr^4\)
diverges and the bounded \(u/r\) argument applies, including
sequences with \(u/r\to0\).

For \(a<2\), first take \(\theta<\theta_*\). The single
horizontal saddle is at \(i a\tan(u/2)\), where \(0\le u<b\)
and \(\theta=u-a\tan(u/2)\). Put \(h=b-u\).
Since \(1-a/(1+\cos u)\asymp b(b-u)\) on \(0\le u\le b\),
the angle deficit is \(\theta_* -\theta\asymp b h^2\).
The magnitude of the negative curvature is
\[
 d=\frac1{2a}-\frac1{2(1+\cos u)}\asymp b h.
\]
The horizontal loss near that saddle is at least
\(c(d x^2+x^4)\), using the same negative second derivative
in \(Y=x^2\) as in (1). The difference from the quadratic
approximation is \(O(u|x|^3+|x|^4)\) because the cubic
coefficient is \(O(u)\). The relative Gaussian replacement
error is therefore
\[
 O\!\left(\frac{u}{\sqrt{n d^3}}+\frac1{n d^2}\right)
\]
by integration against \(e^{-cndx^2}\).
The fixed-distance loss and Gaussian tail again bound the
remaining contour. Outside the bounded Airy window on this
side, pass to a subsequence with \(h/(b\lambda)\to\infty\).
Then \(nbh^3\to\infty\); the displayed error tends to zero
by \(u\le b\), \(d\asymp b h\), and \(h\le b\).
The integral is nonzero.

At \(\theta=\theta_*\) the Airy result applies. Above the
edge, \(u>b\) and the two separated saddles exist. If \(u/r\)
is bounded then \(r\gtrsim u\ge b\) and \(nr^4\to\infty\).
Otherwise pass to a subsequence with \(r/u\to0\).
By (4), \(u\sim b\) and \(u-b\asymp r^2/b\). The exact
angle difference is
\[
 \theta-\theta_*=(u-b)-(\sin u-\sin b)
                  +(1-r/\sinh r)\sin u\asymp b r^2 .
\]
Here the first term group is nonnegative and \(O(b r^2)\),
and the last term is comparable to \(b r^2\). If \(nu r^3\)
is bounded, the scaled Airy distance
\((\theta-\theta_*)/(b^3\lambda^2)\) is bounded as well.
If it tends to infinity, the cubic packet comparison applies.
The two alternatives exhaust subsequences. This excludes
common roots in the entire remaining critical transition
and proves the theorem. \(\square\)
