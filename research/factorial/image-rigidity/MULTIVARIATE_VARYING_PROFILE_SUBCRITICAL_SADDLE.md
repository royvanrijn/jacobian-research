# A uniform positive-axis saddle for varying subcritical profiles

Use the endpoint-normalized ordinary factorial coefficients and the
functions \(L_n,G,H,W\) from the
[varying-profile phase estimates](MULTIVARIATE_VARYING_PROFILE_PHASE.md).
In particular define
\[
 P_n(z)=\sum_{k=0}^n\binom nk \exp(L_n(k/n))z^k.
\]
All varying exponent rows have both entries at least one.

**Theorem.** Fix \(W_0<\infty\), \(0<\delta<1\), and a compact
positive scale interval \([c,C]\). Allow the profile and its
dimension to vary with \(n\), subject to \(W\le W_0\) and
\[
 \sup_{0<u<1}u(1-u)G''(u)\le1-\delta.                       \tag{1}
\]
There are \(\theta_0>0\) and \(n_0\), depending only on the four
bounds, such that \(P_n(z)\ne0\) for
\(n\ge n_0,\ c\le |z|\le C,\ |\arg z|\le\theta_0\).
In particular any consecutive common root with small argument
must avoid this scale interval after each of the two endpoint
normalizations. Condition (1) holds whenever \(W\le2(1-\delta)\).

For the last assertion, each term in \(G''\) has the form
\(w_i/(1+(2u-1)r_i)\) with \(|r_i|<1\). Therefore
\[
 \frac{u(1-u)}{1+(2u-1)r_i}
 \le \frac{1+|2u-1|}{4}\le\frac12.
\]

## The uniformly nondegenerate saddle

Write \(z=e^\omega\), \(\omega=w+\mathrm i\theta\), and define
\[
 \phi(u,\omega)
 =-u\log u-(1-u)\log(1-u)+G(u)+\omega u .
\]
On the real interval at \(\theta=0\), (1) gives
\(\phi_{uu}(u,w)\le-\delta/[u(1-u)]\le-4\delta\).
There is a unique maximum \(u_0\), which stays uniformly in the
interior of \((0,1)\). Indeed \(G'\) is uniformly bounded by
the formula for \(G''\) and \(G(0)=G(1)=0\), whereas the entropy
derivative tends to \(+\infty,-\infty\) at the two endpoints.
Write \(M=\phi(u_0,w)\). Strong concavity gives
\[
 \phi(u,w)\le M-2\delta|u-u_0|^2 \quad (0\le u\le1).        \tag{2}
\]

On a fixed complex rectangle around all such \(u_0\), staying
away from the vertical lines through \(0,1\), the phases and
their derivatives are uniformly bounded. This follows from the
varying-profile estimates and from the uniform lower bound on
\(\Re(1+(2u-1)r_i)\). Choose a small radius \(r>0\) so that
\(\Re\phi_{uu}\le-2\delta\) for \(|u-u_0|\le 3r\), uniformly
in all profiles and \(w\). The second derivative is independent
of \(\theta\). For small \(|\theta|\) there is precisely one
stationary point \(u_*\) in \(|u-u_0|<r/4\), with
\(|u_*-u_0|=O(|\theta|)\) uniformly. To verify this, write
\(\phi_u=i\theta+\phi_{uu}(u_0,w)(u-u_0)+O(|u-u_0|^2)\),
where the linear coefficient has magnitude at least \(4\delta\).
Rouché's theorem first gives one zero on the fixed disk and
then puts it in a disk of radius \(O(|\theta|)\).
It follows that \(|\Re\phi(u_*,\omega)-M|=O(|\theta|)\).

## The lattice sum and its contour

In the interior rectangle interpolate the summands by
\[
 F_n(u,\omega)=
 \frac{\Gamma(n+1)}{\Gamma(nu+1)\Gamma(n(1-u)+1)}
                   \exp(L_n(u)+n\omega u).
\]
It is analytic. Uniform Stirling on this rectangle and the
first phase estimate give
\[
 F_n(u,\omega)=n^{-1/2}A(u)e^{n\phi(u,\omega)}
                 \bigl(1+O(n^{-1})\bigr),\qquad
 A(u)=\frac{e^{H(u)}}{\sqrt{2\pi u(1-u)}} .                 \tag{3}
\]
The square root is positive on the real interval.
The analytic amplitude and its derivative are uniformly
bounded, and its modulus is uniformly bounded away from zero
on the chosen rectangles. [Uniform Stirling expansion](https://dlmf.nist.gov/5.11.E1)

Choose half-grid points \(a_n,b_n\), of the form \((k+1/2)/n\),
within \(1/n\) of \(u_0-r,u_0+r\). For all sufficiently large
orders, the terms with \(k/n\) outside this interval have
total modulus \(O(n e^{n(M-\eta)})\) for some uniform
\(\eta>0\). Indeed
\(\binom nk\le\exp(n[-u\log u-(1-u)\log(1-u)])\)
at \(u=k/n\), while \(H\) and the error in \(L_n-nG\)
are uniformly bounded. Apply (2). The factor \(e^{\mathrm i n\theta u}\)
has modulus one at real \(u\).

Complete a rectangle with these vertical sides and horizontal
sides at imaginary parts \(\pm h\), for a fixed small \(h>0\).
For real \(v\in[a_n,b_n]\) and \(|y|\le h\),
\[
 \Re\phi(v+\mathrm i y,\omega)
 \le \phi(v,w)+O(h^2)+|\theta|h .
\]
The first-order Taylor term at real \(v,w\) is purely imaginary.
Choose \(h\) and then \(\theta_0\) so that this increase is
smaller than both \(\eta/2\) and \(\pi h\). The vertical sides
retain a fixed phase gap by (2), after decreasing \(\eta\) if
needed. We also take \(\theta_0\) small enough that
\(|\Im u_*|<h\) and
\(|\Re\phi(u_*,\omega)-M|\) is less than half of both
the vertical gap and \(\pi h\).

By residues the part of the sum inside the rectangle is
\[
 \frac{n}{2\mathrm i}\oint F_n(u,\omega)\cot(\pi n u)\,du .
\]
On the half-grid walls the cotangent has modulus at most one.
On the horizontal sides it equals \(\mp\mathrm i+O(e^{-2\pi n h})\),
with the minus sign at height \(h\). Thus the segment
contributions other than the two constant-kernel horizontal
integrals are exponentially small relative to
\(e^{n\Re\phi(u_*,\omega)}\), using (3) and the phase estimates.
The two horizontal integrals enter with the same sign and
each has prefactor \(n/2\) when directed left to right.
Deform both to the horizontal line through the saddle.
The connecting segments along the vertical walls retain
the phase gap. Including the omitted lattice terms, this gives
\[
 P_n(e^\omega)
 =n\int_{a_n+\mathrm i\Im u_*}^{b_n+\mathrm i\Im u_*}
                  F_n(u,\omega)\,du
   +O\!\left(e^{n(\Re\phi(u_*,\omega)-\eta_1)}\right)
\]
for a uniform \(\eta_1>0\), absorbing polynomial factors
into a smaller exponential gap.

On this line \(\Re\phi_{uu}\le-2\delta\) and
\(\phi_u(u_*,\omega)=0\). Put \(u=u_*+x/\sqrt n\).
The two endpoints have distances at least a fixed positive
constant from \(u_*\) before rescaling. The real part of the
scaled exponent difference is at most \(-\delta x^2\).
Taylor expansion gives the Gaussian with quadratic coefficient
\(\phi_{uu}(u_*,\omega)/2\), with a remainder
\(O(n^{-1/2}|x|^3)\). Interpolating between the quadratic and the
exact exponent, both with real part at most \(-c_1x^2\) for a
uniform \(c_1>0\), bounds the integral of the error by
\(O(n^{-1/2})\). The same estimate follows for replacing the
amplitude by its saddle value from the uniform derivative bound.
Equation (3) therefore gives
\[
 P_n(e^\omega)=
 e^{n\phi(u_*,\omega)}
 \frac{A(u_*)\sqrt{2\pi}}{\sqrt{-\phi_{uu}(u_*,\omega)}}
 \bigl(1+O(n^{-1/2})\bigr).
\]
Use the square root on the right half-plane. The amplitude
and curvature bounds make the leading coefficient uniformly
bounded away from zero. This proves the theorem.

The proof concerns a positive-axis neighborhood on compact normalized
scales. Interference between saddles at other angles and the loss of
the uniform curvature margin in (1) require separate estimates.
