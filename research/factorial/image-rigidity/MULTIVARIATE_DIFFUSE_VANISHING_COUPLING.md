# Uniform eventual diffuse coprimality at vanishing coupling

Write
\[
 B_n(z;s)=\sum_{j=0}^n\binom nj e^{-s j(n-j)/(2n)}z^j .
\]

**Theorem.** There is an integer \(n_0\) such that
\(B_n(z;s)\) and \(B_{n+1}(z;s)\) are coprime for all
\(n\ge n_0\) and every \(s>0\). The threshold is unevaluated.

By the [critical transition theorem](MULTIVARIATE_DIFFUSE_CRITICAL_TRANSITION.md),
which is uniform on \(s\ge s_->0\), it suffices to exclude
a sequence of common roots with \(n\to\infty\) and \(s\to0\).
The argument below controls this regime after rescaling the
angular field and the saddle coordinate by \(\sqrt{s/2}\).

## Rescaled phase and arguments away from the edge

Let \(a=s/2\), \(\varepsilon=\sqrt a\). The roots are on the
unit circle by the circle theorem used in the
[regular-arc estimate](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#localization-of-possible-common-roots-to-the-edge).
By conjugation consider \(0\le\theta\le\pi\), and set
\[
 x=\frac{\pi-\theta}{\varepsilon},\qquad
 g_a(Z)=\frac{2}{\varepsilon}\sinh(\varepsilon Z/2).
\]
The centered Gaussian integral, with \(t=\varepsilon W\),
expresses \(e^{-in\theta/2}B_n(e^{i\theta};s)\), up to the
positive factor \(e^{-na/4}\varepsilon^n\sqrt{n/(4\pi)}\),
as the integral of
\[
 \left[e^{-W^2/4}i g_a(W-ix)\right]^n
\]
over real \(W\). The function \(g_a(Z)\) is analytic in \(a\)
and entire in \(Z\), with \(g_0(Z)=Z\). This limit gives the
Gaussian representation of a scaled Hermite polynomial.

Shift to \(W=z+iY\) and put \(V=x-Y\). The phase is
\[
 F(z)=-(z+iY)^2/4+\log(i g_a(z-iV)).
 \tag{1}
\]
For \(0<\varepsilon V<2\pi\) the argument of the logarithm
has positive real part on the horizontal. This remains true
in the analogous open strip of imaginary displacements.
The symmetry \(F(-\bar z)=\overline{F(z)}\) makes the
normalized polynomial twice the real part of the left
integral. At \(V=0\) the logarithm on the left is taken with
argument \(-\pi/2\) for the log factor, and the integer-order
integral retains the same symmetry. At bounded \(Y,V\) and
bounded imaginary displacements the tail has a uniform
Gaussian majorant for the base in (1). Indeed
\(|g_a(Z)|\le |Z|e^{\varepsilon|\Re Z|/2}\).
On the real horizontal with bounded \(V\),
\[
 |g_a(X-iV)|^2=\frac{2}{a}(\cosh(\varepsilon X)-\cos(\varepsilon V))
       \ge X^2,
\]
so the logarithm has at most quadratic growth in the phase
on the far tails. These estimates also justify deformation
on any of the horizontal strips used below.

The edge is \(\theta_*=b-\sin b\), where \(b=\arccos(a-1)\).
Define
\[
 V_c=\frac{\pi-b}{\varepsilon},\qquad Y_c=\frac{\sin b}{\varepsilon},
 \qquad x_*=V_c+Y_c .
\]
Both \(V_c\) and \(Y_c\) tend to \(\sqrt 2\).
Below the angular edge \(\theta<\theta_*\) one has \(x>x_*\);
the single stable saddle satisfies
\[
 V=\frac{\pi-u}{\varepsilon}>V_c,qquad
 Y=\varepsilon\cot(\varepsilon V/2),qquad x=V+Y.
 \tag{2}
\]
These are the stationary formulas
\(0\le u<b\), \(\theta=u-a\tan(u/2)\) in scaled coordinates.
Here \(0\le Y\le2/V\), and \(V_c\ge\sqrt2\). For sufficiently
large \(x\) this implies \(V\ge\pi\) and therefore
\(1+\cos u=1-\cos(\varepsilon V)\ge2a\), using
\(1-\cos t\ge2t^2/\pi^2\) on \([0,\pi]\).
For bounded \(x\) away from \(x_*\) on this side, (2) tends to
\(x=V+2/V\) and forces \(V\) to stay separated from \(\sqrt2\).
In either case \(a/(1+\cos u)\) is bounded strictly below one.

On the horizontal, the derivative of the real phase with
respect to \(X=\Re z>0\) is
\[
 \frac X2\left(\frac a{P_{\cos u}(\varepsilon X)}-1\right),
 \qquad P_C(r)=r(\cosh r+C)/\sinh r .
\]
The increasing property of \(P_C\) thus bounds the loss
from zero by \(-c X^2\) for a uniform \(c>0\).
The stationary curvature is bounded away from zero and
the phase derivatives through order three are uniformly
bounded on a fixed disk about zero. This holds also for
unbounded \(V\): the log argument at zero has modulus
\((2/\varepsilon)\sin(\varepsilon V/2)\ge2V/\pi\), and the
derivatives of \(g_a(z-iV)\) of positive order are uniformly
bounded on a fixed disk for \(0\le\varepsilon V\le\pi\).
Gaussian comparison gives a nonzero leading term with
relative error \(O(n^{-1/2})\). Hence the individual
polynomial has no root here.

Above the angular edge, use the two stationary points
from the regular-arc construction. Write
\[
 r=\varepsilon R,quad u=\pi-\varepsilon V,\quad
 Y=\frac{r}{\sinh r}\frac{\sin(\varepsilon V)}{\varepsilon},
 \quad x=V+Y,quad
 \frac{1-\cos(\varepsilon V)}{a}
     =\frac{1-\cosh(\varepsilon R)}a
                  +\frac{\sinh(\varepsilon R)}{\varepsilon R}.
 \tag{3}
\]
The separated saddles in (1) are at \(z=\pm R\).
The inequality \(a\ge r\tanh(r/2)\) gives bounded \(R\);
also \(0\le V\le V_c\). At \(a=0\), (3) becomes
\(R^2+V^2=2\), \(x=2V\). Thus if \(x\) stays separated from
\(x_*\) on this side, \(R\) has a positive lower bound.
The real horizontal phase tends there to a constant plus
\[
 -X^2/4+\tfrac12\log(X^2+V^2).
\]
It has a unique strict maximum on the left at \(-R\),
with strictly negative curvature \(-R^2/2\).
Uniform convergence of the bases gives strict modulus
loss off a fixed saddle neighborhood on bounded pieces,
including at \(V=0,X=0\) where the limiting base vanishes.
The Gaussian tail bound and local analytic convergence
near \(-R\) complete the uniform Gaussian estimate.
The imaginary part of the left saddle phase converges to
\[
 RV/2-\arctan(R/V),
\]
with endpoint value \(-\pi/2\) when \(V=0\).
This lies in \([-\pi/2,0)\) and stays separated from zero
for \(R\) bounded below: at \(V>0\) the first term equals
\((R/V)/(1+(R/V)^2)\). At fixed \(a,\theta\) the left
packet's change in argument at consecutive orders is
\(\Im F(-R)+O(n^{-1/2})\), since the nonzero Gaussian
amplitude is the same. It is separated from integer
multiples of \(\pi\), so the two normalized polynomials
cannot both vanish.

## Uniform estimates through the rescaled edge

It remains to treat \(x-x_*\to0\). The edge phase, centered
at \(iY_c\), with \(\eta=x_*-x\), is
\[
 F_a(z,\eta)=-(z+iY_c)^2/4
                    +\log(i g_a(z-iV_c+i\eta)).
\]
It is analytic on a fixed neighborhood of \((z,\eta)=(0,0)\)
including \(a=0\): \(V_c=\arccos(1-a)/\sqrt a\) and
\(Y_c=\sqrt{2-a}\) extend analytically there. For \(V_c\), this
follows by inverting the analytic function
\((1-\cos\sqrt{y})=y/2+O(y^2)\) near zero and taking the
positive square root of \(y/a\). At \(\eta=0\) the first
two \(z\)-derivatives vanish. The mixed derivative, cubic
coefficient, and quartic coefficient at the origin are
\[
 (F_a)_{z\eta}=i/2,\qquad
 -i A_a=-iY_c/12,\qquad -(3-a)/48 ,
 \tag{4}
\]
respectively. These follow by differentiating the log cosh
phase before the rescaling, with one factor of
\(\varepsilon\) for each derivative in \(z\) or \(\eta\).
The edge real phase has strict loss on the old horizontal
outside every fixed small saddle neighborhood, uniformly
in small \(a\). On bounded horizontal sets this follows by
continuity from \(-X^2/4+\tfrac12\log(X^2+2)\), and at
infinity the Gaussian tail bound applies.

For the bounded Airy window \(\eta=n^{-2/3}v\), take the
rays of arguments \(-5\pi/6,-\pi/6\) to a small fixed
radius, with circular arcs joining the old horizontal.
By (4), the real loss on the rays at zero field displacement
is at least \(c|z|^3\). On the arcs the cubic term is
negative away from the real endpoints, and the quartic
term is strictly negative close to those endpoints.
This gives uniform strict loss on both arcs, as in the
[fixed subcritical Airy proof](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#the-bounded-airy-window-in-the-subcritical-range).
Set \(\lambda=n^{-1/3}\), \(z=\lambda\xi\).
After division by the nonzero centered exponential and
the length factor \(\lambda\), the ray exponent is
\[
 \lambda^{-3}(F_a(\lambda\xi,\lambda^2v)
              -F_a(0,\lambda^2v)).
\]
It converges uniformly on bounded \(\xi,v\) to
\(iv\xi/2-i A_a\xi^3\) and its real part is at most
\(-c|\xi|^3+C|v||\xi|\) on the rays. The first
\(\lambda\)- and \(v\)-derivatives of this exponent
are bounded by \(C(1+|\xi|^4)\) on the scaled ray segments.
Indeed the difference from the limit is a sum of
\(\lambda(v\xi^2 R_2+\xi^4 R_4)\) and
\(\lambda^2(v^2\xi R_1+v\xi^3 R_3)\), with the remainders
and their derivatives uniformly bounded in the fixed
analytic neighborhood. The moving endpoint terms, the
old horizontal outside this neighborhood and the arcs are
exponentially small after normalization and differentiation,
by their strict loss and the Gaussian tail. The normalized
integral and its \(v\)-derivative therefore converge uniformly
to the Airy integral \(2\pi c_a^{-1}\mathrm{Ai}(-v/(2c_a))\)
and its derivative, where \(c_a=(3A_a)^{1/3}\).
Its \(\lambda\)-derivative is bounded. These estimates give a
bounded \(\lambda\)-derivative on each simple real zero branch
\(v_k(a,\lambda)\) in a fixed bounded window, by the implicit
function theorem. The scaled zeros are positive and uniformly
separated. The physical angle
\(\theta_*+\varepsilon\lambda^2 v_k(a,\lambda)\) is strictly
increasing with \(\lambda\), and distinct branches cannot
match at consecutive orders because scaled coordinates
at the same physical angle differ by the factor
\((1+1/n)^{2/3}\). This treats bounded \(v\), allowing a
larger convergence interval to contain both scaled values.

Below the edge in angle with unbounded Airy distance, let
\(h=V-V_c>0\) in (2). The derivative in \(V\) of the right
side for \(x\) vanishes at \(V_c\), and its second derivative
has a positive limit there. Therefore \(x-x_*\asymp h^2\).
The scaled negative curvature
\[
 d=\tfrac12-a/(2(1-\cos(\varepsilon V)))\asymp h
]
at the horizontal saddle is positive. The real loss is
at least \(c(d X^2+X^4)\) locally, from the negative quartic
term in (4) and evenness of the real phase on the horizontal.
The cubic approximation error to the quadratic phase is
\(O(|X|^3)\); the relative Gaussian error is consequently
\(O((nh^3)^{-1/2})\), tending to zero. The retained far
pieces have the same strict loss and Gaussian tail as at
the edge. Thus the polynomial does not vanish on this side.

Above the edge in angle, (3) gives \(V,Y\) as analytic
functions of \(R^2\) near zero, uniformly including \(a=0\).
The limiting derivative of \(x\) in \(R^2\) at zero is
\(-1/\sqrt2\), since at \(a=0\), \(x=2\sqrt{2-R^2}\).
Thus \(\eta=x_*-x\asymp R^2\). For \(m=nR^3\) bounded,
the Airy estimate applies. For \(m\to\infty\) use the
[diagonal deformation for a cubic packet](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#intermediate-arguments-above-the-subcritical-edge)
in the rescaled \(z\) coordinate. Its uniform hypotheses
follow from (3) and (4): if \(F\) is (1) with the actual
stationary parameters, put
\[
 M(z)=i A_a(3R^2z-z^3),\qquad E=F-M .
\]
Then \(E'(-R)=0\), \(E''(z)=O(R^2+|z|^2)\) on a fixed
small disk, and the left stationary phase satisfies
\(\Im F(-R)=-2 A_a R^3+O(R^4)\).
Here \(F(0)\) is real, and integration of the bound on
\(E''\) controls \(E(-R)-E(0)\). The same bound supplies
the diagonal and endpoint losses of that deformation
relative to \(F(-R)\), with \(R\) in place of its saddle
separation. Its fixed-radius connector has uniform strict
loss by (4), and the outer horizontal has the strict
loss and Gaussian tail just established. Thus the scaled
left packet, after its stationary exponential and real
prefactors are removed, is \(C_{a,R}+O(m^{-1/2})\) with
\(C_{a,R}\) bounded away from zero. Its \(m\)-derivative
is \(O(m^{-1/2})\): the derivative inserts the scaled
exponent, which vanishes quadratically at the saddle,
is bounded on the final segment with linear loss, and
has polynomial growth on the tails. Hence the argument
change at consecutive orders is
\(-2 A_a R^3+O(R^4+R^3m^{-1/2})\), strictly between
\(-\pi\) and zero. Both left packets cannot be purely
imaginary. This completes the exclusion at the rescaled
edge and proves the theorem. \(\square\)
