# A computable cutoff at vanishing shared curvature

**Theorem.** There are computable positive rational \(W_0\) and
integer \(N_0\) such that every pair of consecutive ordinary
factorial moments with terminal \(N\ge N_0\) and
\[
               0<W=\sum_i (b_i-a_i)^2/((a_i+b_i)/2)\le W_0
\]
is coprime. The total-degree difference and the number of rows
are arbitrary. The constants are specified by the finite bounds
below, and have not been numerically optimized or evaluated.

This makes the
[vanishing-curvature comparison](MULTIVARIATE_VANISHING_CURVATURE_EVENTUAL.md)
effective. It uses the exact midpoint densities, including their
complex order dependence, and a separate quantitative Airy
localization at the only limiting fold.

## Exact midpoint normalization and effective density bounds

Delete equal rows, write \(m_i=(a_i+b_i)/2\), \(d_i=b_i-a_i\),
\(\Delta=\sum_i d_i\), and put \(\varepsilon=\sqrt{W/2}\).
Take \(\varepsilon\le1/2\), so all unequal rows are shared.
For real \(j\), let \(X_i\) be independent unit-rate Gamma
variables of shapes \(jm_i+1\). Set
\[
 U_j=
  \frac{\sum_i d_i(\log X_i-\log(jm_i))}{\varepsilon},
 \qquad g_\varepsilon(z)=\frac{2\sinh(\varepsilon z/2)}{\varepsilon}.
\]
Use \(t=-\exp(-\sum_i d_i\log(jm_i)+\varepsilon\xi_j)\).
The normalized exact moment at integer \(j\) is
\(\mathbb E_j[g_\varepsilon(U_j+\xi_j)^j]\).
For \(n=N-1\), linked fields satisfy
\[
                   \xi_N-\xi_n=\Delta\log(N/n)/\varepsilon\in\mathbb R .
                                                               \tag{1}
\]
At a root one can take \(|\Im\xi_j|\le\pi/\varepsilon\).
Also \(|\Re\xi_j|\le8\varepsilon\).
Indeed \(V=\sum d_i^2/\min(a_i,b_i)\le2W\), and the
[shared annulus](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md#shared-rows-and-large-degree-imbalance)
gives \(|\log|G_jt||\le V/2\).
For the endpoint scale \(\log G_j\), the integral of
\(\psi(1+jx)\) from \(a_i\) to \(b_i\) differs from
\(d_i\log(jm_i)\) by at most
\(d_i^2/(4\min(a_i,b_i))+|d_i|/(j\min(a_i,b_i))\).
Use \(0\le\psi(1+y)-\log y\le1/y\) and integrate
\(|x-m_i|/\min(a_i,b_i)\) on the row interval.
Thus \(\varepsilon|\Re\xi_j|\le(3/4+1/j)V\),
which gives the stated bound.

The following bounds on the analytic density hold on every
fixed bounded \(s\)-disk. All the constants, and the required
upper bound on \(\varepsilon\) and lower bound on \(n\), can
be computed from the disk. One may take a complex order disk
\(|j-n|\le\eta n\) with fixed computable \(\eta>0\).
There are analytic \(A(s)\) and \(a_j(s)\) such that
\[
  p_j(s)=\sqrt j\, e^{j A(s)}a_j(s),\qquad
  A(s)=-s^2/4+O(\varepsilon),\qquad
  a_j(s)=a_0(s)+O(1/n),\quad
  a_0(s)=1/\sqrt{4\pi}+O(\varepsilon).
                                                               \tag{2}
\]
Here \(p_j\) agrees with the density of \(U_j\) at positive
real \(j\). Bounds for any fixed number of derivatives in \(s\)
follow by Cauchy on larger disks.

For details, set \(e_i=d_i/\varepsilon\), \(u_i=e_i/m_i\).
Integrality gives
\[
 \max |u_i|\le2\varepsilon,\quad \sum|u_i|\le2\varepsilon,
 \quad\sum m_i u_i^2=2,\quad \sum1/m_i\le W .
\]
The Gamma moment function is
\[
 T_j(h)=\prod_i
      \frac{\Gamma(jm_i+1+e_i h)}{\Gamma(jm_i+1)}
                         (jm_i)^{-e_i h}.
\]
On each fixed \(t\)-disk its sectorial Stirling expansion reads
\[
 \log T_j(jt)=jK(t)+J_1(t)+O(W/n),\qquad
 \begin{aligned}
 K(t)&=\sum_i m_i[(1+u_i t)\log(1+u_i t)-u_i t],\\
 J_1(t)&=\tfrac12\sum_i\log(1+u_i t).
 \end{aligned}
\]
For example, the error in \(K(t)-t^2\) follows by summing its
geometric log series using
\(\sum_i m_i|u_i|^k\le2(2\varepsilon)^{k-2}\), \(k\ge2\).
The log bound for \(J_1\) uses \(\sum_i|u_i|\).
The Stirling remainder is bounded by the explicit first
omitted term and the sector factor, summed over
\(\sum_i1/m_i\). [DLMF Stirling remainder](https://dlmf.nist.gov/5.11.ii)

For the inverse Fourier integral
\[
                 p_j(s)=\frac1{2\pi i}\int_{i\mathbb R}
                                  e^{-hs}T_j(h)\,dh ,
\]
the global vertical-line estimate also has computable constants.
On \(h=n(\alpha+i v)\) for \(\alpha\) in any fixed bounded
interval, require \(\eta<1/8\) and
\(2\varepsilon|\alpha|\le1/4\). The real parts of the Gamma
arguments without the \(+1\) are then at least \(nm_i/2\),
while their moduli are at most \(2nm_i(1+\varepsilon |v|)\).
The second \(v\)-derivative of the real leading Stirling
exponent is therefore at most
\[
             -\frac{n}{8(1+\varepsilon^2v^2)} .
\]
This follows directly from the row expression
\(-(n e_i)^2\Re(jm_i+n e_i\alpha)/
 |jm_i+n e_i(\alpha+i v)|^2\).
Initial values and first derivatives at \(v=0\) are \(O(n)\)
with computable bounds from \(K,K'\) on a fixed \(t\)-disk.
The half-log factors have real parts bounded above in sum by
\(C\varepsilon(1+|v|)\), and the remainder has absolute value
\(O(W/n)\). On a fixed bounded \(s\)-disk this proves
\[
     |e^{-hs}T_j(h)|\le
          \exp\!\left(C n(1+|v|)
                -c n v^2/(1+\varepsilon |v|)\right).       \tag{3}
\]
For this integration of the second derivative, restrict to
\(0\le |v'|\le\min(|v|/2,1/\varepsilon)\), where the
quadratic kernel is at least \(|v|/2\) and the denominator
at most two. One can first choose a fixed large \(v_0\)
and then decrease \(\varepsilon\) so the loss in (3)
is at least any specified multiple of \(n|v|\) for
\(|v|\ge v_0\), since \(|v|/(1+\varepsilon |v|)\) increases.
This proves holomorphy and allows vertical shifts with the
same effective tail estimate.

Rouché applied to \(K'(t)=s\) on fixed small circles
about \(s/2\) gives the analytic solution \(t(s)=s/2+
O(\varepsilon)\). Set \(A=K(t(s))-s t(s)\).
Shift to the line through \(h=jt(s)\).
Its local quadratic exponent in \(h=jt(s)+inv\) is
\(-n^2K''(t(s))v^2/(2j)\), with uniformly negative real
quadratic coefficient. On any bounded part outside a fixed
small saddle interval the same loss follows by the explicit
\(C^k\) closeness of \(K\) to \(t^2\); the rest is (3).
Symmetric integration removes the odd first-amplitude and
cubic corrections. The remaining \(v^2,nv^4,n^2v^6\)
terms are integrated against the strict Gaussian loss;
the Stirling error is \(O(W/n)\).
This gives (2) with
\(a_0=\exp(J_1(t(s)))/\sqrt{2\pi K''(t(s))}\).
These arguments give computable errors by fixed Taylor
bounds and elementary Gaussian moments.

For a fixed real \(c>0\), the estimates also give
\(\mathbb E_j e^{jc|U_j|}\le 2e^{C_c j}\).
For bounded fields \(|\xi|\le F\), the inequality
\[
            |g_\varepsilon(s+\xi)|
                   \le (|s|+F)e^{\varepsilon(|s|+F)/2}
\]
therefore makes the true tails beyond \(|s|\ge R\)
exponentially small with any prescribed rate, by taking
\(c\ge2\), \(\varepsilon\le1/2\) and \(R\) sufficiently large.
The cutoff \(R\) is fixed before applying (2) on a larger
complex disk. The moment bound here is obtained from the
fixed \(t=\pm c\) in the Gamma expansion, so its constant
does not depend on that larger disk.

## Uniform regular comparison off a fixed fold neighborhood

First, for sufficiently large fixed \(Y\), there are no
single-order roots at \(|\Im\xi|\ge Y\).
On the allowed strip \(|\Im\xi|\le\pi/\varepsilon\) and
\(|\Re\xi|\le8\varepsilon\),
\[
  |g_\varepsilon(\xi)|\ge2|\Im\xi|/\pi,\qquad
  \frac{g_\varepsilon(s+\xi)}{g_\varepsilon(\xi)}
   =\cosh(\varepsilon s/2)+
                      \coth(\varepsilon\xi/2)\sinh(\varepsilon s/2).
\]
On fixed complex \(s\)-disks the last expression differs
from one by at most \(C_R(\varepsilon^2+1/Y)\).
For real \(s\) its modulus is bounded by
\(\exp((\varepsilon/2+C/Y)|s|)\).
Choose a fixed real tail cutoff using the preceding
exponential moment bound. Then choose \(Y\) large and
\(\varepsilon\) small so that the divided compact phase
\(A(s)+\log(g_\varepsilon(s+\xi)/g_\varepsilon(\xi))\)
is close with derivatives to \(-s^2/4\).
It has the one continued Gaussian maximum near zero.
Deform a fixed small interval through that saddle, while
the tails and finite endpoint connectors retain strict
loss by the uniform closeness.
The same local Gaussian comparison as in (2) gives a
nonzero contribution and excludes roots uniformly.

It remains to consider \(0\le\Im\xi\le Y\) by conjugation.
The model at \(\xi=iy\) is
\(\phi_0(s,iy)=-s^2/4+\log(s+iy)\).
Put \(y_c=2\sqrt2\).
Fix any rational distance from \(y_c\).
For \(y<y_c\) the two saddles on the line \(\Im s=-y/2\)
have real coordinates
\(\pm R_y\), \(R_y=(2-y^2/4)^{1/2}\).
On that line the real phase is
\(-x^2/4+\tfrac12\log(x^2+y^2/4)\) up to its constant,
so each half has one strict quadratic maximum.
At \(y=0\) the central base vanishes, still with strict
loss in modulus. On the compact truncated path the
strict-loss checks can use \(e^{\Re A}|g_\varepsilon|\),
which is continuous even at a vanishing base.
The Gaussian computations need the logarithm only near
each saddle, where the base stays away from zero.
For \(y>y_c\), shift to
\(\Im s=-(y-\sqrt{y^2-8})/2\). The real phase there has
the unique nondegenerate maximum at the imaginary saddle.
The finite connectors and outer tails retain strict loss
by taking the real endpoint radius large before making
\(\varepsilon\) small. These model maxima and gaps are
explicit functions of \(y\), so compact interval bounds
verify them effectively for \(y\) away from \(y_c\).
Their strict inequalities and (2) give the continued
contours and the actual Gaussian estimates.
The large-\(y\) region above \(y_c\) is thus zero-free.

On the two-saddle compact arc below \(y_c\), write
\(D_\varepsilon(\xi)\) for the right-minus-left stationary
difference. Its derivative tends to \(R_y\) on the
imaginary model arc, and the limiting imaginary phase
difference lies in \([-\pi,0)\), reaching zero only at
\(y_c\). Indeed the stationary derivatives are
\(1/(\pm R_y+iy/2)\), so \(d\Im D/dy=R_y\),
and the difference starts at \(-\pi\) and ends at zero.
The strict lower bounds on \(\Re D_\varepsilon'\) and on
\(|e^{D_\varepsilon}-1|\) are thus computable uniformly
away from \(y_c\), using \(|\Re\xi|\le8\varepsilon\).
The leading relative amplitude ratio is analytic,
nonzero, bounded above and away from zero, and independent
of order for a fixed profile by (2).
At two exact zeros the real saddle difference at each
order is \(O(1/n)\). Since the imaginary fields agree,
monotonicity on the real joining segment gives
\(|\xi_N-\xi_n|\le C/n\).
Equation (1) excludes nonzero integer \(\Delta\) for
sufficiently small \(\varepsilon\).
For \(\Delta=0\) the fields agree exactly, and dividing
the two complex cancellation equations gives
\(e^{D_\varepsilon}=1+O(1/n)\), contradicting the
computed phase gap.
 
## A uniform cubic chart with effective constants

At \(\xi_c=2\sqrt2\,i,\ s_c=-i\sqrt2\) the model has
\[
   \phi_s=\phi_{ss}=0,\qquad \phi_{sss}=i/\sqrt2,\qquad
   \phi_{s\xi}=1/2 .
\]
The limiting horizontal line through \(s_c\) has the
exact real loss \(-x^2/4+\tfrac12\log(1+x^2/2)\) at
\(s=s_c+x\). It is strictly negative for every real
\(x\ne0\).
On the chosen finite retained interval the log branch is
analytic for \(\xi\) in a fixed small disk, since \(s+\xi\)
stays in the upper half-plane at a definite distance from zero.
The same is true of \(g_\varepsilon(s+\xi)\) by its uniform
closeness to \(s+\xi\). On this disk put
\(\phi_\varepsilon(s,\xi)=A(s)+\log g_\varepsilon(s+\xi)\).
It is close with all fixed derivatives to the model,
by the effective density and elementary Taylor bounds.

The two equations
\((\partial_s\phi_\varepsilon,\partial_s^2\phi_\varepsilon)=0\)
have one actual fold near the model fold. To obtain a
constructive bound, use the inverse of the limiting two-by-two
Jacobian and a Newton fixed-point map on a fixed small
polydisk. Taylor bounds first make the derivative defect
less than \(1/4\) for the model and then less than \(1/2\)
for the perturbation. A sufficiently small \(\varepsilon\)
makes the displacement of the center less than one quarter
of the radius. This gives a computably small fold
displacement and an analytic cubic unit at the actual fold.

Here is an effective local construction of the exact cubic
normal form. At the actual fold absorb its analytic cubic
unit by taking the cube root of the phase increment divided
by the cube of the centered position. Its value at the
center is bounded away from zero by \(\phi_{sss}\ne0\);
on a small explicit disk the unit differs little from its
center value. The resulting coordinate \(v\) has derivative
close to the model \(-i/\sqrt2\). Its inverse on a fixed disk
follows by Rouché against its linear term. At the fold the
phase is \(C_*+v^3/3\).
On a fixed smaller \(v\)-disk use the norm
\(\|f\|_R=\sum_k |f_k|R^k\) on its Taylor coefficients.
For \(\xi\) in a small parameter disk, the change \(G(v,\xi)\)
of the phase from \(C_*+v^3/3\) has norm at most
\(C|\xi-\xi_{\rm fold}|\) by Cauchy on a larger disk.

To solve for the exact coordinate \(u=v+h\), constant
\(C=C_*+d\), and Airy parameter \(\lambda\), use
\[
 G=v^2h+d-\lambda v+
                 \bigl(vh^2+h^3/3-\lambda h\bigr).
\]
The linear inverse removes the first two Taylor terms
and divides the remainder by \(v^2\). With triple norm
\(R^2\|h\|_R+|d|+R|\lambda|\), this inverse has norm one.
The nonlinear part on a ball of radius \(b\) has norm at
most \(2b^2/R^3+b^3/(3R^6)\) and Lipschitz constant at
most \(4b/R^3+b^2/R^6\).
Taking \(b=R^3/100\) and \(\|G\|_R\le b/2\) gives a
contraction, analytic in the parameter. On a smaller disk
the derivative of \(h\) is small, giving an analytic inverse
coordinate. This proves the exact form
\[
               C(\xi)+u^3/3-\lambda(\xi)u .               \tag{4}
\]
At the model center, comparison of the linear response in
\(v\) gives \(\lambda'(\xi_c)=-i/\sqrt2\).
Its derivative at the actual fold is close to this value
by the fold and density bounds. Cauchy estimates on the
parameter disk give computable smaller convex disks on
which \(-\Im\lambda'\) is bounded below and \(\lambda\)
has an inverse on a fixed disk about zero. The inverse
derivative there has strictly positive imaginary part,
bounded uniformly below, by the same closeness.
These inverse assertions follow also by Rouché against
the leading linear term using the derivative bound.

At the model fold, writing \(h=s-s_c\), the increment
is \(ih^3/(6\sqrt2)-h^4/16+O(h^5)\).
Thus the cubic coordinate is
\(-ih/\sqrt2+h^2/8+O(h^3)\).
The two small real endpoints \(h=\pm r\) enter the right
Airy decay sectors about \(\pm\pi/3\).
Choose a fixed small rational \(r\) using the convergent
log series in \(h\). Connect each endpoint in its decay
sector to a common-radius point on the corresponding
central ray. These paths stay in a fixed small coordinate
disk and have strict cubic real loss. This can be seen
by a circular arc followed by a radial segment inside the
open sector, with both radii positive.
The chart, endpoints, and losses persist for the actual
phase and sufficiently small physical parameter disk.
Hence the local contour is the elementary Ai class.
Outside this local chart the truncated horizontal contour
has uniform strict loss by the explicit model loss.
The end connectors back to the real contour and the
true real tails have the same exponential control,
by choosing the fixed real cutoff sufficiently large.
Only the finite horizontal contour need be interpolated
in the complex order \(j\). Restrict its order disk so
the strict losses survive multiplication by \(j\), using
the bound for the imaginary part on the finite paths.
All these strict losses and choices are computable.
At integer orders the omitted tails and connectors
therefore give a normalized \(O(e^{-c n})\) residual,
relative to \(e^{j C(\xi)}\) up to fixed polynomial factors.

## Effective Airy residual and linked roots

Normalize the analytic finite contour integral in the fold
chart by its common exponential and its power of \(j\).
Write \(z=j^{2/3}\lambda(\xi)\), with the branch continued
from real \(j>0\). It takes the form
\[
 I(j,\xi)=P(j,\xi)\operatorname{Ai}(z)
    +j^{-1/3}Q(j,\xi)\operatorname{Ai}'(z)+E(j,\xi),\qquad
 |E|\le \frac{C}{n}
    \bigl(|\operatorname{Ai}(z)|+
           (1+|z|)^{-1/2}|\operatorname{Ai}'(z)|\bigr).
                                                               \tag{5}
\]
Here \(P\) is bounded away from zero, \(Q\) is bounded, and
they are analytic with bounded field derivatives.
All constants are computable.

To prove (5), divide the holomorphic local amplitude
\(a_j(s)\,ds/du\) by \(u^2-\lambda\).
Using a fixed Cauchy circle \(|v|=r_1\) inside the chart,
the constant and linear remainder coefficients are
the integrals of \(a(v) v/(v^2-\lambda)\) and
\(a(v)/(v^2-\lambda)\), each divided by \(2\pi i\).
The divided amplitude is the integral of
\(a(v)/[(v^2-\lambda)(v-u)]\) on that circle.
These formulas give explicit bounds on a smaller disk
when \(|\lambda|\) is small, and \(P\) is bounded near
the nonzero center amplitude. The constant and linear
terms give the displayed Airy factors.
Integration by parts against \(u^2-\lambda\) gives
\(1/j\) times an integral of the derivative of the
bounded divided amplitude, plus fixed-radius boundary
terms. Apply the
[effective elementary Airy norm](AIRY_EFFECTIVE_ELEMENTARY_JOINT_NORM.md)
after \(v=j^{1/3}u\). The order disk rotates the two
rays by a small angle, and its radius condition holds
by decreasing \(|\lambda|\) inside a fixed physical
coordinate disk. Formula (1) of that estimate gives
the joint norm in (5).
The lower bound (2) there absorbs the exponentially
small boundary and exterior terms into (5), shrinking
the fixed fold field disk if necessary.
This also bounds the actual integer-order residual by
\(O(e^{-c n})\) times the joint norm, since
\(|z|\le |\lambda|\,j^{2/3}\).

We use a quantitative local Airy consequence.
If \(|\operatorname{Ai}(z)/\operatorname{Ai}'(z)|\le q\)
and \(q^2(|z|+2q)\le1/32\), there is an Airy zero
within \(2q\) of \(z\).
In fact \(Y(v)=\operatorname{Ai}(z+v)/\operatorname{Ai}'(z)\)
obeys
\[
 Y(v)=Y(0)+v+\int_0^v(v-t)(z+t)Y(t)\,dt .
\]
On \(|v|\le2q\) this bounds \(\sup|Y|\) by \(4q\)
and its difference from the affine part by \(q/4\).
Rouché on the circle proves the claim.
All Airy zeros are negative real and simple (simplicity
also follows from the differential equation).
[DLMF Airy zero law](https://dlmf.nist.gov/9.9.i)

At a true integer-order root, (5) and the exponentially
small residual give
\[
                 |\operatorname{Ai}(z)|
                     \le C_* n^{-1/3}|\operatorname{Ai}'(z)|
\]
with computable \(C_*\), after absorbing the term
\(O(1/n)|\operatorname{Ai}(z)|\).
If \(\operatorname{Ai}'(z)=0\), this is impossible.
Choose a small bound on \(|\lambda|\) after fixing
\(C_*\), so the local Airy test applies uniformly.
Thus \(z=a_k+O(n^{-1/3})\) for a negative real Airy
zero \(a_k\), with effective constant, and
\(|\Im\lambda|\le C/n\).
On a convex fold field disk \(-\Im\lambda'\) is
bounded below. Equation (1) therefore bounds the
horizontal difference of two common-root fields by
\(C/n\), which excludes all nonzero integer \(\Delta\)
after decreasing \(\varepsilon\).

For \(\Delta=0\) the fields agree exactly.
Then the two scaled Airy coordinates differ by
\(O(|\lambda|n^{-1/3})\), and their Airy labels
must agree. For details, at \(a=a_k\) normalize
\(Y(v)=\operatorname{Ai}(a+v)/\operatorname{Ai}'(a)\).
Its equation and data \(Y(0)=0,Y'(0)=1\) give
\[
 Y(v)=v[1+O((1+|a|)|v|^2)],\qquad
 Y'(v)=1+O((1+|a|)|v|^2)
\]
whenever the product in the error is at most a sufficiently
small absolute constant. This follows from the same integral
equation by bounding the norm and integrating once more
for the derivative. It excludes any other Airy zero in
a disk of radius a fixed multiple of \((1+|a|)^{-1/2}\).
Here \(|a|n^{-2/3}\) is at most the chosen small
bound on \(|\lambda|\) plus \(O(1/n)\).
Decreasing that bound makes the excluded disk larger
than the two localization errors and their coordinate
difference. This proves common labeling effectively.

Fix that label and make the analytic substitution
\(\xi=\lambda^{-1}(j^{-2/3}z)\) on the complex order
disk and \(|z-a_k|\le R_1 n^{-1/3}\), with a large
fixed \(R_1\). For small enough \(|\lambda|\) this
stays in the previously constructed chart even at
complex orders. Dividing (5) by
\(P(j,\xi)\operatorname{Ai}'(a_k)\) yields
\[
 Y(z-a_k)+j^{-1/3}(Q/P)Y'(z-a_k)+O(1/n)
\]
uniformly there, by the ODE bounds.
Choose \(R_1\) before the final decrease of the
\(\lambda\)-disk, so the middle term is less than
half the lower boundary bound for \(Y\).
Rouché gives a unique simple zero
\(z=a_k+\chi(j)\), analytic in \(j\), with
\(|\chi|\le R_1 n^{-1/3}\).
Cauchy's estimate on a smaller order disk gives
\(|\chi'(j)|\le C n^{-4/3}\) for real \(n\le j\le N\).

On a smaller fixed multiple of this \(z\)-disk the
derivative of the normalized left side is within \(1/4\)
of one, by further shrinking the bound on \(|a_k|n^{-2/3}\)
and increasing \(n\).
Indeed the derivative of \(Y\) is controlled by its
integral equation. In the derivative of the middle term,
\(\partial_z(Q/P)=O(n^{-2/3})\) and
\[
 n^{-1/3}|Y''(z-a_k)|
      =O((1+|a_k|)n^{-2/3})
\]
by the Airy equation and the bound for \(Y\).
Cauchy on a slightly larger disk bounds the last
derivative by \(O(n^{-2/3})\).
One can keep both the actual and the interpolating zero
in this smaller disk by performing the same Rouché
test on an intermediate radius and choosing \(R_1\)
larger at the outset.
The integer-order error after normalization is
\(O(e^{-c n})\), by the joint lower bound (2) of the
Airy estimate at \(a_k\), decreasing the fold disk if
needed. The derivative bound and integration along the
straight segment to the analytic zero show that the
actual integer-order root differs from it by
\(O(e^{-c n})\) in the \(z\)-coordinate.

Finally, the branch in the physical scaled field is
\(\xi(j)=\lambda^{-1}(j^{-2/3}(a_k+\chi(j)))\), with
\[
     \xi'(j)=(\lambda^{-1})'
                    \bigl(-\tfrac23 a_k j^{-5/3}+O(n^{-2})\bigr).
\]
The inverse derivative has a uniformly positive imaginary
part. The negative Airy zeros have a computable positive
lower bound on their absolute values, for example from
\(\operatorname{Ai}(0)>0\) and its local series.
Hence \(\Im \xi'(j)\) has a strict fixed sign and modulus
bounded below by \(c|a_k|n^{-5/3}\) for all sufficiently
large \(n\). The actual roots' exponential errors cannot
remove this difference between the two imaginary
coordinates. This contradicts (1) and excludes the fold.

## Choosing the thresholds

First fix the large-field bound \(Y\), the finite model
fold contour, and its analytic chart.
Choose the loss, inversion and Airy constants on those
fixed disks, then shrink the admitted fold field disk
to meet the Airy product and error conditions.
Choose a covered fixed interval about \(y_c\) so that
\(|\Re\xi|\le8\varepsilon\) puts all such fields in that
fold disk once \(\varepsilon\) is sufficiently small.
The compact regular comparison then applies on the
remaining subintervals of \([0,Y]\).
One may use larger density disks or real tail cutoffs
for these finite contours. Their constants were obtained
on fixed model domains, before the final decreases in
\(\varepsilon\) and increases in \(n\).
Finally choose a rational upper bound on
\(W=2\varepsilon^2\) and an integer terminal satisfying
all the strict inequalities just given.
All row estimates depend on the profile only through
the displayed total bounds. This produces the claimed
uniform computable constants. Rational Bézout gives
coprimality over every characteristic-zero field.
