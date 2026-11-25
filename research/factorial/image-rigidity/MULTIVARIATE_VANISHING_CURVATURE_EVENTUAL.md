# Uniform eventual coprimality at vanishing shared curvature

Let \(a,b\) be distinct ordinary exponent columns of arbitrary
finite dimension. Delete their equal rows, put \(d_i=b_i-a_i\),
\(m_i=(a_i+b_i)/2\), and write
\[
         W=\sum_i d_i^2/m_i,\qquad \Delta=\sum_i d_i .
\]
**Theorem.** There exist absolute constants \(W_0>0,N_0\)
such that every consecutive pair of ordinary factorial moment
polynomials at terminal \(N\ge N_0\) with \(W\le W_0\) is
coprime, over every characteristic-zero field. The total
degree difference is arbitrary. The constants are unevaluated.

It suffices to exclude a sequence of counterexamples with
\(N\to\infty\) and \(W\to0\). Such a sequence would exist if
no pair of constants worked. In this regime every unequal
row is eventually shared, since a private row of degree at
least one contributes at least two to \(W\). Write \(n=N-1\)
and \(\varepsilon=\sqrt{W/2}\). The profile and \(\varepsilon\)
are fixed during each comparison of the two orders.

## Midpoint density on a complex order disk

For each real positive order parameter \(j\), let \(X_i\)
have independent Gamma distributions of shapes \(j m_i+1\)
and unit rate. Use the exact order-dependent center
\[
 \mu_j^0=\sum_i d_i\log(jm_i),\qquad
 U_j=\frac{\sum_i d_i\log X_i-\mu_j^0}{\varepsilon}.
\]
For the original scalar coefficient \(t\), put
\[
 g_\varepsilon(z)=2\sinh(\varepsilon z/2)/\varepsilon,\qquad
 t=-\exp(-\mu_j^0+\varepsilon\xi_j).
\]
After removal of a nonvanishing factor the ordinary moment
at integer order \(j\) is
\[
                     {\mathbb E}_j[g_\varepsilon(U_j+\xi_j)^j].
                                                               \tag{1}
\]
Indeed the weight \((X^a X^b)^{j/2}\) is the midpoint
Gamma tilt, and \(X^a+tX^b\) under this factor gives a
hyperbolic sine times a nonzero scalar. The same logarithm
of \(-t\) is used at both orders, so
\[
                 \xi_N-\xi_n=\Delta\log(N/n)/\varepsilon
                 \quad\hbox{is real}.                    \tag{2}
\]
The [root annulus](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md#shared-rows-and-large-degree-imbalance)
gives, at any root and for small \(W\),
\[
       |\Re \xi_j|=O(\varepsilon),\qquad
       |\Im \xi_j|\le\pi/\varepsilon .                     \tag{3}
\]
For the first estimate, \(V=\sum d_i^2/\min(a_i,b_i)=O(W)\)
and the log center \(\log G_j\) of the root annulus
differs from \(\mu_j^0\) by \(O(W)\).
To see the center comparison, subtract \(\Delta\log j\)
from each. The log center from the factorial endpoints
equals \(\sum_i(F(b_i)-F(a_i))+O(V/j)\), with
\(F(x)=x\log x-x\), by the integral or Stirling comparison
in the annulus proof. The difference between this leading
part and \(\sum_i d_i\log m_i\) is \(O(W)\), by integrating
\(F'(x)-\log m_i\) over the row interval. The second estimate
in (3) follows by choosing \(\arg(-t)\in[-\pi,\pi]\).

Here is an analytic density estimate uniform in the number of
coordinates. In this section \(s\) denotes the scalar density
coordinate. Put \(e_i=d_i/\varepsilon\) and \(u_i=e_i/m_i\).
The integer gaps give
\[
 \max_i |u_i|\le 2\varepsilon,\qquad
 \sum_i |u_i|\le 2\varepsilon,\qquad
 \sum_i m_i u_i^2=2,\qquad \sum_i 1/m_i\le W.              \tag{4}
\]
Choose a fixed sufficiently small \(\eta>0\) and consider
complex \(j\) with \(|j-n|\le\eta n\). Set
\[
 H_j(h)=\prod_i
   \frac{\Gamma(jm_i+1+e_i h)}{\Gamma(jm_i+1)}
             (jm_i)^{-e_i h}.
\]
For real \(j\), this is the moment generating function of
\(U_j\). Choose a branch real at positive \(j\) for
\(\log(jm_i)\). As an analytic function of \(j,s\) on fixed
bounded \(s\)-disks, define
\[
       p_j(s)=\frac1{2\pi i}\int_{i\mathbb R} e^{-h s}H_j(h)\,dh.
                                                               \tag{5}
\]
For real \(j\), this is the density by Fourier inversion.
We give the required integrability and saddle estimates.

For bounded \(t\), sectorial Stirling applied at \(h=jt\) gives
\[
 \log H_j(jt)=j K(t)+J_1(t)+O(W/n),\qquad
 \begin{aligned}
 K(t)&=\sum_i m_i[(1+u_i t)\log(1+u_i t)-u_i t],\\
 J_1(t)&=\tfrac12\sum_i\log(1+u_i t).
 \end{aligned}                                             \tag{6}
\]
The branch is chosen to give value zero at \(t=0\).
Here \(K(t)\to t^2\) and \(J_1(t)\to0\) uniformly on every
fixed compact set, including derivatives. This follows by
Taylor expansion and (4). The Stirling bound follows by
summing \(O(1/(nm_i))\). All arguments remain in fixed
right-half-plane sectors on bounded \(t\)-sets.

For a global estimate on the line \(h=n(\alpha+i v)\),
with real bounded \(\alpha\), every Gamma argument has real
part at least \(nm_i/2\) for sufficiently small
\(\varepsilon,\eta\). Uniform Stirling in the right half-plane
gives the leading logarithm in row \(i\) as
\[
 (jm_i+e_i h)\bigl(\log(jm_i+e_i h)-\log(jm_i)\bigr)-e_i h.
\]
The second \(v\) derivative of its real part is
\[
       -(n e_i)^2
        \frac{\Re(jm_i+n e_i\alpha)}
             {|jm_i+n e_i(\alpha+i v)|^2}.
\]
Summing and using (4) bounds it above by
\(-c n/(1+\varepsilon^2v^2)\) with \(c>0\) absolute
for bounded \(\alpha\).
The leading sum and its first \(v\) derivative at zero are
\(O(n)\) by (6) on bounded disks. The half-log terms from
Stirling have real part at most \(C\varepsilon(1+|v|)\)
by (4) and \(\log(1+x)\le x\). The summed Stirling remainder
is \(O(W/n)\). Integrating the second derivative twice
therefore gives
\[
 |e^{-h s}H_j(h)|
 \le \exp\left(C n(1+|v|)
             -\frac{c_1 n v^2}{1+\varepsilon|v|}\right)
                         \qquad (|s|\le C_s).             \tag{7}
\]
For the double integral of \((1+\varepsilon^2 v^2)^{-1}\),
use the subinterval of length
\(\min(|v|/2,1/\varepsilon)\). Constants in (7) may depend
on the fixed \(\alpha,s\) bounds. Since \(v/(1+\varepsilon v)\)
increases for \(v>0\), (7) gives integrable exponential
decay uniformly beyond a sufficiently large fixed \(v\),
provided the sequence is taken sufficiently far.
This proves holomorphy of (5) in both variables and permits
shifts between the indicated vertical lines. The same bounds
on horizontal connectors justify those shifts.

On a fixed sufficiently large \(s\)-disk the equation
\(K'(t(s))=s\) has a holomorphic solution with
\(t(s)\to s/2\), uniformly on the disk, by (6) and Rouché
on small circles about \(s/2\). Put
\[
 A(s)=K(t(s))-s t(s),\qquad
 a_0(s)=\frac{\exp(J_1(t(s)))}{\sqrt{2\pi K''(t(s))}} .
\]
Then \(A(s)\to-s^2/4\) and \(a_0(s)\to1/\sqrt{4\pi}\),
including derivatives on smaller disks.
Shift (5) to the vertical line through \(h_0=j t(s)\).
On bounded \(v\) in \(h=h_0+i n v\), the leading exponent
minus \(j A(s)\) is stationary at \(v=0\) and has quadratic
term \(-n^2 K''(t(s))v^2/(2j)\). Its real part is strictly
negative quadratically on a fixed small neighborhood,
uniformly on the order disk. On the remaining bounded
parts the same strict loss follows from \(K\to t^2\).
The large-\(|v|\) part is exponentially small relative to
the stationary value by (7), after increasing the fixed
threshold. Taylor expansion of the local Gaussian integral
now gives
\[
      p_j(s)=\sqrt j\, e^{j A(s)} a_j(s),\qquad
      a_j(s)=a_0(s)+O(n^{-1}).                             \tag{8}
\]
Indeed the local cubic and linear amplitude corrections at
relative order \(n^{-1/2}\) are odd in the centered \(v\)
coordinate. The bounded higher Taylor derivatives from (6)
give an \(O(n^{-1})\) integrated correction against the
Gaussian majorant. The analytic square root is positive
on the real order ray. The estimate is uniform for complex
\(j\) in the disk and bounded complex \(s\), and \(a_j\) is
analytic there by (5). In particular \(a_j\) and its
\(s\)-derivatives are uniformly bounded, and \(a_j\) stays
away from zero on each fixed compact \(s\)-set.

For real orders, (6) also gives
\({\mathbb E}_j e^{j c |U_j|}\le e^{C_c j}\) at every fixed
\(c>0\), by using \(H_j(j c)+H_j(-j c)\).
Consequently in (1), for bounded \(\xi_j\), the tails
\(|s|>R\) on the real line can be made exponentially small
with any prescribed constant in the exponent by choosing
\(R\) sufficiently large and fixed. Indeed
\[
 |g_\varepsilon(s+\xi_j)|
    \le (|s|+C)\exp(\varepsilon(|s|+C)/2)
\]
and multiplication by \(e^{-c|s|}\) bounds the tail
by the stated exponential moment. This establishes a
finite exact contour reduction at integer orders, with
arbitrarily small absolute exponential tails. Bounded
complex deformations of the retained interval use (8).

## Root localization and the regular arc

By conjugation take \(\Im\xi_j\ge0\); the even-order moment
is positive at every real physical parameter, so a common
root is nonreal. We first exclude individual roots with
\(|\Im\xi_j|\to\infty\). By (3),
\(|g_\varepsilon(\xi_j)|\ge 2|\Im\xi_j|/\pi\) and
\(|\cosh(\varepsilon\xi_j/2)|\) stays bounded.
The quotient \(g_\varepsilon(s+\xi_j)/g_\varepsilon(\xi_j)\)
tends to one uniformly on bounded complex \(s\)-sets, by
the hyperbolic addition formula. For real \(s\) its modulus
is at most
\(\exp((\varepsilon/2+C/|\Im\xi_j|)|s|)\), by the same
formula and \(|\sinh x|\le |x|e^{|x|}\).
The fixed-rate exponential moment bound thus controls
the real tail after division by \(g_\varepsilon(\xi_j)^j\).
On a bounded real interval the divided phase
\[
             A(s)+\log(g_\varepsilon(s+\xi_j)/g_\varepsilon(\xi_j))
\]
converges analytically to \(-s^2/4\). A small local
deformation through its simple saddle and the strict outer
loss give a nonzero Gaussian packet. The deformation stays
on the fixed bounded \(s\)-disk, where (8) holds. This
excludes an escaping zero.

At bounded fields the limiting phase is
\[
                   -s^2/4+\log(s+\xi),\qquad \xi=i y,
                    \quad y\ge0,                          \tag{9}
\]
by (3) and (8); local logarithms are chosen near the
saddles. On the retained interval the integrand
\(e^{jA(s)}g_\varepsilon(s+\xi)^j a_j(s)\) at integer
orders is analytic even at zeros of the last base.
If \(y>2\sqrt2\), set
\(Y=(y-\sqrt{y^2-8})/2\) and shift the interval to
\(\Im s=-Y\). Its limiting horizontal real phase is a
constant plus
\[
             -x^2/4+\tfrac12\log(x^2+(y-Y)^2).
\]
Since \(y-Y>\sqrt2\), this has a single nondegenerate
maximum at zero, which is also the complex saddle.
One can make the endpoint connectors exponentially
small by choosing the real endpoints sufficiently far
out before taking the sequence limit, using \(A\to-s^2/4\)
on the resulting finite rectangle. The outer real tails
are already exponentially small. The same Gaussian
comparison and contour shift apply to the nearby
phases, so no root limits to this field.

For \(0\le y<2\sqrt2\), put \(R_y=\sqrt{2-y^2/4}\).
On the line \(\Im s=-y/2\) the two simple saddles of (9)
are \(s_\pm=\pm R_y-iy/2\). The real horizontal phase
is a constant plus
\(-x^2/4+\tfrac12\log(x^2+y^2/4)\) and is uniquely
maximal on each half at its saddle, with negative
second derivative there. The center is strictly below
both maxima. At \(y=0\) the limiting base vanishes at
the center, which still gives strict loss in modulus.
For nearby actual parameters, split the contour at
the center and deform each half through its respective
continued saddle. Only a bounded neighborhood is
needed, and the midpoint and endpoint connectors retain
strict phase gaps. The real outer tails and both connectors
are exponentially small relative to either packet.
Thus the two actual local integrals have Gaussian
expansions with nonzero leading amplitudes, uniformly
on a small parameter disk about \(iy\). Their phase
difference \(D_\varepsilon(\xi)\), taken as right minus
left, is independent of order at a given profile.
The leading amplitudes likewise are independent of order
by (8), and the leading relative error is \(O(1/n)\).
This follows by the local Gaussian expansion on the two
deformed contours and the uniform fixed-neighborhood
derivative bounds in (8). The logarithm near the negative
saddle at \(y=0\) is taken with limiting imaginary part
\(\pi\), compatible with its integer-order factor.

The limiting \(D\) obeys \(D_\xi=R_y>0\) on the arc:
the stationary-value derivatives are
\(1/(\pm R_y+iy/2)\). Along the arc \(D\) goes from
\(-i\pi\) to zero, since \(d\Im D/dy=R_y\) and
the two branches coincide at the fold. These phase
branches match along the arc by continuation.
At a proposed common root both \(\xi_n,\xi_N\)
therefore limit to the same \(iy\) in the closed
interval, by their common imaginary part and (3).

Suppose the limit is regular, \(y<2\sqrt2\).
Modulus cancellation of the two packets requires
\(\Re D_\varepsilon(\xi_j)=O(1/n)\) at either order.
Since \(\Re D_\varepsilon'\) is uniformly positive near
this limit, (2) implies \(|\xi_N-\xi_n|=O(1/n)\).
A nonzero integer \(\Delta\) would make this shift at
least \(\log(N/n)/\varepsilon\), a contradiction.
Thus \(\Delta=0\) and the fields are identical.
The ratios of the leading Gaussian amplitudes at the
two orders agree up to \(1+o(1)\) by (8).
Dividing the two cancellation equations would give
\(e^{D_\varepsilon(\xi_n)}\to1\), contrary to
\(-\pi\le\Im D<0\). This covers the regular arc including
its real negative boundary.

## The Hermite fold and its adjacent zero motion

The remaining limit is \(\xi_c=2\sqrt2\,i\), with
\(s_c=-i\sqrt2\). The limiting phase (9) has exact
horizontal real loss
\(-x^2/4+\tfrac12\log(1+x^2/2)\) on that line.
It is strictly negative away from zero. Thus the finite
horizontal interval at \(\Im s=-\sqrt2\), with end
connectors to the real interval, has a fixed strict
loss outside any small central neighborhood of the fold.
Its initial shift and connectors are justified just as
above. The loss persists at the nearby phases and for
complex orders in a fixed smaller order disk, after
factoring out the centered exponential. On the retained
horizontal interval \(s+\xi\) stays away from zero on
a fixed parameter disk about \(\xi_c\), and \(g_\varepsilon\)
has no zero there for sufficiently small \(\varepsilon\).
Thus \(\log g_\varepsilon(s+\xi)\) admits one continuous
analytic branch on this interval and its neighborhood.
Together with (8), this gives an analytic interpolation
of the exact retained horizontal integral for complex
\(j\). At integer orders the difference between that
integral and (1), divided by its common Gaussian density
factor, is exponentially small relative to the fold
center. It consists of the real tails and the end
connectors, with the uniform losses just verified.

At the limiting fold set \(h=s-s_c\).
The increment of (9) is
\(i h^3/(6\sqrt2)-h^4/16+O(h^5)\).
The analytic cubic coordinate can therefore be chosen as
\(u=-ih/\sqrt2+h^2/8+O(h^3)\).
Both small real endpoints \(h=\pm r\) go into the two
right Airy decay sectors about the rays of arguments
\(\pm\pi/3\). The common analytic cubic normal form
\[
                 C_\varepsilon(\xi)+u^3/3-\lambda_\varepsilon(\xi)u
                                                               \tag{10}
\]
for the actual phases, with uniform coordinate disks,
follows from the
[two-stationary-point cubic construction](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#localizing-the-exact-contour-at-a-fold).
Its hypotheses hold at the limiting fold
\(\Phi_s=\Phi_{ss}=0\), \(\Phi_{sss}=i/\sqrt2\ne0\), and
\(\Phi_{s\xi}=1/2\). The charts and derivatives converge.
In particular \(\lambda_0'(\xi_c)=-i/\sqrt2\), by
matching the linear \(h\) coefficient. The two sector
endpoints persist, hence the class is the elementary
Airy contour. Along the discarded finite pieces the
strict loss established above persists.

We give the consecutive comparison using the
[varying-index Airy argument](MULTIVARIATE_SHARED_SUBCRITICAL_THERMAL_FOLD.md#a-fixed-cubic-chart-for-the-thermal-contour).
With \(z=j^{2/3}\lambda_\varepsilon(\xi)\) and after
removal of the common exponential and saddle prefactor,
the analytic interpolation satisfies
\[
 \begin{split}
 I(j,\xi)&=P(j,\xi)\mathrm{Ai}(z)
        +j^{-1/3}Q(j,\xi)\mathrm{Ai}'(z)+E(j,\xi),\\
 |E(j,\xi)|&\le \frac{C}{|j|}
     \left(|\mathrm{Ai}(z)|+
       (1+|z|)^{-1/2}|\mathrm{Ai}'(z)|\right).             \tag{11}
 \end{split}
\]
Here \(P\) is bounded away from zero, \(Q\) is bounded,
and both are analytic with bounded \(\xi\)-derivatives.
Indeed divide the analytic amplitude \(a_j(s)\,ds/du\)
by \(u^2-\lambda_\varepsilon(\xi)\). The constant and
linear terms give the displayed Airy factors. The
remainder has a factor \(j^{-1}\) by integration by
parts. The
[joint Airy contour estimate](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#a-uniform-joint-norm-for-the-actual-cubic-class)
bounds the remaining integral in the displayed norm,
uniformly for the complex order disk. Its finite connectors
and exterior pieces are exponentially small in that norm
on a smaller fold disk, by their strict loss. Rotation of
the scaled cubic coordinate for complex \(j\) preserves
the two open decay sectors on a sufficiently small fixed
order disk. The integer-order tails and end connectors
omitted from the interpolation have the same exponential
control, since the fold disk can be reduced until the
Airy joint norm costs less than their fixed phase loss.

At an actual root, (11) yields
\(|\mathrm{Ai}(z)|\le C n^{-1/3}|\mathrm{Ai}'(z)|\).
The [Airy continuation argument](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md#the-elementary-airy-class-and-its-transverse-band)
then puts \(z=a_k+O(n^{-1/3})\) for a negative real
Airy zero \(a_k\). For its product guard,
\(n^{-2/3}|z|=O(|\lambda_\varepsilon(\xi)|)\to0\)
at the fold. In particular
\[
        |\Im\lambda_\varepsilon(\xi)|=O(n^{-1}),\qquad
        |a_k| n^{-2/3}\longrightarrow0 .                 \tag{12}
\]
Because \(\Im\lambda_\varepsilon'\) stays bounded away
from zero on the segment between the two fields, (2) and
(12) give \(|\xi_N-\xi_n|=O(n^{-1})\). This excludes every
nonzero integer \(\Delta\), again by \(\varepsilon\to0\).

For \(\Delta=0\) the fields agree exactly. The scaled
coordinates \(n^{2/3}\lambda_\varepsilon(\xi_n)\) and
\(N^{2/3}\lambda_\varepsilon(\xi_n)\) differ by
\(O(n^{-1/3}|\lambda_\varepsilon(\xi_n)|)\), and hence
the selected Airy zeros have the same index. To see the
required separation, the function
\(Y(v)=\mathrm{Ai}(a_k+v)/\mathrm{Ai}'(a_k)\) satisfies
\(Y''=(a_k+v)Y,\ Y(0)=0,\ Y'(0)=1\). Its integral equation
gives \(Y(v)=v[1+O((1+|a_k|)|v|^2)]\) on a small
fixed multiple of the disk of radius
\((1+|a_k|)^{-1/2}\). No other Airy zero lies there.
By (12), this radius is much larger than \(n^{-1/3}\).

Now fix that index \(k\) and make the analytic substitution
\(\xi=\lambda_\varepsilon^{-1}(j^{-2/3}z)\)
on the complex order disk \(|j-n|\le\eta n\) and
\(|z-a_k|\le R_1 n^{-1/3}\), for large fixed \(R_1\).
This remains in the fold disk by (12).
The same Airy differential equation gives on this disk
\(Y(v)=v(1+o(1))\) and \(Y'(v)=1+o(1)\).
After division by \(P(j,\xi)\mathrm{Ai}'(a_k)\), the
second term of (11) is \(O(n^{-1/3})\) and its error is
\(O(n^{-1})\). Rouché therefore gives a unique simple
interpolating root \(z=a_k+\chi_k(j)\), where
\(\chi_k(j)=O(n^{-1/3})\) is analytic in \(j\).
Cauchy's estimate, on a smaller order disk, gives
\(\chi_k'(j)=O(n^{-4/3})\).
The derivative of the normalized interpolation in \(z\)
is \(1+o(1)\) on a smaller fixed multiple of the
\(n^{-1/3}\) disk. In fact \(Y'(v)=1+o(1)\), and the
derivative of the \(j^{-1/3}(Q/P)Y'(v)\) term is \(o(1)\)
by the analytic bounds on \(Q/P\) and
\(Y''=(a_k+v)Y\), using (12).
Cauchy's estimate for the normalized \(O(n^{-1})\)
error gives derivative \(O(n^{-2/3})\).
It follows that each actual integer-order root
differs from this interpolating root in \(z\) by
\(O(e^{-c n})\). The normalization by
\(\mathrm{Ai}'(a_k)\) preserves the exponential loss
by the joint bound used in (11), since
\(|a_k|n^{-2/3}\to0\).

In the physical scaled field the interpolating root is
\(\xi_k(j)=\lambda_\varepsilon^{-1}(j^{-2/3}(a_k+\chi_k(j)))\).
Differentiating gives
\[
 \xi_k'(j)=(\lambda_\varepsilon^{-1})'
               \left(-\tfrac23 a_k j^{-5/3}+O(n^{-2})\right)
                  \qquad (n\le j\le N).
\]
The inverse derivative tends to \(i\sqrt2\).
Thus its imaginary part has fixed sign and magnitude
at least \(c |a_k| n^{-5/3}\) throughout the order interval.
The negative Airy zeros are bounded away from zero.
The exponentially small errors cannot remove this
difference between the physical imaginary coordinates.
It contradicts (2) and completes the fold exclusion.

All original moment coefficients are integers.
Complex coprimality gives a rational Bézout identity
and hence the same conclusion over every
characteristic-zero field.
