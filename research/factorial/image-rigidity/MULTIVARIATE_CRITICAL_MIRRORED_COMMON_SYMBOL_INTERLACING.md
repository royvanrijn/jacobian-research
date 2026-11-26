# Eventual interlacing in the mirrored moving critical window

For a fully shared mirrored profile omit equal rows and write
\(W=\sum_i d_i^2/(\ell_i+d_i/2)\) over both orientations, where
\(d_i\) is a row's positive integer gap and \(\ell_i\ge1\) is
its smaller exponent. Let \(L_j\) be the endpoint-normalized
log-factorial symbol per degree at order \(j\).

**Theorem.** For every fixed \(K\ge0\) there is \(n_0(K)\) such
that if \(n\ge n_0(K)\), \(N=n+1\), and
\[
                         |W-4|\le K/\sqrt n,
\]
then the polynomials
\[
 \sum_{k=0}^n\binom nk
       e^{n[(1-\lambda)L_n(k/n)+\lambda L_N(k/n)]}z^k,\qquad
 \sum_{k=0}^N\binom Nk e^{N L_N(k/N)}z^k
\]
have simple unit-circle roots strictly interlacing in argument
on the circle cut at \(+1\), uniformly for \(0\le\lambda\le1\).
The smaller exponents may be real. The common value \(W=4\)
is included, as are the actual pair and the frozen upper symbol.
The terminal threshold is unevaluated.

## Uniform phase and circle localization

Put \(x=q-1/2\), \(G=L_\infty\), and
\[
 H(x)=G(1/2+x)-(1/2+x)\log(1/2+x)
                         -(1/2-x)\log(1/2-x).
\]
Write \(m_i=\ell_i+d_i/2\), \(r_i=d_i/(2m_i)\), and
\(W_i=d_i^2/m_i\). If \(W\le5\), then \(r_i\le1/2\).
For gaps one and two this follows from \(\ell_i\ge1\).
A larger gap with \(r_i>1/2\) would give a pair contribution
\(4d_i r_i>6\). The phase is symmetric and has the expansion
\[
 H(x)=H(0)+\tfrac12(W-4)x^2-\kappa x^4
                          +h_6 x^6+h_8 x^8+O(x^{10}),
 \qquad
 \kappa=\frac{32-8\sum_i W_i r_i^2}{24}\in[11/12,4/3].
\]
The Taylor bounds on a fixed disk about zero are uniform.
Indeed \(G''(q)\) is the sum of the row terms
\(W_i/(1\pm r_i(2q-1))\), so all its derivatives on smaller
disks are bounded using \(\sum W_i\le5\).

On a fixed small real interval the expansion gives
\[
                 H(x)-H(0)\le C|W-4|x^2-c x^4.
\]
Outside any such interval there is a uniform strict loss from
\(H(0)\) for large \(n\). In fact the second derivative of the
row contribution to \(\rho(q)=q(1-q)G''(q)\) is
\[
 -\frac{2d_i^2 a_i b_i}{((1-q)a_i+qb_i)^3}\le -4W_i/9 .
\]
This follows from \(|r_i|\le1/2\) and
\(2(1-r_i^2)/(1+|r_i|)^3\ge4/9\). Reflection gives
\(\rho'(1/2)=0\), and
\(H''(x)=(\rho(1/2+x)-1)/(1/4-x^2)\).
At each fixed sufficiently small positive \(x\), \(H'\) and
\(H-H(0)\) are strictly negative as \(n\to\infty\), by the
uniform quartic expansion. The bound on \(\rho''\) makes
\(H''<0\) at every larger \(x<1/2\), for large \(n\).
Reflection gives the other endpoint.

Actual mirrored polynomials have simple roots on the circle
by the [Gegenbauer factorization](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md).
For the interpolated symbol the
[uniform same-grid freezing argument](MULTIVARIATE_SUBCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md#uniform-freezing-of-the-lower-grid)
also gives simple unit-circle roots for large \(n\). That
estimate uses a bounded upper value of \(W\) and a positive
lower bound on \(L=\min_i\ell_i/d_i\), as in the
[supercritical application](MULTIVARIATE_SUPERCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md#away-from-the-positive-axis).
Here \(L\ge1/2\) follows from \(r_i\le1/2\). Its error bounds
are uniform in \(\lambda\), since the symmetric jump measure
is just multiplied by \(\lambda\).

## Bounded critical fields

Use \(m=n,N\), \(t_m=(4-W)\sqrt m\), and \(\theta=y m^{-3/4}\).
For the interpolated symbol at real or complex order \(m\) near
\(n\), write \(L_{m,\lambda}=(1-\lambda)L_m+\lambda L_N\),
where \(L_m\) is defined by the same analytic Gamma formula.
The upper order is exactly \(L_N\) for every \(\lambda\).
On a fixed disk about \(x=0\), analytic Stirling gives even
functions \(A,B,A_0,A_1\) with uniform analytic bounds, with
\(A\) nonzero and \(A(0)>0\), such that
\[
 \begin{split}
 L_m(q)&=G(q)+A_0(q)/m+A_1(q)/m^2+O(|m|^{-3}),\\
 \frac{\Gamma(m+1)e^{mL_m(q)}}
      {\Gamma(mq+1)\Gamma(m(1-q)+1)}
   &=m^{-1/2}e^{mH(x)}A(x)(1+B(x)/m+O(|m|^{-2})).
 \end{split}
\]
Evenness in the first line means evenness as a function of
\(x\). These expansions hold also for \(|m-n|\le\eta n\)
with a fixed small \(\eta>0\). Each row's sector margin
follows from \(r_i\le1/2\), and the sums in the remainder are
bounded by \(\sum_i 1/m_i\le W\). The endpoint chords in \(A_0\)
and \(A_1\) are uniformly bounded on the disk by integrating
their second derivatives. The binomial expansion only uses
arguments bounded away from zero on this disk. On the full
real interval the factorial multiplier contributes
\(m(L_m-G)=O(1)\); the binomial coefficient is bounded by
its exponential entropy estimate. The same bound applies
to the interpolated factorial multiplier.

At \(m=n\), the lower interpolation adds the factor
\[
 \exp[n\lambda(L_N-L_n)(q)]
                       =\exp[-\lambda A_0(q)/N+O(n^{-2})].
\]
After division by its positive value at \(q=1/2\), its
difference from \(1\) is \(O(x^2/n+n^{-2})\) on the disk.
Thus it has no contribution through order \(n^{-1}\) in the
critical expansion below.

Center the polynomial value at \(e^{i\theta}\) by
\(e^{-im\theta/2}\), and divide by the positive factor
\(m^{1/4} e^{mH(0)} A(0)\) and by the interpolation factor
at the center. The result is
\[
 Q_\kappa(t_m,y)+m^{-1/2}Q_1(t_m,y)+m^{-1}Q_2(t_m,y)
                              +O(m^{-3/2}),                 \tag{1}
\]
where
\[
 Q_\kappa(t,y)=\int_{\mathbb R}e^{-\kappa v^4-tv^2/2+iyv}\,dv .
\]
The functions \(Q_1,Q_2\) for a fixed profile are independent
of \(m,\lambda\). The estimates are uniform for bounded real
\(t_m\) and bounded complex \(y\). To get (1), scale \(v=m^{1/4}x\)
in the exact lattice sum. Its spacing is \(m^{-3/4}\).
On \(|v|\le m^{1/64}\) expand the sextic and octic phase terms
and the even analytic amplitude through order \(1/m\).
The remaining integrand is bounded by \(m^{-3/2}\) times
a fixed polynomial in \(|v|\) times
\(\exp(-c v^4+C v^2+C|v|)\).
The real phase losses bound the original discarded tails
by \(O(\exp(-c' m^{1/16}))\). The expanded integrands have
the same quartic envelope there. For each expanded integrand
the full lattice Riemann sum differs from the integral by
\(O(m^{-9/4})\): Poisson summation at mesh \(h=m^{-3/4}\)
bounds the nonzero frequencies by
\(C h^3\|g'''\|_1\), independently of the lattice shift.
All these derivative norms are uniformly bounded on the
stated compact parameter sets. This proves (1), including
the lattice parity change.

The positive zeros \(\beta_j(\kappa,t)\) of \(Q_\kappa(t,\cdot)\)
are simple and infinite in number. Their exact inequality
\[
                      3\beta_j-2t\,\partial_t\beta_j>0       \tag{2}
\]
follows by rescaling \(u=(4\kappa)^{1/4}v\) in the
[quartic root-speed theorem](MULTIVARIATE_DIFFUSE_MOVING_CRITICAL_WINDOW.md#a-quartic-root-speed-inequality).
The rescaling of the Fourier variable is constant in \(t\)
and the rescaling of \(t\) is linear, so it preserves (2).
For fixed \(j\) the root, its derivative, the gaps to its
neighbors, and the positive expression in (2) have uniform
bounds, with positive lower bounds where needed, on the
compact \((\kappa,t)\) range. These follow by simplicity,
analytic dependence, and compactness.

Taylor expansion of the simple roots in (1) now gives the
ordered small positive roots in scaled coordinates as
\[
 y_{m,j}=\beta_j(\kappa,t_m)+m^{-1/2}b_{1j}(t_m)
                        +m^{-1}b_{2j}(t_m)+O(m^{-3/2})
\]
for each fixed finite group of labels. The coefficients
\(b_{1j},b_{2j}\) depend on the profile but are uniformly
bounded with their \(t\)-derivatives, using (1), the analytic
quartic integrals, and the lower derivative bounds at the
simple roots. Rouché in small fixed scaled disks justifies
the root expansion and identification. Outside those disks
on a bounded real scaled interval, uniform convergence
excludes any further roots. Equivalently this identification
can be made locally in \((\kappa,t)\) and then in finitely
many compact neighborhoods.

Since \(t_N-t_n=t_n/(2n)+O(n^{-2})\), the same-label
physical angles satisfy
\[
 \theta_{n,j}-\theta_{N,j}
 =\frac{n^{-7/4}}4
      \left(3\beta_j(\kappa,t_n)-2t_n\partial_t\beta_j(\kappa,t_n)\right)
                       +O(n^{-9/4})>0.                     \tag{3}
\]
Their different-label gaps have order \(n^{-3/4}\) within
each fixed group. Thus the interlacing starts with the
degree-\(N\) zero and holds through every fixed group of
positive labels, uniformly in the moving critical window.

## Small angles outside the bounded critical field

Let \(r=\theta^{1/3}\) for real positive \(\theta\), and
let \(\Lambda=n r^4\). This section uses
\[
                       \Lambda\ge\Lambda_0,\qquad
                       0<\theta\le\theta_0,
\]
where \(\Lambda_0\) is chosen large depending on \(K\) and
\(\theta_0\) is fixed small. The coefficients in the
complex-order bounds are uniform in these ranges.

First convert the centered lattice sum to a fixed interval
integral. With
\[
 a_m(x)=
  \frac{\sqrt m\,\Gamma(m+1)
       \exp(m L_{m,\lambda}(1/2+x)-mH(x))}
       {\Gamma(m(1/2+x)+1)\Gamma(m(1/2-x)+1)}
\]
on a fixed disk, \(a_m\) and its reciprocal and its first
few \(x\)-derivatives are uniformly bounded on smaller
disks for \(|m-n|\le\eta n\). In fact complex-order Stirling
bounds \(m(L_N-L_m)=(m/N-1)A_0+O(1/n)\) there.
At the two integer orders, after scaling by
\(m^{-1/2}e^{-mH(0)}e^{-im\theta/2}\), the polynomial equals
\[
 I_m(\theta)+\epsilon_m(\theta),\qquad
 I_m(\theta)=\int_{-\delta}^{\delta}
      a_m(x)e^{m(H(x)-H(0)+i\theta x)}\,dx,                  \tag{4}
\]
with \(\epsilon_m,\epsilon_m'=O(e^{-c n})\).
Here \(\delta>0\) is fixed small and the remainder estimate
holds on a complex neighborhood of the small real angular
interval. To verify (4), use the residue formula for the
mesh \(1/m\) with the kernel \(\cot(\pi m(x+1/2))\), using
vertical sides at the half-mesh points nearest \(\pm\delta\)
and fixed height \(h>0\). On the horizontal sides the kernel
is \(\mp i+O(e^{-2\pi m h})\), and on the vertical sides
its modulus is at most one. The exterior real lattice
tails and the corner neighborhoods have a fixed real
phase loss by the preceding estimates. Choose \(h\) small.
On the horizontal panels
\(\Re(H(x\pm i h)-H(x))\le Ch^2\), uniformly.
The maximum of \(H(x)-H(0)\) on the real segment is
\(O(K^2/n)\) by the quartic bound, and \(i\theta(x\pm i h)\)
costs at most \(|\Re\theta|h+\delta|\Im\theta|\).
The horizontal kernel error then has strict exponential
loss for sufficiently small \(\theta_0\) and
\(|\Im\theta|\). Replacing the constant-kernel panels
by the real interval changes the integral by the
strict-loss vertical sides; moving the endpoints by at
most one mesh has the same bound. Cauchy's estimate on a
smaller fixed complex neighborhood gives the derivative
bound for \(\epsilon_m\). For real \(\theta\), the error is real
by reflection of the centered polynomial and the integral.

On the deformed contour scale \(x=r v\), and put
\[
 \Psi(v)=r^{-4}(H(rv)-H(0))+iv,\qquad
                 e=(W-4)/r^2 .
\]
The parameter \(|e|\le K/\sqrt{\Lambda}\) is small. On bounded
complex \(v\)-disks the phase and its derivatives equal
\(-\kappa v^4+iv+ev^2/2+O(r^2)\) and their corresponding
derivatives. Consequently there are two continued simple
saddles \(v_+=a+ib\), \(v_-=-a+ib\), close to the quartic
saddles \((4\kappa)^{-1/3}e^{i\pi/6}\) and
\((4\kappa)^{-1/3}e^{5i\pi/6}\). This uses the implicit
function theorem uniformly for the compact \(\kappa\) range,
and reflection gives the common imaginary height.
Let \(\Phi=r^4\Psi(v_+)\). For real \(\theta\) the right
saddle satisfies
\[
               \Im\Phi\asymp r^4,\qquad a\asymp1,\qquad
               \partial_\theta\Phi=i r v_+ .               \tag{5}
\]
The sign in the first estimate comes from the critical
quartic value \(3\kappa v_+^4\) at the unperturbed saddle.
The derivative uses the stationary equation with the
profile held fixed.

Deform the real interval in (4) to height \(rb\) with
vertical endpoint connectors. Let \(J_m\) be the right
half integral, from \(irb\) along the horizontal line
to \(\delta+irb\) and then to \(\delta\). For real \(m,\theta\)
the left half is its conjugate. The centered exact
polynomial is thus \(2\Re\widehat J_m\) where
\(\widehat J_m=J_m+\epsilon_m/2\).

The right horizontal segment has strict Gaussian loss
from the saddle, uniformly in relative complex orders.
For the quartic model the real phase on the horizontal
line \(v=u+ib\) is
\(-\kappa(u^4-6b^2u^2+b^4)-b\). On \(u\ge0\) it has
one nondegenerate maximum at the right saddle. On
bounded \(u\)-intervals, strict loss away from the saddle
and strict negative real Hessian near it persist under
the small \(r^2,e\) perturbation. For \(u\) above a fixed
large constant and at most \(\delta/r\), the real
difference from the saddle is at most \(-c u^4\).
Indeed the quartic dominates \(e v^2/2+iv\), and the
analytic remainder \(O(r^2|v|^6)\) is at most
\(C(\delta+rb)^2|v|^4\). Its imaginary difference is
\(O(u^4)\). The vertical connector has strict loss on the
unscaled phase. In all three parts the bound on the
modulus of the phase difference is a fixed multiple of
the real loss, using stationarity on the local segment.
Therefore replacing \(n\) by \(|m-n|\le\eta n\) preserves
the strict losses for small fixed \(\eta\).

Gaussian integration on the horizontal segment gives
\[
 J_m=e^{m\Phi} r(mr^4)^{-1/2} C_m(1+R_m),\qquad
 C_m=a_m(rv_+)\sqrt{\frac{2\pi}{-\Psi''(v_+)}},\qquad
 |R_m|\le C/\Lambda .                                      \tag{6}
\]
Use the square root continued from the Gaussian integral
on the horizontal direction, where the real part of
\(-\Psi''\) is positive. The other square root uses the
positive branch at real \(m>0\). To bound the remainder,
on a symmetric scaled \(v\)-segment of radius
\(\Lambda^{-5/12}\) expand the analytic amplitude to
first order and the phase to third order about \(v_+\).
The odd first-amplitude and cubic-phase terms integrate
to zero. The fourth-order phase term, squared cubic term,
mixed cubic and linear terms, and second-amplitude term
give relative \(O(1/\Lambda)\) by Gaussian moments.
Beyond this segment the quadratic loss is
\(\exp(-c\Lambda^{1/6})\) until the fixed strict region,
and the quartic loss controls the integration up to
\(\delta/r\). The endpoint loss bounds its contribution
by a relative \(O(1/\Lambda)\) as well. These estimates
hold throughout the relative complex-order disk. Both
\(C_m\) and \(C_m^{-1}\) are uniformly bounded there.
Thus Cauchy's estimates give
\(\partial_m\log C_m=O(1/n)\) and
\(\partial_m\log(1+R_m)=O(1/(n\Lambda))\) on the segment
from \(n\) to \(N\). The integer errors \(\epsilon_m\) have
exponentially small size relative to \(J_m\) by decreasing
\(\theta_0\) in (5)–(6).

Use continuous argument branches \(f_m=\arg\widehat J_m\)
at the two integer orders, paired by continuation in
complex order in (6) and by the exponentially small
integer corrections. Equations (5)–(6) give
\[
                   f_N-f_n=\Im\Phi+O(1/n).
\]
For sufficiently large \(\Lambda_0\) this lies strictly in
\((0,\pi)\). The two phases are also increasing. Indeed
differentiating \(J_m\) in \(\theta\) inserts \(i m x\) in
the integral and a boundary term at \(irb\). The same
Gaussian calculation gives a relative complex first
moment \(x=rv_++O(r/\Lambda)\): the leading centered
linear term is odd, and the first amplitude and cubic
phase corrections each give \(O(r/\Lambda)\).
The boundary term divided by \(J_m\) is \(O(1/\theta)\)
because \(|(rb)'|\le C r/\theta\), its relative scaled
height loss is \(e^{-c\Lambda}\), and the packet
prefactor is of order \(r/\sqrt{\Lambda}\).
The bound on \((rb)'\) follows by differentiating the
stationary equation and using \(|H''(rv_+)|\asymp r^2\).
The differentiated integer error is still exponentially
small relative to the packet. It follows that
\[
                    f_m'=m r a+O(1/\theta)>0
\]
by the same large cutoff. In this angular range the zeros
are exactly those at which \(f_m\) is an odd half-integer
multiple of \(\pi\). At two consecutive lower zeros,
\(0<f_N-f_n<\pi\) and the increasing phases give precisely
one upper zero between them.

## Regular angles and combination of labels

For any fixed closed angular range separated from \(+1\),
take a retained-and-diffuse limit of mirrored profiles
in the moving window. It has \(W=4\). The critical
imaginary-axis divider in the
[mirrored contour argument](MULTIVARIATE_SUPERCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md#the-mirrored-contour-hypotheses)
still has \(Y(v)>0\) for \(v>0\). The displayed formula
there for \((v^2+1/4)Y'\) starts at zero and strictly
increases. Its summands are increasing since every
retained \(a_i>1/2\); a purely diffuse limit has positive
diffuse mass and a strictly increasing diffuse term.

Here is the application of the regular witnesses at that
limit. The proper two-sheet map, the negative-real-gap
inverses, and the divider yield distinct simple left and
right saddles at every interior point on the imaginary
field line. On this line their heights agree by reflection.
It is the unique equal-height point at each imaginary
height, since the real derivative of their phase
difference is strictly signed by their separation in
real part. The two traversals are activated at the
negative-real-ray boundary by the two vertical contour
pieces. The
[finite sublevel-contour proof](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#finite-contour-deformation-for-the-limiting-profile)
propagates them up the connected equal-height line in
the open strip. Its escape and cut checks apply here:
there are retained gaps of both signs if any are retained,
and otherwise the diffuse curvature is positive. The
construction only uses simple saddles and strict-loss
pieces locally, and the negative boundary anchor. In
particular the upper boundary argument for two distinct
real maxima is unnecessary on this closed angular
range. At the negative boundary the same vertical
gap panels supply the witness across the ray.
The compactness and stable strict contour neighborhoods
are as in the [regular cover](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md#finite-contour-witnesses-from-actual-limiting-contours).
This argument also works with real smaller exponents by
the integer-gap and real-parameter compactness checks in
the mirrored contour argument cited above.

In each of finitely many such charts, both integer orders
have the same leading saddle amplitudes. The frozen
multiplier is \(1+O(1/n)\) uniformly on each chart's
contour, by analytic Stirling on its fixed cut-free
neighborhood, so the same leading amplitudes are used
for the interpolated lower order. The exact centered Rice kernel, up to a real scalar, is
\[
 {\cal K}_m(q)=\frac{e^{mL_{m,\lambda}(q)}}{\prod_{k=0}^m(mq-k)}
                e^{im(\theta-\pi)q-im\theta/2}.
\]
For real \(\theta\) it satisfies
\({\cal K}_m(1-\bar q)=-\overline{{\cal K}_m(q)}\).
Indeed the denominator contributes \((-1)^{m+1}\)
and the centered exponential contributes \((-1)^m\)
relative to their conjugates. Reflection followed by
reversal preserves the contour orientation around the
coefficient poles. With the Rice prefactor \(1/i\),
the reflected and reversed period of an oriented
segment is the conjugate of its original period.
Average the two exact winding contours to pair their
saddle traversals and exterior pieces under this
reflection. The nonzero saddle coefficients survive
the averaging, since the original and reflected
representations have identical periods and their
leading coefficients at either saddle are fixed by
uniqueness on the nearby corresponding strict-dominance
side. The leading amplitudes are therefore conjugate
in the centered two-saddle sum. The Rice Stirling
factor \((-1)^m\) cancels the prefactor sign
\((-1)^m\); the remaining real scalar has no effect
on the zeros. In a negative boundary chart the same
identities use the joined real-gap kernels and the
branch \(\theta-\pi\) through zero.
Consequently the leading zero equation on real angles
locally takes the form
\[
                     m d(\theta)+b(\theta)
                             \ \in\ \pi(\mathbb Z+1/2)+O(1/m),
                                                                  \tag{7}
\]
where \(d\) is the imaginary part of the centered phase
of the right saddle. An additive common real
normalization or a phase branch modulo \(\pi\) does not
change the equation. Both \(d\) and \(b\) are shared by
the two orders in the chart.

On \(0<\theta<\pi\), \(d'(\theta)=\Re(q_R-1/2)>0\).
At the critical limit the small-angle quartic saddle
gives \(d(\theta)>0\) near zero and \(d(0)=0\) by
continuity. At \(\theta=\pi\) the right real gap gives
\(d(\pi)=\pi/2\) using the continued entropy logarithm.
Indeed \(1-q_R\) approaches the negative ray from below,
and its entropy term contributes \(\pi(1-q_R)\) to the
imaginary phase; the angle term contributes
\(\pi(q_R-1/2)\). Thus \(d\) is strictly between \(0\)
and \(\pi\) on the fixed chart ranges, also when a chart
extends slightly beyond the negative boundary. The
strict phase margins and positive derivative persist
for sufficiently nearby profiles using their actual
stationary values. No approximation of an exponential
by its limiting phase is needed.

For clarity, the zero description (7) gives exactly one
root near each of the leading phase sites in any compact
subinterval of a chart. On a complex angular disk of
radius \(\epsilon/m\) centered at a leading site, the
centered two-exponential leading sum has a simple zero.
After division by one exponential and a nonzero analytic
amplitude it is bounded below on the boundary by a
positive constant depending only on sufficiently small
fixed \(\epsilon\) and the chart's lower bound for \(d'\).
The relative \(O(1/m)\) contour remainder is smaller
there, so Rouché gives one zero. The real zero equation
places it within \(O(m^{-2})\) of its site and excludes
other real roots between the site disks. The phase
difference between the two leading equations is \(d\),
bounded away from both \(0\) and \(\pi\). Therefore each
gap between adjacent lower zeros contained in the
interior of a chart contains precisely one upper zero.
The chart-sized overlaps account for all regular gaps
because the local root spacing is \(O(1/n)\). These
bounds are uniform by taking the finite compact cover.

To join the three ranges, choose a fixed \(M\) so large
that \(n^{3/4}\theta\ge M\) implies
\(\Lambda\ge\Lambda_0\). Choose enough bounded-field
labels that their last member has scaled root larger
than \(2M\) throughout the compact \((\kappa,t)\) range.
Such a finite choice exists by the infinitude of the
quartic zeros, continuity of each label, and compactness.
The bounded-field label comparison (3) then overlaps
the small-angle phase comparison. That comparison in
turn overlaps the regular charts on a fixed positive
angular interval below \(\theta_0\). Positive phase
derivatives give \(O(1/n)\) spacing on this last
overlap. The negative-ray chart includes the possible
root at \(-1\). These three overlapping ranges give
one upper zero in each consecutive lower-root gap
on the positive angular half, including gaps through
the negative boundary chart. Reflection gives the
negative half. The first and last zeros on the circle
cut at \(+1\) belong to the upper order by (3) and
reflection. The \(n+1\) upper zeros are thereby accounted
for and the strict interlacing follows.
