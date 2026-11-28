# Geometry forced by three exact Gaussian score lines

For a three-cell planar maximizer, the three full adjusted
score equalities impose restrictions on the strict-winner
components and on the directions of transfers out of its
exact simplex germ. These restrictions apply with unequal
stationary prices and without agreement at infinity.

For a compact discrepancy from the centered simplex, high
angular harmonics give another restriction that uses no
stability inequality. A pair score difference is nonconstant
on every centered circle. Consequently a stationary
compact discrepancy with consecutive circular slices is
the simplex even without maximality.

## Each strict-winner cell is disconnected away from the simplex

Let \(A\) be an actual balanced planar maximizer at
\(0<\rho<1\), and suppose its adjusted scores
\(V_i=T_\rho1_{A_i}-\beta_i\) tie on each full
reference pair line. The lines are those of a centered
regular simplex. The central tie makes the three gradients
distinct. Since each gradient difference is normal to its
exact pair line, the cyclic sum of the differences forces
one common nonzero proportionality factor to those of the
regular triangle. The pair factorization at the origin
makes the strict winners equal to one of the two simplex
orientations. In what follows \(T\) denotes that orientation.
[Exact germ and pair factorization](GAUSSIAN_FINITE_SIMPLEX_GERM_PROPAGATION.md#consecutive-circular-arcs-force-the-planar-simplex)

Let \(C_i\) be the connected component of the strict-winner
set for label \(i\) containing its germ near the origin.
The two exact pair lines involving \(i\) miss this
strict-winner set. They separate the plane into open
sectors, and \(C_i\) therefore lies completely within \(T_i\).
If the mass of \(A_i\setminus T_i\) is zero, balance implies
\(A_i=T_i\) almost everywhere, and one-cell rigidity
forces the entire partition to be \(T\). Consequently
every nonsimplex maximizer in this class has all three
strict-winner cells disconnected.
[One simplex cell determines the partition](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition)

This argument is uniform on any compact correlation
band \(J\subset(0,1)\). Indeed, the uniform one-cell
separation from a centered simplex wedge gives a
\(\delta_J>0\) such that
\[
 \gamma(A_i\setminus C_i)
   \ \ge\ \gamma(A_i\setminus T_i)
   \ =\ \gamma(T_i\setminus A_i)
   \ \ge\ \delta_J ,\qquad i=1,2,3,                         \tag{EL1}
\]
where the symmetric-difference separation has been
divided by two. Thus every reference sector also
contains a positive mass of other labels. No assertion
about boundedness of the additional components follows.
[Uniform one-cell separation on compact bands](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition)

## Both cyclic transfer orientations occur

Set \(\sigma^2=1-\rho^2\), \(\mu=\gamma_{\sigma^2 I_2}\),
and write \(h_i=1_{A_i}-1_{T_i}\). Let \(e_i\) be the
unit outward ray between reference labels \(i,i+1\),
using cyclic subscripts. The exact adjusted line gives
a constant \(\beta_i-\beta_{i+1}\) for the unadjusted
pair score on that line. The conditional convolution
factors in the coordinates along and perpendicular
to the line. Injectivity of the longitudinal Gaussian
convolution shows that the transverse
\(\gamma_{\sigma^2}\)-average of \(h_i-h_{i+1}\) is
\(\beta_i-\beta_{i+1}\) on almost every longitudinal
slice. The reference contrast has transverse average
zero by reflection. Hence for any \(g\) integrable
against \(\gamma_{\sigma^2}\),
\[
 \int_{\mathbb R^2}\sum_i
    g(e_i\cdot x)\bigl(h_i(x)-h_{i+1}(x)\bigr)d\mu(x)
 =\left(\int g\,d\gamma_{\sigma^2}\right)
               \sum_i (\beta_i-\beta_{i+1})=0.              \tag{EL2}
\]
The individual equalities retain the pressure terms.
Only their cyclic sum eliminates them.
[Weighted Gaussian slice deconvolution](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)

At a changed point in \(T_i\), call the transfer
forward if its new label is \(i+1\), reverse if it
is \(i-1\). Take \(g(t)=t_+^p+2^p t_-^p\) in
(EL2), with \(p>1\). Write \(s=e_i\cdot x\),
\(t=e_{i-1}\cdot x\), and \(u=e_{i+1}\cdot x\).
On the interior of \(T_i\), \(s+t+u=0\) and
\(u<\min(s,t)\). Forward transfers contribute
\(g(t)+g(u)-2g(s)>0\) and reverse transfers
contribute \(2g(t)-g(s)-g(u)<0\). For the strict
inequalities, when \(s,t\ge0\) use
\(g(u)=2^p(s+t)^p>2\max(g(s),g(t))\).
If \(s>0>t\), use \(s=|t|+|u|\) and \(|t|<|u|\):
strict convexity gives \(g(t)+g(u)>2s^p\), and
\(g(t)<\min(g(u),g(s))\). The other order follows
by exchanging \(s,t\). Thus (EL2) forces positive
mass of both transfer orientations whenever \(A\)
differs from \(T\).

The two transfer masses have a uniform positive
lower bound for actual nonsimplex maximizers on a
compact band. Otherwise fix \(T\) by common
isometries and take a strong maximizer limit.
The prices normalized by \(\sum_i\beta_i=0\)
are \(U_i(0)-1/3\), so the exact adjusted lines
pass to the limit. The uniform central germ
radius keeps the orientation fixed. The limiting
partition is nonsimplex by (EL1), but by strong
convergence at least one orientation has zero
mass. This contradicts (EL2).
[Strong compactness and uniform local growth](POSITIVE_CORRELATION_CONTINUATION.md)
[The compact version of the cyclic test, CD3–CD4](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#both-cyclic-orientations-occur-in-the-defect)

## High angular harmonics exclude centered circles of pair-score equality

Let \(A\) now be any ordinary partition agreeing almost
everywhere with \(T\) outside the disk \(B_S\), and fix a
pair \(i,j\). Here the partition is not assumed to
satisfy any extremality condition. Use \(\eta>0\) for
an observation radius and \(u>0\) for a source radius.
For \(f=1_{A_i}-1_{A_j}\), write
\[
 \widehat f_m(u)=\frac1{2\pi}\int_0^{2\pi}
      f(u\cos\theta,u\sin\theta)e^{-im\theta}\,d\theta.
\]
Put \(a=\rho\eta/\sigma^2>0\). The \(m\)-th angular
coefficient of its score on the observation circle is
\[
 \widehat{T_\rho f}_m(\eta)=
    \frac{e^{-\rho^2\eta^2/(2\sigma^2)}}{\sigma^2}
       \int_0^\infty \widehat f_m(u)\,
                  I_m(au)e^{-u^2/(2\sigma^2)}u\,du,          \tag{EL3}
\]
where the angular integral of
\(e^{z\cos\theta}e^{-im\theta}\) divided by \(2\pi\)
is the modified Bessel function \(I_m(z)\).
This follows directly from the conditional Gaussian
density in polar coordinates.

For \(g=1_{T_i}-1_{T_j}\), the source coefficient is
independent of \(u\). For \(m>0\),
\[
 \widehat g_m(u)=\frac{\sin(m\pi/3)}{\pi m}
      (e^{-im\theta_i}-e^{-im\theta_j}),
 \qquad |\widehat g_m(u)|=\frac{3}{2\pi m}
                            \quad(m\not\equiv0\bmod3),       \tag{EL4}
\]
where \(\theta_i,\theta_j\) are the sector centers.
The size identity uses their angular difference \(2\pi/3\).
For \(z\ge0\), expanding the angular exponential as
\(e^{(z/2)e^{i\theta}}e^{(z/2)e^{-i\theta}}\) gives
\[
 \frac{(z/2)^m}{m!}\ \le I_m(z)\ \le\
                  \frac{(z/2)^m}{m!}e^{z^2/4}.
\]
Write \(q_\eta=e^{-\rho^2\eta^2/(2\sigma^2)}/\sigma^2\).
Using (EL4) and the lower bound in (EL3), for
\(m\not\equiv0\bmod3\) the reference contribution has
absolute value at least
\[
 q_\eta \frac{3}{2\pi m}\frac{(a/2)^m}{m!}\,
             \sigma^2(\sqrt2\,\sigma)^m\Gamma(1+m/2).
\]
The contribution of \(f-g\), supported in \(B_S\)
with \(|f-g|\le2\), has absolute value at most
\[
 q_\eta\,C_{a,S}\frac{(a/2)^m}{m!}S^m
\]
by the upper bound for \(I_m\) on \(0\le u\le S\).
Take \(S\ge1\) if necessary. The ratio of the reference
lower bound to the error upper bound tends to infinity
by the gamma function growth. Thus every sufficiently
high harmonic with \(m\not\equiv0\bmod3\) survives.
In particular \(T_\rho f\) is nonconstant on every
centered circle. Subtracting any stationary price
difference leaves these coefficients unchanged.

## All far-field score orderings agree with the reference

For a stationary compact discrepancy, all pair contrasts
eventually have the reference signs, including comparisons
between the two labels that lose to the spectator. Only
the corresponding reference line remains an equality
outside a sufficiently large disk. The radius can be
uniform for discrepancies supported in a fixed disk and
correlations in a compact band.

Fix the \(i,j\) pair and write \(v=x\cdot n\) for the
signed distance to its reference line, with \(n\) pointing
into \(T_i\). First take \(v>0\). In this half-plane
the angular distance to the reference cone \(T_i\) is
at most \(\pi/3\), so
\(\operatorname{dist}(x,T_i)^2\le3|x|^2/4\).
Reflect the reference \(j\)-integral into \(T_i\).
The integrand for the difference is the conditional
Gaussian density at \(y\in T_i\) multiplied by
\(1-\exp(-2\rho v(y\cdot n)/\sigma^2)\).
Retain only a fixed-radius ball obtained by adding a
fixed interior ball of \(T_i\) to the orthogonal
projection of \(\rho x\) onto its closed cone. Such
a ball lies in \(T_i\) by convexity and has \(y\cdot n\)
bounded below by a positive constant. On it the
multiplier is bounded below by \(c\min(v,1)\), and
\[
 |y-\rho x|\le \rho\operatorname{dist}(x,T_i)+C .
\]
Integration gives, with constants allowed to depend
on the correlation band,
\[
 U_i^T(x)-U_j^T(x)
 \ge c\min(v,1)\exp\!\left(-\frac{3\rho^2}{8\sigma^2}|x|^2
                               -C|x|-C\right).            \tag{EL5}
\]

The compact contribution \(T_\rho(f-g)\) and its
normal derivative are bounded by
\(C\exp(-\rho^2|x|^2/(2\sigma^2)+C|x|)\) for \(|x|\ge1\),
by the conditional kernel formula on the bounded source
disk. As in the next section, exterior stationarity
makes the pair pressures equal and both contrasts
zero on the full reference line. If \(0<v<1\), integrate
the derivative from the projection of \(x\) onto that
line. The segment stays within unit distance of \(x\).
This gives
\[
 |T_\rho(f-g)(x)|\le C\min(v,1)
               \exp\!\left(-\frac{\rho^2}{2\sigma^2}|x|^2
                                 +C|x|\right).
\]
The extra Gaussian rate makes this error strictly
smaller than (EL5) outside a large disk. For \(v<0\),
interchange \(i,j\). This proves the claimed ordering
without a lower bound on the distance to the line.

## Circular propagation for stationary compact discrepancies

Assume in addition that \(A\) satisfies its adjusted
score winner rule and that almost all of its circular
slices consist of one arc per label, allowing empty arcs.
Exterior agreement gives adjusted equality along each
reference ray, and real analyticity extends it to the
complete reference line. The unadjusted difference there
tends to zero on the forward ray by the exterior
binary-halfspace limit, so the pressures are equal.
The normal derivative of the unadjusted pair contrast
along this ray has a strictly positive limit in the
inward normal orientation. This follows by differentiating
the Gaussian conditional density: its expectation is
\(\rho/\sigma\) times the expectation of the normal
Gaussian coordinate against the source contrast, which
converges to its absolute first moment as the source
distribution sees the two reference halfspaces.
Thus the normal derivative restricted to the complete
line is analytic and not identically zero. Its zero set is
discrete.
[Exact lines for a stationary compact discrepancy](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)

By (EL3)–(EL4) there are no centered circles of identically
tied pair contrasts. On each circle the finite angular
tie sets are null. Consequently the angular strict-winner
partition is continuous in radius in angular \(L^1\), by
dominated convergence of the pointwise strict winners
away from those tie sets.

Let \(I\subset(0,\infty)\) be the set of radii where all
three widths lie strictly between zero and \(\pi\).
It is open by continuity and contains a final interval
by exterior agreement. At almost every radius in \(I\)
the three arcs tile the circle, and the angular derivative
of each contrast at the two intersections with its
reference line is nonzero. For the \(i,j\) pair, their
spectator arc of width less than \(\pi\) misses a
neighborhood of at least one of those antipodal points.
There the two pair scores have opposite orders on the
two sides of their transverse equality, and the only
possible labels nearby are \(i,j\). So the chosen
point is the unique endpoint joining those two arcs.
This holds for every pair. Each arc width is therefore
a positive integer multiple of \(\pi/3\); because every
width is strictly less than \(\pi\) and their sum is
\(2\pi\), they all equal \(2\pi/3\). The continuous
widths thus equal \(2\pi/3\) everywhere on \(I\), making
\(I\) closed in \((0,\infty)\) and hence equal to it.

At almost every radius, the three endpoints are now
spaced \(2\pi/3\) apart along the prescribed lines, and
the incident pair labels identify the arcs as those of
\(T\) or its antipode. Angular continuity and exterior
agreement force \(T\) at every radius. This classifies
the stationary compact branch under the arc-slice
hypothesis without using maximality or an interface
reciprocal-gradient estimate.
[The endpoint propagation argument for actual maximizers](THREE_CELL_ANGULAR_SYNCHRONIZATION.md#circular-connectivity-propagates-the-simplex-germ)
