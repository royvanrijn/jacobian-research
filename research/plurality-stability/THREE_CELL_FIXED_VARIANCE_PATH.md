# A fixed-variance path for compact simplex discrepancies

Let \(A\) be an actual balanced planar Gaussian noise-stability
maximizer at \(0<r<1\) that agrees with a centered regular simplex
\(T\) outside a bounded set. Put \(\sigma^2=1-r^2\) and
\(B=\sigma^{-1}A\). The path
\[
      A(t)=\sqrt{1-t^2}\,B,\qquad 0<t<1,
\]
keeps the conditional Gaussian variance fixed in the source
coordinates of \(B\). Its normalized deficit from \(T\) is entire
in \(t\). For a nonsimplex \(A\) it is positive just above zero
and at \(t=1\), where it gives the strict Euclidean Gaussian
junction deficit. At a tie at \(r\), its first two derivatives
can be compared exactly with the fixed-partition derivatives.

## Entire deficit and its endpoint values

Write \(\gamma=\gamma_2\), \(u=1_T\), \(f=1_B\) for the vector
indicators and set \(h=f-u\). Define the relative integral
\[
 F(z)=\iint e^{z x\cdot y}\bigl[u(x)\cdot u(y)-f(x)\cdot f(y)
                         \bigr]\,d\gamma(x)d\gamma(y).
                                                               \tag{FV1}
\]
Expand its bracket as \(-u(x)\cdot h(y)-h(x)\cdot u(y)
-h(x)\cdot h(y)\). In every term at least one source is in the
bounded support of \(h\). On any bounded complex \(z\)-disk the
integrand and each fixed derivative are dominated by integrable
Gaussian densities times a polynomial and an exponential at
most linear in the unbounded source coordinate. Thus \(F\) is
entire. Changing both observations by the scale \(\sqrt{1-t^2}\)
in the planar Gaussian pair gives
\[
 S_t(T)-S_t(A(t))=(1-t^2)F(t),\qquad 0<t<1 .               \tag{FV2}
\]
The marginal masses of \(A(t)\) can depend on \(t\).

The exact shared-line slice identities give equal \(\gamma\)
masses for \(B\) and a first-moment configuration \(m_i(B)=c w_i\),
where \(w_i\) is the unit center of \(T_i\). For \(T\) let
\(m_i(T)=c_T w_i\) with \(c_T>0\). At a nonsimplex \(B\) the strict
pointwise upper and lower centroid comparisons give
\(|c|<c_T\), since equality would give \(T\) or its antipode and
the exterior excludes the antipode. Consequently
\[
           F(0)=0,\qquad F'(0)=3(c_T^2-c^2)>0.             \tag{FV3}
\]
[Weighted projections and first moments, PG10–PG12](THREE_CELL_PRESSURE_GRADIENT_TAIL.md#moment-constraints-on-the-compact-discrepancy)

At \(z=1\) the product of the Gaussian densities with \(e^{x\cdot y}\)
equals \(G(x-y)/(2\pi)\), with \(G\) the unit-variance Euclidean
Gaussian convolution kernel. The bracket in (FV1) is the
ordinary disagreement of \(B\) minus that of \(T\). The strict
compact-defect comparison for the Euclidean tripod thus gives
\[
              F(1)=\frac{\mathcal J(B)}{2\pi}>0             \tag{FV4}
\]
when \(B\ne T\). These endpoint comparisons allow interior zeros.
[Strict comparison for compact Euclidean tripod changes](EUCLIDEAN_THREE_CELL_JUNCTION.md#strict-comparison-by-fiber-rearrangement)

If the defect of \(B\) is supported in a disk of radius \(R\), the
number of zeros in \(0<t<1\), with multiplicity, is at most
\[
 \frac{\log(5e^{4R^2}/F'(0))}{\log 2}.                       \tag{FV5}
\]
On \(|z|=2\), each mixed term of the expansion of (FV1) is bounded
in absolute value by \(\int e^{2R|y|}\,d\gamma(y)\le4e^{4R^2}\).
The defect-defect term is at most \(2e^{2R^2}\), since the dot
product of two transfer vectors has magnitude at most two.
Thus the entire function \(F(z)/z\) has modulus at most
\(5e^{4R^2}\) on that circle and equals \(F'(0)\) at the center.
Jensen's formula gives (FV5).

The bound can be made uniform over actual nonsimplex maximizers
on a compact positive correlation band with the defect of \(A\)
in a common bounded disk. Under a convergent sequence of
correlations, strong maximizer compactness gives subsequential
label convergence in the fixed plane. The common exterior
passes to the limit. Local measurable growth separates the
nonsimplex sequence from the simplex orbit, so the limit is
also nonsimplex. Rescaling by the convergent \(\sigma\) preserves
label convergence, by local continuity of dilations in \(L^1\)
and the common bounded defect support after rescaling. Hence
\(F'(0)\) is bounded away from zero by (FV3), and the radius
in (FV5) is uniformly bounded.
[Strong planar compactness](SMALL_POSITIVE_CORRELATION_INTERVAL.md#a-rigorous-continuation-boundary-in-correlation)
[Local measurable growth on compact bands](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md)

## Fixed and rescaled derivatives at a tie

Write \(D(t)=S_t(T)-S_t(A)\) for the same shape at all
correlations. On a regular interface of \(A\) let \(N\) be
its oriented normal, \(b=x\cdot N\) its dilation speed, \(a\)
the positive normal score coefficient and \(\nu\) the
Gaussian surface measure. Put \(H=\int a b^2\,d\nu\) and
\(q=Q_0(b)\) for the full ordinary second variation form.
If \(D(r)=0\), then
\[
\begin{split}
 (1-r^2)F'(r)&=D'(r),\\
 (1-r^2)F''(r)-4rF'(r)
   &=D''(r)-\frac{4H}{1-r^2}-\frac{2r^2q}{(1-r^2)^2}.
\end{split}                                                     \tag{FV6}
\]
All pair pressures of \(A\) vanish by the shared exterior
line identity. The unpriced first variation under dilation
is therefore zero at \(r\). If \(G(t,\lambda)=S_t(\lambda A)\),
then \(G_\lambda(r,1)=0\),
\(G_{\lambda\lambda}(r,1)=-2q\), and
\(G_{t\lambda}(r,1)=2\int g b\,d\nu\), with \(g\) the
fixed-shape correlation derivative of the pair score.
For the pure dilation variation, the first indicator derivative
is the oriented signed surface transfer with speed \(b\).
Integration by parts in the second derivative supplies
\(-\int a b^2\,d\nu\) in half the Hessian, since the score
differences themselves vanish on the wall. The interaction
of the first transfer with itself gives the other term of
\(-q\). One may use a smooth compact radial field equal to
\(x\) on the defect and tangent to the unchanged exterior
rays. The mixed derivative is the derivative of the first
surface pairing. The dilation identity (G6) gives
\[
               \int g b\,d\nu=-H/r-rq/(1-r^2).
\]
Along \(A(t)=\lambda(t)A\), where
\(\lambda(t)=\sqrt{1-t^2}/\sigma\), one has
\(\lambda'(r)=-r/(1-r^2)\). The chain rule and (FV2) now
give both identities.
[Zero pair pressures for a compact discrepancy](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)
[Translation and dilation variations, G1–G6](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

The curvature decrement on the right of (FV6) is strictly
positive for every nonsimplex \(A\). Indeed (G1) gives
\(rq+(1-r)H\ge0\), so
\[
   \frac{4H}{1-r^2}+\frac{2r^2q}{(1-r^2)^2}
                \ge \frac{2(2+r)H}{(1-r)(1+r)^2}>0 .
                                                               \tag{FV7}
\]
Here \(H>0\) because vanishing of the normal dilation
speed gives distributional invariance under positive
dilations and the prescribed exterior would then identify
the whole partition with \(T\). Equations (FV6)–(FV7)
give an exact path comparison at either sign of a slope.
They impose no balanced optimality comparison for
\(A(t)\) away from \(r\).

## Finitely many stationary parameters in an affine orbit

Suppose the original compact discrepancy \(A\) is nonsimplex.
Up to centered orthogonal transformations, its invertible affine
images satisfy the stationary interface equations at only finitely
many pairs \((\lambda,s)\), where the image is \(\lambda A\),
\(\lambda>0\) and \(s\in(0,1)\). The stationary equation here is a
constant pair-score difference on each actual pair interface. The
affine image may have different masses from \(A\).

First let the image be \(LA+p\). Each image sector of the
reference simplex has opening less than \(\pi\). On the far reverse
direction of the supporting line of any actual \(i,j\) exterior
ray, the conditional source lies increasingly deep in the third
sector. The analytic pair contrast is constant on the entire
supporting line by stationarity on the ray. Its limit on the reverse
tail is zero, so its pair pressure vanishes. On the forward tail
the local scores are those of two opposite normal halfspaces
separated at offset \(b\), with an error tending to zero. The pair
contrast on the line tends to \(2\Phi((s-1)b/\sqrt{1-s^2})-1\)
when the normal is chosen into \(i\). Hence \(b=0\). This applies
to each pair, so all image lines pass through the origin and
\(p=0\).

Write \(C=LT\) for the centered conical exterior. The zero
pair-score line gives the weighted slice identity (CD1) for the
image partition. On a sufficiently distant positive or negative
longitudinal slice \(u\), the affine image agrees with \(C\) for
all transverse coordinates up to null sets. Homogeneity of its
contrast \(f=1_{C_i}-1_{C_j}\) and rescaling the normal coordinate
then give
\[
      \int_{\mathbb R} f(\operatorname{sgn} u,v)
                  e^{-u^2v^2/[2(1-s^2)]}\,dv=0
\]
for almost every such \(u\). The expression is analytic as a
function of the Laplace parameter \(u^2/[2(1-s^2)]\) in the
right half-plane. Substituting \(v^2\) on each of the two normal
half-lines shows that the even part of \(f\) in \(v\) vanishes.
For uniqueness one may evaluate the analytic identity on a
vertical line of positive real part and apply Fourier uniqueness
to the exponentially weighted integrable function of \(v^2\).
The contrast distinguishes the three labels, so reflection in the
supporting line exchanges the two adjacent sectors. Equality of all
three pairs of incident sector angles makes \(C\) regular. An
invertible linear map sending the three regular rays to three
regular rays is a similarity. Indeed the unit ray vectors sum to
zero in both configurations, which forces equal lengths of the
three images of the original unit vectors. This proves the
affine reduction to \(\lambda A\) up to centered rotations and
reflections.

For a stationary \(\lambda A\) at \(s\), apply the zero-pressure
slice identity in the source coordinates of \(A\). Its transverse
variance is
\[
                            \tau=(1-s^2)/\lambda^2 .
\]
There are only finitely many admissible \(\tau\). Indeed the
nonzero even slice contrast and finite sign-change argument of
(CD2) apply unchanged to these variances, with Laplace parameter
\(1/(2\tau)\) in place of \(1/[2(1-r^2)]\). The slice is selected
from the original fixed \(A\), whose maximizing scores give
the finite step structure for almost every longitudinal coordinate.
[Nonzero even slice contrasts and the finite sign-change bound](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#a-fixed-compact-discrepancy-has-finitely-many-stationary-correlations)

Fix one admissible \(\tau\). For a regular actual \(i,j\) interface
point \(x\) of \(A\), the pair score of \(\lambda A\) at \(\lambda x\)
is
\[
 G_{\tau,x}(s)=\frac1{2\pi\tau}\int_{\mathbb R^2}
   (1_{A_i}-1_{A_j})(y)
       \exp\!\left(-\frac{(y-sx)\cdot(y-sx)}{2\tau}\right)dy .
                                                               \tag{FV8}
\]
The dot product is extended bilinearly for complex \(s\), so this
is entire by Gaussian domination on compact parameter sets. If
infinitely many distinct \(s\) were stationary at this variance,
intersect the full surface-measure sets where their pair-score
equalities hold. Choose \(x\) in this intersection and off its
reference pair line. Such a regular piece exists, since if all
actual interfaces were on their reference lines, the indicators
would have zero distributional derivatives in each open chamber
of those lines. Each chamber extends to infinity, where the
prescribed assignment is \(T\).

At this \(x\), infinitely many distinct \(s\in(0,1)\) would give
zeros of (FV8), so it would vanish identically. But the line
\(\mathbb R x\) has a direction in the interior of \(T_i\) or
\(T_j\). Its two directions cannot both belong to the closed
third sector, and the line differs from the \(i,j\) line.
On the corresponding real \(s\)-tail the Gaussian source of
(FV8) lies with probability tending to one in the corresponding
exterior sector. The contrast tends to \(1\) or \(-1\), a
contradiction. Thus there are finitely many \(s\) for each
\(\tau\), and \(\lambda=\sqrt{(1-s^2)/\tau}\) is determined.
