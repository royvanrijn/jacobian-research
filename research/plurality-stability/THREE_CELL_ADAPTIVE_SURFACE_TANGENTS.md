# Perimeter growth and conical tangents of adaptive surface minimizers

Fix \(0<r<1\) and let \(A\) be a balanced planar
minimizer of the adaptive surface functional. At every
finite point \(p\), its total interface length in
\(B(p,t)\) is \(O_p(t)\) as \(t\downarrow0\).
At points of the interface support it is also bounded
below linearly in \(t\), uniformly on compact sets.
Every sequence of dilations about \(p\) has a subsequence
converging locally in label measure to an ordinary
conical partition. There is a weighted density
\(\theta(p)\) that is positive on the interface support.
Every tangent at these points is nonconstant. The
argument does not classify nonconstant cones.
[Adaptive surface variational problem](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md)
[Stationary mixture and reciprocal tensions (ASV2–ASV4)](THREE_CELL_ADAPTIVE_SURFACE_STATIONARITY.md)

## Linear growth at a fixed point

Let \(\mathcal K\) denote the compact set of
attaining dual tests of \(A\), with their weights
\(w^G_e=|T_{\sqrt r}G_i-T_{\sqrt r}G_j|\) on the pair
\(e=(i,j)\). Denote the total tagged Gaussian measure
by \(\nu_{\rm sum}\). On any fixed small disk, the
weights are bounded and Lipschitz uniformly in the
unit-ball tests, and their evaluations depend
continuously on the tests in Gaussian \(L^2\).

Fix \(p\). First consider the tags for which some
\(G\in\mathcal K\) has \(w^G_e(p)=0\). For such a tag
use this one test in the reciprocal bound (ASV4).
Its local Lipschitz bound gives
\[
 \frac{\nu_e(B(p,t))}{t}
    \le C\int_{B(p,t)}\frac{d\nu_e}{w^G_e}
                      \longrightarrow0 .                    \tag{AT1}
\]
The integral tends to zero because it is locally
finite and \(\nu_e\) has no atoms. On every other
tag, \(w^G_e(p)>0\) for all \(G\in\mathcal K\).
Compactness gives a positive lower bound independent
of \(G\) on those tags, and the same lower bound
up to a constant factor holds on a neighborhood
of \(p\) for tests sufficiently close to
\(\mathcal K\) in Gaussian \(L^2\).

For almost every small \(t\), fill \(B(p,t)\) with
one label. The added interface length is at most
\(2\pi t\) by slicing on the boundary circle. Restore
balance by two exterior flows, outside a fixed disk
about \(p\) with all labels still of positive mass
in its complement. Surjectivity of the exterior
mass derivatives and the estimates (AR1)–(AR2)
apply with arbitrary local labels. Their proofs
use the copy of the direction fields on unchanged
interfaces and the bounded kernel on the compact
change region. The mass correction is \(O(t^2)\).
More specifically, if \(B\) is a local competitor,
\(B^*\) its balanced correction, and \(G\) an
attaining test at \(B^*\), the comparison on the
local change region \(U\) is
\[
 \sum_e\int_U w^G_e\,d\nu_e^A
    \le\sum_e\int_U w^G_e\,d\nu_e^B + C|\delta|,              \tag{AT2}
\]
where \(\delta\) is the vector of mass errors.
The changed circle is included when \(U\) is the
closed filled disk.
[Exterior mass correction and stability of nearby dual tests (AR1)](THREE_CELL_ADAPTIVE_SURFACE_TWO_PHASE_REGULARITY.md#nearby-dual-tests-and-an-exterior-mass-correction)

For the filled competitor, the old local perimeter
tends to zero by non-atomicity, while the added
perimeter and correction also tend to zero.
Thus (AR1) brings its attaining tests arbitrarily
close to \(\mathcal K\). Apply the lower bounds for
the positive tags in (AT2). Their old interface
measure is \(O_p(t)\), since the right side is
\(O_p(t)\). The other tags are \(o(t)\) by (AT1).
Monotonicity of the unweighted interface measure
in the ball radius extends the resulting estimate
from almost every radius to every sufficiently
small radius. Gaussian arclength and ordinary
arclength are comparable near \(p\).

The same comparison gives an estimate for general
competitors changed compactly inside \(B(p,t)\).
If the competitor has greater local perimeter,
the estimate below holds already. Otherwise both
local perimeters are \(O_p(t)\) by the linear
growth bound, so (AR1) again applies. Retaining
the lower bounds on the positive tags and adding
the other tags using (AT1) gives
\[
 P(A;B(p,t))\le C_p P(B;B(p,t))+o_p(t)+C_p t^2 .             \tag{AT3}
\]
Here \(P\) is ordinary interface length, counted
once. The \(o_p(t)\) term depends on the fixed
reference partition and point, through (AT1).
This estimate is centered at \(p\) and has an
additive error even when the mass errors vanish.

## Weighted density and radial varifold limits

Choose the stationary test law \(\pi\) and pressures
from (ASV2). Write \(\bar w_e=\int w^G_e\,d\pi(G)\).
Let \(V\) be the one-dimensional varifold of the
tagged interfaces with weight
\(\mu=\sum_e\bar w_e\,\nu_e\) and their approximate
tangent lines. Its first variation against a compact
smooth ambient field is the integral of the
tangential divergence. From (ASV1)–(ASV2) and the
Gaussian boundary flux formula, the first variation
is integration of a locally bounded vector field
against the raw tagged measures. On pair \(e=(i,j)\)
this field is
\((\beta_i-\beta_j)N_e-\nabla\bar w_e+x\bar w_e\).
Here \(\nabla\bar w_e\) on the interface can be
evaluated by averaging the gradients of the
attaining weights. These gradients exist almost
everywhere for each test and are uniformly bounded
on compacts. Fubini and the bound justify the
average.

The reciprocal bound also gives a local \(L^2(\mu)\)
density for the first variation. In fact Jensen's
inequality and (ASV4) give
\(\int_D\bar w_e^{-1}\,d\nu_e<\infty\) on every
compact \(D\). Dividing the locally bounded raw
vector density on each tag by \(\bar w_e\) therefore
gives an \(L^2(\mu)\) vector \(H\), with sign chosen
so that \(\delta V(X)=-\int H\cdot X\,d\mu\).
This argument allows a tangential component of
\(H\); integer multiplicity is not assumed.

Here is the radial calculation used at an arbitrary
center \(p\). Write \(d=|x-p|\), let \(\psi\) be a
nonnegative smooth nonincreasing cutoff supported
in \([0,1)\), constant near zero and bounded by one,
and put \(I(t)=\int\psi(d/t)\,d\mu\). Let
\(X_t(x)=\psi(d/t)(x-p)\). If \(P_\tau\) denotes
orthogonal projection onto the marked tangent line,
differentiation gives
\[
 tI'(t)-I(t)=-\delta V(X_t)
     +\int\frac d t[-\psi'(d/t)]
       \bigl(1-|P_\tau((x-p)/d)|^2\bigr)\,dV .                \tag{AT4}
\]
The integrand involving \(d^{-1}\) is taken as zero
at the center. The second term is nonnegative.
In particular, the fixed-point raw upper bound
gives \(|\delta V(X_t)|\le C_p t^2\), so
\(I(t)/t+C_p t\) is nondecreasing for small \(t\).

There is also a bound uniform in the center on
compact sets. With \(\Lambda\) the \(L^2(\mu)\) norm
of \(H\) on a compact neighborhood containing the
balls, Cauchy–Schwarz yields
\(|\delta V(X_t)|\le t\Lambda\sqrt{I(t)}\).
By (AT4), \(Q(t)=I(t)/t\) thus satisfies
\(Q'(t)\ge-\Lambda\sqrt{Q(t)/t}\).
Apply this to \(\sqrt{Q+\eta}\) and let
\(\eta\downarrow0\) after integrating. It follows
that \(\sqrt{I(t)/t}+\Lambda\sqrt t\) is
nondecreasing. Approximating the indicator of the
open ball in either monotonicity statement gives a
well-defined finite limit
\[
                      \theta(p)=\lim_{t\downarrow0}
                                   \frac{\mu(B(p,t))}{t}.
                                                                \tag{AT5}
\]
Finiteness follows also from the raw perimeter
upper bound and locally bounded weights.

In particular, suppose \(p\) is in the support
of a tag with \(\bar w_e(p)>0\). Continuity gives
a positive lower bound on this tag's weight near
\(p\). Its reduced-boundary points have ordinary
one-dimensional density one almost everywhere
for its arclength measure. Thus for every small
\(t\) one can choose \(q\in B(p,t/2)\) with
weighted density bounded below by a fixed positive
constant, using the positive Gaussian factor and
the weight on that tag. The uniform square-root
monotonicity estimate at \(q\) makes
\(\mu(B(q,t/2))\) bounded below by a positive
constant times \(t\), for small \(t\).
Since this ball is in \(B(p,t)\), it follows that
\(\theta(p)>0\).

For a sequence \(t_n\downarrow0\), translate by
\(-p\), dilate by \(1/t_n\), and multiply the
pushforward measure of \(V\) by \(1/t_n\).
These rescaled varifolds have locally bounded mass
by the upper bound. Every weak subsequential limit
\(V_0\) is stationary: the rescaled first variation
of a test \(Y\) is the original first variation of
\(Y((x-p)/t_n)\), which is \(O_Y(t_n)\) by (ASV2)
and the raw perimeter upper bound. The mass of
\(V_0\) on each centered ball \(B_R\) is
\(\theta(p)R\), first at almost every radius by
weak convergence and then at all radii by
continuity of this linear expression. Hence
\(I_0(R)/R\) is constant for the radial cutoffs in
(AT4), and stationarity makes the tilt term zero.
The tangent line is radial \(V_0\)-almost everywhere.
In fact \(V_0\) is a cone: stationarity tested
against \(Y(x)=\chi(|x|)h(x/|x|)x/|x|\), supported
away from the origin, gives zero integral of
\(\chi'(|x|)h(x/|x|)\). Disintegrating in polar
coordinates gives Lebesgue radial measure times a
fixed angular measure. There is no atom at the
origin by the linear ball mass.

## Homogeneity and nonconstancy of cell limits

Let \(A^{(n)}=(A-p)/t_n\). Its perimeters on every
fixed compact set are bounded by the raw upper
bound. Local \(BV\) compactness gives a subsequence
converging strongly in local label measure to an
ordinary partition \(A^0\). Extract locally weak
limits of the rescaled Euclidean pair arclength
measures marked by their oriented normals.
For the tags of (AT1) these measures vanish in
the limit. For every remaining tag, continuity of
\(\phi_2\bar w_e\) at \(p\) identifies the contribution
to \(V_0\) as the normal-marked limit with the
positive constant weight \(\phi_2(p)\bar w_e(p)\),
with the normal mark replaced by its orthogonal
tangent. Thus each surviving normal is orthogonal
to the radial direction, since the full varifold
has radial tangents and all these weights are
positive. Passing the distributional cell
gradients to the limit gives
\[
                 y\cdot D1_{A_i^0}=0
                     \quad\text{on }\mathbb R^2\setminus\{0\}.
                                                                \tag{AT6}
\]
In polar coordinates, test this equality on
annuli with radial vector fields
\(\chi(r,\omega)\omega/r\). Their divergence is
\(\partial_r\chi/r\). The polar area factor
cancels \(1/r\), showing that each indicator
has zero weak derivative in \(r\) on annuli.
The indicators are therefore invariant under
positive dilations, up to null sets.

If \(\theta(p)>0\), a pure tangent would contradict
minimality. Indeed, suppose along \(t_n\) the
rescaled partition converges locally to a single
label \(i\). In the annulus of radii \(t_n,2t_n\)
the Euclidean volume of labels other than \(i\) is
\(o(t_n^2)\). By coarea one can choose a radius
\(\ell_n\in[t_n,2t_n]\) where the trace of these
labels has arclength \(o(t_n)\), and where the
finite-perimeter gluing formula applies. Filling
this ball with \(i\) adds only \(o(t_n)\) of
interface and changes the volume by \(o(t_n^2)\).
Apply (AT2) with the exterior correction. As
before, all positive tags have weights bounded
below for the filled competitor's attaining test.
Their former total perimeter is bounded below by
a positive multiple of \(t_n\), because
\(\mu(B(p,\ell_n))/\ell_n\to\theta(p)>0\), the
other tags have \(o(t_n)\) measure, and all
averaged weights are bounded above. This
contradicts (AT2).

Conversely, if \(\theta(p)=0\), every tag with
positive weight at \(p\) contributes \(o(t)\)
raw perimeter, by (AT5) and continuity of that
weight. The other tags already have \(o(t)\)
perimeter by (AT1). All tangent partitions
therefore have zero perimeter and are constant.
Moreover the relative isoperimetric inequality
in the ball makes each of the three volume
fractions tend uniformly to either zero or one
as the perimeter divided by radius tends to
zero. For all sufficiently small radii there is
a unique label of fraction greater than \(9/10\).
Continuity of the fractions in the radius
prevents this label from changing. Hence the
pure density label is unique at a null weighted
density point.

## A linear lower bound on the support

Let \(E\) denote the support of the total interface
measure. This is also the support of \(\mu\), since
the averaged weights are positive almost everywhere
on their tagged interfaces. On each fixed compact
set there are \(c>0\) and \(r_0>0\) such that
\[
        P(A;B(p,t))\ \ge\ c t
                \qquad(p\in E,\ 0<t<r_0).                   \tag{AT7}
\]
The proof chooses a circle at a Lebesgue point
where the coarea count of reduced interfaces is zero.

Fix a compact set of centers and a slightly larger
compact neighborhood. All the first-variation
density bounds are uniform in this neighborhood.
The reciprocal measure
\(\sum_e\bar w_e^{-1}\,d\nu_e\) there is finite and
has no atoms. Its mass on balls therefore tends
to zero uniformly in their centers as their radii
tend to zero.

Use a fixed threshold \(0<b<1/4\). Suppose
\(p\in E\) and \(P(A;B(p,t))<b t/16\).
Choose a reduced-interface point \(q\in B(p,t/4)\)
where the ordinary interface length has lower
one-dimensional density at least one. Such points
have full interface measure. Put \(\ell=t/4\) and
\[
 a=\sup\{s\in(0,\ell):P(A;B(q,s))\ge b s\}.
\]
The lower density at \(q\) makes \(a>0\).
For \(s\in[\ell/4,\ell)\) the original bound at
\(p\) gives \(P(A;B(q,s))<b\ell/4\le b s\).
Hence \(a\le\ell/4\). Left continuity of the mass
of open balls gives \(P(A;B(q,a))\ge b a\),
and \(P(A;B(q,2a))<2b a\) by the definition of
the supremum.

Let \(N(s)\) count reduced-interface points on
\(\partial B(q,s)\), summing over tags. The
rectifiable coarea formula for the distance to
\(q\) gives
\[
           \int_a^{2a}N(s)\,ds\le P(A;B(q,2a))<2b a .
\]
The integer count vanishes on a set of positive
measure in this interval. Choose \(R\) in that
set which is a Lebesgue point with \(N(R)=0\)
and at which the sphere has zero interface measure.
Then \(P(A;B(q,R))\ge b a\ge b R/2\).

Use smooth approximations to the radial field
\((x-q)1_{B(q,R)}(x)\) in the weighted first
variation. Concretely, multiply \(x-q\) by a
cutoff changing from one to zero on the radial
interval \([R-\eta,R+\eta]\), with derivative
bounded by \(C/\eta\). The term from this
derivative in the tangential divergence tends
to zero. In absolute value it is bounded by
\[
        \frac{C(R+\eta)}{\eta}
                 \int_{R-\eta}^{R+\eta}N(s)\,ds .
\]
To see this, the tangential derivative term
contains the square of the tangent projection
of the radial unit vector. This is bounded by
the absolute projection, the coarea Jacobian
of the distance function. The remaining
Gaussian factors and tensions are bounded on
the fixed compact neighborhood. The bound tends
to zero by the choice of the Lebesgue point.
The non-derivative part of the divergence tends
to \(\mu(B(q,R))\), since the sphere has zero
interface measure. The raw first-variation
density bound therefore gives
\[
                   \mu(B(q,R))
                        \le C R\,\nu_{\rm sum}(B(q,R)).
\]
By Cauchy–Schwarz applied to the tagged measures,
\[
 \sum_e\int_{B(q,R)}\frac{d\nu_e}{\bar w_e}
     \ \ge\ \frac{\nu_{\rm sum}(B(q,R))^2}{\mu(B(q,R))}
     \ \ge\ \frac{\nu_{\rm sum}(B(q,R))}{C R}
     \ \ge\ c_0>0 .
\]
The last lower bound uses the perimeter bound at
\(R\) and the positive minimum of the Gaussian
density on the compact neighborhood. The ball
\(B(q,R)\) lies in \(B(p,t)\). The lower bound
contradicts the uniform small-ball mass of the
reciprocal measure when \(t\) is sufficiently
small, proving (AT7).

If \(\theta(p)=0\), the argument after (AT6) gave
\(P(A;B(p,t))/t\to0\). Thus (AT7) rules out every
zero-density point on \(E\). In particular there
is no support point at which all the locally
supported tags have zero averaged tension.
