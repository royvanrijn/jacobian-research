# A finite simplex region forces a stationary Gaussian winner globally

Fix \(k\ge 2\) and a compact correlation band \(J=[a,b]\subset(0,1)\).
Let \(T=(T_i)\) be a centered regular simplex partition with active
space \(E\) of dimension \(k-1\), in any ambient dimension \(d\ge k-1\).
Write \(P\) for orthogonal projection onto \(E\), and let \(H_{ij}\) be
the reference plane separating the scores of labels \(i,j\). Put
\(C_R=\{x:|Px|\le R\}\). Consider a partition \(A\) that obeys the
winner rule for its own adjusted scores \(V_i=T_\rho 1_{A_i}-\beta_i\)
at \(\rho\in J\).

There is an explicit finite \(R=R(k,J)\) such that if \(A=T\) almost
everywhere on \(C_R\), then \(A=T\) globally. This holds for a
stationary winner and so does not require global maximality. The
radius is independent of the ambient dimension. If \(A\) depends
only on \(Px\), the same conclusion holds when \(A\) is sufficiently
close to \(T\) in Gaussian measure on \(C_R\) and every adjusted
pair equality holds on its full reference plane.

For three-cell planar maximizers with the three exact adjusted
reference score lines, this rules out all partitions into
consecutive circular arcs almost everywhere, including with
noncompact defects, unequal prices or arbitrary conditional
circular masses. Every nonsimplex maximizer in this class has
slices uniformly far from consecutive arcs on a fixed annulus
on every compact correlation band. No correlation-derivative
condition is needed.

## A comparison in a growing cylinder

The reference score contrast has a global lower bound on the
winning sector. For \(x\in T_i\), write
\(t=\operatorname{dist}(x,H_{ij})\), \(\sigma=\sqrt{1-\rho^2}\),
\(m=\rho Px/\sigma\), and let \(n\) be the inward normal to
\(H_{ij}\). Let
\[
 \ell=2\int_{T_i\cap E}(z\cdot n)\,d\gamma_E(z)>0 .
\]
By simplex symmetry \(\ell\) depends only on \(k\). Reflecting
the reference \(j\)-score into \(T_i\), then retaining
\(m+(T_i\cap E)\subset T_i\cap E\) in the nonnegative difference,
bounds the contrast below by
\(\int_{T_i\cap E}(1-e^{-2m_N z_N})\,d\gamma_E(z)\), where
\(m_N=\rho t/\sigma\) and \(z_N=z\cdot n\). Use
\(1-e^{-u}\ge u/(1+u)\) followed by Cauchy–Schwarz and
\(\int_{T_i\cap E}z_N^2\,d\gamma_E\le1/2\). This gives
\[
 U_i^T(x)-U_j^T(x)
 \ \ge\ \frac{\ell^2m_N}{\ell+2m_N}
 \ \ge\ \vartheta\min(t,1),
 \qquad \vartheta=\frac{a\ell^2}{\ell+2}.                  \tag{FG1}
\]
Here \(U_i^T=T_\rho 1_{T_i}\) and the final comparison uses
\(m_N\ge at\).
[The saturated reference score margin](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)

If the labels of \(A\) and \(T\) agree on \(C_R\), the
winner inequalities imply \(V_i=V_j\) on an open
relative subset of every reference facet in \(C_R\), by
approaching from its two incident sectors. Gaussian
convolution is real analytic, hence \(V_i=V_j\) on all
of \(H_{ij}\). This is the only use of exact agreement
in obtaining the score equalities. Define
\[
 f_{ij}=(1_{A_i}-1_{A_j})-(1_{T_i}-1_{T_j}),\qquad
 D_{ij}=V_i-V_j-(U_i^T-U_j^T).
\]
Thus \(D_{ij}=T_\rho f_{ij}-(\beta_i-\beta_j)\) vanishes
on \(H_{ij}\), and
\(\beta_i-\beta_j=T_\rho f_{ij}(0)\).

Choose \(\lambda=(1+b)/(2b)>1\), and put
\(h=1-b\lambda=(1-b)/2\), \(\sigma_b=\sqrt{1-b^2}\).
If \(x\in C_{\lambda R}\) and
\(Y=\rho x+\sigma Z\) with \(Z\) standard Gaussian, the event
\(Y\notin C_R\) implies \(|PZ|\ge hR/\sigma\ge hR\). Write
\(p_R=\mathbb P_{\gamma_E}(|PZ|\ge hR)\). As \(f_{ij}\)
is zero on \(C_R\) and \(|f_{ij}|\le2\), on \(C_{\lambda R}\)
\[
 |T_\rho f_{ij}|\le 2p_R,\qquad
 \left|\partial_n T_\rho f_{ij}\right|
 =\left|\frac{\rho}{\sigma}
     \mathbb E[(Z\cdot n) f_{ij}(\rho x+\sigma Z)]\right|
 \le \frac{2b}{\sigma_b}\sqrt{p_R}.                       \tag{FG2}
\]
The derivative identity is the derivative of the conditional
Gaussian density. No probability tail in inactive coordinates
is used. In particular \(|D_{ij}|\le4p_R\). Integrate its
normal derivative from the projection of \(x\) onto
\(H_{ij}\), which stays in the same cylinder along the
whole normal segment. The two bounds show
\[
 |D_{ij}(x)|\le
 \max\{4p_R, (2b/\sigma_b)\sqrt{p_R}\}\min(t,1).            \tag{FG3}
\]
For example \(p_R\le 2^{(k-1)/2}e^{-h^2R^2/4}\) by the
Gaussian exponential moment of \(|PZ|^2/4\). Take \(R_0\ge1\)
large enough that the coefficient in (FG3) is at most
\(\vartheta/2\) at \(R_0\). By (FG1) the adjusted winner is
exactly the reference sector throughout \(C_{\lambda R_0}\)
up to ties of measure zero. The estimates improve with
\(R\). Iterate with radii \(\lambda^mR_0\) to obtain global
agreement.

For the approximate version assume \(A\) depends on \(Px\)
and that \(V_i-V_j=0\) on \(H_{ij}\) from the outset.
At the first step \(R=R_0\), separate the part of \(f_{ij}\)
in \(B_{R_0}\cap E\), where its \(L^1(\gamma_E)\) norm is
bounded by twice the changed-label mass \(\delta\).
On \((B_{\lambda R_0}\cap E)\times(B_{R_0}\cap E)\), the
Mehler kernel and its normal derivative are bounded
uniformly on \(J\). These contributions to the error
in (FG3), including the value at zero, are at most
\(C_{k,J,R_0}\delta\min(t,1)\). Choosing
\(\delta\le\vartheta/(4C_{k,J,R_0})\) forces exact equality
on the larger cylinder by (FG1). The subsequent steps use
(FG2) with zero core discrepancy. This gives a computable
positive threshold for the changed mass from the kernel
and the explicit margin.

This argument uses stationarity through the winner rule and
the exact pair-plane identities. An energy tail bound for
arbitrary distant competitors does not give these identities.

## Consecutive circular arcs force the planar simplex

Let \(k=3\), \(d=2\), and assume actual balanced maximality
and exact reference lines for the adjusted scores. At the
origin all three adjusted scores tie and the three maximal
gradients are distinct. The gradient differences are normal
to the corresponding reference lines. Since their cyclic
sum vanishes and the normals have a unique linear
dependence, these differences have a common nonzero
proportionality factor to those of the reference simplex.
The ordinary pair factorization at the origin therefore
gives the centered reference germ or its antipode, even
with unequal prices.
[Maximal-tie separation at the origin](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#uniform-separation-at-every-maximal-tie)
[The local factorization of the pair contrasts](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#exact-germ-and-uniform-radius)

The circular propagation proof applies with \(V_i-V_j\)
as the pair contrasts. The analytic angular winner and
endpoint arguments depend only on these adjusted differences.
At a possible first exceptional circle on which a contrast
vanishes identically, the two germ labels still meet on
their reference ray all the way to that circle from below.
The reciprocal-gradient integrability on this actual wall
prevents a zero of the pair gradient at its limit. Vanishing
on both the circle and the reference ray would force that
zero. This uses the global regular-wall estimate and uses
central maximal-tie separation only at the origin.
[Propagation of consecutive circular slices](THREE_CELL_ANGULAR_SYNCHRONIZATION.md#circular-connectivity-propagates-the-simplex-germ)
[The reciprocal-gradient bound on active walls, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Consequently, if each label is an arc up to angular null
sets for almost every radius below \(R_0(3,J)\), then the
partition agrees with its germ on that disk. The cylinder
comparison with \(E=\mathbb R^2\), applied to that germ,
forces it on the whole plane. Both orientations have the
same exact reference lines. This proves the exclusion
without assumptions on the exterior.

To make the conclusion uniform on \(J\), let
\(d_{\rm arc}(A;s)\) be the summed angular \(L^1\) distance
of the three indicators on the circle of radius \(s\) to a
partition into one arc per label. Empty and full arcs are
allowed. The family of such triples is compact in angular
\(L^1\), by compactness of the arc center and length
parameters and closedness of the partition identities.
Distance is measurable by taking the infimum over a
countable dense set. Its integral below \(R_0\) is
\(1\)-Lipschitz in summed planar Gaussian \(L^1\).

If that integral tended to zero along a sequence of
nonsimplex exact-line maximizers on \(J\), fix the reference
lines by common isometries and take a strongly convergent
subsequence of actual maximizers. The adjusted equalities
pass to the limit on each line. In fact pressures normalized
by \(\sum_i\beta_i=0\) are
\(\beta_i=U_i(0)-1/3\), and local score convergence follows
from the Mehler formula. The limit has arc slices
almost everywhere on \(B_{R_0}\), hence is the reference
simplex or its antipode globally. But the approximate
cylinder comparison now forces the original maximizers
to equal that same simplex for all sufficiently large indices.
This contradicts their choice. Thus the integral of
\(d_{\rm arc}\) over \(0<s<R_0\) has a positive uniform
lower bound. As \(0\le d_{\rm arc}\le2\), one can discard
a sufficiently small central disk and obtain an annulus
and a set of radii of uniformly positive Gaussian radial
measure on which the angular distance is uniformly
positive.
[Strong compactness of positive-correlation maximizers](POSITIVE_CORRELATION_CONTINUATION.md)

The argument leaves open whether arbitrary maximizing
partitions have exact adjusted reference score planes,
and whether the necessary nonconsecutive slices can
occur at a maximizer with those planes.

## Spherically convex slices propagate an exact germ for all cell counts

Let \(k\ge3\), \(n=k-1\), and let \(A\) be an actual balanced
maximizer in \(\mathbb R^n\) that agrees with \(T\) on some
neighborhood of the origin. Suppose that for almost every
\(0<s<R\) each cell slice \(s^{-1}A_i\cap S^{n-1}\) agrees
up to spherical null sets with the section of a closed convex
cone at the origin. Then \(A=T\) almost everywhere on \(B_R\).
In particular, the spherical convexity hypothesis only needs
to hold below \(R_0(k,J)\) to force global equality. No
lower solid-angle bound or common fan for these slices is
assumed.

The germ makes all adjusted pair differences vanish on their
full reference planes by analyticity, as in the first
comparison. A radius is exceptional if some pair contrast
vanishes identically on its sphere. These radii are locally
finite away from zero: the restriction of each contrast to a
fixed ray in an incident germ sector is analytic and nonzero.
There are no exceptional radii sufficiently close to zero.
At every nonexceptional radius the angular score ties have
spherical measure zero by real analyticity. Thus the
angular winner indicators vary continuously in \(L^1\)
at such a radius.

We first verify how near-reference convex slices are pinned.
Suppose that sections of closed convex cones \(C_m\) converge
to \(T_i\cap S^{n-1}\) in angular \(L^1\). Then for every
\(\eta>0\) all their directions eventually belong to the
\(\eta\)-neighborhood of this reference sector. If a
subsequence has unit directions \(z_m\in C_m\) tending to a
point outside the sector, choose \(n\) linearly independent
unit vectors \(v_j\) strictly inside \(T_i\). Angular
convergence gives \(v_{j,m}\in C_m\) tending to \(v_j\).
The convex hull of \(0,v_{1,m},\ldots,v_{n,m}\) contains a
fixed Euclidean ball \(B\) for all large \(m\), since the
limiting simplex has interior. For fixed small
\(u>0\), the ball \((1-u)z_m+uB\) lies in \(C_m\) and
its radial projection contains a fixed positive-measure
spherical cap outside \(T_i\). This contradicts angular
convergence.

For each pair \(i,j\), fix a direction in the relative
interior of its reference facet where the analytic
normal derivative of \(V_i-V_j\) along that ray is not
identically zero. Such a ray exists since the gradient
is nonzero almost everywhere on the actual germ wall
by the interface variation of [Heilman–Tarter,
Lemma 4.8](https://arxiv.org/pdf/2011.05583). The
derivative is then nonzero at this ray at all radii
outside a discrete set. On an angular convex slice
sufficiently close to the reference, the preceding cone
argument excludes all spectator labels on a neighborhood
of this fixed direction, since it is away from their
reference cones. At a radius with the nonzero pair
derivative, the exact pair equality on \(H_{ij}\) gives
opposite signs across it. The angular derivative is
nonzero there because the pair normal lies in the
sphere's tangent plane. Both adjacent cones therefore
contain their own sides locally. A separating hyperplane
for these two convex cones with disjoint interiors is
through the origin and contains this relative open patch
of \(H_{ij}\cap S^{n-1}\). The patch spans \(H_{ij}\), so
the separator is exactly \(H_{ij}\) and has the reference
orientation. Doing this for every pair puts every angular
cone inside its reference sector. As the sectors and
the cones both partition the sphere almost everywhere,
the slice is exactly the reference.

On the interval from zero to the first exceptional radius
below \(R\), or up to \(R\) when there is none, the radii
at which the slice is sufficiently close to the reference
form a relatively open set containing the germ. The previous
argument gives exact equality there at almost every radius,
and hence everywhere there by angular continuity. The set
is therefore relatively closed as well, proving agreement
throughout that interval.

An exceptional sphere at \(s_0<R\) cannot stop this
propagation. Here it suffices to use the active flat
wall immediately below the sphere, so no claim about
gradient separation at dormant ties is needed. To give
the required local estimate, lift the actual maximizer
cylindrically to \(\mathbb R^{n+1}\), where it is still optimal
by [Heilman–Tarter, Theorem 1.9](https://arxiv.org/pdf/2011.05583).
On a compactly supported smooth patch of a reference
pair wall before \(s_0\), use speed \(g(x)\zeta(z)\), with
\(\zeta\) an odd cutoff of the extra coordinate \(z\).
These speeds have zero volume flux. In fact shifting
the normal threshold from zero to \(t g(x)\zeta(z)\)
gives exactly balanced masses for sufficiently small \(t\)
by oddness and symmetry of the normal Gaussian. Lemma 4.8
and the product Mehler kernel, followed by
\(\zeta(z)\to z\) in \(L^2(\gamma_1)\), give
\[
 \int_{\Sigma_{ij}} a_{ij}g^2\,d\nu
 \ \ge\ 2\rho\iint_{\Sigma_{ij}^2}
       M_\rho(x,y)g(x)g(y)\,d\nu(x)d\nu(y),
 \qquad a_{ij}=|\nabla(V_i-V_j)| ,
\]
where \(d\nu\) is Gaussian surface density on the wall.
Here all tests can stay within the open reference facet.
On a bounded wall region, \(M_\rho\) has a positive lower
bound. Apply the inequality to nonnegative smooth
approximations of \((a_{ij}+\varepsilon)^{-1}\) restricted
to compact subpatches, then exhaust the region and let
\(\varepsilon\downarrow0\). Since \(a_{ij}(a_{ij}+
\varepsilon)^{-2}\le(a_{ij}+\varepsilon)^{-1}\), this gives
finite integral of \(a_{ij}^{-1}\) on the wall up to
the sphere. The argument of (G1) and (G20) works on
this flat patch in every active dimension.
[The extra-coordinate Hermite test (G1)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Now take the pair with contrast identically zero on
the exceptional sphere and a compact angular patch of
its relative open reference facet. At the intersection
with the sphere the contrast vanishes both on
\(H_{ij}\) and on the sphere, whose tangent spaces
span the ambient space, so its gradient is zero.
Along radial lines in the patch smoothness bounds
\(a_{ij}(s\omega)\) above by \(C(s_0-s)\) on the actual
wall \(s<s_0\). Gaussian surface density and the polar
factor \(s^{n-2}\) are bounded below near the positive
radius \(s_0\). The integral of \(a_{ij}^{-1}\) therefore
diverges, contradicting the preceding estimate. This
proves propagation up to \(R\).

The same argument starts from the exterior. A maximizing
compact discrepancy from the centered simplex has exact
adjusted score equalities on the full reference planes,
by continuation from the exterior facets. If almost all
its spherical slices are convex cone sections, the
near-reference pinning propagates inward from the exterior
on the component of nonexceptional radii adjoining
infinity. In fact the compact score perturbation tends
to zero along every ray, so the reference-plane equalities
give equal prices by taking a limit along these planes.
Each pair contrast tends to one along an interior ray
of its first reference cell. Hence no exceptional radii
occur sufficiently far out. An exceptional radius, if
there is one, has a largest value by local finiteness.
The active reference wall approaches that sphere from
above, giving the
same logarithmic divergence with \(s-s_0\) in the
denominator. Thus there is no exceptional radius and
the maximizer is the simplex. This inward version uses
exterior agreement in place of the central germ.
