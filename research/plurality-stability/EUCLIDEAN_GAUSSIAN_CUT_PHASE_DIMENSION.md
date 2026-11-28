# Phase counts for Euclidean Gaussian-cut minimizers in low dimension

Let \(A\) be a finite-label partition of \(\mathbb R^d\)
that minimizes the unit Euclidean Gaussian cut against
all compact changes, with free volumes. Suppose it
depends only on its orthogonal projection to an
\(m\)-dimensional affine space. The number of its
nonempty phases is at most two when \(m=1\), at most
three when \(m=2\), and at most four when \(m=3\).
There is no homogeneity assumption.

The bridge to local perimeter is dimension-independent.
After dilation about any fixed origin and passage to a
subsequence, the limit minimizes perimeter against
compact changes and every original phase meets the
same origin in positive proportion at all radii.
Consequently any bound on the number of phases incident
to one point of a perimeter minimizer bounds the number
of nonempty phases of a Euclidean Gaussian-cut minimizer
in that dimension.

## A cylindrical minimizer gives a minimizing base

Write \(G^{(n)}_\varepsilon\) for the density of
\(N(0,\varepsilon^2 I_n)\). In a cylinder, let \(x\)
be the base variable and \(z\in\mathbb R^s\) the
unused variable. Suppose a base replacement \(B\)
differs from \(A\) only on a bounded set \(D\) and
has strictly negative relative ordered cut
\(\mathcal J(B,A)\). In the full space, put that
replacement over \(Q_L=[-L,L]^s\) and keep the old
partition elsewhere. The changed set is bounded.

For pairs with both transverse coordinates in \(Q_L\)
the relative cut is exactly
\((|Q_L|-\eta_L)\mathcal J(B,A)\), where
\[
 \eta_L=\int_{Q_L\times Q_L^c}G^{(s)}_1(z-z')\,dz\,dz'.
\]
For crossing transverse pairs the cut difference is
bounded in absolute value by \(2|D|\eta_L\), since the
base coordinate of the changed endpoint must lie in
\(D\). For each displacement \(h\), the set of
\(z\in Q_L\) with \(z+h\notin Q_L\) has volume at most
\(|h|\mathcal P(Q_L)\), by the directional translation
bound for an indicator. Integrating against \(G^{(s)}_1\)
gives \(\eta_L=O(L^{s-1})\). Thus the full cut difference
is negative for large \(L\), a contradiction. Invariance
of the Euclidean kernel under translations permits moving
an affine base through the origin. It therefore suffices
to prove the phase counts in dimension \(m\).

## Every phase occupies a fixed fraction at large radii

Use the unit Gaussian \(G=G^{(m)}_1\) in the base and
write \(k\) for the number of nonempty phases.
A minimizing partition assigns almost every point to
a largest convolution score \(P_i=G*1_{A_i}\).
Indeed, a transfer on a bounded set with strictly
positive score gap would decrease the cut at first
order in its volume, while the quadratic interaction
error is bounded by a constant times the square of
that volume. Hence for each phase there is a point
\(x_i\) at which \(P_i(x_i)\ge1/k\).
[Euclidean score condition](EUCLIDEAN_THREE_CELL_JUNCTION.md#homogeneous-minimizing-junctions-and-stationary-compact-defects)

Fix \(i\), abbreviate \(x_i=x\) and
\(V(r)=|A_i\cap B_r(x)|\), and put \(E_r=A_i\cap B_r(x)\).
For \(l\ne i\), recolor \(E_r\) into label \(l\).
Half of the relative ordered cut is
\[
 L(r)-\int_{E_r\times A_l}G(u-v)\,du\,dv,\qquad
 L(r)=\int_{E_r\times(A_i\setminus B_r(x))}
                       G(u-v)\,du\,dv .
\]
Minimality and summation over \(l\ne i\) imply
\[
 P_G(E_r):=\int_{E_r\times E_r^c}G(u-v)\,du\,dv
 \le kL(r).                                                     \tag{ED1}
\]

For every bounded measurable set \(E\) of positive
volume \(V\), unit-scale mollification \(v=\psi*1_E\)
with a fixed smooth unit-ball kernel satisfies
\(\|v-1_E\|_1+\|\nabla v\|_1\le C P_G(E)\).
Both terms are controlled by the average of
\(\|1_E(\cdot-z)-1_E\|_1\) over the unit ball,
using \(\int\nabla\psi=0\) for the derivative term
and the positive lower bound for \(G\) on that ball.
If \(\|v-1_E\|_1\ge V/4\) this gives a lower bound
linear in \(V\). Otherwise \(v\ge1/2\) on a subset
of \(E\) of volume at least \(V/2\). The \(BV\)
Sobolev inequality bounds \(\|v\|_{m/(m-1)}\) by
a constant times \(\|\nabla v\|_1\), with \(L^\infty\)
when \(m=1\). In either case
\[
 P_G(E)\ge c_m\min(V,V^{(m-1)/m}).                              \tag{ED2}
\]
This is the dimension-\(m\) form of the phase-density
estimate for the Gaussian near-one problem, here with
constant weight and no volume correction.
[Phase-density argument, P6](THREE_CELL_NEAR_ONE.md#uniform-phase-density-and-two-label-blow-ups-on-open-rays)

Integrating \(L(r)\) for \(R/2\le r\le R\) bounds
the integral by \(C_m V(R)\): each pair contributing
has \(u\in A_i\cap B_R(x)\), and the interval of
eligible radii has length at most \(|u-v|\). The
Gaussian first moment is finite. Combining with
(ED1)–(ED2) yields
\[
 V(R)\ge c_{m,k}R
       \min(V(R/2),V(R/2)^{(m-1)/m}).
\]
Choose \(R_0\) so that the unit Gaussian tail outside
\(B_{R_0}\) is less than \(1/(2k)\). The score at \(x\)
implies \(V(R_0)\ge [2k\|G\|_\infty]^{-1}\).
For \(R\ge2R_0\) the minimum is bounded below by
a positive constant times \(V(R/2)^{(m-1)/m}\).
At dyadic radii, \(V(R)/R^m\) is therefore bounded
below by a positive constant times its preceding
value to the power \((m-1)/m\). Iteration gives
a uniform lower bound, and monotonicity of \(V\)
covers the intervening radii. For every label \(i\),
\[
           |A_i\cap B_R(x_i)|\ge c_{m,k}R^m
                  \quad(R\ge R_0).                             \tag{ED3}
\]

## Dilation limits minimize perimeter with free volumes

Let \(R_j\to\infty\), put \(\delta_j=1/R_j\), and
write \(A^j(u)=A(R_j u)\). The partition \(A^j\)
minimizes the \(G^{(m)}_{\delta_j}\) cut against all
compact replacements, by a change of variables.
For a partition \(F\), write \(J_j(F;B_r)\) for its
ordered disagreement integral on pairs with at least
one endpoint in \(B_r\), divided by \(\delta_j\).
Filling \(B_r\) with a single label leaves only crossing
pairs. The directional translation bound for the ball
then gives \(J_j(A^j;B_r)\le C_m\mathcal P(B_r)\).

The local difference-quotient argument gives
precompactness in \(L^1_{\mathrm{loc}}\).
In detail, mollify each indicator at scale \(\delta_j\)
using a smooth kernel supported in \(B_{\delta_j}\).
On \(B_r\), the \(L^1\) error is bounded by
\(C\delta_j J_j(A^j;B_r)\) and the mollified gradient
has \(L^1\) norm at most \(C J_j(A^j;B_r)\).
Both bounds follow by writing the increments inside
the convolution, using the zero integral of its
derivative and the lower bound for \(G^{(m)}_{\delta_j}\)
on \(B_{\delta_j}\). Local \(BV\) compactness and
diagonal extraction give an ordinary partition \(U\).
The directional lower bound is
\[
 \liminf_j J_j(A^j;B_r)\ \ge\
           \sqrt{2/\pi}\,\mathcal P(U;B_r),                    \tag{ED4}
\]
where the interface perimeter \(\mathcal P\) counts
each wall once. For each fixed displacement
\(\delta_j z\), use the characterization of
\(|z\cdot D1_{U_i}|(B_r)\) by smooth tests on \(B_r\).
Testing the difference quotients and passing to the
local \(L^1\) limit gives the directional lower bound.
Fatou and the Gaussian average then give (ED4),
using \(\frac12\sum_i|D1_{U_i}|\) for the interface
measure. This argument is the constant-weight case
of the intermediate-scale limit.
[Intermediate-scale perimeter convergence, P7](THREE_CELL_NEAR_ONE.md#intermediate-scales-and-localization-of-triple-ties)

For the minimizing property, take a locally \(BV\)
replacement \(F=U\) outside \(B_b\), and fix
\(b<a<c\). Glue \(F\) inside \(B_t\) to \(A^j\)
outside \(B_t\), for \(a<t<c\). The possible increase
of its cut over \(J_j(F;B_t)\) uses only crossing
pairs where the outer endpoint's old label differs
from \(F\). Averaging in \(t\), the length of the
eligible interval of radii is at most the pair
distance. The contribution with outer endpoint in
\(B_{c+1}\setminus B_a\) after the division by
\(\delta_j\) is thus bounded by a constant times
the \(L^1\) label mismatch there, which tends to
zero because \(F=U\) on this annulus. The remaining
contribution is a Gaussian tail, since the pair
distance there is at least one. Some choice of
\(t\) has error \(o(1)\), and minimality gives
\(\limsup_j J_j(A^j;B_a)\le\limsup_j J_j(F;B_c)\).

For the fixed \(F\), the directional translation
bound gives
\(\limsup_j J_j(F;B_c)\le\sqrt{2/\pi}\,
\mathcal P(F;\overline{B_c})\). Indeed the difference
at displacement \(\delta_j z\) counts endpoints in
\(B_c\) and is bounded using the directional
variation in \(B_{c+2\delta_j|z|}\).
For \(|z|\le1/\delta_j\) the bound is dominated
by \(|z|\) times the total variation in \(B_{c+2}\).
The remaining \(z\) contribute a Gaussian tail.
Let \(c\downarrow a\) through radii with zero sphere
perimeter and choose \(a>b\) with this same property.
Together with (ED4) this proves the local
perimeter comparison for unrestricted replacement
volumes. In particular the limit satisfies the
local perimeter minimizing-cluster inequality
with zero volume-error parameter.

## All phases meet the origin of a perimeter limit

For each fixed \(r>0\) and large \(j\), the ball
\(B_{rR_j/2}(x_i)\) is contained in \(B_{rR_j}(0)\).
Applying (ED3) there and passing to the local \(L^1\)
limit shows
\[
               |U_i\cap B_r(0)|\ge c_{m,k}(r/2)^m
                 \quad\text{for every }r>0.                   \tag{ED5}
\]
Thus every original label has infinite volume in \(U\)
and is present at the same origin at every radius.

In dimension one, locally finite perimeter makes
the boundary locally a finite set, so at most two
labels can be present at a single point. In dimension
two the limit qualifies as a locally isoperimetric
partition with compact volume-preserving competitors,
since it minimizes over all compact competitors.
Novaga–Paolini–Tortorelli's Theorem 4.2 applies to an
arbitrary finite number of phases and gives at most
three infinite-area phases.
[Novaga–Paolini–Tortorelli, Theorem 4.2](https://arxiv.org/html/2312.13709#S4)

In dimension three, the regularity theorem of Taylor
gives a local interface model of a plane, a triple
junction, or a tetrahedral junction at each boundary
point of a perimeter-minimizing cluster. In its
application to locally finite-perimeter clusters
by Leonardi–Maggi, Theorem 2.1 and (2.12), the
local comparison (1.1) is precisely the comparison
proved above, with \(\Lambda=0\) and any fixed
\(r_0>0\). Although Theorem 2.1 lists chamber volumes
in its constants, its Appendix A argument is local here.
For the current construction, truncate the cluster
outside a ball containing the support of the
Lipschitz deformation and its image, assigning
the exterior to one label. The truncated cluster
has finite total perimeter. The comparison in
(A.19)–(A.21) agrees with that cluster outside the
deformation support, and can be pasted to the
original cluster there. Every infinite chamber
retains positive volume after this compact change.
This replaces the role of the volume estimate at
(A.22) and admits the competitor in (A.23).
The perimeter comparison in (A.25) is local to the
deformation support, giving the Lipschitz set
comparison (2.10) with zero error. Taylor's models
therefore apply. A small coordinate ball in the
chart (2.12) has at most four complementary sectors,
each belonging to one label. Relation (ED5) bounds
the number of original phases by four.
[Taylor, 1976](https://annals.math.princeton.edu/1976/103-3/p08)
[Leonardi–Maggi, Theorem 2.1, (2.12) and Appendix A](https://arxiv.org/html/1505.06709#S2)

For a near-one microscopic Gaussian maximizer limit
at a full winning tie, (JD6) gives a free-volume
Euclidean compact minimizer with every label present.
The result proves that such a limit cannot depend
on at most two coordinates when \(k\ge4\), or on
at most three coordinates when \(k\ge5\). This
concerns the number of effective coordinates. It
leaves the higher-dimensional junction shapes
undetermined.
[Microscopic minimizers at full winning ties, JD6](GAUSSIAN_NEAR_ONE_JUNCTION_DEGREE_AND_COMPONENTS.md#microscopic-minimizing-partitions-at-the-ties)

## Conical consequences and the four-phase perimeter limit

For every \(m\), a Euclidean Gaussian-cut minimizer
with \(k\) nonempty phases produces a conical local
perimeter minimizer in \(\mathbb R^m\) with all \(k\)
phases. For \(m=3\) and \(k=4\), every large-scale
perimeter limit \(U\) of the Gaussian-cut minimizer is
itself a regular tetrahedral cone centered at the
dilation origin.

To prove these claims, dilate the local perimeter
minimizer \(U\) from (ED4)–(ED5) around the origin.
Each rescaling remains a local perimeter minimizer.
Filling a ball by one label bounds the perimeter in
that ball by the area of its sphere. The bound
can first be taken at radii carrying no perimeter
on their spheres and then at all radii by continuity
from above. Hence the rescalings are locally \(BV\)
precompact. The limit minimizes perimeter and the
perimeters converge on balls whose spheres have
zero limiting perimeter. For the upper bound and
the comparison, paste a fixed competitor for the
limit inside an intermediate sphere to the rescaled
partition outside it. The excess interface area on
that sphere is bounded by its trace-label mismatch.
Its average over an annulus avoiding the competitor's
change support tends to zero by local \(L^1\)
convergence and polar integration. Choosing good
pasting radii and using perimeter lower semicontinuity
gives the comparison. Using the limit itself as
competitor also gives perimeter convergence.

For \(m\ge2\) the rectifiable interface varifold
of \(U\), with each interface counted once, is
stationary by compact smooth deformations. Write
\(\mu\) for its interface measure and \(x^\perp\)
for the projection of \(x\) normal to its approximate
tangent plane. The stationary radial first variation gives
\[
 \frac{\mu(B_b)}{b^{m-1}}-\frac{\mu(B_a)}{a^{m-1}}
    =\int_{B_b\setminus B_a}
           \frac{|x^\perp|^2}{|x|^{m+1}}\,d\mu(x)               \tag{ED6}
\]
at radii with zero sphere measure. By the ball comparison
the normalized area is bounded, so it has limits at
zero and at infinity. Perimeter convergence under
rescaling shows that every tangent at the origin
and every blow-down of \(U\) have constant normalized
area, given by the respective endpoint limit in (ED6).
Equality in (ED6) makes the radial direction tangent
at almost every interface point. For every label this
implies \(x\cdot D1_{U_i}=0\) in the corresponding
constant-density partition. Differentiating test-function
pairings under radial dilation, this distributional
identity gives invariance of each label under all
positive dilations.
[Allard, stationary varifold monotonicity, §5.1](https://doi.org/10.2307/1970868)

Thus every tangent and blow-down of \(U\) is a conical
perimeter minimizer. Relation (ED5) is preserved by
rescalings and local \(L^1\) limits, so all \(k\)
phases survive. For \(m=1\) a tangent is directly
a constant label or a half-line partition, since
locally finite perimeter gives finitely many switches.
This proves the general conical consequence.

When \(m=3\) and \(k=4\), (ED5) and Taylor's local
classification force the origin tangent of \(U\)
to be the regular tetrahedral cone. A blow-down of
\(U\) is also such a cone, since its four phases meet
the origin and it is homogeneous. Both endpoint
limits of the normalized area in (ED6) thus equal
the tetrahedral area density. The normalized area
of \(U\) is constant and the same equality argument
makes \(U\) conical. Its tangent identifies \(U\)
as the regular tetrahedral cone.

## A uniform lower bound on independent score directions

Fix \(k\ge4\) and a physical ball \(B_M\) in the reduced
dimension \(d=k-1\). At a full winning tie \(x_\rho\in B_M\)
of an actual balanced maximizer, let
\[
 Q_{\rho,i}(u)=T_\rho1_{A_i}(x_\rho+eu),\qquad
 H_\rho=\sum_i\int_{B_1}\nabla Q_{\rho,i}(u)
                           \otimes\nabla Q_{\rho,i}(u)\,du,
 \qquad e=\sqrt{2(1-\rho)}.
\]
Order the eigenvalues of \(H_\rho\) from largest to smallest.
There is \(c_{k,M}>0\) such that for all sufficiently large \(\rho<1\),
\[
 \lambda_3(H_\rho)\ge c_{k,M}\quad(k=4),\qquad
 \lambda_4(H_\rho)\ge c_{k,M}\quad(k\ge5).                       \tag{ED7}
\]

Indeed, for any sequence of such ties as \(\rho\uparrow1\),
(JD6) gives strong local convergence of the rescaled indicators
to \(1_{U_i}\), where \(U\) is a compact-replacement Gaussian-cut
minimizer with all \(k\) labels. In noise-length variables the
conditional Gaussian has mean \(\rho u+(\rho-1)x_\rho/e\)
and variance \((1+\rho)/2\). The means converge uniformly on
bounded \(u\) to \(u\) and the variance converges to one.
Convolution and its first derivatives therefore converge
uniformly on \(B_1\), by local \(L^1\) convergence and uniform
Gaussian derivative tails. In particular \(H_\rho\) converges
to the same integral for \(Q_i=G_1*1_{U_i}\).

A vector in the kernel of the limiting matrix annihilates
all directional score derivatives on \(B_1\), by their continuity
and the nonnegative squares in the quadratic form. Analyticity
extends this annihilation to the whole space. Thus all the
smoothed scores are invariant under translations in the common
kernel. Gaussian convolution is injective on bounded functions:
its Fourier multiplier is smooth and nonzero, so it can be
divided out on compactly supported tests for the Fourier
transforms. Applying this to each translated indicator
difference gives invariance of each indicator, modulo null
sets, under translations in that kernel. Equivalently the
indicators depend only on the orthogonal projection to its
complement, as follows by mollifying the translation-invariant
indicators and passing to local \(L^1\). A failure of (ED7)
would give a limiting matrix of rank at most two for \(k=4\),
or at most three for \(k\ge5\), violating the phase bound.
