# Restrictions on a competing three-cell maximizer

Fix $0<\rho<1$. After dimension reduction to the active plane, every
global maximizer has three pairwise interfaces of positive surface measure.
A balanced maximizer with convex cells, or one invariant under dilations
about the origin, is the regular three-sector partition. Thus a nonsimplex
maximizer at the first continuation endpoint must escape both classes.
The last section gives a counterexample to a circular rearrangement
relaxation that drops disjointness while retaining all circular masses.

## The known global interval and its continuation endpoint

Heilman's July 2026 preprint proves the dimension-independent balanced
inequality for $0<\rho\le2/5$. The precise statement is Theorem A.1,
(134), for every measurable three-cell partition in $d\ge2$ of masses
$1/3$. The proof uses the earlier Theorem 1.5 of arXiv:2306.03312v3 for
$\rho<1/10$, and on $[1/10,2/5]$ it uses planar dimension reduction,
aligned circular arcs and the strict mean-zero norm certificate of
Lemma A.4. The latest proof therefore reaches beyond the older explicit
interval. [Heilman, *Sharp Hardness for MAX-3-CUT and Quantum MAX-CUT*, arXiv:2608.00333v1](https://arxiv.org/pdf/2608.00333v1).

At $2/5$ the proof bounds the excess over the simplex by
$(h_+-2)\sum_i\langle p_i-1/3,M_{J_+}(p_i-1/3)\rangle$, in the radial
measure of Appendix A. The weight $J_+$ is positive at every positive
radius and $h_+<2$. Thus equality forces $p_i(r)=1/3$ almost everywhere.
The circular comparison below gives equality rigidity on this class.
It follows that the first transition correlation, if below one, is
strictly greater than $2/5$: the
[three-cell continuation implication](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)
would require a nonsimplex tied maximizer at a first transition.

The preprint's Theorem 1.3 already gives the corresponding discrete
interval. In the precise positive-correlation formulation, fix $q=3$,
$\rho\in(0,2/5]$, and $\epsilon>0$. There is a $\tau>0$ such that for
any $m$ and balanced $f:\{1,2,3\}^m\to\{1,2,3\}$, if all influences of
the three indicator coordinates of $f$ are at most $\tau$, then
$S_\rho(f)\le\lim_{n\to\infty}S_\rho(\mathrm{PLUR}_{n,3})+\epsilon$.
The noise retains each vote with probability $\rho$ and otherwise
resamples it uniformly and independently. The invariance implication
uses Isaksson–Mossel Theorems 1.10 and 7.1 and their plurality limit in
Lemma 7.3. Extension of this discrete conclusion to every positive
correlation requires the corresponding Gaussian bound. [Isaksson–Mossel, arXiv:0903.3362](https://arxiv.org/pdf/0903.3362).

## Translation equations on an unconstrained space of speeds

Let $A$ be a three-cell maximizer in $\mathbb R^2$, with positive prescribed
masses. Its cylinder in $\mathbb R^3$ is a maximizer by Heilman–Tarter's
dimension independence theorem. We use their Lemmas 2.4, 3.2, 4.8 and
5.4 for the reduced interfaces, the stationary scores, the second
variation and the normal translation identity. Theorem 1.9 supplies the
reduction for actual maximizers from any dimension $d\ge2$. The results
require global maximality and fixed positive masses at positive
correlation. [Heilman–Tarter, Forum of Mathematics, Sigma 9 (2021), e65](https://doi.org/10.1017/fms.2021.56).

Orient each existing interface $e=(i,j)$ from cell $i$ into cell $j$.
Write $N_e$ for the resulting unit normal, $d_e=e_i-e_j$ for its label
vector, $d\nu_e=\phi_2\,d\mathcal H^1|_{\Sigma_e}$, and

\[
 U_i=T_\rho 1_{A_i},\qquad
 a_e=|\nabla(U_i-U_j)|,
 \qquad K(e,x;e',y)=(d_e\cdot d_{e'})M_\rho(x,y).
\]

The gradient satisfies $\nabla(U_i-U_j)=-a_eN_e$ on each regular
interface and $a_e>0$ almost everywhere. Integrals in the next formulas
are over the disjoint union of oriented interfaces with measure $\nu$.
For a compactly supported smooth speed $g$ on the regular interiors,
test the cylinder with speed $g(x)\zeta(z)$, where $\zeta$ is smooth, odd
and compactly supported in the extra coordinate. This has zero first
variation of every mass. It is realized locally by smooth vector fields.
Smooth mass corrections at second order exist by the implicit function
theorem using independent bumps on a spanning tree of the interface
graph. That graph is connected, since a union of cells having no outgoing
interface would have distributional gradient zero and thus constant
indicator. The stationary score terms cancel in the exactly constrained
second variation. The cylinder's stability gives

\[
 \|\zeta\|_{L^2(\gamma_1)}^2\int a g^2d\nu
   \ \ge\ \langle\zeta,T_\rho\zeta\rangle_{\gamma_1}
                 \iint K(x,y)g(x)g(y)d\nu(x)d\nu(y).
\]

Since odd cutoffs of $z$ tend to the linear OU eigenfunction in $L^2$,
this proves

\[
 Q^{(1)}(g)=\int a g^2d\nu-\rho\iint K(x,y)g(x)g(y)d\nu(x)d\nu(y)\ge0.
 \tag{G1}
\]

There is no mass constraint on $g$. This is the first Hermite block in
the extra coordinate. For translations $s_v(x)=v\cdot N_e(x)$,

\[
 \rho K s_v=a s_v.
 \tag{G2}
\]

Here and below the operator associated to the kernel integrates against
$\nu$. The identity follows by Gaussian integration by parts in the
convolution of $1_{A_i}-1_{A_j}$ and the gradient formula above.

For the use of bounded, possibly noncompact speeds, integrability can
be deduced from the same variation. Let $L=\Delta-x\cdot\nabla$ and
$0\le\chi_R\uparrow1$ be smooth radial cutoffs. Gaussian divergence on
the cells, with the pairwise stationary normal sign, gives

\[
 \int\chi_R a\,d\nu=
 \sum_i\int_{A_i}[-\chi_R LU_i-\nabla\chi_R\cdot\nabla U_i]d\gamma_2.
\]

The first two derivatives of $U_i$ are bounded at fixed $\rho$. The
cutoff errors vanish and the remaining integrands have integrable
linear growth. Thus

\[
 \int a\,d\nu=\sum_i\int_{A_i}(-L)U_i\,d\gamma_2=\rho S'_\rho(A)<\infty.
 \tag{G3}
\]

The last step uses the OU spectral identity $-LT_\rho=\rho\partial_\rho
T_\rho$. Now take compact nonnegative speeds tending to one on a single
edge in (G1). Fatou gives
$2\rho\iint_{\Sigma_e\times\Sigma_e}M_\rho\,d\nu_e d\nu_e\le\int_{\Sigma_e}a_e\,d\nu_e$.
Mixed Mehler integrals are finite by Cauchy–Schwarz for its positive
quadratic form on compact surface measures followed by exhaustion.
Consequently $|K|$ is integrable on the product space, and (G1) extends
to bounded measurable speeds by bounded smooth approximation on the
regular patches and dominated convergence. These estimates also justify
the integral in (G2) almost everywhere. To obtain the translation
formula there, first apply the convolution divergence identity with a
compact cutoff in $y$, then let the cutoff tend to one. The volume error
vanishes by the Gaussian decay of the conditional density, and the
boundary integrals converge absolutely by the Mehler integrability.

One consequence is a matrix inequality valid at every maximizer. For
any measurable subset $D$ of its oriented interface union and any $v$,
apply (G1) to $(1-2\mathbf1_D)s_v$ and subtract the zero form on $s_v$.
After canceling the local term the result is

\[
 \int_{D\times D^c} K(x,y)(v\cdot N_e(x))(v\cdot N_{e'}(y))
                         d\nu(x)d\nu(y)\ge0.
 \tag{G4}
\]

Equivalently the symmetric matrix integral of
$K(x,y)(N_e(x)N_{e'}(y)^T+N_{e'}(y)N_e(x)^T)/2$ is positive semidefinite
for every $D$. The sign of an individual integrand may vary. A second
consequence concerns a tie $S_r(A)=S_r(T)$ at the first loss of global
optimality. Then $S_\rho(T)-S_\rho(A)\ge0$ for $0<\rho<r$ implies
$S'_r(A)\ge S'_r(T)$. By (G3), the tied partition must have at least the
simplex's total score-weighted interface measure at $r$.

There is also an exact sign requirement for dilation at the actual
correlation. Put $b(x)=x\cdot N_e(x)$ and $g_e=\partial_\rho(U_i-U_j)$
on $e=(i,j)$. Differentiating the conditional convolution when the cells
are dilated gives for the interface kernel

\[
 Kb=-x\cdot\nabla(U_i-U_j)-\frac{1-\rho^2}{\rho^2}\Delta(U_i-U_j)
    =\frac{ab}{\rho^2}+\frac{1-\rho^2}{\rho}g.
\]

In the first equality the conditional mean and standard deviation are
divided by the dilation factor. In the second use $-L U=\rho\partial_\rho U$.
This is the dilation identity in Heilman, Lemma 5.5. For its integration
in (G1), polynomial-weighted cutoff divergence as in (G3) gives
$\int a(1+|x|^m)d\nu<\infty$ for every fixed $m\ge0$. The same
single-edge argument with speeds bounded by $1+|x|^m$, followed by
Mehler Cauchy–Schwarz, gives the needed absolute product moments.
Derivatives in $\rho$ of the convolution are polynomially bounded at
fixed positive correlation. We obtain

\[
      \int b g\,d\nu\le -\frac{1}{\rho(1+\rho)}\int a b^2\,d\nu.
 \tag{G6}
\]

In particular, a proof of the opposite nonnegative sign at a first tie
would force conicality. The global value slope from (G3) controls
$\int a$, whereas (G6) involves the correlation derivative of each
score on its own interface. No implication between the two signs is
available here. [Heilman, Lemma 5.5](https://arxiv.org/html/2209.11216v1).

## Complete adjacency

If a pair of cells had no common interface of positive measure, the
connected graph would have exactly two edges. Orient both outward from
the common cell. Then $d_e\cdot d_{e'}$ is either 2 or 1, and $K$ is
strictly positive. By (G1) and (G2), $s_v$ cannot change sign on sets of
positive surface measure. Otherwise its absolute value would give a
strictly smaller form than $Q^{(1)}(s_v)=0$, since the local terms agree
and $K>0$ on the pairs with unlike signs.

It follows that all the oriented normals equal a fixed unit vector almost
everywhere. In fact two distinct unit normals admit a vector whose dot
products with them have opposite strict signs. Apply this to two points
of the essential support of the distribution of normals, using a finite
positive measure equivalent to $\nu$ if needed. The indicator of the
common cell has zero distributional derivative in the tangent direction
and a nonpositive distributional derivative along that fixed normal.
Hence it depends only on the normal coordinate and is a nonincreasing
binary function of it. It is a halfspace, since its mass lies in $(0,1)$.
The other two cells both lie in its connected open complementary
halfspace. Their missing interface would force each indicator to have
zero distributional gradient inside that halfspace. This contradicts
their two positive masses.

The local Gaussian perimeter analogue of positive adjacency appears in
Milman–Neeman, Theorem 12.4. Its hypotheses and operator concern the
Gaussian perimeter, whereas the argument here uses the nonlocal
cylinder form (G1). [Milman–Neeman, Theorem 12.4](https://arxiv.org/html/1805.10961v3).

## Convex and conical maximizers

**Convex cells.** Suppose each cell has a convex representative in the
plane. For each pair choose the separating line containing their
positive-length shared interface. Each cell is the intersection of its
two pairwise closed halfplanes, modulo null sets. Indeed the intersections
contain the cells, have pairwise disjoint interiors, and the cells cover
the plane. Use affine functions $L_{12},L_{23},L_{31}$ positive on the
side of their first index. Their normals span the plane. If they were
parallel, the cells would be intervals in a single coordinate and the
outer intervals could not be adjacent. The image of the plane under
$(L_{12},L_{23},L_{31})$ is an affine plane in $\mathbb R^3$ missing both
open orthants with all coordinates positive or all negative. Such points
would have no winning cell. Its defining normal has all coordinates of
one weak sign, since mixed signs would let it meet both orthants. Its
offset must be zero, since a nonzero offset would let it meet one of the
two orthants. Thus all dividing lines pass through a common point $w$.

To check centering even at a possible T-junction, take a line containing
an interface ray of cells $i,j$. Their stationary score difference is
constant on the whole line by real analyticity of Gaussian convolution.
Let $p$ be the conditional probability of the $i$ halfplane when the
conditional mean is $\rho x$ and the observation $x$ runs along the line.
This constant lies in $(0,1)$ and equals $1/2$ exactly when the line
passes through zero. Along the shared ray, the other boundary constraints
recede to infinity, so the score difference tends to $2p-1$. When the
three dividing lines are distinct, each cell angle is strictly below
$\pi$ and the reverse ray is in the interior of the third wedge. The
opposite limit is zero. Therefore every line passes through zero and
$w=0$. If two dividing lines coincide, one cell is a halfplane, with the
other two splitting its complement at $w$. Choose $i$ to be the
halfplane cell and $j$ the cell adjacent on one ray of its boundary.
On the reverse ray the other cell fills the opposite side, and the
limit of $U_i-U_j$ is $p$. Constancy would give $2p-1=p$, impossible.
These limits follow by dominated convergence after writing each
evaluation as Gaussian convolution of shifted halfplane constraints.
The analytic continuation and opposite-ray limit for nondegenerate
translated simplicial cones agree with Heilman–Mossel–Neeman, Lemma 4.1
and Corollary 5.1. The argument here also handles the degenerate line
configuration allowed by convexity. [Heilman–Mossel–Neeman, arXiv:1403.0885](https://arxiv.org/abs/1403.0885).

The equal Gaussian masses now force all three wedge angles to be
$2\pi/3$.

**Dilation invariance about zero.** Suppose instead that the labels
depend only on angle, up to null sets. An actual maximizer assigns labels
by the largest $U_i$ minus a constant price for each label, as in the
[score-assignment argument](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality).
No pairwise score difference minus prices is identically zero. This
follows from injectivity of Gaussian convolution and the disjoint cells
of positive mass. For each pair choose a direction where this analytic
function of radius is not identically zero. Outside its discrete zero
set, the pairwise angular score difference on the corresponding circle
is not identically zero and has finitely many zeros. Choose a radius
good for all pairs and for the almost-everywhere angular assignment.
Then the conical partition has finitely many rays between open sectors
of constant label.

At a ray separating $i,j$, rotate its supporting line to the horizontal
axis. The score difference is constant on this line by analyticity. At
the origin it equals zero, since isotropic Gaussian scaling preserves
the masses of cones. Write $f=1_{A_i}-1_{A_j}$ and
$\sigma=\sqrt{1-\rho^2}$. Gaussian convolution along the line is zero,
so one-dimensional Gaussian deconvolution yields

\[
 \int_{\mathbb R} f(x,y)e^{-y^2/(2\sigma^2)}dy=0
                              \quad\text{for almost every }x.
\]

Deconvolution is valid for bounded functions by injectivity of
convolution with a Gaussian on tempered distributions. Next fix either
sign of $x$ and substitute $y=|x|t$. Homogeneity and continuity in the
positive Laplace parameter give

\[
 \int_0^\infty [f(\pm1,t)+f(\pm1,-t)]e^{-s t^2}dt=0
                                      \quad\text{for all }s>0.
\]

With the variable $t^2$, uniqueness of the Laplace transform implies the
even part is zero almost everywhere. For example, multiply by one fixed
exponential to get a finite measure and apply uniqueness of Fourier
transforms after analytic continuation in $s$. Thus reflection in the
line interchanges labels $i,j$ and fixes the third throughout the plane.

The sectors on either side of each ray have equal angular widths because
this reflection preserves the set of rays. Hence all widths are equal.
The labels on either side of a sector cannot coincide. If there were a
run $i,j,i$, reflection at the $j/i$ interface would force the next
sector to be $j$, and iteration in both directions would use only two
labels. Therefore the labels repeat in the order $1,2,3$ up to
permutation. If there are at least two repetitions, a nonidentity planar
rotation fixes each cell and every cell first moment vanishes. Every
translation is then balanced to first order. For some $v$, $s_v$ is
nonzero on an interface, so the ordinary translation form
$\int a s_v^2-\iint K s_v s_v=(1-1/\rho)\int a s_v^2$ is strictly
negative. Realization by cutoffs and balancing corrections follows the
translation argument in Heilman–Tarter, Lemma 5.5 and Section 6. This
contradicts planar stability. There is exactly one repetition.

The dilation-invariant classification in Heilman's hyperstability
argument uses additional hypotheses, including a mixed correlation
variation. The proof above uses positive-correlation maximality once
conicality is given. Conicality of an arbitrary maximizer remains the
missing step. [Heilman, arXiv:2209.11216v1, Definition 1.3, Theorem 1.5 and Lemmas 8.2–8.3](https://arxiv.org/html/2209.11216v1).

## Circular balance and compatibility of first angular moments

For a planar partition let $p_i(r)$ be the fraction of the circle of
radius $r$ in label $i$, and let $u_i(r)$ be its first angular vector
moment, with angular probability measure $d\theta/(2\pi)$.
Write $b(p)=\sin(\pi p)/\pi$. Maximizing the projection of the first
angular moment in its own direction gives $|u_i(r)|\le b(p_i(r))$ by
choosing the largest cosine values. Equality for $0<p_i(r)<1$ forces
a centered arc in that direction. The partition constraint also gives
$\sum_i u_i(r)=0$.

The conditional angular kernel, without its radial factor, is
$\exp(z\cos(\theta-\varphi))$, $z=\rho rs/(1-\rho^2)>0$. Split it into
$z e^{-z}\cos(\theta-\varphi)$ and the remainder with kernel function
$t\mapsto e^{zt}-ze^{-z}t$, nonnegative and nondecreasing on $[-1,1]$.
Apply the circular rearrangement theorem to the remainder. Its inputs
are the two slice indicators and a nonnegative kernel symmetric and
decreasing in circular distance. It upper bounds each same-label
angular pairing by that of aligned arcs with the corresponding widths.
[Friedberg–Luttinger, *Rearrangement Inequality on a Circular Domain*, Theorem I](https://digital.library.unt.edu/ark:/67531/metadc1025197/m2/1/high_res_d/4446751.pdf).

In particular for $p_i(r)=1/3$ at almost every radius and all labels,
put $b_*=b(1/3)$. The exact lower bound from the separated first mode is

\[
 S_\rho(T)-S_\rho(A)\ge\iint_0^\infty
 \frac{rs}{1-\rho^2}e^{-(r^2+s^2)/(2(1-\rho^2))}z e^{-z}
     \sum_i[b_*^2-u_i(r)\cdot u_i(s)]\,dr\,ds.
 \tag{G7}
\]

Each summand is half the sum of the two squared-length deficits and
$|u_i(r)-u_i(s)|^2$. It vanishes almost everywhere only if all the
$u_i$ are constant in radius with norm $b_*$. The slices are then
fixed arcs and disjointness forces the regular three-sector partition.
This proves strict global comparison in the circularly balanced class
for every positive correlation and supplies the equality step at $2/5$.

For varying masses define $E(p,q)$ as the maximum possible aligned
first moment sum $\sum_i b(p_i)b(q_i)$ minus the maximum of
$\sum_i u_i\cdot v_i$ subject to $\sum u_i=\sum v_i=0$,
$|u_i|\le b(p_i)$ and $|v_i|\le b(q_i)$ in $\mathbb R^2$. This compact
optimization gives a nonnegative deficit. The same kernel split bounds
each actual pair of slices by the aligned-arc envelope minus
$ze^{-z} E(p(r),p(s))$. For interior $p,q$, equality $E=0$ forces each
moment to saturate its length and $u_i,v_i$ to point in the same
direction. The closed triples are nondegenerate triangles, since the
numbers $\sin(\pi p_i)$ are the side lengths up to a common factor
of a triangle with angles $\pi p_i$. The common directions imply that
the side lengths are proportional, hence $p=q$ by uniqueness of the
ordered angle profile. Conversely $p=q$ permits a closing saturated
triple with $u=v$. This gives a strict radial coupling term for
different interior profiles. Its dominance over the relaxed excess
near one remains open.

## Independent circular arcs can exceed the simplex value

Consider the relaxation in which each label on each circle is a centered
arc, with centers $0,2\pi/3,4\pi/3$. Impose that the three widths sum to
$2\pi$ at every radius and that each label has mass $1/3$ overall.
Allow overlaps and gaps. This is the class used when each label is
rearranged independently and only its circular measure is retained.
Heilman's proof of Theorem A.1 explicitly first passes to independently
aligned arcs and then estimates that envelope for $\rho\le2/5$.
The following counterexample precludes extending that particular
relaxed comparison to all correlations. [Heilman, arXiv:2608.00333v1, proof of Theorem A.1](https://arxiv.org/pdf/2608.00333v1).

Write the relaxed cells as $A_i(t)$ with half-width
$\pi/3+t c_i h(r)/r$, where $c=(1,-1,0)$,
$h\in C_c^\infty(0,\infty)$ and $\int h\,d\nu=0$. Here $\nu$ is the
three-ray measure in (S14) of the [spectral note](SIMPLEX_INTERFACE_SYMMETRY.md).
The pointwise width sum and all three Gaussian masses are exactly
preserved. The outward normal speed of each endpoint ray of cell $i$
is $c_i h(r)$. For each $r$, the old endpoint scores agree by symmetry.
Their outward derivatives are $-a_\rho(r)/2$ by reflection between the
two adjoining labels. Expand the exact stability difference into its
linear score pairing and transfer interaction. In polar coordinates
$r$ is fixed during a boundary change. The linear score terms of order
$t$ cancel since $\sum_i c_i=0$, and the order-$t^2$ coefficient from
the two moving endpoints of cell $i$ is $-c_i^2\int a_\rho h^2d\nu$.
Each cell's transfer measure has two rays with the same outward speed.
The interaction coefficient is therefore
$c_i^2\iint(2M_1+2M_{-1/2})hh\,d\nu d\nu$. Compact support away from
zero makes the expansions valid by dominated convergence. Thus the
total order-$t^2$ coefficient is

\[
 \|c\|^2\left[\iint M_{-1/2}(r,s)h(r)h(s)d\nu(r)d\nu(s)
                   -\rho\mathfrak q_{\rho,\mathrm{standard}}(h)\right].
 \tag{G5}
\]

This differs from the actual standard-sector quadratic form by one
cross-ray kernel. To see the effect as $\rho\uparrow1$, first take
$h(r)=\exp(-r^2/2)-1/\sqrt2$, which is bounded, Lipschitz, and of zero
$\nu$-mean, with $h(0)\ne0$. By (S17), the physical bracket
$\rho\mathfrak q_{\rho,\mathrm{standard}}(h)$ equals the $K_+$
squared-difference energy minus $(1-\rho)\int(a_\rho/\rho)h^2d\nu$.
The latter term is $O(\sqrt{1-\rho^2})$. The same-ray part of the
squared-difference energy has this order as well. Indeed the Lipschitz
bound and extension of the nonnegative integral to the full line give
an upper bound by a constant times
$\sigma^{-1}\mathbb E(X-Y)^2=2(1-\rho)/\sigma$ for correlated scalar
standard Gaussians. In the different-ray part, scale $r=\sigma u$ and
$s=\sigma v$, where $\sigma=\sqrt{1-\rho^2}$. The joint surface density
becomes $(2\pi)^{-2}\exp[-(u^2+v^2+\rho uv)/2]du\,dv$, so the
squared difference is $O(\sigma^2)$ after integration. But the first
term in (G5), with no squared difference, tends to
$h(0)^2/(6\pi\sqrt3)>0$ by the same rescaling and the Gaussian quadrant
formula. Dominated convergence applies on the positive quadrant.

For each sufficiently large fixed $\rho<1$, approximate this $h$ in
$L^2(\nu)$ by smooth functions compactly supported in $(0,\infty)$,
subtracting an appropriate multiple of a fixed compact bump to retain
zero mean. The quadratic forms at that $\rho$ are bounded, so the
strictly positive coefficient in (G5) survives. Take small nonzero $t$.
Thus the relaxed class strictly exceeds the simplex value for all
$\rho$ sufficiently near one. A global comparison via circular
rearrangement must keep track of the compatibility between the labels.

These are written analytic arguments with internal review. Complete
adjacency, failure of the relaxed comparison, and rigidity for convex or
conical maximizers still allow nonsimplicial maximizers with curved
interfaces and all three pairs adjacent.
