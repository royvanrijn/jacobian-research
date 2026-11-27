# Restrictions on a competing three-cell maximizer

Fix $0<\rho<1$. Every global maximizer with three positive prescribed
masses has a two-dimensional active plane, full moment rank and three
pairwise interfaces of positive surface measure.
A balanced maximizer with convex cells, or one invariant under dilations
about the origin, is the regular three-sector partition. Thus a nonsimplex
maximizer at the first continuation endpoint must escape both classes.
The last section gives a counterexample to a circular rearrangement
relaxation that drops disjointness while retaining all circular masses.

## The known global intervals and their continuation endpoints

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
The [direct interval extension](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#proposition-dg4-global-comparison-through-rho-equal-to-2150)
proves the same inequality through $21/50$. Its strict radial-mass
remainder (DG20), with $J_+(t)>0$ for $t>0$, forces circular balance
throughout $[2/5,21/50]$ at equality. Dimension reduction applies to
each equal-value optimizer and the circular comparison gives
simplex rigidity on the active plane. Thus the first transition
correlation is **strictly greater than $21/50$**: the
[three-cell continuation implication](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)
would require a nonsimplex tied maximizer at any first transition.

The [exact comparison near one](THREE_CELL_NEAR_ONE_EXACT.md)
gives a second interval of global optimality with equality rigidity,
ending at $1$. If strict suboptimality occurs in between, each
connected interval of strict suboptimality has a nonsimplex tied
maximizer at both endpoints. The same near-one comparison yields the
balanced low-influence plurality inequality on that interval.

The preprint's Theorem 1.3 gives the discrete interval through $2/5$.
The same implication with the extended Gaussian bound now gives the
following statement for $\rho\in(0,21/50]$. Fix $\epsilon>0$. There is
a $\tau>0$ such that for
any $m$ and balanced $f:\{1,2,3\}^m\to\{1,2,3\}$, if all influences of
the three indicator coordinates of $f$ are at most $\tau$, then
$S_\rho(f)\le\lim_{n\to\infty}S_\rho(\mathrm{PLUR}_{n,3})+\epsilon$.
The noise retains each vote with probability $\rho$ and otherwise
resamples it uniformly and independently. Theorem 7.1 of
Isaksson–Mossel supplies an ordinary Gaussian label function in
dimension $2m$ with the same means and noise stability within
$\epsilon$. Its Lemma 7.3 identifies the plurality limit with the
planar simplex value. As in their proof of Theorem 1.10, the hypothesis
can use influences truncated at a suitable degree instead.
Extension to every positive correlation requires the corresponding
Gaussian bound. [Isaksson–Mossel, arXiv:0903.3362](https://arxiv.org/pdf/0903.3362).

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

## A first-crossing obstruction at the level of Hermite weights

At a tangential first tie, allowing the competitor to move with
correlation gives a second-order condition. Let $A$ tie the simplex $T$
at the first endpoint $r$, put $D(\rho)=S_\rho(T)-S_\rho(A)$, and suppose
$D'(r)=0$. For a smooth compactly supported balanced speed $h$ on the
regular interfaces of $A$, let
$Q_0(h)=\int a h^2\,d\nu-\langle h,Kh\rangle$ at $r$ and use the
interface function $g_e=\partial_\rho(U_i-U_j)|_r$ from (G6). Then

\[
                  \left|\int gh\,d\nu\right|^2
                      \le \tfrac12 D''(r)Q_0(h).
 \tag{G10}
\]

To see this, realize $h$ by an exactly balanced smooth deformation
$A_t$ using two mass-correcting bumps as for (G1). Its fixed-correlation
first variation vanishes at $r$ and its second derivative is
$-2Q_0(h)$. The mixed derivative is $2\int gh\,d\nu$ by differentiating
the first variation in $\rho$ at fixed $A$. Differentiation is valid
on the compactly supported flows by the smooth Gaussian joint density.
For each fixed $c\in\mathbb R$ apply global comparison to $A_{c\epsilon}$
at $r-\epsilon$. Its deficit, divided by $\epsilon^2$, tends to
$D''(r)/2+2c\int gh\,d\nu+c^2Q_0(h)$. Nonnegativity for every $c$
proves (G10), including when $Q_0(h)=0$.

Suppose also that $\int_{A_i}|x|^2d\gamma_2=2/3$ for all $i$.
Then the dilation speed $b=x\cdot N_e$ is balanced by Gaussian divergence.
Cut it off and subtract fixed compact mass bumps, whose coefficients
tend to zero. The polynomial-weighted bounds used in (G6) justify
passing to the limit in (G10). Put $A_b=\int a b^2d\nu$ and $q=Q_0(b)$.
The dilation identity yields
$\int gb\,d\nu=-A_b/r-rq/(1-r^2)$. For $A_b>0$, (G10) forces $q>0$
and gives

\[
 \qquad D''(r)\ge
   2\big(A_b/r+rq/(1-r^2)\big)^2/q
   \ \ge\ \frac{8 A_b}{1-r^2}.
 \tag{G11}
\]

In particular $D'(r)=D''(r)=0$ together with equal second radial
moments forces $b=0$ on the interfaces. This gives dilation
invariance by distributional differentiation along radial rays, and
hence the simplex by the conical rigidity below. The extra radial
moment hypothesis and the slope equality remain open at a general tie.

An exact condition also holds at every tie with the simplex, including
a transverse first tie. Write $f_i=1_{A_i}-1_{T_i}$. For any measurable
$E$ with $\int_E f_i\,d\gamma_2=0$ for every $i$, put $f_E=1_E f$.
Patching the labels from $A$ into $T$ on $E$ or on $E^c$ produces two
balanced ordinary partitions. Their summed deficit at any $\rho$ is

\[
 S_\rho(T)-S_\rho(A)
       +2\sum_i\langle (f_E)_i,T_\rho(f_{E^c})_i\rangle .
 \tag{G12}
\]

This follows by expanding both quadratic energies around $T$.
Consequently the cross term is nonnegative at a tied global optimum,
with equality exactly when both patches are optimal. The inequality
holds for every rotation and relabeling of $T$. If the cut is saturated
at a first tie, its cross term has derivative at most
$(S'_r(A)-S'_r(T))/2$, by differentiating (G12) from the left.

### Limits of a coefficient-only test

Fix any proposed crossing $s\in(0,1)$. Nonnegative Hermite masses,
the exact zeroth coefficient, Parseval's sum, the smaller first moment
of the competitor and a strictly larger leading perimeter-angle
coefficient do not by themselves preclude a transverse first crossing.
The following example uses formal Hermite masses and makes no
assertion that they belong to a measurable partition.

Write $W_m=W_m(T)$ for the actual nonnegative Hermite masses of the
simplex, so $S_\rho(T)=\sum_m W_m\rho^m$, $W_0=1/3$ and
$\sum_m W_m=1$. Both $W_1$ and $W_3$ are positive. For the latter,
integrate the degree-three Hermite polynomial $x_1^3-3x_1$ over the
wedge centered on $e_1$. With $R$ the planar Gaussian radius, the
radial relation $\mathbb E R^3=3\mathbb E R$ and angular integrals
$\int_{-\pi/3}^{\pi/3}\cos^3\theta\,d\theta=3\sqrt3/4$ and
$\int_{-\pi/3}^{\pi/3}\cos\theta\,d\theta=\sqrt3$ give a nonzero
inner product.

First decrease $W_1$ by $\epsilon s$, increase $W_2$ by
$\epsilon(1+s)$, and decrease $W_3$ by $\epsilon$, keeping every
other mass fixed. For small $\epsilon>0$ this yields a nonnegative
sequence with the same sum and formal deficit
$\epsilon\rho(1-\rho)(s-\rho)$. To impose strict positivity near one,
put

\[
 H(\rho)=\sum_{m\ge2}W_m(\rho-\rho^m)
        =1-S_\rho(T)-(1-W_0)(1-\rho),
 \qquad J(\rho)=\frac{H(\rho)}{\rho(1-\rho)}.
\]

Add to the deficit the term $\eta[H(\rho)-J(s)\rho(1-\rho)]$.
Its coefficients add nonnegative amounts $\eta W_m$ to the proposed
competitor masses for $m\ge3$. The coefficient adjustments in degrees
one and two preserve nonnegativity and leave $W_1$ strictly smaller
if $\eta>0$ is small enough relative to $\epsilon$. The sum and
zeroth mass stay fixed. For $0\le\rho\le s$ the function
$J(\rho)=\sum_{m\ge2}W_m(1+\rho+\cdots+\rho^{m-2})$ has bounded
derivative, by comparison with the power series of a geometric sum.
Choose $\eta\sup_{[0,s]}|J'|<\epsilon$. The complete deficit

\[
 \rho(1-\rho)\,[\epsilon(s-\rho)+\eta(J(\rho)-J(s))]
 \tag{G9}
\]

is positive for $0<\rho<s$, zero at $s$ and has strictly negative
derivative there. At the other endpoint its leading coefficient in
$\arccos\rho$ is $3\eta/(2\pi)>0$ by the
[simplex expansion](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex).
These formal masses therefore produce another crossing between $s$
and one. A first-crossing proof must use geometric information on
realizable partitions or stronger relations among their Hermite masses.

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

## One active dimension is impossible

Every three-cell global maximizer with positive masses in $d\ge2$
uses two active directions, and its Gaussian cell moments span that
plane. Suppose first that the cells depend on just one coordinate.
The score differences minus prices are analytic on the line with
discrete zeros, by Gaussian convolution injectivity and the distinct
indicators, so the transitions are locally discrete. Choose a
transition at zero if one exists. Otherwise choose the smallest
positive transition $x$, reflecting the coordinate if necessary.
Let the labels immediately to its left and right be $i,j$. Testing
(G1) at a single-transition speed requires

\[
 (T_\rho(1_{A_j}-1_{A_i}))'(x)\ge2\rho\phi_\sigma(x-\rho x),
 \qquad \sigma=\sqrt{1-\rho^2},
 \tag{G8}
\]

where $\phi_\sigma$ is the standard-deviation-$\sigma$ scalar normal
density. The diagonal interaction coefficient is 2 times the joint
diagonal density divided by the marginal density, giving the right
side after the Hermite factor $\rho$.
Compare $g=1_{A_j}-1_{A_i}$ to the step $s(y)=\operatorname{sgn}(y-x)$,
whose smoothed derivative is exactly the right side. The difference is

\[
 \frac{\rho}{\sigma^2}\int(g(y)-s(y))(y-\rho x)\phi_\sigma(y-\rho x)dy.
\]

On $y>x$ the first factor is nonpositive and the second positive.
On $y<\rho x$ their signs are reversed. In the possible interval
$(\rho x,x)$ the left label is $i$ throughout, so the first factor
vanishes. The third label has positive measure outside that interval,
giving strict negativity. This contradicts (G8).

If the centered cell moments spanned a smaller space than the active
plane, translation in a nonzero direction orthogonal to all moments
would have zero first-order mass change by Gaussian divergence.
Heilman–Tarter's translation test then forces all normals to be
orthogonal to that direction, since otherwise the constrained second
variation is strictly positive. Their cutoff and volume correction
argument applies as in their proof of Theorem 1.9. Distributional
differentiation would make the cells invariant in that direction,
which was just excluded.

## Each cell uses both active directions

In the active plane of a three-cell global maximizer with positive
prescribed masses, no individual cell is invariant under a nonzero
translation direction. In particular, none of the cells is a halfplane,
a parallel slab, or a union of parallel slabs.

Suppose $A_3$ were invariant along a unit vector $v$. The normal
translation speed $s_v$ vanishes on interfaces $31$ and $32$, and is
nonzero by full active dimension. On its support in $\Sigma_{12}$ the
interaction kernel is $2M_\rho>0$. Thus (G1)–(G2) applied to $|s_v|$
force $s_v$ to have only one sign. Orient $12$ and if necessary reverse
$v$ so that $D_v1_{A_1}$ is a nonnegative nonzero measure. Then
$D_v1_{A_2}=-D_v1_{A_1}$ and $D_v1_{A_3}=0$.
Convolution with a strictly positive Gaussian gives
$\partial_v U_1>0$, $\partial_v U_2<0$ and $\partial_v U_3=0$
throughout the plane. This follows by convergent Gaussian integration
by parts. Equivalently, slice along $v$ and use that the first two
indicators are oppositely monotone with a nonzero measure of
nonconstant slices.

The transverse indicator of $A_3$ has locally finite variation by
the interface regularity of actual maximizers. Its reduced boundary
consists of parallel lines, locally finite in the transverse coordinate.
By complete adjacency there is a positive-length interface of $A_3$
with at least one other label on one such line. However, the relevant
stationary pair score is strictly monotone on the line, whereas
stationarity makes it constant on that interface. This contradiction
establishes the claim.

## A strict independent pair-splitting bound

Let $A$ be an actual balanced three-cell global maximizer, and let
$p_{kk}=\Pr[X,Y\in A_k]$ at its correlation $\rho$. Put
$h_\rho=1-\arccos(\rho)/\pi$ and $s=S_\rho(A)$. Then

\[
 s-2\frac{s-h_\rho/3}{1+h_\rho}
            <p_{kk}<\frac{s-h_\rho/3}{1+h_\rho}\qquad(k=1,2,3).
 \tag{G13}
\]

Work in the cylinder of the active plane. Keep the cylinder of $A_k$
and split the cylinder of its complement by the sign of the
independent third coordinate. Each new cell has mass $1/3$. The
contribution of the split pair is $h_\rho(1/3+p_{kk})$: the
probability of both active observations lying outside $A_k$ is
$1/3+p_{kk}$, and the independent sign labels agree with probability
$h_\rho$ by the planar Gaussian angle formula. Optimality bounds
this new contribution by $s-p_{kk}$. Equality would make the new
partition a global maximizer, so dimension reduction would give
a common invariant direction. Since the unchanged cell uses both
active directions, this direction would have to be parallel to the
third axis, but the two new cells depend on its sign. This proves
the strict upper bound, and summing it for the other two cells
gives the lower bound. At a first tie one can use the known simplex
value for $s$. The pointwise derivative condition in (G3) involves
different data from these three diagonal entries.

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

The [polar deficit decomposition](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md)
gives the equal-circle comparison and a separate obstruction to dropping
the angular Cauchy–Schwarz remainder. Here the first angular moment
supplies equality rigidity and a bound retaining compatibility across
different radii.

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

## Angular mode deficits for compatible circular arcs

Suppose almost every centered circle is divided into three
contiguous single arcs, one per label, allowing empty arcs and
either cyclic order. For each nonzero angular frequency $n$ not
divisible by three, the contribution to the noise stability of
this partition is at most the corresponding simplex contribution
at every positive correlation. The exact mode deficit retains
the relative phases at different radii.

Write $c_n(r)$ for the vector of indicator Fourier coefficients
on the circle of radius $r$, using angular probability measure.
If $\theta_1,\theta_2,\theta_3$ are the circular endpoints there,
put $z_j(r)=e^{-in\theta_j(r)}$. Each component of $c_n(r)$ is a
signed difference of the corresponding two endpoints divided by
$2\pi i n$. Hence

\[
 \|c_n(r)\|^2
       =\frac{9-|\sum_j z_j(r)|^2}{4\pi^2 n^2}.
\]

This also holds for coincident endpoints. The simplex attains the
upper bound $9/(4\pi^2 n^2)$ if $3\nmid n$, and its mode vanishes
if $3\mid n$. For $3\nmid n$ the exact deficit between radii $r,s$
is

\[
 \frac{9}{4\pi^2 n^2}-\operatorname{Re}\langle c_n(r),c_n(s)\rangle
   =\frac{|\sum_j z_j(r)|^2+|\sum_j z_j(s)|^2}{8\pi^2 n^2}
               +\frac12\|c_n(r)-c_n(s)\|^2 .
 \tag{G14}
\]

In the polar Mehler kernel each frequency has a nonnegative
weight. To see this expand
$\exp(\kappa\cos(\theta-\varphi))$ as the product of the
exponentials of $\kappa e^{i(\theta-\varphi)}/2$ and its conjugate.
The coefficient at $n$ is a sum of nonnegative terms for
$\kappa>0$. The same expansion underlies the arc calculation
in [Heilman, (136)](https://arxiv.org/pdf/2608.00333v1).
Thus the modes still requiring a joint bound on this class
are the radial mode and the nonzero multiples of three.
The arc hypothesis is not known for a general tied maximizer.

An unweighted prefix domination for the nonzero frequencies
already fails on a single circle. Use endpoints $-u,u,\pi$
and vary $u=\pi/3+h$. The squared coefficient norm is
$[9-(2\cos(nu)+(-1)^n)^2]/(4\pi^2 n^2)$. In the block
$n=3j-2,3j-1,3j$, the two good deficits have quadratic
coefficients $3/(4\pi^2)$ each and cubic coefficients
$\sqrt3(3j-2)/(4\pi^2)$ and $-\sqrt3(3j-1)/(4\pi^2)$.
The last mode contributes $-3h^2/(2\pi^2)+O(h^4)$ to the
deficit. Thus for every fixed $J\ge1$ the sum of the first
$3J$ positive-frequency norm deficits is
$-J\sqrt3 h^3/(4\pi^2)+O_J(h^4)$, negative for small $h>0$.
This is a slice calculation. It leaves the radial balance and
the specific Mehler weights essential in a full comparison.

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
