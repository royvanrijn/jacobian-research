# Compact changes of the Euclidean Gaussian junction

The Euclidean test in [the near-one analysis](THREE_CELL_NEAR_ONE.md#a-compact-euclidean-junction-test)
uses the planar $120^\circ$ three-sector partition $Y$ with the
unit-variance Gaussian convolution kernel. Every nontrivial compact
change holding one cell fixed increases its cut, as does every compact
change with transfers in only one cyclic orientation. In contrast,
the convex total-variation relaxation has a strict fractional
improvement with the same exterior labels. This distinction matters
when using calibration arguments for the junction.

## Finite transfers with one phase fixed

Write $G$ for the density of $N(0,I_2)$ and $\mathcal J(B)$ for the
relative cut in (P5), with ordered pairs. Fix the third phase $Y_k$,
put $\Omega=\mathbb R^2\setminus Y_k$, and allow arbitrary compact
transfers between the remaining labels $i,j$. There are no volume
constraints in this comparison. The relative cut has contributions
only from $\Omega\times\Omega$. Use the signed coordinate $t$ normal
to the $ij$ line, positive on the $Y_i$ side, so on $\Omega$ the
original binary indicator is $\eta=1_{\{t>0\}}$. For the changed
binary indicator $u=1_{B_i}$, the pointwise inequalities

\[
 |u(x)-u(y)|\ge s(x,y)(u(x)-u(y)),\qquad
 |\eta(x)-\eta(y)|=s(x,y)(\eta(x)-\eta(y)),
 \quad s(x,y)=\operatorname{sgn}(t_x-t_y)
\]

yield

\[
 \mathcal J(B)\ge2\int_\Omega (u-\eta)(x)D(x)dx,
 \qquad D(x)=\int_\Omega G(x-y)s(x,y)dy.
 \tag{J1}
\]

Absolute integrability after subtraction follows from the compact
change and $G\in L^1$. The whole-plane divergence of $s$ for the
translation-invariant kernel is zero, so $D$ is the negative of
the same integral on $Y_k$. In tangential coordinates, the normal
slices of $Y_k$ are centered intervals $[-a,a]$ or empty. For
$t=t_x>0$ the scalar integral on a nonempty slice, after removing
the positive tangential density factor, is

\[
 \int_0^\infty \phi(z)
   \big[1_{\{|t-z|\le a\}}-1_{\{|t+z|\le a\}}\big]dz.
\]

The first interval contains the second and their difference has
positive length when $a>0$. Consequently $D(x)<0$ when $t_x>0$
and $D(x)>0$ when $t_x<0$. This makes the right side of (J1) strictly
positive on any positive-measure transfer. This is a binary sign
calibration with the divergence contributed by the removed wedge.
For the halfspace version over the entire space, see Pagliari's
calibration criterion and sign construction. [Pagliari, arXiv:1905.00623, Theorem 2.3 and Lemma 2.4](https://arxiv.org/pdf/1905.00623).

## Cyclic changes with no symmetry assumption

Suppose compact transfers use only the arrows $1\to2$, $2\to3$ and
$3\to1$, or only their reverses. They may have different areas and
shapes. Let $f_{ij}$ be the indicator of the transfer set for the
arrow $i\to j$, and let $d_{ij}=e_j-e_i$ in label space. For a
competitor $B$ with vector indicator difference $h=1_B-1_Y$, the
relative cut is $2L(h)-Q(h)$, where $L(h)$ is the linear regret
against the scores $G*1_{Y_i}$ and
$Q(h)=\sum_i\langle h_i,G*h_i\rangle$. This follows by expanding
the difference of the two agreement indicators and using the same
total number of labels at each point. For the allowed transfers
$h=\sum_{i\to j}d_{ij} f_{ij}$. Each distinct pair of allowed arrows
has label inner product $-1$. Since $G>0$, the cross terms in $Q$
are nonpositive. The regrets add, and therefore

\[
 \mathcal J(B)\ge \sum_{i\to j}\mathcal J(B^{ij}),
 \tag{J2}
\]

where $B^{ij}$ makes only that one transfer from $Y$. Each term is
nonnegative by (J1), and it is positive for a nonzero transfer.
All mixed transfers with no two distinct arrows sharing a source
or sharing a destination are covered. In a three-vertex graph such
a set is contained in an oriented three-cycle or consists of a
reciprocal pair. The latter holds the unused phase fixed. This
excludes both orientations of arbitrary finite cyclic transfers
at the Euclidean junction. A negative ordinary perturbation would
need a phase sending to two distinct labels or receiving from two.

## Homogeneous minimizing junctions and stationary compact defects

An origin-homogeneous Euclidean compact-perturbation minimizer with
three positive angular masses is $Y$ up to orthogonal maps and labels.
Indeed a minimizer assigns each label to a largest convolution score.
If this fails on a set of positive measure, a sufficiently small
change inside a bounded subset with strictly positive score gap
improves the cut by its linear contribution of order the changed area,
while the interaction error is bounded in absolute value by its square.
Each cone approaches the origin and $G*1_{B_i}(0)$ is its angular
fraction. Continuity of the scores therefore makes all three fractions
equal. Analyticity on a generic circle gives finitely many switching
rays, just as in the [Gaussian conical argument](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-and-conical-maximizers).
The pair score difference is zero on each switching ray. Analyticity
extends it to the line, and the same two deconvolutions from that
argument apply with the Euclidean Gaussian of variance one. They
show that each line reflection swaps the incident labels and fixes
the third. Thus the angular widths are equal and the label sequence
is $1,2,3$ repeated $m$ times.

Homogeneity and dilation imply compact minimality at every kernel
length $\delta>0$. If $m\ge2$, group the rays on a circle of
radius $R$ starting with a $1/2$ boundary, followed by a $2/3$
boundary and a $3/1$ boundary. In each group replace the three
interior radial segments by chords connecting the first to the
second and the second to the third boundary point. Assign labels
2 and 3 to the respective caps and label 1 to the remainder of
the disk. The caps are disjoint modulo null sets, and labels match
across the circle. If $\alpha=2\pi/(3m)\le\pi/3$, the interface length
in each group becomes $4R\sin(\alpha/2)\le2R$ in place of $3R$.
For these two fixed piecewise linear partitions the cut difference
divided by $\delta$ tends to $\sqrt{2/\pi}$ times the negative
perimeter difference. The fixed BV lower and upper bounds in
(P7) apply inside any larger ball, with constant midpoint weight;
the relative energy from outside that ball has an exponentially
small tail. This contradicts minimality at all $\delta$. Hence $m=1$.

There is also an area bound for a stationary compact defect of $Y$.
Suppose $B=Y$ outside a compact set and the labels of $B$ satisfy
the winner condition for their Euclidean convolution scores. If the
area $m$ of changed labels is below a fixed absolute constant,
then $B=Y$ almost everywhere. In fact the incident scores tie on
each open exterior ray by continuity and the winner condition.
Analyticity extends these ties along the pair lines. The perturbation
$H_{ij}=G*[(1_{B_i}-1_{Y_i})-(1_{B_j}-1_{Y_j})]$ vanishes on its pair
line and has $C^1$ norm at most $2m\|G\|_{C^1}$. Thus for $x\in Y_i$,
$|H_{ij}(x)|\le C m\min(t_{ij}(x),1)$, where $t_{ij}$ is inward
distance. The [Euclidean score margin](THREE_CELL_NEAR_ONE.md#a-compact-euclidean-junction-test)
is at least $c_*\min(t_{ij},1)$. For $m<c_*/C$, every score winner
in the interior of $Y_i$ is $i$, proving the claim. Its bound is
uniform in the support radius of the defect. General microscopic
limits are not known to be homogeneous or compact perturbations.

## Nonzero gradients at winning score ties

Let $A$ minimize the Euclidean cut against all compact changes, and
write $P_i=G*1_{A_i}$. Whenever $P_i(x)=P_j(x)=\max_k P_k(x)$ for
$i\ne j$, one has $\nabla(P_i-P_j)(x)\ne0$. In particular the interface
at a strictly two-label tie is an analytic curve. At a triple tie the
three gradients are pairwise distinct.

Translate the candidate tie to $0$. We use $\phi$ for the scalar
standard Gaussian density and $G(0)=1/(2\pi)$. If all three scores and
all three gradients agree at $0$, their common score is $1/3$ and
their gradients vanish because $\sum_iP_i=1$. For every unit vector
$v$ and every $i$,

\[
 D^2P_i(0)[v,v]
  =\mathbb E[(Z_v^2-1)1_{A_i}(Z)]\le
   \mathbb E[(Z_v^2-1)_+]=2\phi(1),
\]

where $Z$ is standard Gaussian in the plane and $Z_v=Z\cdot v$.
Taylor expansion gives $\max_iP_i(x)-1/3\le(\phi(1)+o(1))|x|^2$.
Recolor $B_r$ entirely to one label and average the three choices.
The average of twice the regret is at most
$(\phi(1)+o(1))\pi r^4$. Write $m_i=|A_i\cap B_r|$ and $m=\pi r^2$.
In the interaction term $Q$ of (J2), replacing $G$ by $G(0)$ makes an
$o(r^4)$ error, since both variables lie in $B_r$ and the total areas
are $O(r^2)$. The average of $Q$ is therefore
$G(0)(m^2/3+\sum_i m_i^2)+o(r^4)\ge\pi r^4/3+o(r^4)$.
As $\phi(1)<1/3$, some recoloring improves the cut.

It remains to test a tie of labels 1 and 2 with equal gradients.
If their common score is strictly above label 3, use any half-disk
of $B_r$. If all scores agree, the first argument lets us assume the
common gradient differs from that of label 3. Use the half-disk on
which the common first-order term of 1 and 2 exceeds that of 3. Let
$E$ be the part of the chosen half-disk occupied by 1 and 2. By the
winner condition and the Taylor remainder, the omitted part lies in a
strip of width $O(r^2)$ against the diameter, hence has area $O(r^3)$.

Put $H=D^2(P_1-P_2)(0)$, and for the unit direction $v_\theta$ write
$q(\theta)=H[v_\theta,v_\theta]=a+b\cos(2\theta-\alpha)$ with $b\ge0$.
If $q$ is nonzero, let $p$ be the fraction of the half-circle on which
$q>0$. Up to $o(r^2)$ in area, the masses of labels 1 and 2 in $E$
are $p\pi r^2/2$ and $(1-p)\pi r^2/2$. This follows after rescaling
the disk by $r$, since the sign of $P_1-P_2$ converges to the sign
of $q$ except at its finitely many angular zeros. Recolor $E$ entirely
to label 1 or entirely to label 2, with equal probabilities. Twice
the average regret is $\int_E|P_1-P_2|=(\pi r^4/8)\mathbb E_\theta|q|
+o(r^4)$, where $\theta$ is uniform on a half-circle. With $m_i$
now the areas in $E$, the average interaction is
$G(0)(m_1^2+m_2^2)+o(r^4)$, or
$(\pi r^4/8)(p^2+(1-p)^2)+o(r^4)$. If $q=0$ the regret is $o(r^4)$
and the interaction is at least $\pi r^4/16+o(r^4)$.

For $q$ of one sign one has $p^2+(1-p)^2=1$ and
$\sup_\theta|q(\theta)|\le\mathbb E|Z_v^2-1|=4\phi(1)<1$.
If $q$ changes sign then $|a|<b$. In an eigenbasis of $H$, the
half-difference of the two directional second derivatives satisfies

\[
 b\le \tfrac12\mathbb E|Z_1^2-Z_2^2|=2/\pi.
\]

Here the expectation uses $|1_{A_1}-1_{A_2}|\le1$, and the last
identity follows by rotating $(Z_1,Z_2)$ by $45^\circ$.
With $u=\arcsin(|a|/b)$ the elementary circular integrals give
$p^2+(1-p)^2=1/2+2u^2/\pi^2$ and
$\mathbb E_\theta|q|=(2b/\pi)(\cos u+u\sin u)$. The latter is
strictly smaller: $\pi^2/2+2u^2-4(\cos u+u\sin u)$ is positive at
zero and has derivative $4u(1-\cos u)\ge0$. In all cases the
average relative cut decreases, contradicting minimality.

At a rank-two triple tie the score differences give local analytic
coordinates and identify three incident analytic arcs. Pairwise
distinct gradients may still be collinear, which the test above
does not exclude.

## Conserved stress and the three affine asymptotes

Every planar compact-perturbation minimizer for $G$ with three
nonempty phases has three concurrent $120^\circ$ asymptotic rays.
Outside a compact set the interfaces are analytic normal graphs
over those rays, with offsets and slopes tending to zero.

The density argument (P6) in the [near-one note](THREE_CELL_NEAR_ONE.md#uniform-phase-density-and-two-label-blow-ups-on-open-rays)
applies with $W(u,v)=G(v-u)$ and no volume repair. Its recurrence
holds at all radii, and winner scores are at least $1/3$ at their
centers. Each winning label thus occupies a fixed fraction of
every ball at its center with radius $R\ge R_0$. Likewise the
intermediate-scale proof (P7) applies to any $U(x_j+R_j\,\cdot)$
with $R_j\to\infty$, constant midpoint factor, $\delta_j=1/R_j$,
and no repair. All fixed-center blow-downs are centered tripods:
each phase has a fixed winner point whose density estimate gives a
positive fraction of each phase at the origin of the limit.
The same bound upgrades local convergence to convergence in local
distance of the closed winner sets, by testing balls centered
outside or inside the limiting classes.

Write $q(x,y)$ for disagreement. The stress tensor is

\[
 \mathsf S(w)=2I\int G(z)q(w,w+z)\,dz
 -\int G(z)z\otimes z\int_0^1 q(w-tz,w+(1-t)z)\,dt\,dz .
 \tag{J5}
\]

Pull back both variables by $x\mapsto x+sV(x)$ for $V$ smooth and
compactly supported. The Jacobian derivatives give the first term.
The kernel derivative is $-G(z)z\cdot(V(x+z)-V(x))$; write the
difference of $V$ as the integral of $DV$ along the segment to
obtain the second term. The differentiated density is dominated
by a polynomial times a Gaussian, with one endpoint in a fixed
compact set. Equilibrium gives $\int\mathsf S:DV=0$, so $\mathsf S$
is divergence-free row by row. For the counterclockwise quarter
turn $J$ this gives global potentials

\[
                       D\Psi=\mathsf S J^T,\qquad \nabla F=J\Psi .
 \tag{J6}
\]

The rows of $\mathsf S J^T$ are curl-free, and symmetry of $\mathsf S$
gives $\operatorname{div}\Psi=0$. Integrate smooth mollifications and
fix constants at the origin to obtain the potentials by local uniform
limits. On a single-phase ball $B_d(w)$ the integrands in (J5) vanish
for $|z|<d$, so $|D\Psi(w)|\le Ce^{-c d^2}$. The resulting bound
also controls increments along any segment of such phase depth,
by mollifying in a smaller tube around it.

For each large $r$, choose a labeled centered tripod agreeing with
the winner sets in local distance to $o(r)$ in $B_{4r}$, and let
$p_i(r)$ be a point on its label-$i$ bisector at radius $r$.
The approximations exist uniformly by blow-down compactness.
At scales $r$ and $2r$ the labeled rays are close by overlap.
The segment from $p_i(r)$ to $p_i(2r)$ has phase depth at least
$cr$ by the angular margin. The dyadic chains give
$\Psi(p_i(r))\to L_i$ and $F(p_i(r))-J L_i\cdot p_i(r)\to C_i$:
along each segment the sum of the remaining increments bounds
$|\Psi-L_i|\le Ce^{-cr^2}$. Integrating $\nabla F-J L_i$ gives a
summable bound for the second limit. Two radii within a factor of
two can be joined by similar deep segments, so these limits are
independent of the starting scales.

Take any sequence of winning pair ties $x_j$ tending to infinity.
For sufficiently small fixed $c_0>0$, only the pair can win in
$B_{10c_0|x_j|}(x_j)$, by the tripod approximation at scale $|x_j|$.
For every $R\le r\le c_0|x_j|$, the pair phases approximate
halfplanes through $x_j$ in $B_{4r}(x_j)$, with discrepancy confined
to a slab of width $\varepsilon r$ for $\varepsilon$ arbitrarily
small as $R,j$ grow. Otherwise take a contradicting sequence at
scales tending to infinity. If $r/|x_j|\to0$, the corresponding
blow-down is a two-phase partition by the third-phase exclusion
and (P7). Density puts its line through the center. If $r/|x_j|$
stays bounded below, take an origin blow-down, which already gives
a single line in the comparison ball. In both cases density turns
local convergence into the slab estimate.

After a subsequence, the unit-scale partitions at $x_j$ converge to
a two-phase Euclidean halfplane pair, by the constant-weight version
of the [noise-length compactness argument](THREE_CELL_NEAR_ONE.md#blow-ups-of-moving-maximizers-minimize-the-euclidean-cut)
with no repair, and the two-phase classification cited there.
The limiting scores tie at zero by local strong convergence and
Gaussian tails, so the line is through zero. Let $n$ point toward
label $i$, the other label being $k$. For fixed large $R$ the slab
normal at radius $R$ is close to $n$, and consecutive dyadic slab
normals are close by overlap. On each side, connect the points
$x_j\pm r n_r$ at consecutive radii, starting with $x_j\pm Rn$.
Each segment has phase depth at least $cr$. Stop at a radius in
$[c_0|x_j|/2,c_0|x_j|]$ and connect the last points to
$p_i(|x_j|),p_k(|x_j|)$ inside the corresponding deep convex
sectors of the approximating tripod. The depth of the last segments
is bounded below by a fixed multiple of $|x_j|$.
For $\ell=i$ on the plus side and $\ell=k$ on the minus side, the
same tail estimates along the chains give

\[
 |\Psi(x_j\pm Rn)-L_\ell|
 +|F(x_j\pm Rn)-J L_\ell\cdot(x_j\pm Rn)-C_\ell|
 \le Ce^{-cR^2}
 \tag{J7}
\]

for all sufficiently large $j$ at each fixed large $R$. The constants
are independent of $R$ and of the chosen subsequence.

For the limiting halfplane $H$, its tensor $\mathsf S_H(sn)$ is even
in $s$ and independent of the tangential coordinate. Direct integration
in (J5) gives, for $\tau=J^T n$,

\[
 \int_{\mathbb R}\mathsf S_H(sn)\,ds
 =2\mathbb E|Z\cdot n| I-\mathbb E[(Z\otimes Z)|Z\cdot n|]
 =\sqrt{2/\pi}\,\tau\otimes\tau .
 \tag{J8}
\]

The disagreement interval of $s$ has length $|z\cdot n|$ at each
displacement $z$. The last identity uses independence of Gaussian
coordinates and $\mathbb E|Z\cdot n|^3=2\mathbb E|Z\cdot n|$.
Strong local label convergence gives local $L^1$ convergence of
the translated stresses: for bounded $z$ shift the integration
variable to each endpoint in (J5), then control the larger $z$
by Gaussian tails. The translated $\Psi$ have locally uniformly
convergent subsequences by their Lipschitz bound and (J7).
A limit has derivative $\mathsf S_H J^T$ and endpoint values
$L_k,L_i$ on the normal line by (J7). Thus (J8) gives
$L_i-L_k=\sqrt{2/\pi}\,\tau$. The microscopic normal $n$ is the
same for every far tie of this pair.

Evenness of $\mathsf S_H$ makes the limiting potential on the
normal line its midpoint $(L_i+L_k)/2$ plus an odd function.
Integrate (J6) along that line. The limit of
$F(x_j+Rn)-F(x_j-Rn)$ for fixed $R$ is $R J(L_i+L_k)\cdot n$.
Evaluate the difference also by (J7) and send $R$ to infinity:

\[
               (J L_i-J L_k)\cdot x_j+C_i-C_k\ \longrightarrow\ 0 .
 \tag{J9}
\]

The translated scores at every far tie converge in $C^1$ to their
halfplane counterparts. For any fixed $v$ with $v\cdot n>0$, the
directional derivative of $P_i-P_k$ in direction $v$ is therefore
positive at all sufficiently far ties. If a blow-down ray for the
pair had inward normal $n_*\ne n$, choose such a $v$ with
$v\cdot n_*<0$. A segment parallel to $v$ crossing that ray away
from the junction, rescaled back far out, has no third winner and
goes from positive to negative score difference. This contradicts
the positive derivative at all its zeros. Thus every pair normal
agrees with its blow-down normal. The oriented rays are at
$120^\circ$, and the three lines in (J9) are concurrent since their
normal vectors and offsets telescope in cyclic label order and the
first two normals are independent.

Far from the junction every interface lies in a tube of any fixed
small relative width about its blow-down ray. Normal segments
across the tube have opposite unique winners at the endpoints
and all zeros are strict crossings with the same orientation.
Each segment meets the interface once. The analytic implicit
function theorem gives the graph. The microscopic halfplane
convergence gives slope tending to zero, and (J9) gives offset
tending to zero from the affine asymptote.

## Exponential convergence of the branch graphs

Use coordinates $(s,y)$ with one affine asymptote $y=0,\ s>0$, and let
$g(s)\to0$ be the graph of the corresponding interface between labels
$i$ above and $k$ below. For some fixed $\eta>0$, the partition on
$s$ large and $|y|<4\eta s$ is that binary graph. Extend $g$ to a
bounded function on the whole line and denote the ordinary two-phase
partition given by this extended graph by $U^g$. We prove that
$|g(s)|\le C e^{-\kappa s}$ for some $\kappa>0$. Constants in this
proof may depend on the partition and the choice of extension.

Let $J$ be counterclockwise rotation by $\pi/2$, and use the stress
potentials $D\Psi=SJ^\top$, $\nabla F=J\Psi$ and their phase constants
$L_a,C_a$ from the affine asymptote proof. At $y=\eta s$ and
$y=-\eta s$, respectively, the deep-phase estimates give
\[
 \Psi(s,\pm\eta s)=L_{i,k}+O(e^{-c s^2}),\qquad
 F(s,\pm\eta s)=J L_{i,k}\cdot(s,\pm\eta s)+C_{i,k}+O(e^{-c s^2}).
\]
The horizontal asymptote implies $C_i=C_k$ and
$J(L_i-L_k)\cdot(1,0)=0$.
Define the signed moment of the graph stress by
\[
 M_g(s)=\int_{\mathbb R} y S^g_{11}(s,y)\,dy .
\]
For $|y|\le\eta s$, the true and graph stresses differ by
$O(e^{-c s^2})$, since every endpoint in the stress formula is
within $|z|$ of $(s,y)$ and the two partitions agree there when
$|z|\le c' s$ with $c'>0$ sufficiently small. For $|y|>\eta s$,
the graph stress is $O(e^{-c(|y|-\|g\|_\infty)^2})$, since points at
distance less than $|y|-\|g\|_\infty$ remain on the same side of the
graph. The Gaussian bounds in both estimates include the polynomial
factor $|z|^2$ in the stress. Thus, writing $D=\eta s$ and using
$\partial_y\Psi_1=S_{11}$ and $\partial_y F=\Psi_1$, we obtain
\[
 \begin{aligned}
 M_g(s)
 &= \int_{-D}^D y S_{11}(s,y)\,dy+O(e^{-c s^2})\\
 &= D[\Psi_1(s,D)+\Psi_1(s,-D)]
    -[F(s,D)-F(s,-D)]+O(e^{-c s^2})
  = O(e^{-c s^2}).
 \end{aligned}
 \tag{E1}
\]
The terms linear in $D$ cancel since $(J L)_2=L_1$. Those linear
in $s$ and the constant difference vanish by the two affine
asymptote identities. The stress primitives are Lipschitz on compact
sets, so the integration by parts holds for almost every $s$,
which suffices below.

Write $P_\sigma$ for convolution in one dimension with the centered
Gaussian of standard deviation $\sigma$, with $P_0=I$. Let $H,B$
be independent standard one-dimensional Gaussians, let
$\Phi(t)=\Pr(B\le t)$, and put $c_0=2\Phi'(0)=\sqrt{2/\pi}$.
At a point $(s,g(s))$ on the true interface the incident scores tie.
Replacing the true partition by $U^g$ in the convolution costs
$O(e^{-c s^2})$ by the same agreement tube, since $g$ is bounded.
For the graph, conditional on the horizontal displacement $H=h$,
the difference of the upper and lower scores is
$2\Phi(g(s)-g(s+h))-1$. Taylor expansion gives
\[
 (I-P_1)g(s)=O(e^{-c s^2})+R_1(s),\qquad
 |R_1(s)|\le (K_1*|g|^3)(s),
 \tag{E2}
\]
where $K_1$ is a finite positive measure with all exponential
moments. Indeed $2\Phi(t)-1=c_0t+O(|t|^3)$, and
$|g(s)-g(s+h)|^3\le4(|g(s)|^3+|g(s+h)|^3)$. One may take $K_1$
to be a fixed multiple of the sum of a unit atom and the law of $H$.

The moment in (E1) gives a second expansion. Let
$\mathcal N(a)=\mathbb E|a-B|$. Then $\mathcal N(0)=\mathcal N''(0)=c_0$, $\mathcal N'(0)=0$,
$\mathcal N(a)=c_0+O(a^2)$, and $\mathcal N'(a)=c_0 a+O(|a|^3)$. The normal
Gaussian integration identity is
$\mathbb E[B|a-B|]=-\mathcal N'(a)$, by Gaussian integration by parts.
The contribution from $2I\int G(z)q(w,w+z)\,dz$ to $M_g$ is
\[
 \mathbb E_H\left[(g(s)+g(s+H))\mathcal N(a)+\mathcal N'(a)\right],
 \qquad a=g(s+H)-g(s).
 \tag{E3}
\]
In fact the two label heights in the disagreement indicator are
$g(s)$ and $g(s+H)-B$. The signed first moment of the interval
between any two heights $p,q$ is $(p+q)|p-q|/2$. This proves
(E3) after the factor two in the stress and expectation in $B$.
Its linear part is $2c_0 P_1g(s)$.

For the other stress contribution put $g_-=g(s-tH)$ and
$g_+=g(s+(1-t)H)$. The heights are $g_-+tB$ and
$g_+-(1-t)B$. The same interval computation gives the following term to
subtract from (E3).
\[
 \frac12\mathbb E_H H^2\int_0^1
   \left[(g_-+g_+)\mathcal N(d)-(2t-1)\mathcal N'(d)\right]dt,
 \qquad d=g_+-g_- .
 \tag{E4}
\]
Its linear part is
\[
 c_0\mathbb E_H H^2\int_0^1 [t g_-+(1-t)g_+]dt
   =2c_0\int_0^1 t\,\mathbb E_H[H^2g(s+tH)]\,dt ,
\]
where reflection of $H$ and substitution $t\mapsto1-t$ give the
last equality. Every remainder in (E3) and (E4) is bounded by a
constant times the sum of cubes of its two values of $g$, with
the factor $H^2$ retained in (E4). The measures with shifts $H$
and $\pm tH$, weighted by $1$ or $H^2$ and averaged over
$0\le t\le1$, have all exponential moments. Hence the combined
remainder is bounded by a convolution $K_2*|g|^3$ with another
finite positive measure $K_2$ of that type.

The multiplier of the resulting linear operator $M_g/c_0$ is
\[
 2e^{-\xi^2/2}-
    2\int_0^1 t(1-\xi^2t^2)e^{-\xi^2t^2/2}\,dt
   = \frac{2(1-e^{-\xi^2/2})}{\xi^2}.
\]
To check the identity, set $u=t^2$ and integrate
\((1-\xi^2 u)e^{-\xi^2 u/2}
 =\frac{d}{du}(2u e^{-\xi^2u/2})-e^{-\xi^2u/2}\).
The right side has value one at $\xi=0$ and is the multiplier of
$A=\int_0^1 P_{\sqrt u}\,du$. This computation is an identity of
convolution measures, so it applies to the bounded extension of
$g$. Equations (E1), (E3) and (E4) therefore give
\[
 Ag(s)=O(e^{-c s^2})+R_2(s),\qquad
 |R_2(s)|\le (K_2*|g|^3)(s).
 \tag{E5}
\]

The sum of (E2) and (E5) gives, for $s$ sufficiently large,
\[
 |(\mathcal B g)(s)|\le C e^{-c s^2}+(K*|g|^3)(s),
 \qquad \mathcal B=I-P_1+A,\quad K=K_1+K_2.
 \tag{E6}
\]
The Fourier multiplier of $\mathcal B$ is
$m(\xi)=(1-e^{-\xi^2/2})(1+2/\xi^2)$, with removable value
$m(0)=1$. In every closed strip \(|\operatorname{Im}\xi|\le
\kappa_0<\sqrt{2}\) it is nonzero. Besides the zeros
$\pm i\sqrt{2}$ of the rational factor, the possible nonzero zeros
are solutions of $\xi^2=4\pi i n$, $n\in\mathbb Z\setminus\{0\}$.
The function $m^{-1}-1$ is analytic on the strip and is bounded
there by $C_{\kappa_0}(1+|\operatorname{Re}\xi|)^{-2}$.
Fourier inversion and contour shifting yield a continuous inverse
transform $k$ with $|k(s)|\le C_{\kappa_0} e^{-\kappa_0|s|}$.
The vertical contour contributions tend to zero by that uniform
quadratic bound. It follows by uniqueness of Fourier transforms
of finite measures that $\delta_0+k(s)\,ds$ is a convolution
inverse of the finite measure defining $\mathcal B$.

Here are the absorption details on the whole line, which allow
for an arbitrary bounded extension of $g$. Fix
$0<\kappa<\kappa_0$. Let $b=\delta_0+|k(s)|\,ds$ and for a
finite positive measure $\lambda$ write
$M_\kappa(\lambda)=\int e^{\kappa|h|}\,d\lambda(h)$.
Choose $\varepsilon>0$ with
$\varepsilon^2 M_\kappa(b) M_\kappa(K)<1/2$. Since $g(s)\to0$,
choose $T$ large enough that $|g(s)|\le\varepsilon$ for
$s>T$. On $s\le T$, bound $|g|^3$ by
$G_*^3 1_{\{s\le T\}}$, where $G_*=\|g\|_\infty$.
On the bounded-above range of $s$ where (E6) has not been shown,
$\mathcal B g$ is bounded since $\mathcal B$ is a finite
convolution measure. Thus globally
\[
 |\mathcal B g| \le \mathcal F+\varepsilon^2 K*|g|,
\]
where $\mathcal F$ is the sum of a constant times a left-half-line
indicator, a Gaussian right tail, and
$G_*^3K*1_{\{s\le T\}}$. In particular $e^{\kappa s_+}\mathcal F(s)$
is bounded. For the last summand, use
$e^{\kappa s_+}\le e^{\kappa|h|}e^{\kappa(s-h)_+}$ on
$s-h\le T$. Applying the convolution inverse gives
\[
 |g| \le b*\mathcal F+\varepsilon^2 (b*K)*|g|.
\]
The bounded weights $w_N(s)=e^{\kappa\min(s_+,N)}$ satisfy
$w_N(s)/w_N(s-h)\le e^{\kappa|h|}$.
Consequently
\[
 \|w_N g\|_\infty
 \le M_\kappa(b)\|e^{\kappa s_+}\mathcal F\|_\infty
   +\varepsilon^2 M_\kappa(b)M_\kappa(K)\|w_N g\|_\infty .
\]
Each norm on the left is finite before absorption because $g$
and $w_N$ are bounded. Absorb uniformly and send $N\to\infty$.
This proves the exponential estimate, first almost everywhere and
then everywhere on the graph by continuity.

The same argument applies to all three branches. Outside a fixed
disk the labels agree with the asymptotic tripod away from the
strips between each graph and its asymptote. Their total area
there is bounded by a constant times
$\sum_{\text{branches}}\int_{s_0}^\infty |g(s)|\,ds<\infty$.
Thus the partition differs from the translated tripod by a set of
finite area. The junction core remains unrestricted by this estimate.

## Rigidity from one exact two-sided separation

Let $U$ be an ordinary Euclidean compact-perturbation minimizer
with three phases and $Y$ its affine asymptotic tripod. Assume
one pair line separates its incident labels everywhere: label
$A$ occurs only on the left and label $B$ only on the right.
Then $U=Y$ almost everywhere. The hypothesis concerns only
this single line.

By exponential convergence, $U$ differs in finite area from $Y$.
Minimality extends to any competitor differing from $U$ in finite
area: truncate its change to a ball and use the estimate that two
competitors differing on a set of area $a$ have relative cuts
differing by at most $2a$. The estimate follows from $\int G=1$
and the fact that the integrand can differ only when an endpoint
lies in the change set.

Choose coordinates with pair line $x=0$ and let $D=\{x>0\}$.
On $D$ let $m$ indicate the third phase $M$, and let $n$ indicate
its reflection from the left. The other labels are determined
by the exact separation assumption. For $p,q\in D$ put
$K^-(p,q)=G(p-q)$, $K^+(p,q)=G(p-\sigma q)$ with $\sigma$ reflection
across the line. Write $s=m+n$ and $d=m-n$. In the agreement
quadratic form the terms depending on $d$ are exactly the product
$d(p)d(q)$ against $K^- -K^+/2$, as is seen by expanding the
same-side agreements and the cross-side terms
$K^+(m(p)n(q)+n(p)m(q))$. Replacing $m,n$ by their pointwise
maximum and minimum leaves $s$ fixed and replaces $d$ by $|d|$.
The coefficient is strictly positive. If $d$ had both signs on
positive-area sets, agreement would strictly increase, contrary
to minimality. Thus $d$ has one sign. Now compare to the two
competitors obtained by duplicating $m$ on both sides or $n$ on
both sides. The same-side terms sum to twice the original terms,
while the summed cross-side agreement increases by
$2\iint K^+(p,q)d(p)d(q)\,dp\,dq$. This is strictly positive
unless $d=0$ almost everywhere. All comparisons change finite
area since the two sides have the same third phase in $Y$.
Hence $m=n$.

For a reflection-symmetric oriented profile, the relative cut
from the tripod reduces to the scalar cut with kernel
$K_*=2K^-+K^+$ and potential $1-\Phi(x)$:

\[
 \iint_{D^2} K_*(p,q)
  \big(|m(p)-m(q)|-|m_Y(p)-m_Y(q)|\big)\,dp\,dq
            -2\int_D (1-\Phi(x))(m-m_Y)\,dp .
 \tag{J10}
\]

This follows by expanding the same-side binary agreements and
$2K^+m(p)m(q)$, using $\int_D K^+(p,q)\,dq=1-\Phi(x)$.
Orient the vertical coordinate so that $m_Y=1_{\{g>0\}}$ for
$g(x,y)=y-f_0(x)$ with affine $f_0$. At any $p$ with $g(p)=c$,
the tripod translated vertically by $c$ has a tie there between
$M$ and $B$. Its score equality reads
$\int_D K_*(p,q)1_{\{g(q)>c\}}\,dq=\Phi(x)$, by the two reflected
convolutions. It follows that $\operatorname{sgn}(g(p)-g(q))$
has kernel divergence $1-\Phi(x)$. The inequality in (J1) for
this sign function then bounds the cut part in (J10) below by
the potential part, so the cut difference is nonnegative.
Equality requires preservation of the order of $g$ for almost
every pair with distinct labels, since $K_*>0$. Hence $m$ is a
level threshold of $g$ up to null sets. Any nonzero threshold
shift changes infinite area. Minimality of $U$ also gives the
reverse cut inequality against $Y$, by finite-area comparison,
so equality holds and $m=m_Y$. This proves the claimed rigidity.

## A fractional obstruction to convex total-variation calibration

For a simplex-valued $u$ that agrees with the vector indicator $1_Y$
outside a ball, define the convex relaxed cut by its finite difference

\[
 \mathcal J_{\rm TV}(u)=\frac12\sum_i\iint G(x-y)
  \big( |u_i(x)-u_i(y)|-|1_{Y_i}(x)-1_{Y_i}(y)|\big)dx\,dy.
 \tag{J3}
\]

It agrees with $\mathcal J$ for ordinary partitions. There exists a
compactly changed fractional $u$ with $\mathcal J_{\rm TV}(u)<0$.

To construct it, choose a unit vector $v$ with opposite nonzero
projections on the two inward normals of $Y_1$. Let $\chi_R$ be a
smooth radial cutoff equal to one on $B_R$, zero outside $B_{2R}$,
with $\|\nabla\chi_R\|_\infty\le C/R$. For $R$ large and $|t|\le t_0$
fixed small, $F_t(x)=x+t v\chi_R(x)$ is a diffeomorphism. Denote the
ordinary image partition by $Y(t)$. We claim

\[
                         |\mathcal J(Y(t))|\le C t^2/R .
 \tag{J4}
\]

The first derivative at zero vanishes. For instance, the boundary
first variation on each facet is a normal speed times the difference
of the two convolution scores, which agree there by reflection.
This formula follows by differentiating the two cell indicators
against the smooth convolution or by divergence on the three sectors
for the compactly supported deformation. To bound the second
derivative uniformly, pull back both variables. Integrate on the
fixed set of differing-label pairs of $Y$, with density

\[
 G(x-y+t v(\chi_R(x)-\chi_R(y)))
 (1+t v\cdot\nabla\chi_R(x))(1+t v\cdot\nabla\chi_R(y)).
\]

Set $z=x-y$. From the cutoff Lipschitz bound and Gaussian derivatives,
the second derivative of this density is bounded by
$C R^{-2}(1+|z|^4)\exp(-|z|^2/4)$, and is zero if both endpoints
are outside $B_{2R}$. For fixed $z$, the admissible $x$ lie in
$B_{2R}$ or its translate and within distance $|z|$ of one of the
three supporting lines. Indeed differing labels imply a supporting
line crosses the joining segment. The total area in those two disks
and three strips is at most $C R |z|$. Integration now gives a bound
$C/R$ for the second derivative, proving (J4).

Take $u=(1_{Y(t)}+1_{Y(-t)})/2$ for a fixed $0<t<t_0$. In a fixed
bounded neighborhood, the two deformations are just translations
when $R$ is large. There are two positive-area rectangles $E,F$ in
strips around the two open facets of $Y_1$ away from the junction.
Choose them so that $E$ enters label 1 under one translation and
leaves under the other, while the two orientations are reversed
on $F$. The choice is possible by the opposite normal projections.
On $E\times F$ there is a fixed, strictly positive integrated Jensen
defect for the first component of (J3). Its two absolute differences
are one, while the absolute difference of the midpoint values is
zero. Convexity elsewhere and (J4) therefore give
$\mathcal J_{\rm TV}(u)\le C t^2/R-c(t)$ with $c(t)>0$ independent of
$R$. Take $R$ large. In particular a certificate proving
minimality of $1_Y$ in this full convex relaxation cannot hold.

This fractional energy differs from the soft agreement extension
using $1-u(x)\cdot u(y)$. Weak approximations by hard labels converge
for compact changes in the soft agreement extension, and their
limit need not have the total-variation value (J3). The strict
decrease in (J3) therefore leaves the ordinary compact-change
comparison for $Y$ open. The restricted comparisons (J1)–(J2) and
the relaxation counterexample are analytic.
