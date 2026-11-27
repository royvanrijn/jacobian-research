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
