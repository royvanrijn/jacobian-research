# Balanced reflected swaps: a Hermite obstruction and full-noise loss

There is no finite constant $C$ such that
$W_2(A)-W_2(T)\le C\delta(A)$ for every balanced four-cell Gaussian
partition with $\delta(A)=W_1(T)-W_1(A)>0$. This fails already in dimension
three, arbitrarily close to $T$ in Gaussian measure and in centroid Gram
space. Consequently the proposed uniform higher-Hermite bound also fails.
Every proper reflected swap nevertheless strictly loses full noise stability
at every $0<\rho<1$, by the kernel comparison below. The first-Hermite theorem
and a uniform interval of positive-correlation optimality remain open.

## Normalization and comparison with the literature

Use the normalized probabilists' Hermite convention of the
[current Hermite argument](BALANCED_FOUR_CLASS_POLYNOMIAL_PARTITIONS.md#hermite-bounds).
For any cell put

$$
m_i=\int_{A_i}x\,d\gamma_3,\qquad
Q_i=\int_{A_i}(xx^{\mathsf T}-I)\,d\gamma_3.
$$

Then $W_1=\sum_i\|m_i\|^2$ and $W_2=\frac12\sum_i\|Q_i\|_F^2$.
Write the tetrahedral vertices as the even sign vectors $s_i$, with
$s_1=(1,1,1)$, $s_2=(1,-1,-1)$, $s_3=(-1,1,-1)$ and
$s_4=(-1,-1,1)$, and let $T_i$ maximize $s_i\cdot x$. The existing
tetrahedral moment evaluations give

$$
m_i(T)=c s_i,\qquad Q_i(T)=h(s_i s_i^{\mathsf T}-I),\qquad
c=\frac1{\sqrt\pi}\left(\frac14+\frac{\arcsin(1/3)}{2\pi}\right),
\quad h=\frac1{2\pi\sqrt3}.
\tag{1}
$$

In particular $W_1(T)=12c^2$ and $W_2(T)=12h^2=1/\pi^2$.

The closest checked four-cell positive-correlation theorem is
[Heilman, arXiv:2209.11216, Theorem 1.5](https://arxiv.org/html/2209.11216).
It assumes that an equal-volume maximizing partition is hyperstable.
Definition 1.3 imposes, in addition to the usual variation conditions, a
vanishing mixed variation in the deformation parameter and correlation at
zero. This hypothesis is not supplied by the unfinished first-Hermite
theorem, and the source gives no unrestricted relative estimate of the
form above. The comparison uses the actual definition and theorem, not
the abstract. No novelty claim is made here.

## An exactly balanced family

Use orthonormal coordinates

$$
U=x_1,\qquad N=(x_2+x_3)/\sqrt2,\qquad V=(x_2-x_3)/\sqrt2.
$$

The reflection $g(U,N,V)=(U,-N,V)$ interchanges $T_1,T_2$ and fixes
$T_3,T_4$. For every real $R\ge1$, let

$$
E_R=\{-R-1<U<-R,\ -4R-1<N<-4R,\ |V|<1\}.
\tag{2}
$$

This box lies strictly inside $T_2$: the required inequalities are
$N<0$ and $U>(N+|V|)/\sqrt2$, and throughout the box

$$
U>-R-1>\frac{-4R+1}{\sqrt2}>\frac{N+|V|}{\sqrt2}.
$$

The middle inequality holds at $R=1$ and its margin increases with $R$.
Thus $gE_R\subset T_1$. Exchange just these two boxes:

$$
A_{R,1}=(T_1\setminus gE_R)\cup E_R,\qquad
A_{R,2}=(T_2\setminus E_R)\cup gE_R,\qquad
A_{R,3}=T_3,\quad A_{R,4}=T_4.
\tag{3}
$$

Reflection invariance gives $\gamma_3(E_R)=\gamma_3(gE_R)=p_R$,
so all four cell volumes remain exactly $1/4$. The changed-label set has
measure $2p_R\to0$. All changes occur outside a ball whose radius tends
to infinity.

Set

$$
a_R=-\int_{E_R}N\,d\gamma_3>0,\qquad
b_R=\int_{E_R}UN\,d\gamma_3>0.
$$

Pointwise on $E_R$ we have $-U>R$ and $-N>4R$, hence

$$
b_R>Ra_R,\qquad a_R>4Rp_R,\qquad
a_R\le\phi(4R)\le\phi(4)<c/\sqrt2.
\tag{4}
$$

For the upper bound, enlarge the box to the halfspace $N<-4R$ and use
$\int_{4R}^\infty t\phi(t)\,dt=\phi(4R)$. The last comparison follows,
for example, from $c>1/8$, $\sqrt2<2$, and
$\phi(4)<1/256<1/16<c/\sqrt2$.

## Exact first- and second-moment changes

Let $u,n,v$ be the coordinate unit vectors for $U,N,V$. Reflection and
the $V\mapsto-V$ symmetry of $E_R$ give

$$
\Delta m_1=-2a_R n,\quad \Delta m_2=2a_R n,\qquad
\Delta Q_1=2b_R(un^{\mathsf T}+nu^{\mathsf T}),\quad
\Delta Q_2=-\Delta Q_1.
\tag{5}
$$

The identity terms in $Q$ cancel because the exchanged volumes agree.
The other two cells are unchanged. From (1),
$m_1(T)-m_2(T)=2\sqrt2c\,n$ and
$Q_1(T)-Q_2(T)=2\sqrt2h(un^{\mathsf T}+nu^{\mathsf T})$.
Expanding the squared norms, including their quadratic terms, yields

$$
\delta(A_R)=8\sqrt2c\,a_R-8a_R^2>0,\qquad
W_2(A_R)-W_2(T)=8\sqrt2h\,b_R+8b_R^2>0.
\tag{6}
$$

In particular $\delta(A_R)\to0$, the first moments converge to those of
$T$, and

$$
\frac{W_2(A_R)-W_2(T)}{\delta(A_R)}
\ge\frac hc\frac{b_R}{a_R}>\frac hc R\longrightarrow\infty.
\tag{7}
$$

This is an analytic counterexample family for all $R\ge1$, not an
extrapolation from numerical samples. Taking products with unused Gaussian
coordinates preserves both masses and balance in every dimension $n\ge3$.

The boxes introduce additional interface components. They are admissible
measurable balanced partitions, and can have arbitrarily small changed-label
measure. This does not assert failure within a class restricted to small
normal graphs over the six original facets. One may instead use the radius
$1/4$ ball centered at $(-R-1/2,-4R-1/2,0)$ in $(U,N,V)$ coordinates,
reflected and exchanged as in (3). It lies inside $E_R$, is $V$ symmetric,
and satisfies all the estimates (4)--(7). The resulting new interfaces are
spheres, so the obstruction does not depend on box corners.

## What fixed infinitesimal interface tests see

For comparison, consider a smooth deformation with normal speeds $f_{ij}$
compactly supported in the interiors of the six original facets. The
deformation can be required to preserve all four cell volumes; in
particular its first mass variations must vanish. Orient
$n_{ij}=(s_i-s_j)/(2\sqrt2)$ into $T_i$ and use $f_{ij}$ as the outward
speed of $T_i$, so the boundary moves with velocity $-f_{ij}n_{ij}$.
Primes below are derivatives at the undeformed tetrahedron. On a facet,
the first-variation score differences for both $W_1$ and $W_2$ vanish.
Consequently the acceleration terms drop out, giving

$$
W_1''=2\sum_i\|m_i'\|^2-2d\sum_{i<j}\int_{F_{ij}}f_{ij}^2\phi_3\,d\mathcal H^2,
\qquad d=2\sqrt2c,
\tag{9}
$$

$$
W_2''=\sum_i\|Q_i'\|_F^2
-4\sqrt2h\sum_{i<j}\int_{F_{ij}}U_{ij}f_{ij}^2\phi_3\,d\mathcal H^2.
\tag{10}
$$

Here $U_{ij}$ is the signed shared coordinate of the pair $s_i,s_j$,
positive on the facet interior. To check the local term, in the facet's
orthonormal frame $Q_i-Q_j=2\sqrt2h(un^{\mathsf T}+nu^{\mathsf T})$.
The outward derivative of $x^{\mathsf T}(Q_i-Q_j)x$ is therefore
$-4\sqrt2hU_{ij}$. Formula (9) follows in the same way from the outward
derivative $-2d$ of $2(m_i-m_j)\cdot x$. The remaining terms are obtained
by differentiating the squared norms of the moments.

Thus the local $W_2$ contribution in (10) is nonpositive, including on
remote portions of the original facets. A bound between the complete
quadratic forms, with rotation directions treated correctly, is a separate
question; the sign of this one term does not prove it. The swaps (3) are
small in Gaussian measure but are not a fixed smooth infinitesimal normal
deformation. Such tangent tests cannot establish the unrestricted bound
refuted by (7).

## Consequences for the proposed bridge

Define

$$
H_\rho(A)=\sum_{k\ge2}\rho^{k-2}(W_k(A)-W_k(T)).
$$

For a fixed balanced $A$, Parseval justifies absolute convergence and
$H_\rho(A)\to W_2(A)-W_2(T)$ as $\rho\downarrow0$. For example the
absolute contribution from $k\ge3$ is at most $3\rho/2$. Therefore no
finite $C$ and $\rho_0>0$ can satisfy
$H_\rho(A)\le C\delta(A)$ simultaneously for every balanced $A$ with
$\delta(A)>0$ and every $0<\rho\le\rho_0$: taking the limit at each fixed
$A_R$ would contradict (7).

There is a stronger obstruction to using $\delta$ as a uniform lower-bound
scale. Changing the label only on a set of measure $2p_R$ changes the
agreement indicator only when at least one member of the Gaussian pair
lies in that set. At every fixed $0<\rho<1$,

$$
|S_\rho(T)-S_\rho(A_R)|\le4p_R,\qquad
\delta(A_R)\ge4\sqrt2c\,a_R>16\sqrt2c\,Rp_R.
$$

Consequently

$$
\frac{|S_\rho(T)-S_\rho(A_R)|}{\rho\delta(A_R)}
\le\frac1{4\sqrt2c\,\rho R}\longrightarrow0.
\tag{8}
$$

Even at a single fixed positive correlation there is no positive constant
$\kappa(\rho)$ such that
$S_\rho(T)-S_\rho(A)\ge\kappa(\rho)\rho\delta(A)$ for every balanced
$A$ with positive deficit. In particular the suggested
$\rho\delta(A)(1-C\rho)$ bound cannot hold uniformly on a nonzero
interval. The estimate does not require a sign for the noise-stability
difference. The next theorem establishes that it is strictly positive for
every member of this family.

## Proper reflected swaps strictly lose full noise stability

Let $n\ge3$ and extend the tetrahedral cells cylindrically to $\mathbb R^n$.
Let $g$ be the orthogonal reflection above, acting trivially on unused
coordinates. Up to null boundaries, $gT_1=T_2$, and $T_1$ lies in
$\{x_N>0\}$. For any measurable $B\subset T_1$, put $C=T_1\setminus B$ and

$$
A_1=C\cup gB,\qquad A_2=B\cup gC,\qquad A_3=T_3,\quad A_4=T_4.
$$

**Reflected-swap theorem.** This partition is balanced. For every
$0<\rho<1$,

$$
S_\rho(T)-S_\rho(A)
=4\int_{B\times C}\bigl[K_\rho(x,y)-K_\rho(x,gy)\bigr]
\,d\gamma_n(x)d\gamma_n(y)\ge0. \tag{11}
$$

Equality holds exactly when $\gamma_n(B)=0$ or $\gamma_n(C)=0$.
In particular every proper swap, meaning $0<\gamma_n(B)<1/4$, has
$S_\rho(A)<S_\rho(T)$ for every $0<\rho<1$.

Here the noise kernel is the density of the correlated Gaussian pair
relative to $\gamma_n\otimes\gamma_n$:

$$
K_\rho(x,y)=(1-\rho^2)^{-n/2}
\exp\!\left(
\frac{2\rho\,x\cdot y-\rho^2(\|x\|^2+\|y\|^2)}{2(1-\rho^2)}
\right). \tag{12}
$$

This follows by dividing the joint Gaussian density, with covariance
$\begin{psmallmatrix}I&\rho I\\\rho I&I\end{psmallmatrix}$, by the two
standard Gaussian densities. In particular the kernel is positive,
symmetric, and satisfies $K_\rho(gx,gy)=K_\rho(x,y)$.

**Proof.** Reflection preserves Gaussian measure, so
$\gamma_n(A_1)=\gamma_n(C)+\gamma_n(B)=1/4$, and likewise for $A_2$.
Write $I_\rho(E,F)=\int_{E\times F}K_\rho\,d\gamma_n d\gamma_n$.
The four self-interactions of $B,C,gB,gC$ cancel between the two
partitions, and the contributions of cells 3 and 4 are unchanged. The
remaining difference is

$$
2\bigl[I_\rho(B,C)+I_\rho(gB,gC)
       -I_\rho(C,gB)-I_\rho(B,gC)\bigr]
=4\bigl[I_\rho(B,C)-I_\rho(B,gC)\bigr].
$$

The equality uses reflection invariance and symmetry. Changing variables
in the last integral proves (11), including its factor four. All terms
are finite since the full kernel integrates to one. Since reflection
changes $x\cdot y$ by $2x_Ny_N$ and preserves norms, (12) gives

$$
\frac{K_\rho(x,y)}{K_\rho(x,gy)}
=\exp\!\left(\frac{2\rho x_Ny_N}{1-\rho^2}\right)>1
\quad\text{for almost every }(x,y)\in B\times C. \tag{13}
$$

If both sets have positive measure, their product has positive measure,
so (11) is strictly positive. If either is null, the partition agrees
with $T$ up to null sets or globally interchanges labels 1 and 2, and
the difference is zero. This proves the equality characterization.

Taking $B=gE_R$ gives exactly (3), so all the tail counterexamples are
strictly suboptimal at every positive correlation. Combining (8) and
(11) yields, for each fixed $0<\rho<1$,

$$
0<\frac{S_\rho(T)-S_\rho(A_R)}{\rho\delta(A_R)}
\le\frac1{4\sqrt2c\,\rho R}\longrightarrow0.
$$

Thus strict full-noise loss coexists with failure of any positive uniform
lower bound on this deficit scale. No Hermite truncation is used in the
kernel proof.

The closest checked rearrangement mechanism is the two-point
polarization argument in [Kindler, O'Donnell and Witmer, *Remarks on the
Most Informative Function Conjecture at fixed mean*, July 1, 2015,
Theorem 4.2 and Lemmas 4.4--4.7](https://www.cs.cmu.edu/~odonnell/papers/continuous-most-informative-function.pdf).
That theorem uses uniform sphere measure, a bounded nondecreasing
inner-product kernel and a convex uniformly continuous outer function.
Its reflection identities and sign comparison motivate (11)--(13);
the Gaussian density and balanced four-cell identity are proved here
directly. This comparison makes no novelty claim and supplies no
reduction of arbitrary four-cell competitors to reflected swaps.

## Balanced three-cell cycles and their cross-interactions

Let $B_i\subset T_i$ for $i=1,2,3$ have the same Gaussian measure
$p$, and put $C_i=T_i\setminus B_i$. Indices in this section are
cyclic modulo three. The forward cycle assigns $B_i$ to label $i+1$:
$A_i^+=C_i\cup B_{i-1}$, with $A_4^+=T_4$. The reverse cycle uses
$A_i^-=C_i\cup B_{i+1}$. Both partitions are exactly balanced.
Write $D_\pm=S_\rho(T)-S_\rho(A^\pm)$ and retain the integral
$I_\rho$ from the reflection proof. Cancellation of self-interactions
gives the exact identity

$$
D_+=2\sum_{i=1}^3
  \bigl[I_\rho(B_i,C_i)-I_\rho(B_i,C_{i+1})\bigr].
\tag{CY1}
$$

Equivalently, with $f_i=1_{B_i}$, the Ornstein--Uhlenbeck operator
$T_\rho$, and $\|f\|_\rho^2=\langle f,T_\rho f\rangle$,

$$
D_+=2\sum_i\int_{B_i}T_\rho(1_{T_i}-1_{T_{i+1}})\,d\gamma_n
       -\sum_{i<j}\|f_i-f_j\|_\rho^2.
\tag{CY2}
$$

Indeed expanding $C_i=T_i\setminus B_i$ in (CY1) leaves
$-2\sum_i I_\rho(B_i,B_i)+2\sum_{i<j}I_\rho(B_i,B_j)$,
which is the last term in (CY2). The kernel norm is nonnegative
for $0<\rho<1$ by the Hermite expansion. These cross-interactions
must be controlled along with the linear score loss.

Let $D_{ij}$ be the loss for swapping only $B_i,B_j$ between their
two original cells, leaving all other labels in $T$ unchanged.
No reflection relation between the two sets is assumed in this
definition. Expanding the same six blocks gives

$$
D_++D_-=D_{12}+D_{23}+D_{31}.
\tag{CY3}
$$

For example, the $ij$ summand on the right is
$2[I_\rho(B_i,C_i)+I_\rho(B_j,C_j)
-I_\rho(B_i,C_j)-I_\rho(B_j,C_i)]$.
Thus (CY3) retains both orientations of every interaction;
it does not assert that either cycle has a sign for arbitrary seeds.

**Equivariant cyclic-swap theorem.** Suppose the three seeds are
equivariant under the tetrahedral orthogonal symmetries permuting
labels $1,2,3$ and fixing label $4$: $hB_i=B_{h(i)}$ up to null
sets for every such symmetry $h$. For every $0<\rho<1$,

$$
D_+=D_-=\frac12(D_{12}+D_{23}+D_{31})
 =6\int_{B_1\times C_1}
       [K_\rho(x,y)-K_\rho(x,g_{12}y)]\,d\gamma_n(x)d\gamma_n(y).
\tag{CY4}
$$

Here $g_{12}$ is the tetrahedral reflection interchanging labels
$1,2$ and fixing $3,4$. Both cycles strictly lose full noise
stability exactly when $0<p<1/4$. The null cycle and the complete
global three-label cycle give equality.

For the proof, any transposition in this symmetry group conjugates
the forward cycle to the reverse cycle, so orthogonal invariance
gives $D_+=D_-$. Every pair swap now has exactly the reflected
form of (11), because $g_{ij}B_i=B_j$ and
$g_{ij}C_i=C_j$. All three pair losses agree by the cyclic
symmetry. Equation (CY3) and the reflection identity give (CY4),
including the factor six and both equality cases. This argument
compares each partition directly with $T$ and does not apply a
second reflection theorem to an already modified pair.

This is a nonempty family of measurable partial cycles. To construct
it, choose $B_1\subset T_1$ invariant under the reflection
interchanging labels $2,3$ and fixing $1,4$, and take its orbit under
the order-three symmetry $(123)$. For example
$B_i=T_i\cap\{\|x\|>R\}$ gives proper cycles for every $R>0$.
The argument uses the full kernel, with no Hermite-tail truncation.

The next coupled test is (CY1) for equal-mass seeds without this
permutation equivariance, including order-three orbit seeds that
are not invariant under the cell-stabilizing reflection. For such
seeds (CY3) alone gives no sign for a chosen orientation. A general
cycle-loss inequality, or a reduction of arbitrary competitors to
these cycles, remains open. The checked polarization theorem cited
above is a two-point sphere theorem; it is not imported as a
three-cell Gaussian rearrangement theorem here.

## Effect on the remaining proof route

The [full-rank first-Hermite stability obligation](BALANCED_AFFINE_CENTROID_OPTIMIZATION.md#boundary-exclusion-and-the-remaining-global-obligation)
is unchanged. A fixed-partition Hermite expansion does not transfer its
eventual solution to one interval uniform over all balanced partitions.
The affine replacement used for $W_1$ preserves neither $W_2$ nor noise
stability, so an affine-only test cannot remove this obstruction.

A uniform small-correlation theorem may still hold. Equations (11)--(13)
settle the full-noise comparison for every single reflected swap of a
tetrahedral cell pair, including the remote tails. A remaining full-kernel
approach could seek balance-preserving rearrangements that improve an
arbitrary competitor toward $T$. No such reduction or sequence of
rearrangements has been established. Comparisons for different cell pairs
cannot simply be iterated using (11): its proof assumes the initial pair
is the reflected tetrahedral pair. A quantitative gap proportional to
the unrestricted $\delta$ cannot be the conclusion, by (8).

The other task is the high-value full-rank first-Hermite stability gap,
using the [five-dimensional shape Hessian](BALANCED_AFFINE_CENTROID_OPTIMIZATION.md#five-shape-directions-and-the-averaged-trace).
That finite-dimensional obligation and a uniform full-noise comparison
remain distinct; resolving the former alone does not resolve the latter.
The radius-changing subspace there targets equal centroid norms. On
the full-noise side, (CY4) settles the permutation-equivariant three-cell
cycles, while (CY1)--(CY3) retain the open interactions for general seeds.

## Reproducible controls and assurance

[The checker](check_relative_tail_swaps.py) independently assembles the full
first- and second-moment tensors from product truncated-normal integrals,
then compares their mass changes with (6), checks box containment, and
checks a reflected positive-$U$ control that reverses the linear $W_2$
term. Run from the repository root:

```sh
python3 research/plurality-stability/check_relative_tail_swaps.py
```

The internal run uses Python 3.12.3 and its standard library at
$R=1,2,4,8$. Floating-point controls check normalization and implementation;
the proof for unbounded $R$ is (2)--(8). The full-noise theorem (11)--(13)
is a written kernel argument for arbitrary measurable $B$; the tensor
checker does not verify it. There is no interval certificate, independent
replay, formal verification or external review.
The checker also compares all 21 symmetric six-block interaction
coefficients exactly over the rationals for (CY1), (CY3), and the
permutation-equivariant reduction in (CY4). A separate coefficient
control rejects equality of the two orientations without equivariance.
This finite algebra check is supplementary to the written measurable-set
proof; it performs no Gaussian kernel quadrature and proves no general
cycle sign.
