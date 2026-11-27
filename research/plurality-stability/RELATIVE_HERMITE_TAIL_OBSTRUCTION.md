# Balanced tail swaps obstruct relative Hermite control

There is no finite constant $C$ such that
$W_2(A)-W_2(T)\le C\delta(A)$ for every balanced four-cell Gaussian
partition with $\delta(A)=W_1(T)-W_1(A)>0$. This fails already in dimension
three, arbitrarily close to $T$ in Gaussian measure and in centroid Gram
space. Consequently the proposed uniform higher-Hermite bound also fails.
The first-Hermite theorem and a uniform interval of positive-correlation
optimality remain open; neither is contradicted by this obstruction.

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
interval. This conclusion uses no assumed sign of the noise-stability
difference and gives no partition beating $T$.

## Effect on the remaining proof route

The [full-rank first-Hermite stability obligation](BALANCED_AFFINE_CENTROID_OPTIMIZATION.md#boundary-exclusion-and-the-remaining-global-obligation)
is unchanged. A fixed-partition Hermite expansion does not transfer its
eventual solution to one interval uniform over all balanced partitions.
The affine replacement used for $W_1$ preserves neither $W_2$ nor noise
stability, so an affine-only test cannot remove this obstruction.

A uniform small-correlation theorem may still hold. A next bridge must
either compare the full Ornstein--Uhlenbeck functional on the tails, where
individual Hermite degrees need not approximate it uniformly, or prove a
reduction to a restricted competitor class that excludes the swaps above.
Any estimate on such a class needs its hypotheses and the reduction proved
separately. A quantitative gap proportional to the unrestricted $\delta$
cannot be its conclusion, by (8).

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
the proof for unbounded $R$ is (2)--(8). There is no interval certificate,
independent replay, formal verification or external review.
