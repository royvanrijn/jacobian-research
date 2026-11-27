# Radial mass concentrations at a simplex tie

For a balanced planar three-cell partition, let $p_i(r)$ be its
three circular masses and put $b_i(r)=p_i(r)-1/3$. Write
$d\mu(r)=r e^{-r^2/2}dr$,
$H^2=\sum_i\|b_i\|_{L^2(\mu)}^2$, and
$M=\max_i\|b_i\|_\infty$. A nonsimplex partition reaching at
least the simplex value at any fixed positive correlation needs a
positive fraction of its squared radial mass error inside a finite
centered disk. Along any such aligned sequence of maximizers near
correlation one, $H$ tends to zero and $M/H$ tends to infinity.
Thus the radial mass error in such a sequence cannot be uniformly
distributed in amplitude. The statements concern arbitrary
measurable partitions, with no hypotheses on their circular
label sets or on derivatives of their radial mass functions.

## The radial residual

Use $P_\rho$, $W_\rho$ and the angular convolution $U_\kappa$
from the [polar deficit identity](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#exact-separation-of-radial-masses-and-angular-changes).
Set $\sigma=\sqrt{1-\rho^2}$. For every pair of radii, circular
rearrangement bounds the pairing of each label by that of two
aligned centered arcs of the same lengths. This gives an upper
bound on $S_\rho(A)$ even when the three centered arcs overlap.
No disjointness of these rearranged arcs is used.
The two-set rearrangement input is exactly the one in the proof
of Theorem A.1. [Heilman, arXiv:2608.00333v1](https://arxiv.org/pdf/2608.00333v1)

Let $I_p$ be the centered arc of angular mass $p$, and put
$u_p=1_{I_p}-p$ and $F_\kappa(p)=\langle u_p,U_\kappa u_p\rangle$.
Write $C_{\rm arc}$ for the right side defining $\mathcal C_\rho$
in (DG2), with the centered arc indicators instead of the
original centered indicators. In particular $C_{\rm arc}\ge0$.
Integrating the endpoints of one arc gives
$F_\kappa''(p)=2[k_\kappa(2\pi p)-1]$, where $k_\kappa$
is the normalized kernel. Taylor's formula about $1/3$ yields

\[
 3F_\kappa(1/3)-\sum_i F_\kappa(p_i)
 \ge \sum_i b_i^2
     -\sum_i b_i^2\int_0^1
          2(1-t)k_\kappa(2\pi(1/3+t b_i))\,dt .
\]

The first derivatives cancel since $\sum_i b_i(r)=0$. With
$q=1/3+t b_i$ one has $\min(q,1-q)\ge(1-t)/3$.
The bound $1-\cos(2\pi q)\ge8\min(q,1-q)^2$ gives

\[
 \int_0^1 2(1-t)k_\kappa(2\pi q)\,dt
 \le \frac{9e^\kappa(1-e^{-8\kappa/9})}{8\kappa I_0(\kappa)}
 \le \frac{C}{\sqrt{1+\kappa}} .
\]

The expression at zero is taken by continuity. For the final
bound use $I_0(\kappa)\ge c e^\kappa/\sqrt{1+\kappa}$, by
integrating at angles $|\theta|\le (1+\kappa)^{-1/2}$.
Conditional on radius $r$, put $x=\rho r/\sigma$. Then
$\kappa$ has the distribution $x|x(1,0)+Z|$ for a standard
planar Gaussian $Z$. For $x\ge1$, split at
$|x(1,0)+Z|\ge x/2$. This bounds the row average of
$C/\sqrt{1+\kappa}$ by $C'/(1+x)$, including $0\le x\le1$
after increasing $C'$.

If instead $p_i(r)\in[\delta,1-\delta]$ at a given radius
with $0<\delta\le1/3$, one can bound the Taylor remainder
there by a smaller row estimate. Write $c_\delta=\cos(2\pi\delta)$
and $\Phi$ for the scalar Gaussian distribution function.
The conditional row integral is explicitly

\[
 R_\delta(x)=
 \int P_\rho(r,s)\frac{e^{\kappa c_\delta}}{I_0(\kappa)}\,d\mu(s)
 =e^{-x^2/2}+c_\delta x\sqrt{2\pi}\,
              e^{-(1-c_\delta^2)x^2/2}\Phi(c_\delta x).
 \tag{RC1}
\]

This follows by putting $s=\sigma u$ and completing the square
in $\int_0^\infty u e^{-u^2/2+c_\delta xu}du$. It decays as
$C_\delta(1+x)e^{-c_\delta' x^2}$ for some $c_\delta'>0$.

By the rearrangement comparison and (DG3), taking row averages
in the Taylor bound proves

\[
 S_\rho(T)-S_\rho(A)\ge
 \sum_i\langle b_i,(I-P_\rho)b_i\rangle_\mu+C_{\rm arc}
           -\sum_i\int m_i(r)b_i(r)^2\,d\mu(r),
 \tag{RC2}
\]

where $m_i(r)$ may be chosen as $C'/(1+\rho r/\sigma)$ everywhere,
or $R_\delta(\rho r/\sigma)$ at any radius where the corresponding
width is in $[\delta,1-\delta]$. For the diagonal estimate one
can also keep $\mathcal C_\rho(A)$ of the original partition,
by applying circular rearrangement only to each self pairing.
In (RC2) the rearranged cross term gives an estimate from the
radial masses alone.

The balanced $b_i$ have mean zero and, considered as Gaussian
radial functions, no odd Hermite degrees. Therefore the
Dirichlet term in (RC2) is at least $(1-\rho^2)H^2$. For
each fixed $\rho>0$ choose $R$ such that
$C'/(1+\rho R/\sigma)\le(1-\rho^2)/2$. If $S_\rho(A)\ge S_\rho(T)$,
(RC2) gives

\[
 \sum_i\int_0^R b_i^2\,d\mu\ge c(1-\rho^2)H^2 .
 \tag{RC3}
\]

Here $c>0$ is absolute. In particular, equality of all circular
masses on $[0,R]$ already forces equality on every circle
almost everywhere for such a competitor. The equal-circle
rigidity result [in the global analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#circular-balance-and-compatibility-of-first-angular-moments)
then gives the simplex.

## Uniform comparison without radial amplitude concentration

For every fixed $K>0$ there are $\rho_0<1$ and $H_0>0$ such
that the simplex is uniquely optimal for planar balanced
partitions satisfying $M\le KH$ and $H\le H_0$, at every
$\rho_0<\rho<1$. The conclusion uses the cross-radius term in
(RC2). A bound from the diagonal correction alone would
miss this term.

On a fixed annulus $0<a<b<\infty$, take two radii
$r,s\in[a,b]$ with $|r-s|\le\sigma$. Their $\kappa$ is
comparable to $\sigma^{-2}$. The centered arcs are nested.
For one label, put $d=|p_i(r)-p_i(s)|$. Then
$u_{p_i(r)}-u_{p_i(s)}$ is, up to sign, the indicator of the
two difference strips minus $d$. Each strip has angle length
$\pi d$. In the quadratic form of $U_\kappa$ the two strip
self pairings are nonnegative. On pairs in one strip with
angular separation $O(\sigma)$ the normalized kernel is
at least $c_{a,b}/\sigma$. For this use
$I_0(\kappa)\le C e^\kappa/\sqrt{1+\kappa}$, obtained from
$1-\cos\theta\ge 2\theta^2/\pi^2$ on $[-\pi,\pi]$.
The area of such pairs is at least a constant times
$\min(d^2,\sigma d)$. The zero-mode subtraction costs $d^2$.
Consequently, for $\sigma,M$ small enough depending on the
annulus, the cross quadratic form is bounded below by
$c_{a,b} d^2/\max(\sigma,2M)$, since $d\le2M$.

The same annular radius pairs have $W_\rho$ density at least
$c_{a,b}/\sigma$ against $dr\,ds$. This follows from the
lower bound on $I_0$ and the Gaussian exponent
$-(r-s)^2/(2\sigma^2)-rs/(1+\rho)$. For any measurable $f$ on
the annulus one has

\[
 \int_a^b\int_{\substack{a\le s\le b\\|r-s|\le\sigma}}
           |f(r)-f(s)|^2\,ds\,dr
     \ge c_{a,b}\sigma^3\int_a^b|f-\bar f|^2\,dr ,
\]

where $\bar f$ is the interval mean. Partition $[a,b]$ into
equal intervals of length between $\sigma/4$ and $\sigma/2$
for small $\sigma$. Pairs in one interval give its variance
times twice its length. Pairs in neighboring intervals give
at least the squared difference of means times the squared
length. The discrete interval Poincaré inequality for means,
with constant at most the square of the number of intervals,
then gives the bound. It follows that

\[
 C_{\rm arc}\ge\frac{c_{a,b}\sigma^2}{\max(\sigma,2M)}
             \sum_i\int_a^b|b_i-\overline b_i|^2\,dr .
 \tag{RC4}
\]

Under $M\le KH$, the normalized profiles $b_i/H$ have summed
squared norm one, mean zero and sup norm at most $K$. Choose
the annulus so the $\mu$-mass outside is sufficiently small
depending on $K$. The summed squared norms inside are then
bounded below, and their means on the annulus are small by
the global zero mean. Hence the displayed sum of unweighted
variances is at least $c_K H^2$, by boundedness of the
radial density on that annulus. For $H$ small, all widths lie
in $[1/6,5/6]$. The exponential multiplier (RC1) integrates
against $\mu$ to $O(\sigma^2)$ for $\rho\ge1/2$, so the
adverse term in (RC2) is $O(\sigma^2 M^2)=O_K(\sigma^2 H^2)$.
It is strictly dominated by (RC4) for $\sigma,H$ sufficiently
small. The case $H=0$ is the equal-circle comparison.

For an aligned sequence of actual planar maximizers near one,
the [compactness result](THREE_CELL_NEAR_ONE.md#compactness-and-the-sharp-perimeter-lower-bound)
implies $H\to0$ by contraction onto radial profiles and
boundedness of the indicators. For every nonsimplex member
$H>0$ by equal-circle rigidity. The comparison just proved
therefore forces $M/H\to\infty$.
