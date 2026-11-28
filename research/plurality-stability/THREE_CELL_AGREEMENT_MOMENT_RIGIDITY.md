# Third and fourth agreement moments at a balanced simplex tie

Let \(A\) be a balanced planar three-cell maximizer that ties the
regular centered simplex \(T\) at \(0<r<1\). A single contracted
quartic moment of the agreement event separates any nonsimplex
candidate at a first scalar contact. Additional cubic and quartic
comparisons use the case where distinct scalar correlations at
which \(T\) is maximizing accumulate at \(r\) and \(A\) has the
same first parameter derivatives as \(T\). The latter means its
centered prices are zero and \(M_A=M_T\). It makes the first and
second agreement moments of the Gaussian midpoint and difference
equal to the simplex's by (VV13).
[Joint covariance and volume response](THREE_CELL_VOLUME_VALUE_DERIVATIVES.md#joint-covariance-and-volume-response).

Two selected higher moments suffice to distinguish every nonsimplex
tie satisfying these hypotheses. If the quartic moment also matches,
the remaining cubic moment yields a strict advantage over the local
simplex branch for certain nearby quotas. The quota comparison uses
only scalar correlations and stays in the active plane.

## A cubic and a quartic comparison

Let \((X,Y)\) have standard Gaussian marginals and cross-covariance
\(rI\) in dimension \(d=2\). Put
\[
 Z=\frac{X+Y}{\sqrt{2(1+r)}},\qquad
 W=\frac{X-Y}{\sqrt{2(1-r)}},\qquad
 K=\frac{|Z|^2-d}{2(1+r)}-\frac{|W|^2-d}{2(1-r)}.
\]
For a partition \(F\), let \(E_F\) be its same-label event at \((X,Y)\).
For now retain the hypothesis of matching first parameter derivatives.
Define the unnormalized moment comparisons
\[
 \mu=\mathbb E[ZK1_{E_A}],\qquad
 \kappa=\mathbb E[K^2(1_{E_T}-1_{E_A})].                    \tag{AG1}
\]
Since the two same-label probabilities agree, equality of the
unnormalized moments is equivalent to equality of their conditional
versions. The simplex's cubic moment is zero because rotation by
\(2\pi/3\) preserves its same-label event and \(K\).

Use the interface conventions of the
[maximizer variations](THREE_CELL_GLOBAL_OBSTRUCTIONS.md). Write
\(b(x)=x\cdot N(x)\), \(s_v(x)=v\cdot N(x)\), and
\(g_{ij}=\partial_r[T_r(1_{A_i}-1_{A_j})]\) on its corresponding
interface. Let \(v\) be the translation whose mass flux equals that
of \(b\); it is unique by full centroid rank. Put \(h=b-s_v\), a
balanced speed. With the weighted normal norm \(\|b\|_a\) and the
balanced second variation form \(Q_0\), one has \(\kappa\ge0\) and
\[
 \frac{\|b\|_a^2}{r(1+r)}\le
 \sqrt{\frac{\kappa Q_0(h)}2}
       +\frac{|v|\,|\mu|}{\sqrt{2(1+r)}} .                  \tag{AG2}
\]
In particular \(\kappa=0\) and \(\mu=0\) imply that \(A\) is the
centered regular simplex up to rotation and relabeling. More
quantitatively, on each compact positive scalar correlation band
there is \(c>0\) such that a nonsimplex tie satisfying these
hypotheses has \(\sqrt\kappa+|\mu|\ge c\).

To prove (AG2), first keep the cells fixed and differentiate the
joint Gaussian density in the raw midpoint and difference variables
\((X+Y)/\sqrt2,(X-Y)/\sqrt2\). The logarithmic derivative is \(K\),
and its derivative at fixed raw variables is
\[
 K'=\frac{d/2-|Z|^2}{(1+r)^2}+\frac{d/2-|W|^2}{(1-r)^2}.
\]
The second density derivative divided by the density is \(K^2+K'\).
The \(K'\) terms cancel in the difference between \(T\) and \(A\)
because their same-label probabilities and second agreement moments
match. Thus for \(D(\rho)=S_\rho(T)-S_\rho(A)\), we have
\(D''(r)=\kappa\). Also \(D(r)=D'(r)=0\); nonnegativity along the
scalar sequence of simplex optimality yields \(\kappa\ge0\).

Next translate the partition by \(tw\), with \(w\) fixed. The derivative
in \(t\) of the density factor after moving the translation from the
cells to the joint density is \(-\sqrt{2/(1+r)}\, Z\cdot w\). Its
mixed correlation derivative divided by the unperturbed density is
\(\sqrt{2/(1+r)}(Z\cdot w)((1+r)^{-1}-K)\). Since
\(\mathbb E[Z1_{E_A}]=0\), this gives
\[
 \left.\partial_\rho\partial_t S_\rho(A+tw)\right|_{\rho=r,t=0}
          =-\sqrt{\frac{2}{1+r}}\,w\cdot\mu
          =2\langle g,s_w\rangle .                         \tag{AG3}
\]
The last equality is the first interface variation for translation,
differentiated at fixed cells. It can be obtained by Gaussian
divergence on compact spatial cutoffs. Passing to the limit uses the
polynomially weighted interface bounds of (G3) and (G6), and compact
correlation subintervals provide an integrable Gaussian majorant for
the density derivatives.

The tangential contact inequality (G10) on \(h\) gives
\(|\langle g,h\rangle|\le\sqrt{\kappa Q_0(h)/2}\). The cutoff and
mass-bump approximation is exactly the balanced dilation argument of
(G21). The one-sided dilation inequality (G6) bounds
\(\|b\|_a^2/[r(1+r)]\) by \(-\langle g,b\rangle\). Now split
\(b=h+s_v\) and use (AG3). This proves (AG2). If its right side
vanishes, \(b=0\) on the interfaces almost everywhere, since the
normal score derivative \(a\) is positive. Distributionally each
indicator has zero radial derivative away from the origin, so the
partition is conical about zero. The conical maximizer rigidity
theorem forces the simplex.
[Balanced dilation, compactness and conical rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#balanced-dilation-at-an-endpoint-tie).

For uniformity, (G25) bounds \(\|b\|_a^2\) below over nonsimplex
balanced maximizers on compact positive scalar bands. The inverse
translation mass map is uniformly bounded by full rank and strong
maximizer compactness, and the dilation mass flux is bounded by
Gaussian second moments. Hence \(|v|\) is bounded. The form satisfies
\(0\le Q_0(h)\le\|h\|_a^2\) since the Mehler interaction on the
signed interface measures is nonnegative, by exhaustion and
positivity of the Mehler operator. The last norm is uniformly
bounded by (G24) for the fields \(x-v\). These estimates yield the
asserted lower bound on \(\sqrt\kappa+|\mu|\).
[Uniform conical gap, G24–G25](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-uniform-gap-from-conical-geometry).

When the cubic vector vanishes, the curvature has a sharper lower
bound without any cellwise radial moment hypothesis. To see this,
pair the dilation identity in (G6) with a translation speed \(s_u\).
Self-adjointness of the interface kernel and the translation identity
in (G2) give
\[
 \langle g,s_u\rangle=-\frac{\langle b,s_u\rangle_a}{r(1+r)}.
\]
By (AG3), \(\mu=0\) therefore makes \(b\) weighted-orthogonal to all
translation speeds. For the balancing translation \(w=s_v\), the
quantity \(B_h=(1+r)\|h\|_a^2+\langle w,h\rangle_a\) of (G21)
then satisfies
\[
 B_h=(1+r)\|b\|_a^2+r\|w\|_a^2.
\]
If \(A\) is nonsimplex, \(B_h>0\) by conical rigidity. Applying
(G22) yields the explicit bound
\[
          \kappa\ge\frac{8}{1-r^2}\|b\|_a^2 .                \tag{AG6}
\]
Thus a vanishing cubic response forces a strictly positive scalar
deficit curvature. This holds even when the balancing translation
is nonzero.

## A single quartic comparison without matched first responses

The following bound applies at every tie, with no price or moment
matching assumed. Put \(H=2d-|Z|^2-|W|^2\) and
\(d_A=S'_r(A)-S'_r(T)\). Then
\[
 \mathbb E[HK(1_{E_T}-1_{E_A})]
     \ge 2d_A+\frac{2}{r(1+r)}\|b\|_a^2 .                     \tag{AG5}
\]
In particular, if scalar simplex optimality accumulates from below,
then \(d_A\ge0\). The right side is strictly positive for a
nonsimplex tie, with a uniform positive lower bound on a compact
positive scalar band by (G25). Thus matching all agreement moments
of degree at most four is impossible at any nonsimplex first scalar
contact, and (AG5) identifies one explicit polynomial witness.

In fact \(H\) is the density score for the dilation \(A\mapsto e^t A\).
The change of variables in the joint integral scales both raw pair
coordinates by \(e^t\). Its Jacobian is \(e^{2dt}\), and the Gaussian
exponent is multiplied by \(e^{2t}\). At fixed raw variables the
correlation derivative of \(H\) is
\[
 \partial_r H=\frac{|Z|^2}{1+r}-\frac{|W|^2}{1-r}
             =2K+d\big((1+r)^{-1}-(1-r)^{-1}\big).
\]
Therefore the mixed dilation and correlation derivative is
\(\mathbb E[(HK+\partial_r H)1_{E_A}]=2\langle g,b\rangle\),
using the same first variation and polynomial cutoff bounds as in
(AG3). This mixed derivative is identically zero for the conical
simplex. Subtract the two formulas. The constant in \(\partial_r H\)
cancels at equal stability, and \(\mathbb E[K1_{E_F}]=S'_r(F)\), so
\[
 \mathbb E[HK(1_{E_T}-1_{E_A})]=2d_A-2\langle g,b\rangle.
\]
The dilation inequality (G6) yields (AG5).

## A quota comparison below the tie

Let \(m_*= (1/3,1/3,1/3)\), and let \(m(t,w)\) be the masses of
the translated partition \(A+tw\). Then
\(m(t,w)=m_*+tC_Aw+O(t^2)\), where \(C_Aw\) is the first
translation mass change \((-z_i\cdot w)_i\). Let
\(V_\rho^T(m)\) be the orbit-local simplex value of (LM1), with
quadratic mass coefficient \(\Lambda_\rho\).
For fixed \(w\) and \(t\downarrow0\) one has the exact limit
\[
\begin{split}
 &\lim_{t\downarrow0}\frac{S_{r-t}(A+tw)-V^T_{r-t}(m(t,w))}{t^2}\\
 &\hspace{1em}=-\frac{\kappa}{2}+\sqrt{\frac{2}{1+r}}\,\mu\cdot w
               +(1-r)w^{\mathsf T}M_Aw-\Lambda_r|C_Aw|^2 .   \tag{AG4}
\end{split}
\]
Indeed the fixed-shape deficit at \(r-t\) has quadratic
coefficient \(\kappa/2\). The mixed term for translating the cells
by \(tw\) with correlation shift \(-t\) is the negative of (AG3).
The prices are zero, so the translation Hessian (AM3) gives the
stability coefficient \((1-r)w^{\mathsf T}M_Aw\). Finally (LM1)
gives the quadratic comparison for the simplex at those masses,
uniformly at correlations tending to \(r\).
[Translation Hessian, AM3](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#balanced-translations-with-a-matrix-correlation).
[Measurable local simplex mass profile, LM1](THREE_CELL_LOCAL_MASS_PROFILE.md#exact-expansion-over-measurable-competitors).

If \(A\) is nonsimplex and \(\kappa=0\), then \(\mu\ne0\) by (AG2).
Choose \(w\) to be a sufficiently small positive multiple of \(\mu\).
The linear term in \(w\) in (AG4) dominates the quadratic terms,
giving a strict excess over the orbit-local simplex value at the
translated quotas for all sufficiently small \(t>0\). If the
balanced simplex is the unique global orbit at \(r-t\), varying the
quotas along the segment from \(m_*\) to \(m(t,w)\) also forces
separated global competitors at some intermediate quota. This follows
from strong maximizer compactness and the excluded-annulus argument
used for (LM8): close to the base tie, no global maximizer lies in a
fixed annulus outside the local simplex neighborhood. The sets of
segment parameters with a near and with a distant maximizer are
closed, cover the segment, and are each nonempty by the two endpoints.
[Coexistence by varying quotas](THREE_CELL_LOCAL_MASS_PROFILE.md#consequences-at-a-first-balanced-tie).

These comparisons apply to the higher-order responses of a balanced
tie. For \(\kappa>0\), formula (AG4) gives a strict quota comparison
only when its right side is positive for some \(w\). The moment
bounds here leave open the exclusion of that general scalar tie.
