# Continuation at positive Gaussian correlation

The [fixed-cell normal-flux proof](SMALL_POSITIVE_CORRELATION_INTERVAL.md#extension-to-every-fixed-cell-count)
establishes a qualitative interval of simplex optimality for each fixed
cell count. This note proves a local comparison at an arbitrary positive
correlation **under an explicit interface-kernel hypothesis**, and sharpens
the obstruction to continuing global optimality. For arbitrary fixed cell
count the hypothesis is proved for sufficiently small correlations. The
three-cell case is now proved throughout $(0,1)$ by the separate symmetry
and reflected-kernel argument linked below.
All proofs and source audits here are internal. Independent mathematical
review and formal verification remain unperformed.

The [symmetry-sector analysis](SIMPLEX_INTERFACE_SYMMETRY.md) now gives
the signed $S_k$ facet decomposition. For three cells it proves an
explicit full-interface gap modulo rotations for every positive
correlation, independent of ambient dimension, by reflected OU kernel
comparisons and an elementary bound at correlation $1/2$. The higher-cell
blocks remain open. These local results do not exclude a distant tied
global maximizer; for three cells such a tie is now necessary at any
first global continuation endpoint below one.

Fix $k\ge2$ and $n=k-1$. Let $T=(T_1,\ldots,T_k)$ be the central regular
simplex partition of $\mathbb R^n$, with cell moments $b_i^*$ and common
moment difference length $\ell>0$. Write $B_*$ for its moment matrix and
$F_k=\|B_*\|_F^2$, so $B_*^TB_*=(F_k/n)I_n$. Put
$n_{ij}=(b_i^*-b_j^*)/\ell$, pointing into $T_i$,
$H_{ij}=n_{ij}^\perp$, $F_{ij}=T_i\cap T_j$, and
$d\nu_{ij}=\phi_n d\mathcal H^{n-1}$ on the pair plane. We use the
balanced facet space $\mathcal H_0$, moment map $\mathcal D$, rotation
speeds $\mathcal R_k$, and moment tangent space $\mathcal V$ defined in
(A2)--(A3) of the linked proof. Speeds are antisymmetric, with positive
$f_{ij}$ denoting outward motion of $T_i$.
Write $d_k(A)=\inf_{T'}\sum_i\gamma_n(A_i\triangle T'_i)$, where
$T'$ ranges over the regular-simplex orbit, including labels.

For an ordinary partition $A$, write
$h_{ij}=1_{T_i\cap A_j}$, $e_i=1_{A_i}-1_{T_i}$, and

\[
 U_{\rho,ij}=T_\rho(1_{T_i}-1_{T_j}),\qquad
 L_\rho(A;T)=\sum_{i\ne j}\int U_{\rho,ij}h_{ij}\,d\gamma_n.
\]

For balanced $A$, quadratic expansion gives

\[
 S_\rho(A)-S_\rho(T)=\widetilde Q_\rho(e)-2L_\rho(A;T),\qquad
 \widetilde Q_\rho(e)=\sum_i\langle e_i,(T_\rho-\Pi_0)e_i\rangle.
 \tag{B1}
\]

The closest inspected general result is
[Christof--Wachsmuth, arXiv:1707.07579v2](https://arxiv.org/html/1707.07579v2),
Theorems 4.4 and 6.12. Their explicit bang-bang theorem concerns a bounded
Lipschitz domain, scalar box constraints and a regular switching surface;
it assumes weak-star continuity of the quadratic term. Our weighted
multicell proof below establishes its own compactness and rotational gauge.
The general multicell variation formula in
[Heilman--Tarter, Lemma 4.6](https://doi.org/10.1017/fms.2021.56)
has the same boundary interaction and gradient terms used below; that lemma
does not assume maximality. Its maximizer-specific Lemma 4.8 is not used
to infer stability of the simplex. Neither comparison supplies a global
continuation theorem.

## Relative exterior control below every fixed correlation bound

**Uniform weighted exterior bound.** For every $0<b<1$, choose

\[
 \frac{b}{2(1+b)}<\alpha<\frac14,\qquad
 w_\alpha(x)=(1+|x|^2)e^{\alpha|x|^2}.
 \tag{B2}
\]

For all $0<\rho\le b$, all ordinary measurable $k$-cell partitions
and $R\ge2$, the estimates (C3)--(C5) hold with $w_\alpha$ and a
dimension-dependent tail function $\beta_{n,\alpha}(R)\to0$.
In particular, with exterior restrictions denoted by superscript out,

\[
\begin{split}
 \left(\sum_i\int_{|x|>R}w_\alpha|e_i|\,d\gamma_n\right)^2
 &\le C\beta_{n,\alpha}(R)L_\rho(A;T)/\rho,\\
 0\le\widetilde Q_\rho(e^{\rm out})
 &\le\varepsilon_{k,b}(R)L_\rho^{\rm out}(A;T),\qquad
 \varepsilon_{k,b}(R)\longrightarrow0.
\end{split}
\tag{B3}
\]

These estimates do not require balance. The constants may depend on $k,b$
and $\alpha$; uniformity as $b\uparrow1$ is not asserted.

Here are the details of the extension beyond small correlation. The exact
Mehler kernel satisfies, for $0\le r\le b$,

\[
 M_r(x,y)\le(1-b^2)^{-n/2}
    \exp\!\left(\frac{b}{2(1+b)}(|x|^2+|y|^2)\right).
\]

This follows by bounding $x\cdot y$ by $(|x|^2+|y|^2)/2$ in its
exponent. Its correlation derivative is a quadratic polynomial with
coefficients bounded in terms of $n,b$, times $M_r$. The strict gap in
(B2) absorbs this polynomial, proving

\[
 \left|\frac{M_\rho(x,y)-1}{\rho}\right|
 \le C_{n,b,\alpha}w_\alpha(x)w_\alpha(y).
 \tag{B4}
\]

For any $0\le h\le1$ supported in $t=x\cdot u\ge0$, the one-dimensional
mass-cost argument in (C3) gives

\[
 \left(\int_{|x|>R}w_\alpha h\,d\gamma_n\right)^2
 \le\beta_{n,\alpha}(R)\int\min\{t,1\}h\,d\gamma_n.
 \tag{B5}
\]

An upper bound for $\beta$, up to a fixed multiplicative constant, is
the sum of the tails of
$(1+|x|^2)^2e^{2\alpha|x|^2}$ against $\gamma_n$ outside $R$ and
$(1+|z|^2)^2e^{2\alpha|z|^2}$ against $\gamma_{n-1}$ outside
$\sqrt{R^2-1}$. For $t\ge1$ use Cauchy--Schwarz. For $t<1$,
$w_\alpha(z+tu)\le2e^\alpha(1+|z|^2)e^{\alpha|z|^2}$ and the
normal mass squared is at most $2\phi_1(0)$ times its first normal
moment. Integrating in $z$ proves (B5). Both tails vanish since
$2\alpha<1/2$. The second tail is zero when $n=1$. Removing the
exterior restrictions proves the finite global bound as well.

The all-cell reflection margin (A4) holds for every $0<\rho<1$, and
$\min\{\rho t/\sqrt{1-\rho^2},1\}\ge\rho\min\{t,1\}$.
Sum (B5) over ordered transfers and use (B4) to obtain (B3).
For the second line restrict the transfers themselves to the exterior.
Positivity of $T_\rho-\Pi_0$ supplies the mixed-term estimate (C5).
The available exponential tail margin $1/2-2\alpha$ tends to zero
as $b\uparrow1$; this argument gives no uniform endpoint control.
For example $\alpha=(1+3b)/(8(1+b))$ is admissible and gives the
Gaussian integrability margin $1/2-2\alpha=(1-b)/(4(1+b))$ explicitly.

## The full interface form at a fixed positive correlation

For $0<\rho<1$ define on $F_{ij}$

\[
 a_{\rho,ij}(x)=\partial_{n_{ij}}U_{\rho,ij}(x),\qquad
 A_\rho(f)=\sum_{i<j}\int\frac{a_{\rho,ij}}\rho f_{ij}^2d\nu_{ij}.
\]

Let $\lambda_i^f=\sum_{j\ne i}f_{ij}\nu_{ij}|_{F_{ij}}$ and set

\[
 B_\rho(f)=\frac1\rho\sum_i\iint M_\rho(x,y)
                         d\lambda_i^f(x)d\lambda_i^f(y),\qquad
 \mathfrak q_\rho(f)=A_\rho(f)-B_\rho(f)\quad(f\in\mathcal H_0).
 \tag{B6}
\]

On $\mathcal H_0$ the constant kernel can equally be subtracted from
$M_\rho$, since every $\lambda_i^f$ has total mass zero.

**Interface structure.** At every fixed positive correlation, the local
coefficient is uniformly positive and bounded on all facets:

\[
 \frac\ell{\sqrt{1-\rho^2}}
 \le\frac{a_{\rho,ij}(x)}\rho
 \le\frac{\sqrt{2/\pi}}{\sqrt{1-\rho^2}}.
 \tag{B7}
\]

The interaction $B_\rho$ is a compact quadratic form on
$\mathcal H_0$. The rotation space is contained in the operator kernel
of $\mathfrak q_\rho$, including zero mixed pairing with every balanced
speed. If $\mathfrak q_\rho\ge0$ and its kernel is exactly
$\mathcal R_k$, then

\[
 \mathfrak q_\rho(f)\ge\kappa_\rho
                  \operatorname{dist}_{\mathcal H}(f,\mathcal R_k)^2
 \quad\text{for some }\kappa_\rho>0.
 \tag{B8}
\]

For (B7), reflect in $H_{ij}$ and put
$\sigma=\sqrt{1-\rho^2}$, $m=\rho x/\sigma$. On the facet, $m_N=0$
and $m\in T_i$. Differentiating the reflected Gaussian difference gives

\[
 a_{\rho,ij}(x)=\frac{2\rho}{\sigma}
                     \int_{T_i-m}z_N\phi_n(z)dz.
\]

Cone addition gives $T_i\subset T_i-m\subset\{z_N\ge0\}$.
The lower bound follows from $2b_i^*\cdot n_{ij}=\ell$; enlarging
the last integral to the halfspace gives the upper bound.

For each pair of facets, the square of $M_\rho$ times their two Gaussian
surface densities is bounded by a constant times

\[
 \exp\!\left(-\frac{1-\rho}{2(1+\rho)}(|x|^2+|y|^2)\right).
\]

This is integrable on the supporting planes. Thus every interaction block
is Hilbert--Schmidt, and the finite sum and its compression to
$\mathcal H_0$ are compact. The form is bounded and symmetric.

For smooth speeds compactly supported on facet interiors, use their
normal-strip deformations and the tree-bump balance correction from
(A2). The exact score difference is zero on each pair plane. Expanding
its regret on a strip gives
$t^2\sum_{i<j}\int a_{\rho,ij}f_{ij}^2d\nu_{ij}/2+o(t^2)$.
The interaction is
$t^2\rho B_\rho(f)+o(t^2)$. Second-order balance corrections contribute
no term because the original partition is critical. Equation (B1) gives

\[
 -\frac1{2\rho}\left.\frac{d^2}{dt^2}S_\rho(A(t))\right|_{t=0}
       =\mathfrak q_\rho(f).
 \tag{B9}
\]

The same calculation is the specialization of Heilman--Tarter Lemma 4.6
with zero interface scores; its acceleration and density terms vanish.
It does not use a maximizing hypothesis. Differentiating a rotation and
one arbitrary deformation at this critical partition, rotation invariance
gives zero mixed second derivative. Gaussian domination permits the
unbounded linear rotational velocity; truncation and density extend the
zero mixed pairing to the whole balanced facet space. In particular a
global maximizing simplex necessarily has $\mathfrak q_\rho\ge0$.

To prove (B8) without a spectral assumption about essential eigenvalues,
suppose unit speeds $f_j\perp\mathcal R_k$ have
$\mathfrak q_\rho(f_j)\to0$. Pass weakly to $f$.
Compactness of $B_\rho$ and weak lower semicontinuity of $A_\rho$
give $\mathfrak q_\rho(f)\le0$. Nonnegativity and the kernel hypothesis
force $f=0$. Now $B_\rho(f_j)\to0$, contradicting
$A_\rho(f_j)\ge\ell$. This also shows that the kernel of a nonnegative
$\mathfrak q_\rho$ is finite dimensional: an infinite orthonormal
sequence in it would give the same contradiction. For $k=2$,
$\mathcal H_0=\{0\}$ and these assertions are understood vacuously.

## An interface gap implies a uniform measurable local comparison

**Conditional local comparison.** Fix $r\in(0,1)$ and assume

\[
 \mathfrak q_r\ge0\text{ on }\mathcal H_0,
 \qquad \ker\mathfrak q_r=\mathcal R_k.
 \tag{B10}
\]

There exist $c,c',\eta,\varepsilon>0$, depending on $k,r$, such that
every ordinary balanced partition of $\mathbb R^n$ with $d_k(A)<\eta$
and every $|\rho-r|<\varepsilon$ satisfy

\[
 S_\rho(T)-S_\rho(A)\ge c L_\rho(A;T_A)
                         \ge c'd_k(A)^2.
 \tag{B11}
\]

Here $T_A$ is obtained from a closest regular moment matrix to the moment
matrix of $A$. The neighborhood of $r$ is restricted to $(0,1)$.
The hypothesis (B10) is part of the theorem, not a conclusion at arbitrary
correlation. No smoothness, support restriction or global interface graph
is required of $A$.

First note two facts about the moment gauge. If partitions approach the
simplex orbit in Gaussian $L^1$, their moments approach the regular
moment orbit, by uniform integrability of $|x|$. After an initial metric
alignment to $T$, every closest regular moment matrix tends to $B_*$.
Since $B_*$ has full column rank, its orthogonal orbit has a continuous
local alignment: if $C_j$ is a regular matrix tending to $B_*$, the
orthogonal matrix solving $C_jO_j=B_*$ is obtained from
$O_j=(n/F_k)C_j^TB_*$ and tends to the identity. Thus the nearest-moment
alignment preserves metric convergence to $T$. Differentiating its
nearest-point condition gives exactly $B(A)-B_*\perp\mathcal V$.
These facts require no global first-moment deficit estimate.

On the closed subspace
$\mathcal G=\{f\in\mathcal H_0:\mathcal Df\perp\mathcal V\}$,
(B8) improves to $\mathfrak q_r(f)\ge\kappa_G\|f\|^2$ for some
$\kappa_G>0$. To check this, write $f=g+u$ with $u\in\mathcal R_k$
and $g\perp\mathcal R_k$. The map
$\mathcal D:\mathcal R_k\to\mathcal V$ is an isomorphism, so
$\mathcal Du=-\Pi_{\mathcal V}\mathcal Dg$ bounds $\|u\|$ by
a constant times $\|g\|$. Zero mixed pairing with rotations and (B8)
then prove the claim. When the rotation space is zero no inverse is needed.
If $A_r(f)\le C_r\|f\|^2$, choose a fixed $0<\delta<1/4$ so small
that $2\delta C_r\le\kappa_G$. We have

\[
 B_r(f)\le(1-2\delta)A_r(f)\qquad(f\in\mathcal G).
 \tag{B12}
\]

If the first inequality in (B11) failed with $c=\delta$, there would
be $\rho_j\to r$ and nonsimplex balanced partitions approaching the
orbit with
$S_{\rho_j}(T)-S_{\rho_j}(A^{(j)})<\delta L_{\rho_j}(A^{(j)};T_{A^{(j)}})$.
Apply the exact moment gauge to fix all reference simplices at $T$.
Write $L_j=L_{\rho_j}(A^{(j)};T)>0$ and $s_j^2=L_j/\rho_j$.
Since score differences are bounded and the changed-label measure tends
to zero, $s_j\to0$. Estimate (B3), including its global version, makes
$\mu_{j,i}=e_i^{(j)}\gamma_n/s_j$ uniformly bounded and tight in
weighted variation. The normal densities

\[
 \theta_{j,il}(z)=\frac1{s_j\phi_1(0)}
      \int_0^\infty h_{il}^{(j)}(z+tn_{il})\phi_1(t)dt
\]

are bounded in their pair-plane $L^2$ spaces. As in (C9)--(C10), their
weak limits $\theta_{il}\ge0$ are supported on $F_{il}$, and their
antisymmetric differences $f_{il}=\theta_{li}-\theta_{il}$ give
$\mu_i=\lambda_i^f$. Normal distances bounded away from zero have
normalized mass tending to zero. Plane $L^2$ bounds exclude junction
mass. Weighted tightness passes balance and first moments, hence
$f\in\mathcal G$.

The sharp regret lower bound now has a variable coefficient. On a compact
tangential set and $0\le t\le h$, reflection and smoothness of the OU
scores give

\[
 \frac{U_{\rho_j,il}(z+tn_{il})}{\rho_j}
 \ge\left(\frac{a_{r,il}(z)}r-o_j(1)-\omega(h)\right)t,
 \qquad \omega(h)\longrightarrow0.
\]

Use compact cutoffs supported in facet interiors, where (B7) bounds the
coefficient away from zero, and take $h$ small enough to retain positivity.
Apply the one-dimensional density bound
$\int_0^h t h_{il}^{(j)}\phi_1dt\ge
(\int_0^h h_{il}^{(j)}\phi_1dt)^2/(2\phi_1(0))$.
After dividing by $\rho_js_j^2=L_j$, take $j\to\infty$, then
$h\downarrow0$, and exhaust facet interiors. The weak $L^2$ limits are
unchanged by the normal cutoff, since omitted normalized mass is
$O(s_j)$ for each fixed $h>0$. This proves

\[
 1\ge\frac12\sum_{i\ne l}\int\frac{a_{r,il}}r\theta_{il}^2d\nu_{il}
       \ge\frac12 A_r(f).
 \tag{B13}
\]

There is no loss on junctions, which are plane-null sets. Weighted
tightness, (B4), and local uniform convergence of
$(M_{\rho_j}-1)/\rho_j$ give

\[
 \frac{\widetilde Q_{\rho_j}(e^{(j)})}{L_j}
       \longrightarrow B_r(f).
 \tag{B14}
\]

For this full-kernel limit, on each compact product approximate its
continuous kernel uniformly by finite sums $p(x)q(y)$; weak convergence
and bounded variations then give convergence of the product integrals.
The weighted exterior bound makes the omitted product tails uniformly
small. Balance cancels the constant kernel in the limit.

By the assumed failure and (B1), the left side of (B14) exceeds
$2-\delta$. But (B12)--(B13) give
$B_r(f)\le2(1-2\delta)<2-\delta$, a contradiction.
This proves the first inequality in (B11).
The global version of (B5) and the margin bound imply
$(\sum_i\gamma_n(A_i\triangle (T_A)_i))^2\le C L_\rho(A;T_A)/\rho$.
As $\rho$ stays near the positive $r$, this proves the second inequality.

## An effective local radius from truncated normal flux

The preceding comparison admits an error bound in the active dimension.
Fix $0<b<1$ and a set of correlations $J\subset(0,b]$ on which
$B_\rho(g)\le(1-\theta)A_\rho(g)$ on $\mathcal G$, with a given
$0<\theta<1$. Here $A_\rho,B_\rho$ are the two interface terms in
(B6), and $\mathcal G$ is the balanced space with rotational moment
component zero. A uniform interface gap gives such a $\theta$ by
the argument of (B12). In particular, the explicit gaps (S37) and
(S83) supply the premise on three-cell closed interior bands and on
the all-cell initial interval, respectively.

Choose the nearest regular moment alignment and put $s^2=L_\rho/\rho$.
For $R\ge2$ and $0<h<1$, retain the transfer $i\to j$ only when
$|x|\le R$, $0\le t=n_{ij}\cdot x\le h$, and its projection
$z=x-tn_{ij}$ belongs to $F_{ij}$. Project each retained transfer,
divided by $s$, onto the facet. Let $\theta_{ij}$ be its nonnegative
density with respect to $\nu_{ij}$ and put
$f_{ij}=\theta_{ji}-\theta_{ij}$. The pointwise one-sided density
inequality of (C9)–(C11) applies to these restrictions as well.
Write $\lambda^f$ for the resulting vector of signed facet measures.

The loss of normalized mass and moment in this projection, and
the error in the normalized interaction, are controlled by

\[
 E(R,h,s)=C_{k,b}\sqrt{\beta_{n,\alpha}(R)}
              +C_{k,b,R}\big(\sqrt h+s/h\big).
 \tag{B15}
\]

Indeed the exterior weighted total variation is bounded by
$C_{k,b}\sqrt{\beta_{n,\alpha}(R)}$ by (B3), with a constant independent
of $R$. Inside the ball, transfers with $t>h$ have mass at most
$C_k s^2/h$ by (A4). If $x\in T_i$ and $t\le h$ but $z\notin F_{ij}$,
some remaining defining inequality of $F_{ij}$ is violated at $z$
by at most $C_k h$. Its tangential linear functional is nonzero.
The corresponding strip in the plane has Gaussian surface measure
$O_k(h)$, and the projected thin densities have uniformly bounded
plane $L^2$ norm by (C9). This discards $O_k(\sqrt h)$ normalized
mass. For $k=2$ there are no such strips. Moving the remaining mass
from $x$ to $z$ changes first moments by $O_k(h)$.
These estimates give the mass and moment bound in (B15).

For the interaction bound, $K_\rho=(M_\rho-1)/\rho$ is uniformly
bounded by $C_{k,b}w_\alpha\otimes w_\alpha$ by (B4), and its first
spatial derivatives on compact products are bounded uniformly for
$0<\rho\le b$, by differentiating the Mehler formula, with its
continuous extension at $\rho=0$. The normalized full transfers and
the projected facet measures have uniformly bounded weighted
variation, the latter by the plane $L^2$ bound. Hence in the double
integral, exterior terms cost $O_{k,b}(\sqrt{\beta_{n,\alpha}(R)})$.
The remaining discards cost $O_{k,b,R}(\sqrt h+s/h)$ and the retained
projections cost $O_{k,b,R}(h)$ by the compact derivative bound.
Thus the difference between $\widetilde Q_\rho(e)/L_\rho$ and the
$K_\rho$ interaction of $\lambda^f$ is bounded by (B15). The centered
kernel is used for the projected measures even before balance is restored.

Reflection and the bounded second spatial derivatives of the normalized
score on the compact ball give
$U_{\rho,ij}(z+tn_{ij})/\rho\ge (a_{\rho,ij}(z)/\rho-C_{k,b,R}h)t$.
The lower bound $\ell$ for $a_{\rho,ij}/\rho$ and the one-sided density
inequality therefore give

\[
 A_\rho(f)\le\frac{2}{1-C_{k,b,R}h/\ell}
 \tag{B16}
\]

when the denominator is positive. This uses
$(\theta_{ji}-\theta_{ij})^2\le\theta_{ij}^2+\theta_{ji}^2$.
The exact balances and moment gauge for the full transfers make the
mass and rotational moment defects of $f$ at most (B15).
Remove the mass defect with constant speeds on a fixed spanning
tree of positive-area facets. Each edge constant is its tree flow
divided by its facet area, so the $L^2$ correction is bounded by a
fixed multiple of the defect. Then remove the remaining skew moment
with a rotation speed using the isomorphism
$\mathcal D:\mathcal R_k\to\mathcal V$. Both operations are finite
linear maps with norms bounded by facet Gaussian integrals and
$B_*$. They produce $g\in\mathcal G$ with
$\|g-f\|\le C_k E(R,h,s)$. The bounds for $A_\rho$ and $K_\rho$
on the weighted facet space make the change in the two quadratic
terms at most a fixed multiple of $E$ when $E\le1$.

Consequently $\widetilde Q_\rho(e)/L_\rho$ is bounded by
$(1-\theta)2/(1-C_{k,b,R}h/\ell)+C_{k,b}E(R,h,s)$,
after enlarging the constants. Choose $R$ first to control the
exterior error, then $h$ to control the strips and the denominator,
then an upper bound on $s$. This makes the expression at most
$2-\theta$ uniformly on $J$. These choices only use the specified gap,
Gaussian tail and strip integrals, compact kernel derivative bounds,
and the finite correction maps.

An effective cutoff can be imposed directly on the aligned linear
regret $r_1=F_k-\langle B_*,B_A\rangle_F=\ell\sum_{i\ne j}\int t h_{ij}$.
The upper reflection bound
$U_{\rho,ij}\le 2\Phi(\rho t/\sqrt{1-\rho^2})-1$ from the
[saturated transfer comparison](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md#a-compact-core-test-with-saturated-transfer-cost)
gives $s^2\le (2\phi_1(0)/(\ell\sqrt{1-b^2}))r_1$. For a partition
at summed label distance $d$ from the orbit, Bessel's inequality
and nearest-moment alignment give $\|B_A-B_*\|_F\le\sqrt d$,
so $r_1\le\sqrt{F_k d}$. This yields an effective form of (B11) in active
dimension, with constants that can also be used as $\rho\downarrow0$
if the relative interface gap is uniform there. Dimension reduction
then suffices for local exclusion of actual global maximizers in higher
dimension. No numerical local radius has been evaluated here.

## A finite test for distant competitors on a correlation band

Fix $I=[a,b]\subset(0,1)$ and work in dimension $n$. The compact-core
test (13)–(14) extends to this band by keeping the exact score loss
and finitely many quadratic Hermite terms. For any $b<r<1$, choose a
weight as in (B2) at $r$. Write $e^c=e1_{\{|x|\le R\}}$,
$L_c=L_\rho^c$, $L_t=L_\rho^{\rm out}$, and

\[
 P_N(\rho,e^c)=\sum_{m=1}^N\rho^m\sum_i\|\Pi_m e_i^c\|_2^2 .
\]

The fixed-$r$ Mehler bound and (B5) bound
$\sum_i\langle e_i^c,T_r e_i^c\rangle$ by a constant times
$L_c/\rho$. For $m>N$ we have $\rho^{m-1}\le
r^{m-1}(b/r)^N$. Thus the omitted interaction is at most
$\tau_N L_c$, where $\tau_N=C_{k,r}(b/r)^N$. Combining with (C5)
gives, for every $u>0$,

\[
 S_\rho(T)-S_\rho(A)\ \ge\
 (2-u-\tau_N)L_c-P_N(\rho,e^c)
 +[2-\varepsilon_{k,b}(R)(1+C_0/u)]L_t .
 \tag{B17}
\]

Use the core body of (14), replacing the truncated cost coordinate by
$L_c$ and retaining the Hermite coefficients through degree $N$.
Take $N\ge1$ so the body retains the core first moments. Its support
function is still the integral of the pointwise maximum over staying
and the possible target labels. In the transfer coordinates of (14),
the additional linear features for $i\to j$ are $U_{\rho,ij}(x)$ and
$(\mathbf e_j-\mathbf e_i)H_\nu(x)$ for the orthonormal Hermites
$1\le|\nu|\le N$. The mass error, distance to the positive semidefinite
moment cone, and lower bound on changed-label mass are the same as in
(14). They include every aligned distant competitor.
With the effective linear-regret cutoff $r_*$ one may replace the
changed-label coordinate by the constraint
$-\langle B_*,M_c\rangle_F\ge r_*-\ell\int_{|x|>R}|x|\,d\gamma_n$.
The last term bounds the omitted linear regret because the total
transfer density at a point is at most one.

A strictly positive lower bound for
$(2-u-\tau_N)L_c-P_N$ on this relaxed distant body, together with
$\varepsilon_{k,b}(R)(1+C_0/u)<2$, excludes distant ties and
improvers over the band. This is a finite-dimensional test. A strict
bound can be enclosed by finitely many support directions and verified
integrals on the fixed ball: the feature vectors are bounded there,
and increasingly fine directional nets give outer polytopes
converging to the compact body. The closed tolerance constraints and
continuity of the objective preserve strict margins in sufficiently
accurate outer bodies. A finite correlation mesh also suffices.
For fixed transfers
$L_c=-\sum_i\langle e_i^c,T_\rho1_{T_i}\rangle$, and the OU and
truncated OU operators have Lipschitz constant at most
$\sup_{m\ge1}mb^{m-1}$ on $I$, by their Hermite eigenvalues.
Since $\sum_i\|e_i^c\|^2\le2$, the objective is uniformly Lipschitz
in $\rho$ with an explicit bound. The core constraints are fixed in
$\rho$, so strict mesh margins larger than the interpolation error
bound the whole band.

Alternatively the body can be held fixed for all $\rho$ in the band.
Truncate the score to $T_\rho\Pi_{\le N}(1_{T_i}-1_{T_j})$, where
$\Pi_{\le N}$ is the sum of Hermite projections through degree $N$.
Its remainder is $T_r T_{\rho/r}(I-\Pi_{\le N})(1_{T_i}-1_{T_j})$.
The input to $T_r$ has $L^2$ norm at most $C(\rho/r)^{N+1}$,
and $T_r$ maps $L^2$ to functions with uniformly bounded gradients
on the fixed ball, by the $L^2$ norm of the spatial derivative of its
Mehler kernel there. The remainder vanishes on the whole pair plane
by reflection, which commutes with the projections. Its pointwise
error is therefore bounded by $C_{R,r}t(\rho/r)^{N+1}$. Using
$t\le R\min(t,1)$ on the core and (A4), the score-loss error is at
most $\delta_N L_c$ with $\delta_N=C_{k,R,r}(b/r)^N$.
The truncated score loss $L_c^{(N)}$ is polynomial in $\rho$ and
linear in the same core Hermite features as $P_N$.
For $2-u-\tau_N>0$ one can thus test
$(2-u-\tau_N)L_c^{(N)}/(1+\delta_N)-P_N$ in (B17).
The support function uses polynomial features and constraints
independent of $\rho$. At a fixed $R$ the additional score error
tends to zero as $N\to\infty$.

Strict global uniqueness on $I$ implies the existence of such a
certificate for every fixed positive distant cutoff $\eta$. Otherwise
there are relaxed cores with $R\to\infty$, correlations converging in
$I$, and $2L_c-\widetilde Q_\rho(e^c)\le o(1)$. Extend each core by
$T$ and take a weak fractional limit in the fixed active dimension.
The tolerances imply balance and the positive semidefinite moment gauge
in the limit, and the mismatch probability stays at least $\eta/2$.
Compactness of the OU operators and norm continuity in $\rho$ pass
the objective to the limit. The additional constant Hermite term
is $\|\int e^c\|^2$ and tends to zero. The limit is a fractional
maximizer, which must be ordinary by the strictly positive quadratic
variation on the balanced feasible set. By uniqueness and the
positive semidefinite gauge it is $T$, contradicting the mismatch
lower bound, or the positive linear-regret lower bound when that
constraint is used. Hence the exact core deficit has a uniform positive
lower bound for all sufficiently large $R$. Choose $u$ and $\tau_N$
small enough to retain some of that bound, using $0\le L_c\le1$,
and then $R$ large for the tail absorption. Truncation subtracts
nonnegative interaction and so cannot reduce this margin further.

[De–Mossel–Neeman, Theorem 2.2](https://theoryofcomputing.org/articles/v015a006/v015a006.pdf)
supplies additive approximate optimization via bounded-degree
polynomial threshold functions and approximate dimension reduction.
The test here instead uses the exact active dimension for actual
maximizers and separates a strict comparison away from the orbit
from a local comparison at the simplex. When the interface gap is
available, the effective local radius and (B17) together give a way
to certify equality rigidity on a band. For three cells, this applies
to the intermediate range between the direct comparison and the
qualitative near-one interval. It supplies no numerical certificate
for that range. The tail and degree choices are for $b<1$ and are
not uniform at the perimeter endpoint.

## The first continuation endpoint has a tie or an additional zero mode

The interface hypothesis is available for every fixed $k$ on some small
positive interval. Indeed the all-cell first-moment form
$Q_k(f)=\ell\|f\|^2-\|\mathcal Df\|^2$ has a gap modulo rotations
by (A2)--(A3). Taylor expansion in Hilbert--Schmidt norm on every facet
pair gives

\[
 \|M_\rho-1-\rho x\cdot y\|_{L^2(\nu\otimes\nu)}\le C_k\rho^2
 \qquad(0<\rho\le1/8).
\]

The second correlation derivative is a degree-at-most-four polynomial
times a Mehler kernel. For this range its square is integrable with
Gaussian surface density, including the polynomial factor. On balanced
speeds the constant term vanishes, so
$B_\rho(f)=\|\mathcal Df\|^2+O_k(\rho)\|f\|^2$.
The lower bound (B7) gives
$\mathfrak q_\rho(f)\ge Q_k(f)-C_k\rho\|f\|^2$.
Apply this on $\mathcal R_k^\perp$ and use the exact rotational mixed
kernel to conclude (B10) and a uniform positive gap for all sufficiently
small $\rho$. No sign outside that interval follows from the expansion.

Let $r_k$ be the endpoint of the initial global optimality interval from
the [positive-correlation compactness argument](SMALL_POSITIVE_CORRELATION_INTERVAL.md#a-rigorous-continuation-boundary-in-correlation).
If $r_k<1$, at least one of the following holds:

* There is a nonsimplex ordinary balanced global maximizer tied with $T$
  at correlation $r_k$.
* The nonnegative form $\mathfrak q_{r_k}$ has a nonrotational zero mode
  in $\mathcal H_0$. Its kernel is finite dimensional.

To prove this sharper alternative, the earlier compactness result gives
improving global maximizers at correlations decreasing to $r_k$ from
above and a strong subsequential limit maximizing at $r_k$.
If that limit is nonsimplex, the first alternative holds. Otherwise the
maximizers approach the simplex orbit. Global optimality of $T$ at $r_k$
and (B9) imply $\mathfrak q_{r_k}\ge0$. If its kernel were exactly
rotational, (B11) would forbid those improving maximizers. Thus there is
an additional zero mode. The finite-dimensional assertion follows from
the compact interaction and positive multiplication term proved above.

This necessary alternative does not exclude a global crossing by distant
maximizers, even if the simplex retains a strictly positive interface gap.
Proving (B10) for every positive correlation addresses only the second
alternative. Global variational rigidity must also rule out the first.
At the large-correlation endpoint the weights in (B2) approach the
integrability threshold. A proof using the perimeter limit still needs
deformation-relative estimates for moving maximizers; none are supplied
by the pointwise perimeter limit or by (B3) with fixed $b<1$.

The independent review obligations are the strict exponent window (B2),
the exterior normalization, the global facet bounds (B7), the compact
interaction, the exact rotational mixed kernel, simultaneous metric and
moment alignment, the variable-coefficient factor $1/2$ in (B13), and
the product-measure limit (B14). The local theorem is an implication
under (B10), not a proof of that hypothesis for arbitrary cell count
throughout $(0,1)$. The separate three-cell theorem now supplies (B10)
at every positive correlation and excludes its extra-zero-mode
alternative. For $k\ge4$ both endpoint alternatives remain open.

## Strict continuation past the three-cell bound \(221/500\)

There is an $\varepsilon>0$ such that for every $d\ge2$ the balanced
regular tripod is globally optimal throughout
$0<\rho\le221/500+\varepsilon$. Every maximizer on
$[11/25,221/500+\varepsilon]$ is a cylindrical regular tripod.
The strict extension uses the equality information in the
[finite-amplitude comparison](THREE_CELL_FINITE_AMPLITUDE_COMPARISON.md)
and the three-cell local gap. Its extra length is qualitative.

Work in the active plane, as permitted for actual maximizers by
Heilman–Tarter. On $[11/25,221/500]$, (FA6) excludes a tie when
the scalar radial sum is negative. If the sum is nonnegative,
(FA6)–(FA7) and the strict determinant in (FA17) force $A=G=0$
in the trial-plus-remainder decomposition. Thus any tying partition
has equal masses on almost every circle. The equal-circle argument
(DG2) gives $\mathcal D_\rho\ge0$ and $\mathcal R_\rho=0$.
Equality in (DG3) therefore forces $\mathcal C_\rho=0$. Its angular
quadratic kernel has strictly positive Fourier multipliers at every
positive concentration parameter. The radial pair density is strictly
positive for positive radii, so almost every pair of slices has identical
centered cell indicators. Fubini gives one angular partition independent
of radius.

For that angular partition the nonnegative self-comparison with the
three equal arcs has zero integral in (DG3). Its kernel parameter
$\kappa=\rho r u/(1-\rho^2)$ has full support in $(0,\infty)$,
and the comparison is continuous in $\kappa$. Equality therefore
holds for every $\kappa>0$. Differentiate at zero. The normalized
circular kernel has derivative $\cos(\theta-\varphi)$ there, so the
sum of squared first circular Fourier moments equals that of the three
equal arcs. For a set of given angular mass $1/3$, the norm of
its first Fourier moment is maximized uniquely by an arc, up to
rotation and null sets: choose its moment phase, then integrate
cosine over its upper level set of that measure. Each cell must
attain this separate bound. The three equal arcs tile the circle
and yield the regular tripod.

If nonsimplex maximizers existed at correlations arbitrarily close
to the right of $221/500$, reduce each to the active plane and take
a strongly convergent subsequence by positive-correlation maximizer
compactness. The limit maximizes at $221/500$ and hence is a tripod.
The full three-cell interface gap at this correlation and (B11)
exclude all those nearby nonsimplex maximizers, a contradiction.
This proves the stated extension. In particular the first failure,
if any, lies strictly beyond this endpoint, and the finite band
test has a strict certificate in principle on a short band containing
the endpoint. No ordering relative to the earlier implicit
threshold extension or overlap with the near-one interval is asserted.
