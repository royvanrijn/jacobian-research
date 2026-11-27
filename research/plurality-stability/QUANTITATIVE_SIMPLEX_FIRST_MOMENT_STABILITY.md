# Quantitative stability from the equal-cell Gaussian first-moment theorem

For every fixed number of cells $k\ge2$, the equal-mass first-moment
endpoint is supplied by Yongxian Zhang's external preprint. The argument
below deduces a dimension-independent square-root distance estimate with
explicit constants, without optimizing them. This is an internal written
deduction from that preprint, not a claim made in its frozen edition and
not an independently reviewed or formally verified theorem. The uniform
positive-correlation interval for four cells now has an
[internal normal-flux compactness proof](SMALL_POSITIVE_CORRELATION_INTERVAL.md)
using these inputs; independent mathematical review is pending.

## External endpoint and application audit

The primary source is *Equal Gaussian cells: the sharp simplex first-moment
theorem*, immutable release
[`gaussian-balanced-simplex-all-k-v1`](https://github.com/mxym/math/releases/tag/gaussian-balanced-simplex-all-k-v1),
commit `ebfd90558439afca21babeb7b77c8d8837ea745f`, 8 October 2026.
The exact [retained manuscript](sources/zhang-equal-gaussian-cells-ebfd905.md)
has SHA-256
`8416f741398ceb4207edcc3ff31964883ae14044698043bf17a75d66c1dec832`.
Its frozen author field says “mxym/math research project”; the attribution
to Yongxian Zhang is consistent with the
[signed prescribed-mass companion](https://opencode.terrychan.org/mxym/math/blob/main/research/gaussian-prescribed-mass-centroid-ellipsoid/README.md).
This dependency is an external preprint, not a paper of this repository.

Put $n=k-1$, $P=I-\mathbf1\mathbf1^T/k$,

\[
c_k=\frac{\mathbb E\max_{i\le k}Z_i}{\sqrt n},\qquad F_k=c_k^2,
\qquad Q_*=P/n.
\]

The preprint's Theorem 1 states that $W_1(A)\le F_k$ for every measurable
Gaussian partition of masses $1/k$ in every dimension $d\ge1$.
For $d\ge n$, equality means a central regular-simplex partition extended
cylindrically, up to orthogonal maps, labels and null sets. Below $n$, the
inequality is strict; no sharp fixed-lower-dimensional value is asserted.
The source also includes fractional partitions. Here the distance theorem
is stated for ordinary partitions and $d\ge n$.

The essential imported geometric theorem is
[Milman–Neeman, *The Gaussian Double-Bubble and Multi-Bubble Conjectures*,
Annals of Mathematics 195 (2022), 89–206, Theorem 1.1](https://arxiv.org/html/1805.10961v3).
Its stated range is $2\le k\le N+1$ in ambient dimension $N$.
Zhang applies it to full-rank affine winning partitions in dimension
$N=k-1$, so every cell count is in range. Its perimeter is
$S=\frac12\sum_iP_\gamma(A_i)$, the sum of the ordinary Gaussian
interface areas. The actual statement, cluster definition and perimeter
normalization were inspected; the published geometric proof is imported,
not independently replayed here.

For clarity, the following is the precise part of Zhang's argument used
below. For centered positive semidefinite $Q$, let

\[
\mathcal C(Q)=\inf_\lambda\left\{k^{-1}\sum_i\lambda_i+
  \mathbb E\max_i(Y_i-\lambda_i)\right\},\quad Y\sim N(0,Q).
\]

On the positive cone in $\mathbf1^\perp$, balancing prices in the
zero-sum gauge exist, are unique and smooth. With facet lengths
$\ell_{ij}$, areas $A_{ij}$, weights $w_{ij}=A_{ij}/\ell_{ij}$,
and their Laplacian $L$, Sections 2–4 prove

\[
\mathcal C(Q)=\operatorname{tr}(LQ)=\sum_{i<j}w_{ij}\ell_{ij}^2,
\quad d\mathcal C(Q)[H]=\tfrac12\operatorname{tr}(LH).
\tag{E1}
\]

For every centered trace-one $Q\succeq0$, set
$Q_t=Q_*+t(Q-Q_*)$. Then

\[
F_k-\mathcal C(Q)^2=\int_0^1D(t)t^{-2}\,dt,
\quad D(t)=\frac2n\left[
 \sum_{e<f}w_e(t)w_f(t)(\ell_e(t)-\ell_f(t))^2+S(t)^2-S_*^2\right],
\quad S_*=c_k\sqrt{n/2}.
\tag{E2}
\]

Here $e,f$ denote unordered pairs of distinct labels, counted once,
and $S(t)\ge S_*$ by Milman–Neeman. All weights are positive for $t<1$.
In particular $\mathcal C(Q)\le c_k$, with equality only at $Q_*$.

The local application audit checked the following substantive points.
The price Hessian is $L$ with these ordinary-interface normalizations;
its positivity and the implicit function theorem justify the envelope
derivative in (E1). Flux gives the winning moments $LM$ for score matrix
$M$; truncating at large balls removes the unbounded-boundary terms.
Weighted Cauchy and the perimeter theorem give
$\mathcal C(Q)\operatorname{tr}L/n\ge F_k$. Along $Q_t$,
$tC'=(C-\operatorname{tr}L/n)/2$, so, with $h=C^2-F_k$,
$(h/t)'=-D/t^2$. Symmetry gives $L_*=c_kP$, hence $h'(0)=0$.
Continuity from coupling matrix square roots covers the singular endpoint
$t=1$. Nonnegativity and the endpoint limits justify the improper
integral (E2). Equality forces equal pair lengths at an interior time,
then $Q=Q_*$; assignment equality forces the winning labels. These
steps supply the imported covariance and moment endpoints without any
stationary-shape classification. This is an internal source/application
audit, with no external peer-review or full formalization assertion.

## Uniform covariance control with explicit constants

Write $q=k(k-1)/2$ and let $v_j=\pi^{j/2}/\Gamma(1+j/2)$ be the volume
of the Euclidean unit $j$-ball, with $v_0=1$. Define

\[
\alpha_k=\frac{v_{n-1}}{\sqrt2\,4^{n-1}(2\pi)^{n/2}}
                  e^{-n(k+2)^2},\qquad
\beta_k=\frac{q\alpha_k^2}{4n}>0.
\tag{C1}
\]

For $0\le t\le1/2$, the eigenvalues of $Q_t$ on
$E=\mathbf1^\perp$ lie between $1/(2n)$ and 1. Indeed
$Q_t\succeq(1-t)P/n$, it has trace one, and is positive on $E$.
We next prove the uniform lower bound $w_{ij}(t)\ge\alpha_k$.

Choose the balancing price gauge $\sum_i\lambda_i=0$. A zero price
vector gives $\mathcal C(Q_t)\le\mathbb E\|Y\|\le1$, while the score
with minimum price gives $\mathcal C(Q_t)\ge-\min_i\lambda_i$.
Consequently all $\lambda_i\ge-1$. The simplex
$\lambda_i\ge-1,\sum_i\lambda_i=0$ has vertices with one coordinate
$k-1$ and the rest $-1$, so convexity of the norm gives
$\|\lambda\|\le\sqrt{k(k-1)}\le k$.

For a fixed pair $i,j$, choose $z\in E$ with
$z_i=z_j=(k-2)/k$ and $z_l=-2/k$ for the other labels.
In the hyperplane $y_i-y_j=\lambda_i-\lambda_j$, the ball of radius
$1/4$ about $\lambda+z$ lies in the top $ij$ facet: the other
score gaps are 1 at its center and change by at most $\sqrt2/4$.
Every point of this ball has norm at most $k+2$, since
$\|z\|^2=2(k-2)/k\le2$. The Gaussian density on $E$ is at least
$(2\pi)^{-n/2}e^{-n(k+2)^2}$ there. In fact
$\det_E Q_t\le1$ and $Q_t^{-1}\preceq2nI_E$.

The coarea formula in score space gives

\[
w_{ij}(t)=\frac1{\sqrt2}\int_{\text{top }ij\text{ facet}}
                         \varphi_{Q_t}(y)\,d\mathcal H^{n-1}(y).
\tag{C2}
\]

To see the factor, the linear score difference $y_i-y_j$ has gradient
of norm $\sqrt2$ on $E$. Equivalently integrate its delta function
first in score coordinates, or first in standard Gaussian coordinates,
where it equals $A_{ij}/\ell_{ij}$. The facet ball and its density
lower bound now prove (C2) is at least $\alpha_k$, including $n=1$
with zero-dimensional ball volume one.

Put $H=Q-Q_*$, $d_e=(e_i-e_j)^TH(e_i-e_j)$. The squared lengths
along the path are $\ell_e(t)^2=2/n+td_e$. Since any trace-one positive
covariance has $\ell_e\le\sqrt2$,

\[
(\ell_e(t)-\ell_f(t))^2\ge t^2(d_e-d_f)^2/8.
\]

Moreover $\sum_e d_e=k\operatorname{tr}H=0$. If the symmetric matrix
$D_H$ has entries $d_{ij}$ off diagonal and zero diagonal,
$H=-PD_HP/2$; hence $\sum_e d_e^2\ge2\|H\|_F^2$. Thus

\[
\sum_{e<f}(d_e-d_f)^2=q\sum_e d_e^2\ge2q\|H\|_F^2.
\]

Keeping only $0<t<1/2$ and the edge-variance term in (E2) proves

\[
 F_k-\mathcal C(Q)^2\ge\beta_k\|Q-Q_*\|_F^2.
\tag{C3}
\]

For $k=2$ the trace-one domain is a singleton, so the right side is
zero and all the statements still apply. No uniform regularity near a
singular endpoint or estimate for its prices is used.

## From moments to regular-score regret

Let $A$ be any ordinary balanced partition of $\mathbb R^d$, $d\ge n$,
and let $B$ have rows $b_i=\int_{A_i}x\,d\gamma_d$. Write
$W=\|B\|_F^2$ and $\delta=F_k-W\ge0$. If $W>0$, put
$Q=BB^T/W$. The feasible-label inequality gives
$W\le\mathcal C(BB^T)=\sqrt W\,\mathcal C(Q)$, so

\[
\|Q-Q_*\|_F^2\le\delta/\beta_k.
\tag{M1}
\]

Choose a regular centered $k\times d$ matrix $B_*$ with Gram matrix
$F_kP/n$, aligned with $B$ by singular value decomposition. This is
possible even when $B$ is rank deficient because $d\ge n$: complete
the right singular vectors in its null directions. If $s_1,\ldots,s_n$
are its singular values and $a=\sqrt{F_k/n}$, then

\[
\begin{aligned}
\|B-B_*\|_F^2
 &=\sum_j(s_j-a)^2
 \le\frac n{F_k}\|BB^T-F_kQ_*\|_F^2\\
 &=\frac n{F_k}\left(W^2\|Q-Q_*\|_F^2+\delta^2/n\right)
 \le L_k\delta,\qquad L_k=1+nF_k/\beta_k.
\end{aligned}
\tag{M2}
\]

Here the middle cross term is zero since $Q-Q_*$ is centered and
trace zero; the last bound uses $W\le F_k$ and $0\le\delta\le F_k$.
For $W=0$, any such $B_*$ satisfies (M2), since its squared distance
to $B=0$ equals $F_k=\delta$.

The central winning partition $T$ for scores $b_i^*\cdot x$ has actual
moment matrix $B_*$, by the regular flux identity in the external
source. Its expected winning score is therefore $F_k$. The regret of
the original labels against these regular scores satisfies the exact
identity

\[
r=\int\left(\max_i b_i^*\cdot x-b_{I_A(x)}^*\cdot x\right)d\gamma_d
=F_k-\langle B_*,B\rangle
=\tfrac12(\delta+\|B-B_*\|_F^2)
\le\tfrac12(1+L_k)\delta.
\tag{M3}
\]

This uses the labels of $A$ themselves. It does not replace $A$ by
its own affine winning partition or require price stability for that
replacement.

## Gaussian margin and partition distance

The difference of any two regular scores is a centered normal with
standard deviation $\sqrt{2F_k/n}$. If $g(x)$ is the gap between the
largest and second-largest scores, a union bound gives, for every $u>0$,

\[
\gamma_d(g\le u)\le M_ku,\qquad M_k=q\sqrt{n/(\pi F_k)}.
\tag{G1}
\]

Indeed the gap event is contained in the union of
$\{|(b_i^*-b_j^*)\cdot x|\le u\}$, each of probability at most
$\sqrt{n/(\pi F_k)}u$. Ties are null.
On the changed-label set with $g>u$, the score loss is at least $u$.
If its total measure is $\eta$, (G1) and Markov's inequality give
$\eta\le M_ku+r/u$. Minimizing in $u$, or taking $u\downarrow0$
if $r=0$, yields $\eta\le2\sqrt{M_kr}$.

Consequently the following **partition stability estimate** holds:

\[
\boxed{\quad d_k(A):=\inf_T\sum_i\gamma_d(A_i\triangle T_i)
 \le K_k\sqrt{F_k-W_1(A)},\qquad
 K_k=4\sqrt{\frac{M_k(1+L_k)}2}.\quad}
\tag{G2}
\]

The factor two between total symmetric difference and changed labels is
included. $T$ ranges over orthogonal images and relabelings of the
central regular-simplex partition, extended in unused coordinates.
The proof covers arbitrary measurable competitors, including all remote
tails. The constants depend on $k$ alone, are explicit and very far
from optimal. The exponent is an upper-bound result; optimality of the
exponent and constants is not claimed.

## Consequences for previous four-cell obligations

At $k=4$, the external endpoint resolves the global first-Hermite
inequality, its equality classification and the
[balanced quantization endpoint](BALANCED_GAUSSIAN_QUANTIZATION.md): optimal
distortion is $d-4F_4$, and equality uses the tetrahedral encoder with
reconstructions $4b_i(T)$. It also makes every planar four-cell value
strictly less than $F_4$, including measurable tie subdivisions.

For a balanced affine score tuple $z$, let $A(z)$ be an optimal
balanced assignment and $B(z)$ its actual moment tuple. The existing
dual objective satisfies

\[
J(z)=2\langle z,B(z)\rangle-\|z\|_F^2
    =W_1(A(z))-\|z-B(z)\|_F^2\le F_4.
\tag{H1}
\]

If $J(z)\ge F_4$, equality forces $z=B(z)$ and the tetrahedral
moment tuple. Thus there are no non-tetrahedral high-value centroidal
stationary tuples. The high-value instability, stable equal-radius and
unequal-radius Hessian endpoints are discharged by this empty-domain
argument; no new Hessian sign computation is asserted. Classification
of all lower-valued full-rank stationary tuples remains a separate,
lower-priority open question.

The published
[Heilman, *Stable Gaussian minimal bubbles*, Theorem 1.6](https://link.springer.com/article/10.1007/s00526-026-03408-z)
was inspected together with Problem 1.15 and Conjecture 1.16. Equal masses
give translation $w=0$. The external first-moment endpoint supplies
exactly this balanced instance of its first-moment hypothesis. In the
surface functional and normalization printed in that theorem, its
conclusion is a lower bound by a positive constant $\varepsilon_k$
times $F_k-W_1(A)$.

One can use Theorem 1.5 directly at the balanced mass vector, avoiding
any assumption about other mass vectors. Let $\mathcal P_H$ be the
interface functional printed in its Problem 1.4, and put
$\mathcal B_H=2\mathcal P_H$, the sum of individual surface terms.
That theorem makes the simplex a global minimizer of
$\mathcal P_H+\varepsilon_k\sqrt{\pi/2}\,W_1$. Thus

\[
\mathcal B_H(A)-\mathcal B_H(T)
\ge 2\varepsilon_k\sqrt{\pi/2}\,(F_k-W_1(A))
\ge \varepsilon_k(F_k-W_1(A))
\ge \varepsilon_k d_k(A)^2/K_k^2 .
\tag{P1}
\]

This recovers the needed balanced specialization, with room in the
constant, without claiming its unequal-mass first-moment hypothesis.
For four cells its stated
choice is $\varepsilon_4=1/400$; for $k\ge5$ its positive constant
is existential and dimension independent. This specialization does not
discharge the unequal-mass conjecture. The normalization in this paragraph
is Heilman's printed surface functional, not a conversion to the ordinary
interface sum used in (E2).

## Localization and the remaining full-kernel theorem

For balanced four-cell partitions, write

\[
S_\rho(A)=\tfrac14+\rho W_1(A)+R_\rho(A),\qquad
0\le R_\rho(A)\le\rho^2(3/4-W_1(A))\le3\rho^2/4.
\]

This follows from nonnegative Hermite masses and Parseval, for $0<\rho<1$.
If $S_\rho(A)>S_\rho(T)$, then

\[
0\le F_4-W_1(A)<3\rho/4,\qquad
d_4(A)<K_4\sqrt{3\rho/4}.
\tag{N1}
\]

The established
[Heilman–Tarter, *Three candidate plurality is stablest for small
correlations*, Forum of Mathematics, Sigma 9 (2021), e65,
Theorem 1.9 and Lemma 2.3](https://doi.org/10.1017/fms.2021.56)
give existence of a noise-stability maximizer at each fixed $0<\rho<1$
and, for four cells in $d\ge3$, its cylindrical reduction to three
dimensions after rotation. The exact hypotheses, theorem statement and
existence lemma were inspected. This is a theorem about maximizers; it
does not reduce every nearby measurable competitor to three dimensions
while preserving its objective.

A precise sufficient **nearby full-kernel obligation** is:
there exist $\eta_0>0$ and $\rho_1>0$ such that every balanced
measurable four-cell partition of $\mathbb R^3$ with $d_4(A)\le\eta_0$
satisfies $S_\rho(A)\le S_\rho(T)$ for every
$0<\rho\le\rho_1$. Require $\rho_1<1$. Neither remote tails nor
asymmetric cyclic exchanges are excluded from this statement.

Given this obligation, take

\[
\rho_0=\min\{\rho_1,\;4\eta_0^2/(3K_4^2)\}>0.
\tag{N2}
\]

An improving competitor in any dimension would have an improving
maximizer, hence a three-dimensional improving maximizer. Equation (N1)
puts it in the neighborhood covered by the obligation, a contradiction.
This proves the desired uniform interval conditional on that one
remaining theorem. No quantitative loss proportional to
$\rho(F_4-W_1(A))$ is required or asserted.

The [reflected-swap and equivariant-cycle kernel results](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md)
are available restricted comparisons. Their counterexamples to uniform
relative Hermite control remain valid arbitrarily close to a tetrahedron.
Metric localization supplies no repair of that retired estimate, and
pair comparisons cannot be iterated on already changed cells.

The maximizer-rigidity route is developed in the
[maximizer-rigidity programme](MAXIMIZER_PERTURBATIVE_RIGIDITY.md): retain
global maximality and prove uniform local uniqueness modulo rotations.
Its implication (R1) uses the same localization and dimension reduction.
The new note proves the nonlinear assignment condition, normalized-price
and bounded-interface estimates, an expanding-ball kernel estimate, and a
spectral gap for the limiting first-moment interface form modulo rotations.
Combining those inputs sharpens actual-maximizer localization to
$d_4(A)=O(\rho)$ and $F_4-W_1(A)=O(\rho^2)$ in its (S1)--(S5).
The [normal-flux compactness proof](SMALL_POSITIVE_CORRELATION_INTERVAL.md)
now supplies a relative exterior estimate and the nonlinear rigidity argument,
under the same first-moment inputs. Its interval is qualitative; the new
argument remains internally reviewed only.
The nearby comparison above remains an alternative sufficient obligation;
its inequality alone would not classify possible tied maximizers.

## Literature and assurance limits

The source review inspected Zhang's frozen Sections 1–5, the exact
Milman–Neeman statement and normalization, Heilman's published balanced
specialization, and the Heilman–Tarter statements above. A focused web
screen for Gaussian equal-cell first-moment quantitative partition
stability found no matching theorem among the inspected sources. This
is not a worldwide literature or priority certification. Zhang explicitly
does not assert a metric stability rate in the frozen source.

The proofs (C1)–(G2), (H1) and (N1)–(N2) are analytic written deductions.
No finite sample, floating computation, repository consistency check or
partial upstream Lean replay proves them. Independent mathematical
review, full formal verification and external review remain unperformed.
The next verification task is independent review of the
[small-correlation proof](SMALL_POSITIVE_CORRELATION_INTERVAL.md), especially
its relative tail control, normal-flux lower bound and exact moment gauge.
Comparison at arbitrary positive correlation and the general asymmetric-cycle
problem remain outside that interval theorem.
