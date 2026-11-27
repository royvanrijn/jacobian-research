# Balanced four-class Gaussian polynomial partitions

Every balanced four-cell Gaussian partition with equal-norm centroid vectors,
other than a regular tetrahedral partition up to rotation, cell relabeling and
null sets, has strictly less noise stability than the tetrahedral partition for
all sufficiently small positive correlations. The threshold below depends on the
individual partition. The proof uses the four-variable Gaussian expected-maximum
theorem of Sun, Hu and Lan. A numerical search in a fifteen-parameter family of
balanced polynomial partitions produced 23 distinct fitted candidates across four
positive correlations. All 23 had lower held-out estimates than the tetrahedron.

## Problem and scope

Let $X,Y$ be jointly standard Gaussian vectors with cross covariance $\rho I$.
For a measurable partition $A=(A_i)_{i=1}^q$, put
$S_\rho(A)=\sum_i\Pr(X,Y\in A_i)$. The discrete problem asks whether plurality
maximizes stability at positive correlation for balanced low-influence
functions $[q]^n\to\Delta_q$. [ProofAtlas #468](https://www.proofatlas.ai/open-problems/)

Isaksson and Mossel identify the plurality limit with the simplex stability
and transfer any strict equal-measure Gaussian improvement at $0<\rho<1$ to a
counterexample to the discrete clause. A partition in dimension less than
$q-1$ can be extended with unused Gaussian coordinates.
[Isaksson and Mossel, Theorems 1.10 and 1.11, Lemma 7.3 and Theorem 7.4](https://arxiv.org/pdf/0903.3362)

Heilman proves the Gaussian inequality for three equal-measure cells at
$0<\rho\le 0.1082$. His theorem for three or four cells at positive
correlation assumes hyperstability of the equal-measure maximizer.
[Heilman, arXiv:2306.03312v3, Theorem 1.5](https://arxiv.org/pdf/2306.03312v3)
[Heilman, arXiv:2209.11216, Theorem 1.5 and Definition 1.3](https://arxiv.org/html/2209.11216)

Known counterexamples with prescribed measures use cells of unequal measures.
[Heilman, Mossel and Neeman](https://arxiv.org/abs/1403.0885)

## Hermite bounds

For normalized probabilists' Hermite polynomials $H_\alpha$ in independent
standard Gaussian coordinates, set
$W_k(A)=\sum_i\sum_{|\alpha|=k}(\int_{A_i}H_\alpha\,d\gamma)^2$.
The Gaussian generating function gives
$\mathbb E[H_\alpha(X)H_\beta(Y)]=1_{\alpha=\beta}\rho^{|\alpha|}$.
Orthogonal expansion and Parseval yield, for equal measures and $0<\rho<1$,

$$
S_\rho(A)=\frac1q+\sum_{k\ge 1}\rho^k W_k(A),\qquad
0\le S_\rho(A)-\frac1q-\sum_{k=1}^D\rho^k W_k(A)
\le \rho^{D+1}\left(1-\frac1q-\sum_{k=1}^D W_k(A)\right).
$$

In particular, writing $b_i=\int_{A_i}x\,d\gamma(x)$ gives
$W_1(A)=\sum_i |b_i|^2$. Suppose $q=4$ and each $|b_i|$ has the same
value $r$. When $r>0$ put $u_i=b_i/r$ and write $L(X)$ for the cell label.
Then

$$
4r=\mathbb E\langle u_{L(X)},X\rangle
\le\mathbb E\max_i\langle u_i,X\rangle\le M,
$$

where $M$ is the expected maximum of four centered jointly Gaussian random
variables with unit variance and cross covariances $-1/3$. The last
inequality holds for every positive semidefinite covariance matrix with unit
diagonal, with equality only at these simplex covariances.
[Sun, Hu and Lan, arXiv:2008.04827v2, Theorem 1.1](https://arxiv.org/pdf/2008.04827v2)

A tetrahedral partition $T$
has centroids $Mv_i/4$, where its unit vertices are $v_i$. Indeed, symmetry
makes each centroid a common scalar multiple of its vertex, and the cell
assignment attains the pointwise maximum of the four linear scores.
It follows that $W_1(A)\le W_1(T)=M^2/4$. This inequality is strict for
$r=0$. For $r>0$, equality requires both regular tetrahedral $u_i$ and
assignment to the maximizing score almost everywhere. Each pair of distinct
linear scores ties only on a Gaussian-null hyperplane.

The first two masses of $T$ are

$$
W_1(T)=\frac{12}{\pi}
\left(\frac14+\frac{\arcsin(1/3)}{2\pi}\right)^2,
\qquad W_2(T)=\frac1{\pi^2}.
$$

To compute $M$, use the four variables $(2/\sqrt3)(Z_i-\bar Z)$ with
independent standard normals $Z_i$. Their expected maximum is
$(2/\sqrt3)m_4$, where $m_4=\mathbb E\max_i Z_i$.
Writing $\phi,\Phi$ for the univariate normal density and distribution
function, integration by parts gives
$m_4=4\int x\phi(x)\Phi(x)^3\,dx=12\int\phi(x)^2\Phi(x)^2\,dx$.
Here $\phi^2$ is $1/(2\sqrt\pi)$ times the $N(0,1/2)$ density. If
$V\sim N(0,1/2)$ and $U_1,U_2$ are independent standard normals independent of
$V$, then $\mathbb E\Phi(V)^2=\Pr(U_1\le V,U_2\le V)$.
The differences $V-U_i$ have correlation $1/3$, so the planar Gaussian
quadrant formula gives $1/4+\arcsin(1/3)/(2\pi)$. That formula follows by
measuring the angle between the two defining half-planes for a rotationally
invariant planar Gaussian.

For $W_2$, take the vertices proportional to the four even sign vectors in
three dimensions. In the cell of $(1,1,1)$, coordinate permutation symmetry
and independence of radial and angular parts give
$\int(x_j^2-1)\,d\gamma=0$ for each $j$, with the integral over this cell.
Write $a_j=|x_j|$. Given the magnitudes, the signs producing this cell are
all plus and the pattern with one minus at the index of the smallest
magnitude. Thus

$$
\int_{T_{+++}}x_1 x_2\,d\gamma
=\frac14\mathbb E[a_1a_2\,1_{a_3<\min(a_1,a_2)}]
=\frac1{2\pi\sqrt3}.
$$

The last equality uses
$\mathbb E[a_1\,1_{a_1>t}]=\sqrt{2/\pi}\,e^{-t^2/2}$ and integrates its
square against the half-normal density of $a_3$. Each of the three
off-diagonal Hermite coefficients has this magnitude in each cell.
Consequently $W_2(T)=4\cdot3/(2\pi\sqrt3)^2$.

For any balanced partition $A$ with $W_1(A)<W_1(T)$, the degree-one upper
bound for $A$ and degree-two lower bound for $T$ now give

$$
S_\rho(T)-S_\rho(A)\ \ge\
\rho\bigl(W_1(T)-W_1(A)\bigr)
-\rho^2\bigl(3/4-W_1(A)-W_2(T)\bigr).
$$

This is positive for
$0<\rho<(W_1(T)-W_1(A))/(3/4-W_1(A)-W_2(T))$.
The denominator is at least $W_1(T)-W_1(A)$ by Parseval.
For example, all balanced four-cell partitions whose cells are invariant
under $x\mapsto -x$ have $W_1(A)=0$, which gives a common threshold of
approximately $0.544498$. The displayed bound approaches zero when
$W_1(A)$ approaches $W_1(T)$ from below.

A uniform replacement by a relative higher-Hermite error is obstructed by
[exactly balanced tail swaps](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md).
They have $\delta(A)=W_1(T)-W_1(A)>0$ tending to zero while
$(W_2(A)-W_2(T))/\delta(A)$ is unbounded. Even at a fixed positive
correlation, no positive uniform lower bound proportional to
$\rho\delta(A)$ holds for all balanced partitions. Nevertheless every
proper reflected swap strictly loses full noise stability for every
$0<\rho<1$, by the [full-kernel identity](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md#proper-reflected-swaps-strictly-lose-full-noise-stability).
A uniform small positive-correlation interval over arbitrary balanced
partitions now has an [internal analytic proof](SMALL_POSITIVE_CORRELATION_INTERVAL.md)
using the external first-moment input. Independent mathematical review is
pending; no conclusion for all positive correlations follows.
The [external all-cell first-moment theorem and quantitative partition
estimate](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md) now localize every
improving competitor to a neighborhood of a tetrahedron. The nearby
full-kernel comparison remains the required new theorem.

## Polynomial family and numerical comparison

Index the classes by $s\in\{-1,1\}^3$ with $s_1s_2s_3=1$. For $j=1,2,3$,
use cyclic coordinate indices and put

$$
F_j(x)=a_jx_j+b_jx_{j+1}x_{j+2}
       +c_jx_j(x_j^2-3)+d_jx_j(x_{j+1}^2-1)+e_jx_j(x_{j+2}^2-1),
\qquad p_s(x)=\sum_j s_jF_j(x).
$$

The cell label maximizes $p_s(x)$. Acting on $x$ by an even coordinate sign
change $g$ sends $F_j(x)$ to $g_jF_j(x)$ and permutes labels by $s\mapsto gs$.
There are no identically tied pairs if at least two of the $F_j$ are nonzero.
Indeed, each pair of labels differs in two coordinates and those polynomials
have disjoint monomial supports. The zero set of each nonzero score difference
has Gaussian measure zero. The Gaussian measure is invariant under even sign
changes, so all four cells have measure exactly $1/4$, and their centroids have
equal norms. Setting $F_j=x_j$ gives the tetrahedral baseline.

The program screened 363 coefficient
lists per correlation, including the baseline, two controls, shared coefficient
perturbations and four scales of seeded anisotropic perturbations. It then
used a separate, fixed sample for coordinate refinement from the baseline, both
controls and the three strongest screened lists. Each refinement tested both
signs of each coefficient step at four dyadic step sizes. The final sample was
separate from screening and refinement. Within each sample the candidate
and baseline used the same Gaussian pairs. The reported difference is the
mean of the per-pair difference of agreement indicators, with a paired sample
standard error. [Search program](search_polynomial_partitions.py)

| $\rho$ | Refined lists tested | Best refinement count vs. baseline | Holdout gap for best by refinement | Paired standard error | Largest holdout gap |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.10 | 6 | 9531 vs. 9445 | −0.000936 | 0.000241 | −0.000124 |
| 0.35 | 6 | 12898 vs. 12781 | −0.000830 | 0.000240 | −0.000091 |
| 0.65 | 6 | 17868 vs. 17820 | −0.001408 | 0.000276 | −0.000436 |
| 0.90 | 5 | 24677 vs. 24574 | −0.000777 | 0.000190 | −0.000563 |

The screening, refinement and holdout sizes were $2^{14}$, $2^{15}$ and
$2^{20}$ pairs respectively. As controls, $F=(x_1,x_2,0)$ and
$F=(x_2x_3,x_3x_1,x_1x_2)$ have exact stabilities
$(1+m)^2/4$ and $(1+3m^2)/4$, where $m=(2/\pi)\arcsin\rho$ is the
correlation of two coordinate signs. The first label agrees when the first
two signs each agree. The second label agrees when the three coordinate signs
either all agree or all switch. Both controls were within 1.40 binomial
standard errors of their exact stabilities at each correlation.

From the repository root, run

    python3 research/plurality-stability/search_polynomial_partitions.py

The recorded run used Python 3.12.3, NumPy 1.26.4 and seed 468. Its JSON
output includes the fitted coefficients, agreement counts, wins and losses
for replay. The holdout differences are numerical estimates. Several are
comparable to their standard errors, and the tested coefficients and
correlations form a finite sample.

## Certification condition

A strict improvement is certified when the lower Hermite bound for a specific
candidate exceeds the upper bound for the tetrahedron, using certified
coefficient intervals. For rational polynomial scores that are pairwise
distinct, one can truncate to a Gaussian box and bound each missing coefficient by
Cauchy–Schwarz. Inside the box, interval bounds on polynomial differences
label subboxes whose winning score is certified. Unresolved subboxes can be
charged against a bound on the absolute Hermite polynomial on the box.
The unresolved volume tends to zero under subdivision because every tie
set is closed and Gaussian-null, and each point away from ties has a
strict score margin. On the resolved subboxes the Hermite integrals factor
by coordinate. In positive univariate Hermite degree, the primitive comes
from $(\phi\,\mathrm{He}_{k-1})'=-\phi\,\mathrm{He}_k$, where
$\mathrm{He}_k$ is the monic probabilists' polynomial. In degree zero the
primitive is the Gaussian distribution function.

The Hermite truncation can certify a strictly separated pair when the
truncation degree, box, mesh and precision are increased. Rational parameters
suffice in this family for any hypothetical strict improvement with
pairwise distinct scores. For almost every Gaussian pair both argmax labels
are locally constant under coefficient perturbations, and writing
$Y=\rho X+\sqrt{1-\rho^2}Z$ for independent $X,Z$ shows that the same reasoning
applies to correlation perturbations. Dominated convergence then preserves
the strict gap in a neighborhood. The numerical search supplies no positive
held-out gap for that test.
