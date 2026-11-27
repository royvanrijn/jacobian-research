# A small positive-correlation interval from Gaussian normal-flux compactness

This note gives an internal written proof of the following Gaussian statement,
using the repository's imported equal-cell first-moment theorem and its
quantitative moment-alignment deduction. The proof has not received independent
mathematical review or formal verification. Its most important new steps are
the relative exterior estimate (C4) and the normal-flux lower bound (C11).
No explicit numerical correlation threshold is claimed.

**Three-dimensional rigidity.** There exists $\rho_0>0$ such that, for every
$0<\rho\le\rho_0$, every balanced measurable four-cell partition $A$ of
$\mathbb R^3$ satisfying $S_\rho(A)\ge S_\rho(T)$ is a central regular
tetrahedral partition, up to an orthogonal map, labels and null sets.

**Dimension-independent interval.** There exists the same type of constant
$\rho_0>0$, independent of $d\ge3$, such that every balanced measurable
four-cell partition of $\mathbb R^d$ satisfies
$S_\rho(A)\le S_\rho(T)$ for $0<\rho\le\rho_0$.
Every global maximizer is tetrahedral and cylindrical after rotation.

The second assertion uses the published existence and dimension reduction
theorem for maximizers, as audited in the
[localization proof](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#localization-and-the-remaining-full-kernel-theorem).
Publication, external review and formal verification are separate from this
internal theorem status.

## Inputs and relation to prior work

The exact first-moment input is the external Zhang preprint and the written
deductions (M2), (M3) and (G2) in
[quantitative first-moment stability](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md).
In dimension three these give, with $B$ the cell-moment matrix,

\[
0\le\delta=F_4-W_1(A),\qquad
\inf_{B_*}\|B-B_*\|_F^2\le L_4\delta.
\]

The central tetrahedron corresponding to the closest regular moment matrix
$B_*$ satisfies
$\gamma_3(I_A\ne I_T)\le C\sqrt\delta$ by (M3)--(G2).
Thus moment alignment and partition localization can be imposed together.
The equality case and moment estimate are imported with their existing
internal source-audit limits; this note supplies no independent validation
of that preprint.

We also use the following results from
[the maximizer/interface proof](MAXIMIZER_PERTURBATIVE_RIGIDITY.md):

* The tetrahedral score differences have the global saturated margin (T1).
* The balanced first-moment interface form (H1)--(H3) is nonnegative and
  its kernel consists exactly of rotations.
* Pairwise reflection makes the tetrahedral OU score difference zero on
  the entire pair plane; its normalized gradient converges locally
  uniformly to the regular moment difference as $\rho\downarrow0$.

Only the kernel characterization of the first-moment form is needed below;
the full-noise Hessian estimate and the expanding-domain interface chart
are additional results, not premises of this closing argument.

[Heilman, *Euclidean Partitions Optimizing Noise Stability*,
arXiv:1211.7138v2](https://arxiv.org/html/1211.7138v2), Section 8 and
Conjecture 3, identifies four-cell first-moment stability as a missing input
for a possible perturbative extension. Its stated three-cell iterative
lemma is not imported for four cells here.

The use of measures concentrated on switching interfaces is an established
second-order method. The closest inspected general comparison is
[Christof--Wachsmuth, *No-Gap Second-Order Conditions via a Directional
Curvature Functional*, arXiv:1707.07579v2](https://arxiv.org/html/1707.07579v2),
Theorems 4.4, 6.4 and 6.12 and Assumptions 6.2 and 6.6.
The explicit bang-bang result assumes a bounded Lipschitz domain, a scalar
box constraint, a regular switching level set and weak-star continuity of
the quadratic term. It is not directly a theorem for an unbounded Gaussian
multicell partition with balance and rotational degeneracy. The weighted
tail argument, ordered transfer measures and moment gauge below establish
the needed assertions directly. No general no-gap theorem is used as an
unverified implication, and no priority certification is asserted.

## Notation and the exact comparison

Fix $T$ and let $b_i^*$ be its moment vectors. Put
$\ell=|b_i^*-b_j^*|=\sqrt{2F_4/3}$ and
$n_{ij}=(b_i^*-b_j^*)/\ell$, pointing into $T_i$.
Write $H_{ij}=n_{ij}^\perp$, $F_{ij}=T_i\cap T_j$, and
$d\nu_{ij}=\phi_3\,d\mathcal H^2$ on that plane or its facet.
In coordinates $x=z+tn_{ij}$,
$d\gamma_3=\phi_1(t)\,dt\,d\gamma_2(z)$ and
$d\nu_{ij}=\phi_1(0)d\gamma_2(z)$.

For an ordinary partition $A$, not necessarily balanced in the tail
estimates, set

\[
h_{ij}=1_{T_i\cap A_j}\quad(i\ne j),\qquad
e_i=1_{A_i}-1_{T_i}=\sum_{j\ne i}(h_{ji}-h_{ij}).
\]

Every ordered transfer has $t=n_{ij}\cdot x\ge0$.
Let $U_{\rho,ij}=T_\rho(1_{T_i}-1_{T_j})$ and define the regret

\[
\mathcal L_\rho(A)=\sum_{i\ne j}\int U_{\rho,ij}h_{ij}\,d\gamma_3.
\]

The saturated margin gives, for $0<\rho\le1/16$,

\[
\mathcal L_\rho(A)\ge c_T\rho
 \sum_{i\ne j}\int\min\{n_{ij}\cdot x,1\}h_{ij}\,d\gamma_3.
\tag{C1}
\]

Indeed $\min\{\rho t/\sqrt{1-\rho^2},1\}\ge\rho\min\{t,1\}$.
For balanced $A$, with $\Pi_0f=\int f\,d\gamma_3$,

\[
S_\rho(A)-S_\rho(T)
=\widetilde Q_\rho(e)-2\mathcal L_\rho(A),\qquad
\widetilde Q_\rho(v):=\sum_i\langle v_i,(T_\rho-\Pi_0)v_i\rangle.
\tag{C2}
\]

This follows by quadratic expansion and $\int e_i=0$. The operator
$T_\rho-\Pi_0$ is positive semidefinite by its Hermite eigenvalues.
The subtraction of constants is essential for cut-off pieces, whose
individual cell masses need not vanish.

## Weighted transfer control and the relative exterior estimate

Set $w(x)=(1+|x|^2)e^{|x|^2/16}$.
For every unit vector $n$, every measurable $0\le h\le1$ supported in
$\{t=n\cdot x\ge0\}$, and $R\ge2$, there is a fixed function
$\beta(R)\downarrow0$ after replacing it by a decreasing envelope such
that

\[
\left(\int_{|x|>R}w(x)h(x)\,d\gamma_3\right)^2
\le\beta(R)\int\min\{t,1\}h(x)\,d\gamma_3.
\tag{C3}
\]

One may take an upper bound
$\beta(R)\le C(1+R^6)e^{-3(R^2-1)/8}$; only convergence to zero is used.
There is also a finite constant for the same assertion without the
restriction $|x|>R$.

Here are the details. Split $t\ge1$ and $0\le t<1$.
For the first part, Cauchy--Schwarz gives the square bound
$\int_{|x|>R}w^2\,d\gamma_3$ times its mass, which is already its
$\min\{t,1\}$ cost. For the second, $|x|>R$ implies
$|z|>\sqrt{R^2-1}$ and
$w(z+tn)\le2e^{1/16}(1+|z|^2)e^{|z|^2/16}$.
Put $m(z)=\int_0^1h(z+tn)\phi_1(t)\,dt$.
Since $h\phi_1\le\phi_1(0)$,

\[
m(z)^2\le2\phi_1(0)\int_0^1t h(z+tn)\phi_1(t)\,dt.
\]

For clarity, any density $v(t)\le\phi_1(0)$ of mass $m$ on the
positive halfline satisfies
$\int t v(t)dt\ge\int_0^{m/\phi_1(0)}(m-\phi_1(0)s)ds
=m^2/(2\phi_1(0))$.
Apply Cauchy--Schwarz in $z$. Its other factor is the two-dimensional
tail integral of $(1+|z|^2)^2e^{|z|^2/8}$ against $\gamma_2$.
Both this factor and the three-dimensional $w^2$ tail tend to zero;
radial integration gives the displayed bound on $\beta$.
Adding the two normal ranges changes only the fixed constant.
Omitting their tail restrictions gives the finite global constant.

The normalized Mehler kernel satisfies

\[
L_\rho(x,y):=\frac{M_\rho(x,y)-1}{\rho},\qquad
|L_\rho(x,y)|\le Cw(x)w(y)
\quad(0<\rho\le1/16).
\]

Indeed $L_\rho=\int_0^1\partial_rM_r|_{r=s\rho}\,ds$.
The derivative is a quadratic polynomial with uniformly bounded
coefficients times $M_r$, and
$M_r\le C\exp[(|x|^2+|y|^2)/34]$ for $0\le r\le1/16$.
Also $L_\rho\to x\cdot y$ uniformly on each compact product.

Apply (C3) to the twelve ordered transfers and use (C1).
If $e^{\rm out}=e1_{\{|x|>R\}}$ and $\mathcal L_\rho^{\rm out}$ is
the regret integrated there, this proves

\[
\begin{split}
\left(\sum_i\int_{|x|>R}w|e_i|\,d\gamma_3\right)^2
 &\le C\beta(R)\mathcal L_\rho(A)/\rho,\\
0\le\widetilde Q_\rho(e^{\rm out})
 &\le\varepsilon(R)\mathcal L_\rho^{\rm out},
\qquad \varepsilon(R)\longrightarrow0.
\end{split}
\tag{C4}
\]

For the second line, apply (C3) to the transfers already restricted to
the exterior. The global versions also give
$\widetilde Q_\rho(e^{\rm in})\le C_0\mathcal L_\rho^{\rm in}$.
Positivity of $T_\rho-\Pi_0$ controls the mixed term, so for every $a>0$,

\[
\widetilde Q_\rho(e)
\le\widetilde Q_\rho(e^{\rm in})+a\mathcal L_\rho^{\rm in}
 +\varepsilon(R)(1+C_0/a)\mathcal L_\rho^{\rm out}.
\tag{C5}
\]

These are uniform, deformation-relative exterior bounds for the exact
kernel. They concern regret against tetrahedral OU scores, not the
unrestricted first-moment deficit. The retired relative-Hermite
counterexamples are consistent with them.

## Normalizing a hypothetical improving or tied sequence

Suppose three-dimensional rigidity fails arbitrarily close to zero.
Then there are $\rho_j\downarrow0$ and balanced partitions $A^{(j)}$
outside the tetrahedral orbit with
$S_{\rho_j}(A^{(j)})\ge S_{\rho_j}(T)$.
The absolute Hermite remainder gives
$0\le\delta_j\le3\rho_j/4$.
Choose the closest regular moment matrix and rotate it to the fixed $B_*$.
The imported alignment and regret estimates give

\[
\eta_j:=\gamma_3(I_{A^{(j)}}\ne I_T)\longrightarrow0,
\qquad B_j-B_*\perp\mathcal V,
\quad\mathcal V=\{\text{rows }\Omega b_i^*: \Omega^T=-\Omega\}.
\tag{C6}
\]

The orthogonality is exact: differentiate squared distance to the compact
orthogonal moment orbit at its nearest point. Label permutations are
already realized in that orbit. All objective values and balance are
preserved by this alignment.

Write $\mathcal L_j=\mathcal L_{\rho_j}(A^{(j)})$ and
$s_j^2=\mathcal L_j/\rho_j$.
The saturated margin makes $\mathcal L_j>0$, since the partition is not
$T$ after alignment. Equation (C2) gives
$\widetilde Q_{\rho_j}(e^{(j)})\ge2\mathcal L_j$.
The kernel bound gives

\[
s_j^2\le\frac{\widetilde Q_{\rho_j}(e^{(j)})}{2\rho_j}
\le C\left(\sum_i\int w|e_i^{(j)}|\,d\gamma_3\right)^2
\longrightarrow0.
\tag{C7}
\]

The last limit uses (C6), $|e_i^{(j)}|\le1$ and $w\in L^1(\gamma_3)$.
Define signed measures $\mu_{j,i}=e_i^{(j)}\gamma_3/s_j$.
The global version of (C4) bounds their weighted total variations, and
its exterior version gives the uniform tightness estimate

\[
\sup_j\sum_i\int_{|x|>R}w\,d|\mu_{j,i}|
\le C\sqrt{\beta(R)}\longrightarrow0.
\tag{C8}
\]

This normalization accommodates arbitrarily small or remote deviations;
no lower bound on $s_j$ relative to a power of $\rho_j$ is imposed.

## Ordered normal fluxes and the sharp regret lower bound

Let $\phi_0=\phi_1(0)$. For each ordered pair define a nonnegative
density on its entire pair plane by

\[
\theta_{j,ik}(z)=\frac1{s_j\phi_0}
 \int_0^\infty h^{(j)}_{ik}(z+tn_{ik})\phi_1(t)\,dt.
\tag{C9}
\]

These twelve densities are bounded in their plane $L^2(\nu)$ spaces.
To verify this for the full normal integral, split $t<1$ and $t\ge1$.
The first mass squared is controlled by its first moment as in (C3).
The second mass squared is bounded by a fixed constant times that mass.
Consequently $(\int_0^\infty h\phi_1)^2\le C\int_0^\infty
\min\{t,1\}h\phi_1$.
Now use (C1) and $\mathcal L_j=\rho_js_j^2$.
After a subsequence, write $\theta_{j,ik}\rightharpoonup\theta_{ik}$
in these $L^2$ spaces, with $\theta_{ik}\ge0$.

Only normal distances tending to zero survive this normalization.
For each fixed $b>0$, (C1) gives
$s_j^{-1}\int_{t\ge b}h^{(j)}_{ik}d\gamma_3\le C_bs_j\to0$.
For a compactly supported continuous test function, replacing its value
at $z+tn$ by that at $z$ therefore has vanishing error: first restrict
$t<b$, use uniform continuity and bounded normalized mass, and then let
$b\downarrow0$.
It follows that the ordered transfer measures converge to
$\theta_{ik}\nu_{ik}$ on $H_{ik}$.
This density is supported on $F_{ik}$: any compact subset of
$H_{ik}\setminus F_{ik}$ has positive distance from the closed cone
$T_i$, so a transfer originating in $T_i$ there needs a fixed positive
normal distance and vanishes by the same estimate.
In particular no limiting transfer mass is concentrated on a triple
ray or on the fourfold point, which have plane measure zero.

Set $f_{ik}=\theta_{ki}-\theta_{ik}$ for $i<k$, with antisymmetric
extension. Then

\[
\mu_{j,i}\rightharpoonup
 \mu_i:=\sum_{k\ne i}f_{ik}\nu_{ik}|_{F_{ik}}.
\tag{C10}
\]

The convergence is against compact continuous functions and also passes
total masses and first moments, by (C8). Hence $f\in\mathcal H_0$ and
$\mathcal Df\perp\mathcal V$, using exact balance and (C6).

The precise lower bound, with its required factor one half, is

\[
1\ge\frac\ell2\sum_{i\ne k}\|\theta_{ik}\|_{L^2(\nu)}^2
\ge\frac\ell2\|f\|_{\mathcal H}^2.
\tag{C11}
\]

To prove it, fix a compact tangential cutoff $0\le\chi\le1$ and a
normal cutoff $0<b<1$. Reflection gives $U_{\rho,ik}=0$ on $H_{ik}$.
Local $C^1$ convergence of the normalized tetrahedral scores gives,
uniformly on the resulting compact cylinder,

\[
U_{\rho_j,ik}(z+tn_{ik})\ge
 \rho_j(\ell-o(1))t\qquad(0\le t\le b).
\]

Discard the nonnegative regret outside these cylinders.
For $m_j^b(z)=\int_0^b h^{(j)}_{ik}(z+tn)\phi_1(t)dt$, the
one-dimensional density bound gives
$\int_0^b t h\phi_1dt\ge(m_j^b)^2/(2\phi_0)$.
Divide the resulting regret inequality by
$\rho_js_j^2=\mathcal L_j$.
The densities $m_j^b/(s_j\phi_0)$ have the same weak $L^2$ limit as
the full densities in (C9): the omitted normal mass tends to zero in
$L^1$ after normalization, and the densities are uniformly bounded
in $L^2$. Approximation of an $L^2$ test by bounded tests justifies this
last implication. Weak lower semicontinuity now gives
$1\ge(\ell/2)\sum_{i\ne k}\int\chi\theta_{ik}^2d\nu$.
Let the cutoffs exhaust the planes.
Finally $\theta_{ik},\theta_{ki}\ge0$ implies
$(\theta_{ki}-\theta_{ik})^2\le\theta_{ik}^2+\theta_{ki}^2$,
which proves the second inequality of (C11).

## The interaction limit and exclusion of the sequence

The weighted tightness (C8) and the kernel bound imply

\[
\frac{\widetilde Q_{\rho_j}(e^{(j)})}{\rho_js_j^2}
\longrightarrow\sum_i\left|\int x\,d\mu_i\right|^2
=\|\mathcal Df\|_F^2.
\tag{C12}
\]

Indeed write the left side as the integral of $L_{\rho_j}(x,y)$
against $\mu_{j,i}\otimes\mu_{j,i}$ and sum over cells.
On each compact product the kernel converges uniformly to $x\cdot y$.
The contribution where either variable is outside a ball is at most
a fixed constant times the weighted tail variation in (C8) times the
uniform global weighted variation. It tends to zero uniformly as the
ball grows. The $x\cdot y$ integral is a product of first moments, whose
convergence follows from (C8) and (C10).
This proves (C12) without assuming that normalized transfers have
bounded support or converge in total variation.

By (C2), the left side of (C12) is at least $2$. Together with (C11),

\[
\|\mathcal Df\|_F^2\ge2\ge\ell\|f\|_{\mathcal H}^2,
\qquad Q(f)=\ell\|f\|_{\mathcal H}^2-\|\mathcal Df\|_F^2\le0.
\]

The first-moment interface form is nonnegative on $\mathcal H_0$ and
its kernel is $\mathcal R$. Therefore $f$ is rotational. But
$\mathcal Df\in\mathcal V$ for a rotational speed and
$\mathcal Df\perp\mathcal V$ by the exact moment gauge.
Thus $\mathcal Df=0$, contradicting its squared norm being at least $2$.
The hypothetical sequence cannot exist. This proves the stated uniform
three-dimensional rigidity theorem.

## Passage to arbitrary ambient dimension and proof-audit boundary

Fix the resulting small-correlation interval. If an improving partition
existed in $d\ge3$, existence of a global maximizer and the published
Heilman--Tarter dimension-reduction theorem would give an improving
three-dimensional maximizer. The rigidity theorem just proved excludes
it. Conversely any maximizer in $d\ge3$ reduces to a three-dimensional
maximizer, which must be a tetrahedron. This establishes the announced
dimension-independent Gaussian interval and the maximizer classification.

The argument is qualitative and provides no explicit numerical $\rho_0$.
It neither proves the conjecture at all positive correlations nor proves
a noise-stability deficit proportional to the unrestricted first-moment
deficit. It uses no finite search to pass to the unrestricted measurable
class. The new relative exterior estimate is the reason arbitrarily small
remote transfers cannot disappear from the normalized argument.

The required independent mathematical review should check (C3)--(C4),
the shared moment/metric alignment (C6), normalization (C7), the absence
of limiting junction mass in (C9)--(C10), the factor in (C11), and the
weighted product-kernel limit (C12), together with the exact imported
first-moment and dimension-reduction hypotheses. Repository consistency
checks are not such a review. No submission or release is asserted.
