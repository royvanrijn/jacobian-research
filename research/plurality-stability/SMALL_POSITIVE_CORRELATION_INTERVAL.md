# A small positive-correlation interval from Gaussian normal-flux compactness

This note gives an internal written proof of the following Gaussian statement,
using the repository's imported equal-cell first-moment theorem and its
quantitative moment-alignment deduction. The proof has not received independent
mathematical review or formal verification. Its most important new steps are
the relative exterior estimate (C4) and the normal-flux lower bound (C11).
No explicit numerical correlation threshold is claimed.

The [fixed-cell extension](#extension-to-every-fixed-cell-count) below proves
the same conclusion for **every fixed $k\ge2$**, with a threshold depending
on $k$ but independent of the ambient dimension $d\ge k-1$. It includes the
all-cell interface-kernel argument and the dimension-dependent tail estimates;
the four-cell proof is retained explicitly as its original specialization.

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

## Extension to every fixed cell count

**Fixed-cell small-correlation theorem.** For every integer $k\ge2$ there
exists $\rho_k>0$ such that, for every $d\ge k-1$, every ordinary measurable
partition $A=(A_1,\ldots,A_k)$ of $\mathbb R^d$ with
$\gamma_d(A_i)=1/k$ satisfies

\[
 S_\rho(A)\le S_\rho(T_k)\qquad(0<\rho\le\rho_k).
 \tag{A1}
\]

Equality holds precisely for central regular-simplex partitions and their
cylindrical extensions, up to orthogonal maps, labels and null sets.
Here $\rho_k$ may depend on $k$; there is no claim of one positive interval
uniform over unbounded cell counts. The proof is internal and analytic,
under the same accepted external first-moment input as the four-cell result.
Independent mathematical review remains open.

Put $n=k-1$, $F=F_k$ and $\ell=\sqrt{2F/n}$. Work first in dimension
$n$. Let $b_i^*$ be the regular cell moments, so their matrix has Gram
matrix $FP/n$, and put
$n_{ij}=(b_i^*-b_j^*)/\ell$. Define $H_{ij}$, $F_{ij}$ and
$d\nu_{ij}=\phi_n\,d\mathcal H^{n-1}$ as above. In dimension one the
pair plane is a point, with its zero-dimensional Gaussian surface measure.
The following checks supply each ingredient of (C1)--(C12) for general $k$.

### The all-cell first-moment interface kernel

Use

\[
\begin{split}
\mathcal H&=\bigoplus_{i<j}L^2(F_{ij},\nu_{ij}),\\
\mathcal H_0&=\left\{f:\ \sum_{j\ne i}\int f_{ij}\,d\nu_{ij}=0
                                      \text{ for every }i\right\},\\
(\mathcal Df)_i&=\sum_{j\ne i}\int x f_{ij}(x)\,d\nu_{ij}(x),\qquad
Q_k(f)=\ell\|f\|_{\mathcal H}^2-\|\mathcal Df\|_F^2.
\end{split}
\]

Speeds are extended antisymmetrically. Positive $f_{ij}$ is outward
motion of cell $i$, whose inward normal is $n_{ij}$.
The same shape derivative as (H1) gives $-W_1''/2=Q_k(f)$ on
$\mathcal H_0$. Here is a justification that does not assume a
classification of stationary simplex shapes. Start with smooth speeds
compactly supported on facet interiors. Disjoint tubular neighborhoods
give actual interface motions. Correct their second-order volume errors
by $k-1$ additional facet bumps supported on a spanning tree of the
complete label graph. Their mass derivative is an incidence matrix of
rank $k-1$, so the implicit function theorem gives exactly balanced
deformations without changing the initial speed. For $k=2,n=1$,
$\mathcal H_0=\{0\}$ and the assertions are immediate.

For the other cases the moment derivative is $\mathcal Df$ and
the exact second derivative of $W_1$ gives the displayed form: in the
score term each interface contributes
$-\ell\int f_{ij}^2d\nu_{ij}$ to $W_1''/2$; the derivative of the
moment square contributes $\|\mathcal Df\|_F^2$. The acceleration terms
cancel because the original regular scores agree on each interface.
The imported global first-moment theorem gives $Q_k\ge0$.

Let $\mathcal M$ be the regular moment orbit and $\mathcal V$ its
tangent space at $B_*$, with rows $\Omega b_i^*$ for $\Omega^T=-\Omega$.
Apply the all-cell quantitative moment estimate (M2) to these exactly
balanced deformations and divide by the square of the deformation
parameter. Since distance to the smooth orbit has first-order tangent
distance, this yields

\[
 \operatorname{dist}(\mathcal Df,\mathcal V)^2\le L_kQ_k(f).
 \tag{A2}
\]

Truncate and smooth arbitrary $\mathcal H_0$ speeds on facet interiors,
then correct their converging linear mass errors by the same tree bumps.
This proves density and extends nonnegativity and (A2) to all speeds.
All facet second moments are finite, so $\mathcal D$ is continuous.

Rotations have speeds $f^\Omega_{ij}=-(\Omega x)\cdot n_{ij}$,
belong to $\mathcal H_0$, have zero $Q_k$, and satisfy
$(\mathcal Df^\Omega)_i=\Omega b_i^*$. Gaussian flux or truncation of
the actual rotational deformation justifies these identities on the
unbounded facets. If $Q_k(f)=0$, (A2) supplies a rotational speed $r$
with $\mathcal Dr=\mathcal Df$. Nonnegativity makes the polarization
against $r$ vanish; hence
$0=Q_k(f-r)=\ell\|f-r\|^2$. The kernel is therefore exactly

\[
 \mathcal R_k=\{f^\Omega:\Omega^T=-\Omega\},\qquad
 \dim\mathcal R_k=n(n-1)/2.
 \tag{A3}
\]

The moment vectors span $\mathbb R^n$, proving injectivity of this
rotation map. Finally $Q_k$ on $\mathcal H_0$ is $\ell I$ minus a
finite-rank operator; its nonzero eigenvalues are bounded below away
from zero. Thus it also has a gap modulo $\mathcal R_k$. The closing
argument needs only nonnegativity and (A3), not a numerical gap.

### Saturated margins for every simplex

The transposition $i\leftrightarrow j$ is the orthogonal reflection
$G_{ij}$ in $H_{ij}$; it interchanges $T_i,T_j$ and fixes the remaining
vertices. Write $v_i=b_i^*/|b_i^*|$. Every inward normal of $T_i$
satisfies

\[
 n_{ij}\cdot v_i=\sqrt{\frac{k}{2(k-1)}}\ge1/\sqrt2.
\]

Consequently $D_i=B(v_i,1/4)$ is contained in $T_i$, and
$z\cdot n_{ij}\ge a:=1/\sqrt2-1/4>0$ on $D_i$.
Its Gaussian mass is at least
$p_n=\operatorname{vol}_n(B(0,1/4))(2\pi)^{-n/2}e^{-25/32}>0$.
For $m=\rho x/\sqrt{1-\rho^2}$ the reflection identity (T3), now in
$n$ dimensions, has nonnegative integrand whenever $x\in T_i$.
Restricting it to $m+D_i\subset T_i$ proves

\[
 U_{\rho,ij}(x)\ge
 p_n(1-e^{-2a})\min\{\rho n_{ij}\cdot x/\sqrt{1-\rho^2},1\}.
 \tag{A4}
\]

Indeed the reflected density ratio is
$e^{-2m_N(m_N+z_N)}\le e^{-2am_N}$, and concavity of
$1-e^{-2au}$ on $[0,1]$ gives the stated bound. Thus (C1) holds with
$c_T$ replaced by the fixed positive constant in (A4), for the
$k(k-1)$ ordered transfers. Also $U_{\rho,ij}=0$ on the whole pair
plane, and

\[
 \rho^{-1}\nabla U_{\rho,ij}(x)\longrightarrow b_i^*-b_j^*
 \quad\text{locally uniformly as }\rho\downarrow0.
 \tag{A5}
\]

For example differentiate the Mehler kernel at zero, obtaining
$\partial_\rho M_\rho|_0=x\cdot y$; its first $x$ derivatives are
dominated on compact $x$ sets by an integrable polynomial times a
Gaussian in $y$ for $\rho$ in a neighborhood of zero. This proves (A5)
without a bound on remote normal derivatives.

### Weighted tails in dimension $k-1$

Keep $w(x)=(1+|x|^2)e^{|x|^2/16}$ and $R\ge2$. The proof of (C3)
works with $z\in\mathbb R^{n-1}$. More explicitly an admissible
decreasing bound, up to a constant depending only on $n$, is

\[
\begin{split}
\beta_n(R)=C_n\bigg[&\int_{|x|>R}(1+|x|^2)^2e^{|x|^2/8}d\gamma_n(x)\\
 &+\int_{|z|>\sqrt{R^2-1}}
              (1+|z|^2)^2e^{|z|^2/8}d\gamma_{n-1}(z)\bigg]
 \longrightarrow0.
 \tag{A6}
\end{split}
\]

When $n=1$ the second integral is zero. Split $t\ge1$ and $t<1$,
use Cauchy--Schwarz on the former and the bounded one-dimensional
density inequality on the latter, exactly as in (C3). The global
constant comes from the same integrals without tail restrictions.
This establishes (C3) in all fixed dimensions; the polynomial $R^6$
bound stated earlier is specific to dimension three and is not imported
unchanged.

The $n$-dimensional Mehler derivative satisfies
$|L_\rho(x,y)|\le C_nw(x)w(y)$ for $0<\rho\le1/16$.
Indeed
$M_r\le C_n\exp((|x|^2+|y|^2)/34)$ in this range and its derivative
is a quadratic polynomial with uniformly bounded coefficients times
$M_r$. The exponent $1/16$ in $w$ leaves room to absorb that polynomial.
Summing the ordered transfers now proves (C4)--(C5) with constants
depending on $k$ and exterior factor tending to zero. No constant has
to be uniform in $k$.

### Normal-flux contradiction and dimension reduction

If (A1) failed arbitrarily close to zero in dimension $n$, choose
$\rho_j\downarrow0$ and balanced partitions outside the simplex orbit
with $S_{\rho_j}(A^{(j)})\ge S_{\rho_j}(T_k)$. Nonnegative Hermite
masses and Parseval give

\[
 0\le\delta_j=F-W_1(A^{(j)})\le(1-1/k)\rho_j.
 \tag{A7}
\]

Align the nearest regular moment matrix with $B_*$. The all-cell
(M2)--(M3) and (G2) give changed-label measure tending to zero and
the exact condition $B_j-B_*\perp\mathcal V$ simultaneously.
The regret is strictly positive by (A4). With
$s_j^2=\mathcal L_{\rho_j}(A^{(j)})/\rho_j$, (C2), the global
weighted bound and localization give $s_j\to0$, exactly as in (C7).
The measures $\mu_{j,i}=e_i^{(j)}\gamma_n/s_j$ have bounded weighted
variation and weighted exterior variation tending to zero uniformly,
by (A6).

For every ordered pair $i,l$ set

\[
 \theta_{j,il}(z)=\frac1{s_j\phi_1(0)}
     \int_0^\infty h_{il}^{(j)}(z+tn_{il})\phi_1(t)dt.
\]

They are bounded in $L^2(H_{il},\nu_{il})$ by (C1) and the same
one-dimensional mass-cost inequality as (C9). Along a subsequence they
converge weakly to nonnegative $\theta_{il}$. Transfers at any fixed
positive normal distance have mass $O(s_j^2)$, so the normalized measures
concentrate on their pair planes. The limits are supported on $F_{il}$
because transfers originate in $T_i$. Junctions have dimension at most
$n-2$ and are null for $\nu_{il}$; the weak $L^2$ limit gives them
no mass. For $n=1$ there are no junctions. With
$f_{il}=\theta_{li}-\theta_{il}$ for $i<l$, antisymmetrically extended,
we obtain precisely

\[
 \mu_i=\sum_{l\ne i}f_{il}\nu_{il}|_{F_{il}},\qquad
 f\in\mathcal H_0,\qquad\mathcal Df\perp\mathcal V.
\]

Total masses and first moments pass to the limit by weighted tightness,
so these last two conditions are exact. The plane normalization is
$d\nu_{il}=\phi_1(0)d\gamma_{n-1}$, including $n=1$.
This is why the factor in (C11) is independent of dimension.
On compact normal cylinders (A5) gives
$U_{\rho_j,il}(z+tn)\ge\rho_j(\ell-o(1))t$. The bound
$\int_0^bth\phi_1dt\ge(\int_0^bh\phi_1dt)^2/(2\phi_1(0))$,
weak lower semicontinuity, and exhaustion of tangential cutoffs therefore
give

\[
 1\ge\frac\ell2\sum_{i\ne l}\|\theta_{il}\|_2^2
       \ge\frac\ell2\|f\|_{\mathcal H}^2.
\]

No simultaneous normal cylinders need be disjoint: regret is already a
sum of separately indexed ordered transfers, and each summand is
estimated only once. Weighted tightness and local uniform convergence
$L_{\rho_j}\to x\cdot y$ then give (C12) in dimension $n$. Thus

\[
 \|\mathcal Df\|_F^2\ge2\ge\ell\|f\|_{\mathcal H}^2,
 \qquad Q_k(f)\le0.
 \tag{A8}
\]

Nonnegativity and (A3) make $f$ rotational. Its moment lies in
$\mathcal V$ and is orthogonal to $\mathcal V$, contradicting (A8).
This proves rigidity for every sufficiently small positive correlation
in dimension $n$, with one threshold for this fixed $k$.

For arbitrary $d\ge n$, the fixed-mass existence argument below, applied
in dimension $d$, gives a global maximizer (also the conclusion of
Heilman--Tarter Lemma 2.3). Their Theorem 1.9 reduces **each global
maximizer** to $n$ dimensions. In their notation set $m=k$ and $n_{\rm HT}+1=d$;
their condition $m\le n_{\rm HT}+2$ is exactly $k\le d+1$.
A hypothetical improving competitor would yield an improving maximizer,
which the $n$-dimensional rigidity excludes. Every equality case is now
a global maximizer and reduces in the same way. This proves (A1) and
its equality statement; individual competitors were never dimension
reduced. If desired, cylindrical embedding also gives the upper bound
for $1\le d<n$, with strict inequality by moment rank.

### Audit boundary of the fixed-cell extension

The all-cell input is the exact frozen first-moment/covariance preprint,
not an extrapolation of the four-cell endpoint. Its covariance path uses
Milman--Neeman in dimension $k-1$, within the stated cell-count range.
The all-cell moment estimate (M2) supplies (A2); no unproved Hessian
classification or compact stationary-shape search enters. The exact
Heilman--Tarter statement and proof of its existence lemma were inspected
again for general $k$. The displayed convex set in that lemma's proof
omits the fixed-mass equations. For the present application the argument
below proves existence directly on the mass-constrained set and proves
ordinaryness by a strictly positive quadratic variation. This fills in
the required fixed-mass step without importing an assertion about the
extreme points of a different feasible set.

The new proof obligations for independent review are (A2)--(A3), the
reflection ball in (A4), the dimension-dependent tails (A6), and their
applications to (C9)--(C12). The $k=2$ zero-dimensional facet and trivial
balanced speed space have been treated explicitly. No finite computation
or source-consistency check establishes the analytic theorem, and no
all-positive-correlation conclusion or publication is asserted.

## A rigorous continuation boundary in correlation

For fixed $k$, work in $n=k-1$ dimensions and let $V_k(\rho)$ be the
maximum balanced noise stability. Write
$G_k(\rho)=V_k(\rho)-S_\rho(T_k)\ge0$. The following facts constrain
a continuation argument; they do not rule out a transition.

**Compactness of maximizers and closed optimality.** On every compact
subinterval of $(0,1)$, $V_k$ is Lipschitz, the set $\{G_k=0\}$ is
closed, and any sequence of global maximizers with correlations tending
to an interior value has a subsequence converging in $L^2(\gamma_n)^k$
to an ordinary global maximizer there.

The existence, Lipschitz and compactness assertions hold in every fixed
finite dimension $d\ge1$ by the same proof, replacing $n$ by $d$.
For continuation of the simplex comparison we use $d=n=k-1$.

To check existence and the last assertion with the balance constraints
explicit, use the weakly compact convex set

\[
 \mathcal K=\{(u_i):0\le u_i\le1,\ \sum_i u_i=1,
                                      \ \int u_i\,d\gamma_n=1/k\}.
\]

On fixed finite-dimensional Gaussian space $T_\rho$ is compact on
$L^2$: degree-$m$ Hermite eigenspaces have finite multiplicity and
eigenvalues $\rho^m\to0$. Hence its quadratic functional is weakly
continuous on $\mathcal K$, and a maximum exists. Every maximizer is
ordinary. Otherwise some pair $u_i,u_j$ is bounded away from zero on
a set of positive measure. On that set both are also bounded away from
one. A nonzero bounded function $g$ supported there with zero Gaussian
integral gives two feasible perturbations $u\pm\varepsilon g(e_i-e_j)$.
Their average objective exceeds that of $u$, since
$\langle g,T_\rho g\rangle>0$ for every nonzero $g$, by the strictly
positive Hermite eigenvalues. This contradicts maximality.

The Hermite formula also gives, for $a\le\rho,\tau\le b<1$,

\[
 \|T_\rho-T_\tau\|_{2\to2}
 \le\left(\sup_{m\ge1}m b^{m-1}\right)|\rho-\tau|.
 \tag{A9}
\]

The factor is finite and $\sum_i\|u_i\|_2^2\le1$ on $\mathcal K$.
Taking maxima proves the Lipschitz statement, and therefore closedness.
For correlations $\rho_j\to r\in(0,1)$, pass a maximizing sequence
weakly in $\mathcal K$. Compactness of $T_r$, (A9), and continuity of
$V_k$ make its limit a maximizer at $r$, thus ordinary by the preceding
argument. Both the ordinary sequence and its limit have squared Hilbert
norm one, so weak convergence is strong. This is compactness at positive
correlation; $T_0$ has no strictly positive higher-degree eigenvalues and
the same argument cannot give ordinary strong limits at zero.

Define the endpoint of the initial optimality interval by

\[
 r_k=\sup\{r\in(0,1):G_k(\rho)=0\text{ for every }0<\rho\le r\}.
\]

By (A1), $r_k>0$. If $r_k<1$, continuity proves $G_k(r_k)=0$.
There are correlations $\rho_j>r_k$, $\rho_j\to r_k$, with
$G_k(\rho_j)>0$; choose their maximizers and pass to a strong limit at
$r_k$. Consequently a finite continuation endpoint has at least one of
two mechanisms:

* The limit is a nonsimplex global maximizer tied with $T_k$ at $r_k$.
* The limit is a simplex, and improving maximizers approach its orbit
  as the correlation approaches $r_k$.

The first mechanism is a global crossing that a strict local interface
Hessian at the simplex would not exclude. The second requires a
comparison uniform in both the correlation and arbitrary measurable
nearby partitions, including remote transfers; a smooth fixed-correlation
second variation alone does not supply it. Thus closedness, local
stability, and optimality on the small-correlation interval cannot
justify continuation by themselves. Proving that neither mechanism can
occur, or directly proving a global comparison at every correlation,
remains the unrestricted objective. The large-correlation perimeter
limit must also control moving maximizers and its error relative to
their deformation; a pointwise limit for each fixed partition is
insufficient for that purpose.

The [positive-correlation continuation proof](POSITIVE_CORRELATION_CONTINUATION.md)
extends regret-relative exterior control uniformly to every fixed range
$0<\rho\le b<1$. Under a full-interface form that is nonnegative with
only rotational zero modes at a fixed positive correlation, it proves a
uniform local comparison for arbitrary measurable nearby partitions,
including remote transfers. Consequently a finite first continuation
endpoint has either a nonsimplex tied global maximizer or an additional
nonrotational interface zero mode. Neither possibility is excluded; the
interface-kernel hypothesis throughout $(0,1)$ remains open.
