# Four-cell maximizers: assignment, interface geometry and the limiting gap

The preferred remaining route is uniform rigidity of actual noise-stability
maximizers near the tetrahedron. This is an open theorem. The deductions below
establish its assignment equation, quantitative control of its prices and
interfaces on bounded regions, an estimate on expanding balls, and a spectral
gap for the limiting first-moment interface form modulo rotations. Combining
these steps sharpens maximizer localization to $d_4(A)=O(\rho)$ and
$F_4-W_1(A)=O(\rho^2)$, with original prices $O(\rho^2)$. They do
**not** establish rigidity for the full noise kernel or a positive-correlation
interval. These are internal written arguments using the
[first-moment stability proof](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md),
including its external preprint input; independent mathematical review and
formal verification have not occurred.

Two further deductions control the linear tail issue: the full-noise
interface Hessian at the tetrahedron has a uniform small-correlation gap
modulo rotations on the entire Gaussian facet space, and every maximizer
has one analytic tetrahedral interface chart on a domain containing a ball
of radius $c/\rho$. A global nonlinear chart and a uniform nonlinear
remainder estimate are still missing.

## Closest inspected primary results

[Heilman--Tarter, *Three candidate plurality is stablest for small correlations*,
Forum of Mathematics, Sigma 9 (2021), e65](https://doi.org/10.1017/fms.2021.56),
Lemma 2.3 and Theorem 1.9, supply existence and reduction of a four-cell
maximizer to three dimensions. Their hypotheses are fixed positive correlation,
positive prescribed masses, and ambient dimension at least the number of cells
minus one. This reduces maximizers, rather than arbitrary competitors.

The closest perturbative method is
[Heilman, *Euclidean Partitions Optimizing Noise Stability*,
arXiv:1211.7138v2](https://arxiv.org/html/1211.7138v2), Sections 3--6,
especially Lemmas 3.1 and 6.1 and Theorem 7.1. The argument bootstraps
first-variation information on expanding regions, controlling kernel errors far
from the origin. Crucially its rigidity lemma assumes **three cells**, aligned
pairwise moment directions in one plane, a pointwise discrepancy enclosure,
and the assignment condition (28) for the **correlation derivative** objective
using $LT_\rho$. It is not a four-cell rigidity theorem for $S_\rho$ with
balancing prices. Its statement and the opening proof steps were inspected;
the seven-step estimate is not imported as a proved four-cell implication.
Section 8 explicitly identifies four-cell first-moment stability as a
missing input for a possible extension; it does not execute that extension.

[Heilman, *Hyperstable Sets with Voting and Algorithmic Hardness Applications*,
arXiv:2209.11216v1](https://arxiv.org/html/2209.11216v1), Definition 1.3 and
Theorem 1.5, classifies three- and four-cell maximizers **assuming
hyperstability**. The definition additionally requires vanishing of a mixed
deformation/correlation variation at correlation zero. Maximality at one
positive correlation does not supply this hypothesis. The argument here does
not assume or deduce it. These comparisons identify applicable methods and
missing hypotheses, not a novelty or exhaustive literature certification.

For the full-noise quadratic form below, the same paper's **general**
second-variation formula, Lemma 4.5, is applicable to motions compactly
supported in facet interiors. The tetrahedral score differences vanish on
their facets by reflection, so its acceleration term vanishes. This uses
stationarity, proved directly here, without assuming that the tetrahedron
maximizes noise stability. The maximizer-specific specialization must not
be used to assume the desired conclusion.

## Quantified maximizer target and its consequence

Write $\mathcal T$ for the orbit of the central tetrahedral partition under
orthogonal maps and labels, and
$d_4(A)=\inf_{T\in\mathcal T}\sum_i\gamma_3(A_i\triangle T_i)$.
A sufficient **uniform maximizer-rigidity obligation** is:

> There exist $\eta_*>0$ and $0<\rho_*<1$ such that, for every
> $0<\rho\le\rho_*$, every balanced global maximizer $A$ of $S_\rho$ over
> measurable four-cell partitions of $\mathbb R^3$ with $d_4(A)\le\eta_*$
> belongs to $\mathcal T$, up to null sets.

The neighborhood and correlation interval must work for all these maximizers.
Global maximality is retained; classification of every stationary solution is
an optional stronger objective. Ruling out strictly improving maximizers would
already suffice for the inequality, without classification of possible ties.

The existing localization bound (N1), existence and dimension reduction imply
that this obligation gives the desired inequality in every $d\ge3$ for

\[
0<\rho\le\rho_0,\qquad
\rho_0=\min\{\rho_*,4\eta_*^2/(3K_4^2)\}>0.
\tag{R1}
\]

Indeed an improving competitor gives an improving maximizer, its
three-dimensional reduction has
$d_4(A)<K_4\sqrt{3\rho/4}\le\eta_*$, and rigidity contradicts improvement.
This is a proved implication with an open premise. The earlier
[nearby full-kernel comparison](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#localization-and-the-remaining-full-kernel-theorem)
remains a separate sufficient route. No implication from that inequality to
uniqueness of maximizers is asserted.

## Assignment and the exact information supplied by global maximality

Fix a balanced global maximizer $A$ in $\mathbb R^3$, put $f_i=1_{A_i}$,
$u_i=T_\rho f_i$, and $b_i=\int x f_i\,d\gamma_3$. For every balanced
partition $B$, set $v_i=1_{B_i}-f_i$. Self-adjointness and the semigroup
identity give

\[
\sum_i\langle f_i-1_{B_i},u_i\rangle
\ \ge\ \frac12\sum_i\langle v_i,T_\rho v_i\rangle
=\frac12\sum_i\|T_{\sqrt\rho}v_i\|_2^2\ge0.
\tag{A1}
\]

This follows by expanding $S_\rho(B)-S_\rho(A)\le0$. In particular $A$
maximizes the linear assignment objective with fixed scores $u_i$. The full
inequality retains global-maximality information; the assignment equation
alone is only a necessary condition.

For all sufficiently small $\rho$, the moment rows are pairwise distinct.
In fact $S_\rho(A)\ge S_\rho(T)$ and the absolute Hermite remainder give
$0\le\delta=F_4-W_1(A)\le3\rho/4$. Equation (M2) of the first-moment
proof permits a rotation and labels such that, with $B$ having rows $b_i$,

\[
\|B-B_*\|_F^2\le L_4\delta\le3L_4\rho/4,
\qquad B_*B_*^T=(F_4/3)P.
\tag{A2}
\]

Here $B_*$ is the actual moment matrix of the aligned tetrahedron.
For $i\ne j$, $u_i-u_j$ is real analytic and nonconstant: its first
moment is $\rho(b_i-b_j)\ne0$. Gaussian convolution of a bounded function
is real analytic, as follows by extending its Gaussian kernel to complex
$x$ and using an integrable bound on each compact complex set. Every level
set of this nonconstant analytic difference is null.

Minimize the finite convex price functional

\[
D_u(\lambda)=\tfrac14\sum_i\lambda_i+
  \int\max_i(u_i(x)-\lambda_i)\,d\gamma_3(x),
\qquad \sum_i\lambda_i=0.
\]

It is coercive: $D_u\ge1/4-\min_i\lambda_i$, which tends to infinity
with the norm in this gauge. Its minimizer exists. Ties are null by the
analytic observation, so differentiation under the integral shows that its
winning partition has all four masses $1/4$. It realizes the largest linear
assignment value. By (A1), $A$ also realizes this value; the integrated
nonnegative pointwise assignment regret is zero. Thus

\[
A_i=\{x:u_i(x)-\lambda_i>u_j(x)-\lambda_j\ \text{for all }j\ne i\}
\quad\text{a.e.}
\tag{A3}
\]

No regularity or topology of $A$ was assumed in this derivation. The
non-strict version differs only on null sets.

## Normalized scores and prices

Set

\[
g_i=(u_i-1/4)/\rho=b_i\cdot x+h_i(x),\qquad
\mu_i=\lambda_i/\rho,\qquad q_i=g_i-\mu_i.
\]

Hermite orthogonality and Parseval give the uniform bound

\[
\sum_i\|h_i\|_2^2
=\sum_{r\ge2}\rho^{2r-2}W_r(A)
\le\rho^2(3/4-W_1(A))\le3\rho^2/4.
\tag{B1}
\]

For a full-rank centered moment matrix $B$ close to $B_*$, define

\[
D_B(\nu)=\tfrac14\sum_i\nu_i+
 \mathbb E\max_i(b_i\cdot X-\nu_i),\qquad \nu\in\mathbf1^\perp.
\]

At $(B_*,0)$ its price Hessian is the Laplacian of positive Gaussian
facet weights $A_{ij}/\|b_i^*-b_j^*\|$. It is positive definite on
$\mathbf1^\perp$. To justify this formula, differentiate the winning
masses: motion of a pairwise score level contributes its Gaussian facet
integral divided by the length of its score gradient. Triple intersections
have codimension two and contribute no first-order mass. A full-rank score
map identifies these integrals with polyhedral Gaussian integrals on
$\mathbf1^\perp$; their dependence is smooth near the regular data.
The implicit function theorem therefore supplies a unique local minimizer
$\nu_B$ with $\|\nu_B\|\le C\|B-B_*\|_F$ and uniform local strong
convexity. It is also the unique global minimizer: convexity makes it global,
and any other minimizer would contradict strict convexity on a segment
starting at $\nu_B$.

Let $D_g$ denote the corresponding price functional with scores $g_i$.
By (B1), uniformly in all prices,

\[
|D_g(\nu)-D_B(\nu)|\le\mathbb E\max_i|h_i|
\le\sqrt3\rho/2.
\tag{B2}
\]

The prices $\mu$ in (A3) minimize $D_g$. All such minimizers stay in a
fixed compact set: $\mathbb E g_i=0$ gives $D_g(\nu)\ge-\min_i\nu_i$,
while $D_g(0)\le\mathbb E\|(g_i)_i\|\le\sqrt{F_4+3\rho^2/4}$.
The zero-sum gauge bounds $\|\nu\|\le\sqrt{12}(-\min_i\nu_i)$.
Equations (A2) and (B2) now imply that these minimizers approach zero:
on this fixed compact set $D_g$ converges uniformly to $D_{B_*}$, whose
unique minimizer is zero. They consequently enter the uniform strong
convexity neighborhood. Minimality and (B2) give
$D_B(\mu)-D_B(\nu_B)\le\sqrt3\rho$, so strong convexity proves

\[
\|\mu-\nu_B\|\le C\sqrt\rho,\qquad
\|\mu\|\le C\sqrt\rho,\qquad
\|\lambda\|\le C\rho^{3/2}.
\tag{B3}
\]

Constants here and below are finite and independent of the maximizer;
their optimization or an explicit usable correlation threshold is not claimed.

For every fixed integer $m\ge0$ and $R<\infty$, Taylor expansion of the
Mehler kernel at $\rho=0$ also gives

\[
\max_i\|h_i\|_{C^m(B_R)}\le C_{R,m}\rho,\qquad
\max_i\|q_i-b_i^*\cdot x\|_{C^m(B_R)}\le C_{R,m}\sqrt\rho.
\tag{B4}
\]

For completeness, relative to $\gamma_3(dy)$ the kernel is

\[
M_\rho(x,y)=(1-\rho^2)^{-3/2}
 \exp\!\left(\frac{2\rho x\cdot y-\rho^2(|x|^2+|y|^2)}
                   {2(1-\rho^2)}\right).
\]

Its value and first correlation derivative at zero are $1$ and $x\cdot y$.
For $|x|\le R$, $|\rho|\le1/2$, the second correlation derivative and
each of its first $m$ spatial derivatives are bounded by a fixed polynomial
in $|y|$ times $e^{C_R|y|}$. This is integrable against $\gamma_3$.
Taylor's integral remainder, $0\le f_i\le1$, and differentiation under
the integral prove the first estimate in (B4). Equations (A2) and (B3)
prove the second, including the constant price terms.

## What this proves about bounded interfaces

Let $K$ be a compact subset of the relative interior of a tetrahedral
facet $F_{ij}$, away from its triple rays. The other two regular scores
are separated below $i,j$ on $K$. In a fixed sufficiently small tube over
$K$, (B4) retains this separation and gives

\[
\partial_{n_{ij}}(q_i-q_j)
=\|b_i^*-b_j^*\|+O_K(\sqrt\rho)>0,
\qquad n_{ij}=(b_i^*-b_j^*)/\|b_i^*-b_j^*\|.
\]

The implicit function theorem, with opposite signs at the ends of each
normal segment, gives one interface graph there, with displacement
$O_K(\sqrt\rho)$ in every fixed $C^m$ norm. Assignment (A3) identifies
the two labels on its sides. Thus an additional interface component cannot
hide in such a fixed tube.

One can also locate the fourfold junction without assuming its existence.
The map $H=(q_1-q_4,q_2-q_4,q_3-q_4)$ has derivative uniformly close on
any fixed ball to the invertible matrix with rows $b_i^*-b_4^*$.
Integration along line segments gives injectivity on that ball for small
$\rho$. Writing $H(x)=J_*x+e(x)$, the equation
$x=-J_*^{-1}e(x)$ is a contraction on a ball of radius $C\sqrt\rho$:
(B4) bounds $e(0)$ and $De$. It has a solution there. Hence there is one
fourfold score junction in any fixed ball containing the origin in its
interior, and its location is $O(\sqrt\rho)$.

These statements concern score interfaces on fixed regions. They do not
assert a global normal-graph representation, tail topology, or uniqueness
of the partition solving (A3).

## An exact kernel estimate on expanding balls

The bounded-region result has a useful larger-domain version. In the
normalized Hermite basis, Cauchy--Schwarz for the degree-at-least-two part
of $f_i$ gives, with $r=\rho^2$ and
$V_r(x)=M_r(x,x)=(1-r^2)^{-3/2}e^{r|x|^2/(1+r)}$,

\[
|h_i(x)|\le\frac{\sqrt3}{4\rho}
 \sqrt{V_r(x)-1-r|x|^2},
\tag{E1}
\]

because $\|f_i-1/4\|_2^2=3/16$. Differentiating the Hermite kernel in
its two arguments gives likewise

\[
|\partial_j h_i(x)|\le\frac{\sqrt3}{4}
 \sqrt{V_r(x)\left(\frac1{1-r^2}
                  +\frac{r x_j^2}{(1+r)^2}\right)-1}.
\tag{E2}
\]

Indeed $\partial_{x_j}\partial_{y_j}M_r(x,y)|_{y=x}
=V_r(x)[r/(1-r^2)+r^2x_j^2/(1+r)^2]$; subtracting the degree-one
contribution $r$ and dividing by $\rho^2$ proves (E2). The diagonal
identity and its derivatives converge locally absolutely for $r<1$,
so Cauchy--Schwarz applies also to their differentiated series.

For $0<\rho\le1/2$, $R\ge1$, and $\rho R\le1/2$, these identities give

\[
\begin{split}
\max_i\|q_i-b_i^*\cdot x\|_{L^\infty(B_R)}
 &\le C[\sqrt\rho(1+R)+\rho(1+R^2)],\\
\max_i\|\nabla q_i-b_i^*\|_{L^\infty(B_R)}
 &\le C[\sqrt\rho+\rho(1+R)].
\end{split}
\tag{E3}
\]

For example, put $a=\rho^2|x|^2\le1/4$ and use
$(1-\rho^4)^{-3/2}\le1+C\rho^4$,
$e^a-1-a\le Ca^2$. The radicand in (E1) is at most
$C\rho^4(1+R^4)$, and that in (E2) is at most
$C[\rho^2R^2+\rho^4]$. Combine with (A2), (B3).

Taking $R_\rho=\rho^{-1/4}$, for sufficiently small $\rho$, makes the
score error $O(\rho^{1/4})$ and gradient error $O(\sqrt\rho)$ throughout
an expanding ball. Its omitted Gaussian mass is at most

\[
\gamma_3(|X|>R_\rho)
\le\sqrt{2/\pi}(R_\rho+R_\rho^{-1})e^{-R_\rho^2/2}.
\tag{E4}
\]

This follows by radial integration and
$\int_R^\infty e^{-t^2/2}dt\le R^{-1}e^{-R^2/2}$.
The error is smaller than every fixed power of $\rho$, but it is not
zero. Neither (E3) nor (E4) excludes a branch of maximizers whose
differences are concentrated still farther out.

## The limiting interface Hessian has only rotational zero modes

This statement is for the first-moment functional $W_1$, rather than
$S_\rho$ at positive correlation. Let $\ell=\|b_i^*-b_j^*\|
=\sqrt{2F_4/3}>0$. Orient $n_{ij}$ into $T_i$ and let $f_{ij}$ be
the outward speed of $T_i$ on $F_{ij}$, so velocity is $-f_{ij}n_{ij}$.
Use $f_{ji}=-f_{ij}$ when referring to the same speed from $T_j$.
Define the Hilbert space and its volume-preserving subspace by

\[
\mathcal H=\bigoplus_{i<j}L^2(F_{ij},\phi_3\,d\mathcal H^2),
\quad
\mathcal H_0=\{f:\ \sum_{j\ne i}\int_{F_{ij}}f_{ij}\phi_3=0
                         \text{ for each }i\}.
\]

Let $\mathcal D f$ have rows
$(\mathcal D f)_i=\sum_{j\ne i}\int_{F_{ij}}x f_{ij}\phi_3$.
It is a bounded finite-rank operator: all Gaussian facet second moments
are finite. The [second variation (9)](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md#what-fixed-infinitesimal-interface-tests-see)
is

\[
-\tfrac12 W_1''=Q(f):=\ell\|f\|_{\mathcal H}^2-\|\mathcal D f\|_F^2.
\tag{H1}
\]

For a skew matrix $\Omega$, rotation $e^{t\Omega}$ has speed
$f^\Omega_{ij}=-(\Omega x)\cdot n_{ij}$. Denote this three-dimensional
subspace of $\mathcal H_0$ by $\mathcal R$. There exists $\kappa>0$ such
that

\[
Q(f)\ge\kappa\inf_{r\in\mathcal R}\|f-r\|_{\mathcal H}^2
\qquad(f\in\mathcal H_0).
\tag{H2}
\]

Here is a proof including the zero-mode issue. Begin with smooth speeds
compactly supported away from facet edges and satisfying the four linear
volume conditions. Such speeds have genuine deformations preserving volumes:
start with their disjoint interface motions, then correct their $O(t^2)$
volume errors with three fixed interior facet bumps. Choose the bump pairs
as a spanning tree on the four labels; the volume derivative has rank three,
so the implicit function theorem gives these corrections with zero initial
speed. Thus (H1) applies. The global first-moment theorem makes $Q\ge0$.

For such a deformation $B(t)=B_*+t\mathcal D f+O(t^2)$ and
$F_4-W_1(A(t))=t^2Q(f)+o(t^2)$. Equation (M2) therefore implies

\[
\operatorname{dist}(\mathcal D f,\mathcal V)^2\le L_4Q(f),
\quad
\mathcal V=\{\text{matrices with rows }\Omega b_i^*:
                                      \Omega^T=-\Omega\}.
\tag{H3}
\]

To see the limiting distance, the aligned tetrahedral moment matrices form
the smooth orthogonal orbit through $B_*$. Taking the first-order distance
of $B_*+t\mathcal D f$ to that orbit gives the distance to its tangent
space $\mathcal V$. Label permutations leave this orbit unchanged.

The chosen smooth speeds are dense in $\mathcal H_0$: truncate and smooth
on facet interiors, then use the same three bumps to correct the converging
linear volume errors. Continuity extends $Q\ge0$ and (H3) to
$\mathcal H_0$.

Rotations preserve both volumes and $W_1$, so $Q(f^\Omega)=0$ and
$\mathcal D f^\Omega$ has rows $\Omega b_i^*$. This can also be checked
by the Gaussian facet integrals in (H1); their polynomial integrands permit
truncation of the unbounded rotational deformation with vanishing errors.
If $Q(f)=0$, (H3) supplies $r\in\mathcal R$ with
$\mathcal D r=\mathcal D f$. A nonnegative quadratic form has zero
polarization against every vector in its kernel. Consequently
$Q(f-r)=0$, while $\mathcal D(f-r)=0$ forces
$Q(f-r)=\ell\|f-r\|^2$. Hence $f=r$. The centroid vectors span
$\mathbb R^3$, so the rotational map is injective and the kernel has
exactly three dimensions.

Finally, on $\mathcal H_0$, the operator of $Q$ is $\ell I$ minus a
finite-rank nonnegative operator. Its restriction outside the finite
dimensional range of the latter is $\ell I$. On that range, after
removing its kernel, its least eigenvalue is positive. Taking the smaller
of this value and $\ell$ proves (H2). This is an existential analytic gap,
not a numerical eigenvalue certificate or a positive-correlation Hessian
estimate.

## Sharper localization using the limiting gap

The preceding estimates improve, for actual maximizers, to

\[
\begin{split}
\|B-B_*\|_F+\|\mu\|&\le C\rho,\qquad \|\lambda\|\le C\rho^2,\\
d_4(A)&\le C\rho,\qquad 0\le F_4-W_1(A)\le C\rho^2,\\
\max_i\|q_i-b_i^*\cdot x\|_{C^m(B_R)}&\le C_{R,m}\rho
\quad\text{for fixed }R,m.
\end{split}
\tag{S1}
\]

The constants and small-correlation threshold are uniform over all such
maximizers. The bounds follow in three steps; no tail is dropped from the
assignment equation.

First, compare winning regions for $g_i=b_i\cdot x+h_i$ and $b_i\cdot x$
at the **same** price $\nu$ in a fixed neighborhood of zero. On
$|x|\le R_\rho=\rho^{-1/4}$, (E1) gives
$|h_i(x)|\le C\rho(1+|x|^2)$. Whenever the winning labels differ, some
pair satisfies

\[
|(b_i-b_j)\cdot x-(\nu_i-\nu_j)|
\le C\rho(1+|x|^2).
\tag{S2}
\]

Write $x=z+tn$, where $n$ is the unit pairwise score normal and $z\perp n$.
The score length is bounded below uniformly by (A2), and
$a=(\nu_i-\nu_j)/\|b_i-b_j\|$ is uniformly bounded. Thus (S2) implies
$|t-a|\le C\rho(1+|z|^2+t^2)$. On this ball
$|t-a|\le R_\rho+C$; use $t^2\le2(t-a)^2+2a^2$ and absorb
$2C\rho(R_\rho+C)|t-a|$ for small $\rho$. The resulting strip is

\[
|t-a|\le C\rho(1+|z|^2).
\]

Conditioning on the two Gaussian tangential coordinates bounds both its
Gaussian mass and its integral of $|x|$ by $C\rho$. Indeed its width is
$C\rho(1+|z|^2)$, the one-dimensional Gaussian density is bounded, and
$|t|\le |a|+C\rho(1+|z|^2)$ inside it; all resulting tangential moments
are finite. Sum over six pairs. The omitted mass and first moment outside
$R_\rho$ are smaller than every power of $\rho$ by the same radial
integration as (E4). This proves an $O(\rho)$ difference of winning
masses and winning moments, uniformly over these prices.

In particular $\|\nabla D_g(\nu)-\nabla D_B(\nu)\|\le C\rho$.
The preliminary price bound (B3) puts $\mu$ and $\nu_B$ in this
neighborhood. Since $\nabla D_g(\mu)=\nabla D_B(\nu_B)=0$, uniform
strong convexity gives $\|\mu-\nu_B\|\le C\rho$. Changing the price
of the affine scores from $\mu$ to $\nu_B$ changes their moments by
$O(\rho)$, again by Gaussian strips. If $\mathcal C(B)$ denotes the
moment matrix of the exactly balanced affine winning partition for $B$,
assignment (A3) consequently gives

\[
\|B-\mathcal C(B)\|_F\le C\rho.
\tag{S3}
\]

Second, $B\mapsto\mathcal C(B)$ is smooth near $B_*$, by the full-rank
Gaussian facet formulas and price implicit function theorem used above.
More explicitly, the score map $x\mapsto Bx$ is an isomorphism onto
$\mathbf1^\perp$. In score coordinates, translate by the balancing price
so each winning region is a fixed central cone. Its Gaussian density, the
linear inverse score map, and the price depend smoothly on $B$; derivatives
are dominated by polynomial Gaussian bounds near $B_*$. This justifies
differentiating the moment integrals despite the moving interfaces.
Its linearization is nondegenerate modulo rotations. Here is the precise
connection with (H2). For a centered matrix perturbation $Z$ with rows $z_i$,
let $f_Z$ be the outward facet speed of this affine winning partition.
Differentiating its score equalities gives

\[
\ell(f_Z)_{ij}=(z_i-z_j)\cdot x-(\nu'_i-\nu'_j),
\qquad f_Z\in\mathcal H_0.
\]

For any $w\in\mathcal H_0$, the price terms cancel on summing facets, so
$\ell\langle f_Z,w\rangle=\langle Z,\mathcal D w\rangle_F$.
With $\Pi_0$ denoting orthogonal projection onto $\mathcal H_0$,

\[
\mathcal C'(B_*)Z=\mathcal D f_Z
=\ell^{-1}\mathcal D\Pi_0\mathcal D^*Z.
\tag{S4}
\]

This is self-adjoint on the finite-dimensional space of centered moment
matrices. If $(I-\mathcal C'(B_*))Z=0$, then
$\mathcal D f_Z=Z$ and
$\ell\|f_Z\|^2=\langle Z,\mathcal D f_Z\rangle_F=\|Z\|_F^2$.
Thus $Q(f_Z)=0$, and (H2) makes $f_Z$ rotational and $Z\in\mathcal V$.
Conversely equivariance of $\mathcal C$ and
$\mathcal C(B_*)=B_*$ show that every $Z\in\mathcal V$ lies in this
kernel. Therefore $I-\mathcal C'(B_*)$ is invertible on
$\mathcal V^\perp$.

Choose the closest tetrahedral moment matrix to $B$ and rotate it to $B_*$.
The error $E=B-B_*$ is then orthogonal to $\mathcal V$, by differentiation
of squared distance along the orthogonal orbit. It tends to zero by (A2).
Taylor expansion of $B-\mathcal C(B)$ and the invertibility just proved
give
$\|B-\mathcal C(B)\|_F\ge c\|E\|_F-O(\|E\|_F^2)$.
Equation (S3) now proves $\|E\|_F\le C\rho$. Smooth price dependence
gives $\nu_B=O(\rho)$, hence $\mu=O(\rho)$ and $\lambda=O(\rho^2)$.

Third, on $B_{R_\rho}$ these refined moment and price estimates and (E1)
give $|q_i-b_i^*\cdot x|\le C\rho(1+|x|^2)$. If $A$ and the regular
partition $T$ assign different labels, their nonnegative regular score
regret is at most $C\rho(1+|x|^2)$. Apply the strip argument with $a=0$.
It bounds the changed-label mass by $C\rho$. It also bounds the total
regular-score regret $r$ by $C\rho^2$: conditionally on $z$, integrate
$\ell|t|\phi(t)$ on $|t|\le C\rho(1+|z|^2)$ and use the finite
fourth tangential moment. Outside the expanding ball the regret is at
most a constant times $\int |x|\,d\gamma_3$, restricted to that tail,
which is $o(\rho^2)$. The exact identity (M3) gives

\[
0\le\delta=2r-\|B-B_*\|_F^2\le2r\le C\rho^2.
\tag{S5}
\]

This proves the distance and deficit assertions in (S1); (B4)'s kernel
remainder proves its fixed-ball $C^m$ assertion with the improved moment
and price rates. The fixed-region facet and junction displacements therefore
also improve to $O(\rho)$. This remains a rate theorem, not an exclusion
of nonzero maximizing branches.

## Global score margins and facet normal derivatives

Write $U_{ij}(x)=T_\rho(1_{T_i}-1_{T_j})(x)$,
$\sigma=\sqrt{1-\rho^2}$, and retain the inward pair normal $n_{ij}$.
For every $0<\rho<1$ and every $x\in T_i$, there is a fixed $c_T>0$
such that

\[
c_T\min\{\rho\,n_{ij}\cdot x/\sigma,1\}
\le U_{ij}(x)
\le\min\{\sqrt{2/\pi}\,\rho\,n_{ij}\cdot x/\sigma,1\}.
\tag{T1}
\]

In particular this is a margin for the **full** kernel, with saturation
instead of unbounded linear growth in remote interiors. On every original
facet $F_{ij}$, including its remote portions,

\[
\frac{\rho\ell}{\sigma}
\le a_{\rho,ij}(x):=\partial_{n_{ij}}U_{ij}(x)
\le\frac{\rho}{\sigma}\sqrt{2/\pi}.
\tag{T2}
\]

Here is a direct kernel proof. Let $G$ reflect across $n_{ij}^\perp$;
then $GT_i=T_j$. With $m=\rho x/\sigma$ and $y_N=n_{ij}\cdot y$,
cone invariance under positive scaling gives

\[
U_{ij}(x)=\int_{T_i}[\phi_3(y-m)-\phi_3(Gy-m)]\,dy,
\qquad
\frac{\phi_3(Gy-m)}{\phi_3(y-m)}=e^{-2m_Ny_N}.
\tag{T3}
\]

For $m\in T_i$, $m_N\ge0$, so the integrand is nonnegative. Enlarging
$T_i$ to the whole positive halfspace bounds the integral above by
$2\Phi(m_N)-1$, proving the upper bound in (T1).
Let $v_i=b_i^*/|b_i^*|$. Each inward facet normal of $T_i$ has inner
product $\sqrt{2/3}$ with $v_i$. Thus the ball
$D=B(v_i,1/4)$ lies in $T_i$, has $z_N>1/2$, and satisfies

\[
\gamma_3(D)\ge p_0:=(\pi/48)(2\pi)^{-3/2}e^{-25/32}>0.
\]

Convex cone addition gives $m+D\subset T_i$. Restrict (T3) to
$y=m+z$ for $z\in D$. It is at least
$p_0(1-e^{-m_N})$, which is at least
$p_0(1-e^{-1})\min\{m_N,1\}$. This proves (T1) with
$c_T=p_0(1-e^{-1})$.

For (T2), take $x\in F_{ij}$, so $m_N=0$ and $m\in T_i$.
Differentiation of (T3) in the normal direction, followed by translation,
gives

\[
\partial_{n_{ij}}U_{ij}(x)
=\frac{2\rho}{\sigma}\int_{T_i-m} z_N\phi_3(z)\,dz
\ge\frac{2\rho}{\sigma}\int_{T_i}z_N\phi_3(z)\,dz
=\frac{\rho\ell}{\sigma}.
\]

Indeed $T_i\subset T_i-m$, because $m+T_i\subset T_i$, and the extra
integrand is nonnegative: $m_N=0$ and $T_i-m\subset\{z_N\ge0\}$.
The last equality is $2b_i^*\cdot n_{ij}=\ell$.
Differentiating the original convolution bounds the derivative above by
$(\rho/\sigma)\mathbb E|Z_N|$, proving its upper bound.
All differentiations are justified by Gaussian domination. Reflection
also gives $U_{ij}=0$ on the entire pair plane, so the gradient on that
plane is normal to it.

For any measurable partition $A$, (T1) yields the exact comparison

\[
\begin{split}
S_\rho(A)-S_\rho(T)
 &=\sum_i\langle e_i,T_\rho e_i\rangle-2\mathcal L_\rho(A),
 \qquad e_i=1_{A_i}-1_{T_i},\\
\mathcal L_\rho(A)
 &=\sum_{i\ne j}\int_{T_i\cap A_j}U_{ij}\,d\gamma_3
 \ge c_T\sum_{i\ne j}\int_{T_i\cap A_j}
       \min\{\rho\,n_{ij}\cdot x/\sigma,1\}\,d\gamma_3.
\end{split}
\tag{T4}
\]

This is a saturated regret bound, not a sign for the complete stability
difference. The nonnegative quadratic interaction in (T4) must still be
controlled. Saturation is essential: a remote swap can have first-moment
loss of order $R$ times its mass while its full-kernel score loss stays
bounded by its mass.

## A uniform full-noise interface gap, including linear tail motions

Use the Hilbert spaces $\mathcal H_0$ and $\mathcal R$ from (H1)--(H2).
Write $d\nu=\phi_3\,d\mathcal H^2$ on each facet. Define

\[
\begin{split}
\mathcal B_\rho(f)
 &=\sum_i\sum_{j,k\ne i}\int_{F_{ij}}\int_{F_{ik}}
 M_\rho(x,y)f_{ij}(x)f_{ik}(y)\,d\nu(y)d\nu(x),\\
\mathcal Q_\rho(f)
 &=\rho^{-1}\left[
   \sum_{i<j}\int_{F_{ij}}a_{\rho,ij}f_{ij}^2\,d\nu
   -\mathcal B_\rho(f)\right].
\end{split}
\tag{K1}
\]

For a smooth balanced deformation compactly supported away from the facet
edges, the general second-variation formula reduces to
$-S_\rho''/(2\rho)=\mathcal Q_\rho(f)$.
In fact differentiation of the moments of the indicator functions gives
the double boundary integral $\mathcal B_\rho$; the local boundary term
is $-a_{\rho,ij}f_{ij}^2$. Acceleration and density-variation terms are
multiplied by the tetrahedral score difference on the interface, which
is zero. No optimality assumption on $T$ is used.

There exist $\rho_h>0$ and $c_h>0$, uniform over all facet speeds, such
that

\[
\mathcal Q_\rho(f)\ge c_h
 \inf_{r\in\mathcal R}\|f-r\|_{\mathcal H}^2
\quad(f\in\mathcal H_0,\ 0<\rho\le\rho_h).
\tag{K2}
\]

In particular the estimate includes speeds supported arbitrarily far out
on the original facets. Its kernel is exactly $\mathcal R$.

To prove this, first note that, for every ordered pair of facets,

\[
\|M_\rho-1-\rho x\cdot y\|_{L^2(\nu\otimes\nu)}
\le C\rho^2\qquad (0<\rho\le1/8).
\tag{K3}
\]

The second correlation derivative of $M_\rho$ is a polynomial of degree
at most four in $x,y$, with uniformly bounded coefficients, times
$M_\rho$. For $|\rho|\le1/8$,
$M_\rho\le C\exp[(|x|^2+|y|^2)/12]$; this follows from
$|x\cdot y|\le(|x|^2+|y|^2)/2$ in the exact Mehler formula.
The square of this bound times the two Gaussian surface densities is
integrable on the two supporting planes, including every polynomial
factor just mentioned. Taylor's integral remainder in this $L^2$ space
proves (K3). There are only finitely many facet pairs, so one constant
works for all of them.

The constant term contributes
$\sum_i(\sum_{j\ne i}\int f_{ij}\,d\nu)^2=0$ on $\mathcal H_0$.
The term $\rho x\cdot y$ contributes $\rho\|\mathcal Df\|_F^2$.
The Hilbert--Schmidt bound (K3), Cauchy--Schwarz, and the three facets
incident to each cell consequently give

\[
|\mathcal B_\rho(f)-\rho\|\mathcal Df\|_F^2|
\le C_1\rho^2\|f\|_{\mathcal H}^2.
\tag{K4}
\]

Together with the **global** lower bound (T2), this proves
$\mathcal Q_\rho(f)\ge Q(f)-C_1\rho\|f\|^2$.
This already gives (K2) on $\mathcal R^\perp$ for
$\rho\le\min\{1/8,\kappa/(2C_1)\}$, with $c_h=\kappa/2$.

It remains to justify removing the rotational component. The tetrahedron
is a critical partition for $S_\rho$, since all its interface score
differences vanish. Rotation invariance of $S_\rho$, differentiated in
one rotational and one arbitrary deformation direction at this critical
partition, gives zero mixed Hessian. Thus
$\mathcal Q_\rho(r,f)=0$ for $r\in\mathcal R$ and every smooth
volume-preserving speed $f$. One can differentiate the actual rotation
and the compactly supported deformation together; the first variation
is zero, so no mixed acceleration term remains. Gaussian decay justifies
the unbounded linear rotational speed. Alternatively truncate it and
pass to the same boundary integrals.
By (T2)'s upper bound and (K4), the full form is bounded on
$\mathcal H_0$ for each such $\rho$, uniformly after normalization.
The density argument used for (H2) therefore extends the mixed identity
to all $f\in\mathcal H_0$. Subtract the orthogonal projection onto
$\mathcal R$ to obtain (K2) in its stated form.

The tail issue is consequently resolved for this **linearized form at
the tetrahedron**. This does not put arbitrary nearby measurable
partitions or actual maximizers into its tangent space with a controlled
nonlinear error. In particular it does not assert the retired relative
Hermite estimate, or a value deficit proportional to the unrestricted
first-moment deficit.

## A common interface chart on the scale of inverse correlation

There are fixed $c>0$ and $\rho_c>0$ such that every maximizer under
consideration and $0<\rho\le\rho_c$ has an analytic interface chart
on a domain containing $B_{c/(4\rho)}$. This includes all junction strata,
rather than only compact portions of individual facet interiors.

First, the normalized scores have the following global estimates:

\[
\begin{split}
|q_i(x)-b_i^*\cdot x|&\le C\rho(1+|x|^2),\\
|\nabla q_i(x)-b_i^*|&\le C\rho(1+|x|),\\
\|D^m q_i\|_\infty&\le C_m\rho^{m-1}\quad(m\ge2).
\end{split}
\tag{P1}
\]

To see the derivative bounds directly, differentiate the Gaussian
convolution using integration by parts in its Gaussian noise variable.
For a multi-index $\alpha$ of length $m\ge1$ and the probabilists'
Hermite polynomial $H_\alpha$,

\[
\partial^\alpha g_i(x)
=\rho^{m-1}\sigma^{-m}
 \mathbb E[H_\alpha(Z)1_{A_i}(\rho x+\sigma Z)].
\tag{P2}
\]

This identity and $0\le1_{A_i}\le1$ give the asserted uniform bounds
for $m\ge2$. At $x=0$, the Mehler formula also gives
$|g_i(0)|\le C\rho$ and $|\nabla g_i(0)-b_i|\le C\rho^2$.
Indeed its factors at zero are
$\sigma^{-3}e^{-\rho^2|y|^2/(2\sigma^2)}$ and, for the normalized
gradient, $y\sigma^{-5}e^{-\rho^2|y|^2/(2\sigma^2)}$.
Their differences from $1$ and $y$ are bounded in Gaussian $L^1$ by
$C\rho^2$, uniformly against $0\le1_{A_i}\le1$.
Taylor's theorem along line segments, the uniform Hessian bound
$\|D^2g_i\|_\infty\le C\rho$, and (S1) now prove (P1).
These polynomially weighted estimates remain valid for all $x$; they
need not be small compared with linear scores when $|x|\gg1/\rho$.

Let $J_*$ have rows $b_i^*-b_4^*$ for $i=1,2,3$ and set

\[
\Phi_A(x)=J_*^{-1}(q_1-q_4,q_2-q_4,q_3-q_4)^T(x).
\tag{P3}
\]

The matrix $J_*$ is invertible. By (P1), with $R=c/\rho$,
$\|D\Phi_A-I\|_{L^\infty(B_R)}\le C'(\rho+c)$ and
$|\Phi_A(0)|\le C'\rho$. Choose $c$ and then $\rho_c$ so that the
derivative bound is at most $1/4$. Integration along segments gives

\[
\tfrac34|x-y|\le|\Phi_A(x)-\Phi_A(y)|
\le\tfrac54|x-y|\qquad(x,y\in B_R).
\tag{P4}
\]

Moreover $B_{R/2}\subset\Phi_A(B_R)$ for small $\rho$.
For $|z|\le R/2$, solve
$x=z-[\Phi_A(x)-x]$ on the closed ball of radius $R$.
The right side is a contraction of factor $1/4$ and stays in the ball,
since its norm is at most $R/2+C'\rho+R/4<R$.
It has a unique fixed point, which lies in the open ball.
Consequently

\[
\Omega_A:=\Phi_A^{-1}(B_{R/2})\cap B_R,
\qquad B_{R/4}\subset\Omega_A\subset B_R,
\tag{P5}
\]

and $\Phi_A:\Omega_A\to B_{R/2}$ is an analytic diffeomorphism.
Its inverse has Lipschitz constant at most $4/3$.
The first inclusion in (P5) follows from
$|\Phi_A(x)|\le C'\rho+5|x|/4$.

By construction,
$q_i(x)-q_j(x)=(b_i^*-b_j^*)\cdot\Phi_A(x)$ for every pair.
The assignment equation therefore gives the **single common chart**

\[
A_i\cap\Omega_A
=\Phi_A^{-1}(T_i\cap B_{R/2})\quad\text{a.e.}
\tag{P6}
\]

With the analytic score-winning representatives of the cells, this
identifies every facet, triple ray and the fourfold junction in this
domain. Their pieces have exactly the tetrahedral incidences and
connectedness; additional interior interface branches are excluded
there. The exterior mass satisfies

\[
\gamma_3(\Omega_A^c)
\le\sqrt{2/\pi}\left(\frac{c}{4\rho}+\frac{4\rho}{c}\right)
       \exp[-c^2/(32\rho^2)].
\tag{P7}
\]

The chart is a theorem on an expanding domain. The positive exterior
mass in (P7), however small, is not discarded from any exact comparison.
No global continuation of this chart, or uniqueness of a maximizing
solution, follows from (P7) alone.

## The unresolved uniform step

The proved facts give actual maximizers an assignment equation, bounded
prices, smooth nearby facet interiors on fixed regions, expanding-ball
score control, the sharper maximizer bounds (S1), the uniform full-noise
interface gap (K2), and the common inverse-correlation-scale chart (P6).
They leave the following work in the maximizer-rigidity obligation:

1. Extend (P6) to an admissible global interface class for every localized
   maximizer, or control the exterior directly from the interior deformation.
   All junctions inside $\Omega_A$ are now covered; remote components beyond
   this domain are not excluded.
2. Control the exact kernel on the tails of that description, uniformly as
   $\rho\downarrow0$, with an error measured relative to the deformation
   being excluded. An absolute exponentially small error can still permit
   an exponentially small nonzero solution. Exact reflected-swap losses
   concern restricted competitors and cannot simply be iterated on changed
   cells or used for asymmetric cycles.
3. Use the full-noise gap (K2) in a uniform nonlinear uniqueness argument
   in that class, with balance enforced and rotations removed. Its linear
   transfer from (H2), including linear tail motions, is now proved.
   A uniform nonlinear remainder and admission of every maximizer remain
   unproved; neither is implied by the Hessian calculation.

For example, a valid contraction on a complete space containing every such
maximizer, for the balanced assignment map of (A3) after fixing rotations,
would finish the argument if it fixes the tetrahedron and has contraction
factor strictly below one uniformly for small $\rho$. Neither existence of
that space and chart nor the contraction estimate is proved here. A
second-variation argument may instead use global maximality through (A1)
and exclude every nonzero admissible deformation.

The tail-swap obstruction to a relative Hermite estimate for arbitrary
competitors remains valid. Its members are not offered as maximizing
solutions of (A3), and it does not refute the maximizer-rigidity route.
Conversely, discarding those competitors does not supply the missing
uniform nonlinear tail estimate. These arguments do not establish a
four-cell small-positive-correlation interval.
