# Balanced affine optimization of the first Gaussian Hermite mass

The four-cell first Hermite problem reduces exactly to a three-dimensional
balanced affine assignment. The tetrahedron is a strict local maximizer in
centroid Gram space up to its orthogonal symmetries, with no equal-norm
assumption. A deterministic search from 76 asymmetric initializations
converged numerically to the tetrahedral Gram matrix. The global balanced
first-order comparison remains open in this work.

## Affine reduction

Let $A=(A_i)_{i=1}^4$ be a measurable partition for standard Gaussian measure
$\gamma_n$, with $m_i(A)=\int_{A_i}x\,d\gamma_n$. Write
$W_1(A)=\sum_i\|m_i(A)\|^2$ and call $A$ balanced when each cell has measure
$1/4$. The [tetrahedral value and Hermite tail bounds](BALANCED_FOUR_CLASS_POLYNOMIAL_PARTITIONS.md#hermite-bounds)
give $W_1(T)=0.353204552849\ldots$. A fixed balanced partition with
$\Delta=W_1(A)-W_1(T)>0$ would beat $T$ for every
$0<\rho<\min\{1,\Delta/(3/4-W_1(T))\}$ by the tail bound after degree one.

Heilman, Jagannath and Naor prove the sharp unbalanced upper bound
$9/(8\pi)=0.358098\ldots$ in dimension three. The three-sector equality
example leaves the fourth cell empty. Heilman explains why the spherical
rearrangement used for three candidates runs into this example at four
cells. The biases below enforce the measure constraint.
[Heilman, Jagannath and Naor, Theorems 1.1–1.2](https://www.stevenheilman.org/~heilman/papers/PropR3.pdf)
[Heilman, Section 1.5](https://arxiv.org/pdf/2306.03312v3)

For a tuple $z=(z_i)_{i=1}^4$ define
$F(z)=\sup_A\sum_i z_i\cdot m_i(A)$ over balanced partitions. The elementary
semidiscrete transport dual is

$$
F(z)=\min_{b\in\mathbb R^4}G(z,b),\qquad
G(z,b)=\mathbb E\max_i(z_i\cdot X+b_i)-\frac14\sum_i b_i.
$$

For any balanced $A$, the pointwise maximum gives the upper bound
$\sum_i z_i\cdot m_i(A)\le G(z,b)$. To see equality for some $A$ and $b$,
fix the gauge $\sum_i b_i=0$. There $G(z,b)\ge\max_i b_i\ge\|b\|_\infty/3$,
so a minimum exists by continuity and coercivity. A group of indices with
identical $z_i$ has equal biases at a minimum. Raising any smaller bias to
the largest in the group leaves the pointwise maximum unchanged and
strictly reduces $G$ through the subtracted mean. One may do this before
returning to the mean-zero gauge. With one common bias per group, ties
between groups have Gaussian measure zero, since their slopes differ.
Differentiating with respect to the common bias gives the winning mass
minus the group's multiplicity divided by four. This derivative vanishes
at the minimum. A winning group of multiplicity $k$ can be subdivided
measurably into $k$ sets of measure $1/4$ by cutting it with thresholds of
one coordinate. All labels in the group still attain the maximum. This
gives an exact balanced maximizing assignment, including degenerate ties.

Completing the square and taking the two suprema in either order gives

$$
\sup_{A\,\mathrm{balanced}} W_1(A)=
\sup_{z_1,\dots,z_4}\left(2F(z)-\sum_i\|z_i\|^2\right).
$$

Translation of all four sites by the same vector leaves $F$ unchanged
because $\sum_i m_i=0$. Centering the sites decreases their squared norm.
The centered sites span at most three dimensions, and their linear scores
depend only on the standard Gaussian projection to that span. Thus the
supremum over all dimensions is attained already in dimension three.
More explicitly, $F$ is positively homogeneous and continuous, since
all centroid tuples are bounded, so the last supremum is the square of the
maximum of $F$ on the compact centered Frobenius unit sphere in dimension
three. The affine assignment attains that value when scaled to the
optimizing site norm. A direct way to preserve a putative strict improvement
is to take $z_i=m_i(A)$ and replace $A$ by the balanced maximizing assignment
$B$. Then

$$
W_1(B)-W_1(A)\ge\sum_i\|m_i(B)-m_i(A)\|^2.
$$

## Full-rank stationarity

Let $F_{ij}$ be the interface of affine cells $i,j$, with Gaussian area
$a_{ij}$ obtained by integrating the ambient Gaussian density against
surface measure. For affinely independent sites, each cell is an invertible
affine image of the positive orthant, and each of the six interfaces has
positive area. Gaussian integration by parts gives

$$
m_i=\sum_{j\ne i}a_{ij}\frac{z_i-z_j}{\|z_i-z_j\|}.
$$

Indeed $x\phi_n(x)=-\nabla\phi_n(x)$, and the outward normal on interface
$ij$ points opposite to $z_i-z_j$. One can apply the divergence theorem on
bounded truncations and let their radii tend to infinity using the Gaussian
tail. At a full-rank optimizer put $z_i=m_i$. This choice is justified by
the same completion of squares, which makes any optimizing partition
maximize the balanced linear score associated to its own centroids. Let
$L$ be the weighted graph Laplacian with edge weights
$a_{ij}/\|m_i-m_j\|$. If $M$ has rows $m_i$ then $LM=M$, and the columns
span $\mathbf1^\perp$ because the rank is three. Thus
$L=I-\mathbf1\mathbf1^\mathsf{T}/4$ and every edge ratio equals $1/4$.

## Local tetrahedral maximality

For any measurable four-cell partition, allowing arbitrary cell measures,

$$
W_1(A)\le \mathcal U(m(A)),\qquad
\mathcal U(z)=2\mathbb E\max_i z_i\cdot X-\sum_i\|z_i\|^2.
$$

This follows by taking $z_i=m_i(A)$ in the pointwise maximum inequality.
At the tetrahedral centroid tuple $z^*$ equality holds. We show that
$\mathcal U$ strictly decreases under small centered shape changes away from
the orthogonal orbit of $z^*$.

Write $z_i^*=c s_i$ for $c>0$ and the four even sign vectors $s_i$ in
$\mathbb R^3$. Put $d=2\sqrt{2}c$ for their common pairwise distance.
Write $H(z)=\mathbb E\max_i z_i\cdot X$. This function is smooth near $z^*$.
For each cell the three independent score differences are an invertible
linear coordinate chart in which the cell is a fixed orthant. Changing
variables in its Gaussian integrals gives derivatives dominated locally
by polynomials times a decaying Gaussian. The first derivative of the
expectation with respect to the site tuple is the tuple of cell centroids,
since ties have zero measure. For a perturbation $V=(V_i)$, the second
directional derivative is

$$
D^2 H(z^*)[V,V]=
\sum_{i<j}\frac{1}{d}\int_{F_{ij}}
      ((V_i-V_j)\cdot x)^2\phi_3(x)\,d\mathcal H^2(x).
$$

For the formula, differentiate the first derivative, expressed as the
sum of integrals of $V_i\cdot x$ over the moving cells. The cell's
linear coordinate chart gives a linear deformation of the fixed cell.
Differentiating the Jacobian and the integrand gives the divergence of
the transport vector field times $(V_i\cdot x)\phi_3(x)$. The outward
velocity on face $ij$ for cell $i$ is $((V_i-V_j)\cdot x)/d$, obtained by
differentiating $(z_i-z_j)\cdot x=0$. Adding the terms from the two cells
gives the square. Gaussian decay removes the boundary term at infinity.
Intersections of faces have zero surface measure.

Every centered nearby tuple is $z_i=P z_i^*$ for a linear map $P$ close
to $I$. By orthogonal invariance and polar decomposition it suffices to
use $P=I+S$ with $S$ symmetric. At $z^*$ the first derivative of
$\mathcal U$ vanishes since the score-cell centroids are $z^*$. Set
$V_i=S z_i^*$. The norm of the perturbation is
$\sum_i\|V_i\|^2=(d^2/2)\|S\|_F^2$.

For a pair of sites with labels $p,q$, the sites agree on one coordinate,
say $a$, and disagree on $j,k$. Use the unit normal $n=(z_p^*-z_q^*)/d$,
the shared-sign axis $u=s_p^{(a)} e_a$, and a unit vector $v$ perpendicular
to $n$ and $u$. The facet in $(u,v)$ coordinates is the wedge of
half-angle $\theta=\arctan\sqrt{2}$, equivalently $U\ge |V|/\sqrt{2}$.
For instance, for $(1,1,1)$ and $(1,-1,-1)$ the equality plane is
$x_2+x_3=0$ and the winning condition there is $x_1\ge|x_2|$.
Let $\mu_u,\mu_v$ be the conditional second moments of the two independent
standard Gaussian coordinates on this wedge. The angle is then uniform
on $[-\theta,\theta]$, independent of the radius whose second moment is
two. Hence

$$
\mu_u=1+\frac{\sin 2\theta}{2\theta}<2,\qquad
\mu_v=1-\frac{\sin 2\theta}{2\theta}<1,
\qquad \mathbb E[UV\mid\text{wedge}]=0.
$$

At the tetrahedron $a_{ij}/d=1/4$ by the full-rank centroid identity
above. The Hessian of the expectation is therefore

$$
\frac{d^2}{4}\sum_{\text{faces}}
\bigl(\mu_u ((Sn)\cdot u)^2+\mu_v ((Sn)\cdot v)^2\bigr).
$$

The two faces with shared coordinate $a$ have normals
$(e_j+e_k)/\sqrt{2}$ and $(e_j-e_k)/\sqrt{2}$ up to orientation.
Their contributions to the sum of $((Sn)\cdot u)^2$ add to $S_{aj}^2+S_{ak}^2$.
Each has $((Sn)\cdot v)^2=(S_{jj}-S_{kk})^2/4$, where symmetry of $S$
cancels the off-diagonal entries. Summing over the three choices of $a$
gives $\|S_{\mathrm{off}}\|_F^2$ and
$\frac32\|\operatorname{diag}(S)-\frac13\operatorname{tr}(S)I\|_F^2$,
respectively, with $S_{\mathrm{off}}$ the off-diagonal part. The Hessian
of the expectation divided by $\sum_i\|V_i\|^2$ for $S\ne0$ is consequently
at most $\max\{\mu_u/2,3\mu_v/4\}<1$. Subtracting the Hessian of the squared
norm in $\mathcal U$ gives strict negativity in the symmetric slice.
Taylor's theorem proves the local claim.

This implies $W_1(A)<W_1(T)$ for every partition with centroid Gram matrix
sufficiently near the tetrahedral one and different from it. This statement
holds in arbitrary dimension since the Gram matrix determines all site
scores in law. To see the local implication explicitly, represent its
centroids as $Pz_i^*$ with $P\colon\mathbb R^3\to\mathbb R^n$. Closeness of
Gram matrices implies $P^\mathsf{T}P$ is close to $I$, so the symmetric
factor $(P^\mathsf{T}P)^{1/2}$ is close to $I$. The isometric factor preserves
the score law. If the Gram matrices agree, equality still requires that
the tetrahedral linear score maximize pointwise almost everywhere. Thus
the cells must be tetrahedral up to an orthogonal map and null sets.
The local argument compares the first Hermite coefficient. It gives no
uniform positive correlation interval for all nearby partitions.

## Deterministic search

For each fixed affinely independent site tuple the computation solves for
balanced biases using the graph Laplacian of interface-area ratios as the
mass derivative. The derivative follows by moving each interface a signed
normal distance $(db_i-db_j)/\|z_i-z_j\|$, or by applying the same divergence
calculation to the mass integral. This is the Euclidean facet formula in
the semidiscrete transport Newton method. Its finite-domain convergence
theorems are unnecessary for this calculation on the full Gaussian space.
[Kitagawa, Mérigot and Thibert, equation (1.8) and Theorem 1.5](https://arxiv.org/pdf/1603.05579)

For a face at signed distance $h$ from the origin, the standard Gaussian
density on the plane has normal factor $\phi(h)$ and a two-dimensional
independent Gaussian in orthonormal tangent coordinates. The two remaining
score inequalities give a bivariate Gaussian probability. For upper
thresholds $x,y$ and correlation $r$, the computation uses

$$
\Phi_2(x,y;r)=\Phi(x)\Phi(y)+\frac{1}{2\pi}\int_0^{\arcsin r}
\exp\!\left(-\frac{x^2-2xy\sin\alpha+y^2}{2\cos^2\alpha}\right)\,d\alpha.
$$

One obtains this formula by differentiating the bivariate Gaussian density
with respect to its correlation. That derivative equals its mixed coordinate
derivative, so the derivative of the probability is the density evaluated
at the two upper thresholds. Integrating from zero correlation gives the
display after the sine substitution. [Plackett, 1954](https://doi.org/10.1093/biomet/41.3-4.351)

At zero biases the three-score orthant probabilities equal $1/8$ plus
$1/(4\pi)$ times the sum of the three pairwise arcsines. This follows by
expanding the product of the three positive-sign indicators and evaluating
pairwise sign correlations by the planar Gaussian angle. At general biases
the mass is the zero-bias mass plus the integral of the mass derivative
along the segment to the requested bias vector. The search uses
Gauss–Legendre rules for both integrals. It then updates the sites to the
new centroids, normalized by their Frobenius norm. The normalization does
not affect the optimal assignment because it is a common positive scale.
In exact arithmetic the replacement inequality in the affine reduction
makes this a nondecreasing iteration for the first Hermite mass of the
cells.

The [search program](search_balanced_affine.py) used NumPy seed 468, a
40-node search rule and an 80-node rebalancing and evaluation rule.
It terminated each run when the Frobenius difference between the old and
new site tuples was less than $10^{-7}$, with a cap of 100 iterations.
The 76 initial tuples consisted of 32 generic Gaussian tuples, 24 tuples
compressed in one coordinate, and 20 asymmetric perturbations of three
planar sites and an axial site. The range of centroid norms at the first
iterate was nonconstant for every one of these runs. The spread between
largest and smallest norm varied between $0.00697$ and $0.30712$.

| Initialization | Runs | Iterations through termination, at most | Final normalized Gram discrepancy, at most |
| --- | ---: | ---: | ---: |
| Generic | 32 | 23 | $1.85\cdot10^{-7}$ |
| Compressed | 24 | 24 | $1.56\cdot10^{-7}$ |
| Axial | 20 | 18 | $1.80\cdot10^{-7}$ |

The discrepancy is $\|Z^\mathsf{T}Z-I/3\|_F$, with sites as the rows of
the centered unit Frobenius matrix $Z$. Zero discrepancy identifies the
regular tetrahedral Gram matrix. Estimated final first Hermite masses
agreed with the tetrahedral value within $2.2\cdot10^{-15}$ and the final
estimated cell masses were within $5\cdot10^{-12}$ of $1/4$. A tetrahedral
control reproduced its analytic value and two finite-difference checks
agreed with the off-center mass derivative and the exact local Hessian.
The quadrature and balance residuals are floating-point diagnostics.
This computation supplies no certified improvement or global exclusion.

From the repository root, the recorded Python 3.12.3 and NumPy 1.26.4 run was

```sh
OPENBLAS_NUM_THREADS=1 python3 research/plurality-stability/search_balanced_affine.py
```

If a future affine candidate shows a larger gap, certified cell mass and
centroid intervals can suffice even with approximately balanced biases.
For any partition into four cells of Gaussian masses $p_i$, let
$\delta=\frac12\sum_i|p_i-1/4|$. Atomless measure permits moving precisely
$\delta$ mass from surplus to deficit cells. Each cell changes in one
direction only. The centroid change from a subset of mass $u$ has norm at
most $\sqrt{u}$ by testing against unit directions and using Gaussian
Cauchy–Schwarz. The tuple of centroid changes therefore has Frobenius norm
at most $\sqrt{2\delta}$. The resulting exactly balanced partition has
$W_1\ge(\sqrt{W_1(A)}-\sqrt{2\delta})_+^2$. A strict certified excess in
this lower bound is sufficient for a small-positive-correlation witness.
