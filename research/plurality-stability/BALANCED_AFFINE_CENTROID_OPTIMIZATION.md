# Balanced affine optimization of the first Gaussian Hermite mass

The four-cell first Hermite problem reduces exactly to a three-dimensional
balanced affine assignment. The tetrahedron is a strict local maximizer in
centroid Gram space up to its orthogonal symmetries, with no equal-norm
assumption. If one first moment lies in the closed convex hull of the other
three, then $W_1\le0.3252948978196\ldots<1/3<W_1(T)$, by the exact bound
below. This holds for arbitrary measurable balanced partitions. A planar
centroidal configuration in strictly convex quadrilateral position is unstable
on adding a third Gaussian coordinate. Together these facts exclude every
rank-deficient global maximum in dimension at least three. The global theorem
therefore reduces to the full-rank stationary stability obligation. The
separate value comparison for all planar configurations remains open.
The balanced Hessian now has an explicit five-dimensional shape reduction,
a trace identity and exact tetrahedral eigenvalues. The local theorem below
also has an explicit polar Frobenius radius $10^{-5}$. These results do not
prove instability of every other high-valued full-rank stationary tuple.
The quantitative boundary theorem below excludes high-value full-rank
stationarity whenever $\lambda_{\min}(C)/\operatorname{tr}C\le10^{-11}$.
A deterministic search from 76 asymmetric initializations
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

The equal-mass, four-cell specialization of
[Heilman, *Stable Gaussian minimal bubbles*, Problem 1.15 and
Conjecture 1.16](https://doi.org/10.1007/s00526-026-03408-z)
is exactly this first-moment problem, apart from the positive factor
$\sqrt{\pi/2}$: equal simplex masses give $w=0$. The published statements
of Problem 1.15, Conjecture 1.16 and Theorem 1.6 were inspected. The latter
retains the first-moment conjecture as a hypothesis when $m>3$; the
balanced four-cell conclusion here would supply that specialization,
not the unequal-mass cases or the full positive-correlation simplex
conjecture. The boundary estimate below is derived directly from (3G)
and the stationary interface identities, without importing an estimate
from this reference or making a priority claim.

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

## Balanced second variation

Put $J(z)=2F(z)-\sum_i\|z_i\|^2$. At a full-rank tuple with balanced
biases $b$, write $d_{ij}=\|z_i-z_j\|$, $w_{ij}=a_{ij}/d_{ij}$, and let
$L$ be the Laplacian with these weights. For a site perturbation $V$, define

$$
q_{ij}(x)=(V_i-V_j)\cdot x,\qquad
r_i=\sum_{j\ne i}\frac1{d_{ij}}\int_{F_{ij}}q_{ij}(x)\phi_3(x)\,d\mathcal H^2(x),
\qquad
Q(V)=\sum_{i<j}\frac1{d_{ij}}\int_{F_{ij}}q_{ij}(x)^2\phi_3(x)\,d\mathcal H^2(x).
$$

Here $\sum_i r_i=0$. The exact Hessian after rebalancing is

$$
D^2J(z)[V,V]=2\left(Q(V)-r^\mathsf T L^+r-\sum_i\|V_i\|^2\right),
\tag{BV}
$$

where $L^+$ is the inverse on $\mathbf1^\perp$, extended by zero on
constants. At a full-rank stationary tuple, the six interface identities
give $L=I-\mathbf1\mathbf1^\mathsf T/4$, so the correction is exactly
$\|r\|^2$. Thus Gaussian first and second moments on six interfaces
determine a test for a positive second variation. The correction cannot
be dropped when testing stability of the balanced problem.

To prove (BV), perturb both sites and biases by $(V,c)$. The moving-face
calculation gives

$$
D^2G(z,b)[(V,c),(V,c)]
=\sum_{i<j}\frac1{d_{ij}}\int_{F_{ij}}
  (q_{ij}(x)+c_i-c_j)^2\phi_3(x)\,d\mathcal H^2(x)
=Q(V)+2r\cdot c+c^\mathsf TLc.
$$

The full-rank orthant charts used below also apply with nonzero biases,
so differentiation is justified by Gaussian domination. All six edge
weights are positive. Differentiating the balance equations in the
mean-zero gauge gives $Lc=-r$. The implicit function theorem applies
because $L$ is positive definite on that gauge. Substitution and
differentiation of the squared norm prove (BV).

This is the Gaussian whole-space calculation. The closest checked
transport result is [Kitagawa, Merigot and Thibert, Theorem 1.3 and
equation (1.8)](https://arxiv.org/pdf/1603.05579), which gives the
interface Hessian for bias variations with compact source support.
Their Theorem 1.5 also assumes compact support and does not by itself
give a convergence or global optimality theorem for this Gaussian
problem. The site variation, elimination of the balancing biases, and
the boundary argument below are derived here; no priority claim is made.

## Five shape directions and the averaged trace

Assume here that the tuple is centered, has rank three, is exactly balanced,
and is centroidal: $m_i=z_i$. Regard $Z$ as the $4$ by $3$ matrix of row
sites, and put $C=Z^{\mathsf T}Z$, $W=\operatorname{tr}C=J(Z)$.
All six weights equal $1/4$. The centered site space has dimension nine.
Its rotation directions are $ZA$ for $A^{\mathsf T}=-A$, and its radial
direction is $Z$. The five-dimensional orthogonal complement is

$$
\mathcal S_Z=\{V:\ \mathbf1^{\mathsf T}V=0,\quad
 Z^{\mathsf T}V\text{ is symmetric},\quad
 \langle Z,V\rangle_F=0\}. \tag{S5}
$$

The metric is the ordinary Frobenius metric of site perturbations.
For anisotropic $C$, taking $V=ZS$ with $S$ symmetric does not in
general give the orthogonal complement to rotations: the correct
condition is $CS$ symmetric.

**Shape reduction.** At such a stationary tuple, the balanced Hessian
vanishes on rotations, has radial eigenvalue $-2$, and restricts to (BV)
on $\mathcal S_Z$. Optimizing positive scale gives

$$
\Psi(Z)=\max_{t\ge0}J(tZ)=F(Z)^2/\|Z\|_F^2,\qquad
D^2\Psi[V,V]=D^2J[V,V]+2\langle Z,V\rangle_F^2/W. \tag{OS}
$$

Thus a positive eigenvalue on $\mathcal S_Z$ is exactly a positive
second variation after removing rotations and optimizing scale.
For the proof, positive homogeneity gives $D^2F[Z,\cdot]=0$.
Rotational invariance and the vanishing gradient of $J$ give
$D^2J[ZA,\cdot]=0$. Differentiating $F^2/\|Z\|_F^2$, using
$F=W$ and $\nabla F=Z$, proves (OS).

An explicit orthonormal frame is available. Diagonalize $C$ with
eigenvalues $\lambda_a>0$ and work in those physical coordinates.
Two diagonal directions are
$V=Z\operatorname{diag}(d_a/\sqrt{\lambda_a})$, where the two unit
vectors $d$ span $(\sqrt{\lambda_1},\sqrt{\lambda_2},
\sqrt{\lambda_3})^\perp$. For each $a<b$, the third through fifth
directions are $V=ZT$ with only

$$
T_{ab}=\sqrt{\frac{\lambda_b}{\lambda_a(\lambda_a+\lambda_b)}},
\qquad
T_{ba}=\sqrt{\frac{\lambda_a}{\lambda_b(\lambda_a+\lambda_b)}}
\tag{SF}
$$

nonzero. Direct multiplication verifies all five orthogonality and norm
conditions. Equal eigenvalues cause no singularity.

For each interface let $\mu_{ij}$ and $T_{ij}$ be its conditional
Gaussian first and raw second moments; set $\tau_{ij}=\operatorname{tr}T_{ij}$.
Define the symmetric cell second-moment matrices

$$
D_i=\int_{A_i}(xx^{\mathsf T}-I)\,d\gamma_3(x)
 =\sum_{j\ne i}(z_i-z_j)
      \left(\frac1{d_{ij}}\int_{F_{ij}}x\phi_3(x)\,d\mathcal H^2(x)\right)^{\mathsf T}.
$$

The second equality follows by integrating each derivative of
$x_b\phi_3$. In particular $\sum_iD_i=0$. For $V=ZT$, the bias
derivative in (BV) is $r_i=\langle D_i,T\rangle_F$.

**Shape trace identity.** For an orthonormal basis of $\mathcal S_Z$,
the sum of its five balanced Hessian values is exactly

$$
\operatorname{Tr}_{\mathcal S_Z}D^2J
 =\sum_{i<j}\tau_{ij}
       -2\sum_i\operatorname{tr}(D_i C^{-1}D_i)-16. \tag{ST}
$$

To verify this, first trace (BV) over the full centered nine-dimensional
site space. The functional $V\mapsto(V_i-V_j)\cdot x$ has squared
dual norm $2\|x\|^2$. Since $a_{ij}/d_{ij}=1/4$, the trace of $Q$
is $\frac12\sum_{i<j}\tau_{ij}$. The metric on $T$ is
$\operatorname{tr}(T^{\mathsf T}CT)$, so the squared dual norm of
$T\mapsto\langle D_i,T\rangle_F$ is
$\operatorname{tr}(D_i C^{-1}D_i)$. This is the full trace of the
balancing correction. The norm term has trace nine. Removing the
three zero rotation eigenvalues and the radial eigenvalue $-2$ gives
(ST). In particular, positivity of its right-hand side is sufficient
for an improving perturbation. Its positivity at all relevant
non-tetrahedral stationary tuples has not been proved.

There is also a specific anisotropy direction in this space:

$$
T=C^{-1}-\frac3W I,\qquad V=ZT,\qquad
\|V\|_F^2=\operatorname{tr}C^{-1}-\frac9W. \tag{AD}
$$

It is nonzero unless $C=(W/3)I$. Isotropy, together with centering and
rank three, means that the four sites form a regular tetrahedron;
balance fixes zero biases, and centroidal stationarity fixes the
tetrahedral scale. Thus (AD) is nonzero at every configuration still
under consideration. Its exact test is

$$
\frac12 D^2J[V,V]
 =\frac14\sum_{i<j}(T(z_i-z_j))^{\mathsf T}T_{ij}(T(z_i-z_j))
  -\sum_i\bigl(\operatorname{tr}(TD_i)\bigr)^2
  -\left(\operatorname{tr}C^{-1}-\frac9W\right). \tag{AT}
$$

Positivity of (AT) on every non-tetrahedral stationary tuple with
$W\ge W_1(T)$ would suffice. This is a candidate scalar inequality,
not an established sign or a replacement for the other shape directions.

### Radius-changing directions and the equal-norm bridge

There is a more targeted three-dimensional test. Set
$g_i=\|z_i\|^2$, $G=ZZ^{\mathsf T}$, $P=I-\mathbf1\mathbf1^{\mathsf T}/4$
and $N=\operatorname{diag}(g)-G/4$. For $a\in\mathbb R^4$ put

$$
V(a)=P\operatorname{diag}(a)Z,\qquad g^{\mathsf T}a=0.
\tag{RV}
$$

These directions belong to (S5): centering follows from $P$,
$Z^{\mathsf T}V=\sum_i a_i z_i z_i^{\mathsf T}$ is symmetric, and
$\langle Z,V\rangle_F=g^{\mathsf T}a=0$. Their metric and the
first variations of the four squared radii are

$$
\|V(a)\|_F^2=a^{\mathsf T}Na,\qquad
\left.\frac d{dt}\|z_i+tV_i(a)\|^2\right|_{t=0}=2(Na)_i.
\tag{RM}
$$

The matrix $N$ is positive definite. Indeed $V(a)=0$ would make
all $a_i z_i$ equal. A nonzero common vector would make every site
collinear. A zero common vector forces $a=0$, since a centered
rank-three four-tuple has no zero site. Also $N\mathbf1=g$.
Consequently $a\mapsto2Na$ maps $g^\perp$ bijectively onto
$\mathbf1^\perp$: (RV) realizes every infinitesimal change of
squared radii whose sum is zero. In particular this space has
dimension three. Its orthogonal complement inside (S5) consists
exactly of the two shape directions preserving all four radii,
because $\langle V(a),U\rangle_F=\sum_i a_i z_i\cdot U_i$.
If this vanishes for all $a\perp g$, then
$(z_i\cdot U_i)_i$ is proportional to $g$; its sum is zero by
(S5), so the proportionality constant is zero.

The complete radius Hessian can be assembled without differentiating
the balance solve. For each edge $e=(i,j)$ let $B_e$ be the $4$ by $3$
matrix with rows $z_i$ at $i$, $-z_j$ at $j$, and zero elsewhere, and put

$$
u_e=\frac1{d_{ij}}\int_{F_{ij}}x\phi_3\,d\mathcal H^2,\qquad
M_e=\frac1{d_{ij}}\int_{F_{ij}}xx^{\mathsf T}\phi_3\,d\mathcal H^2,
\qquad R=\sum_e(e_i-e_j)(B_eu_e)^{\mathsf T}.
$$

Here $e_i$ denotes a label coordinate vector. Translation cancels
from $V_i-V_j$, so $q_{ij}=a^{\mathsf T}B_ex$ and $r=Ra$.
At centroidal stationary tuples, (BV) therefore gives

$$
D^2J[V(a),V(a)]=a^{\mathsf T}H_Ra,\qquad
H_R=2\left(\sum_e B_eM_eB_e^{\mathsf T}-R^{\mathsf T}R-N\right).
\tag{RH}
$$

Away from stationarity the same formula uses $R^{\mathsf T}L^+R$.
If the columns of a $4$ by $3$ matrix $E$ span $g^\perp$, the sufficient
test is a positive generalized eigenvalue of
$(E^{\mathsf T}H_RE,E^{\mathsf T}NE)$. The second matrix is positive
definite. At the tetrahedron these three eigenvalues are exactly
$\lambda_{\rm off}$ in (TS). Indeed $z_i=c s_i$ makes the constraint
$\sum_i a_i=0$, and $V=ZS$ with
$S=(4c^2)^{-1}\sum_i a_i z_i z_i^{\mathsf T}$ symmetric and with
zero diagonal. Thus these directions fill the off-diagonal eigenspace.
A sign for unequal radii must use
stationarity and the high-value hypothesis, not just proximity to $T$.

The **equal-radius stability obligation remains open**: every centered
rank-three centroidal exactly balanced tuple with $J\ge W_1(T)$ and
nonpositive balanced second variation in every site direction should
have $g_1=g_2=g_3=g_4$. This is sufficient for the existing high-value
stability gap. The [equal-centroid-norm theorem](BALANCED_FOUR_CLASS_POLYNOMIAL_PARTITIONS.md#hermite-bounds)
already gives $W_1\le W_1(T)$ in that case, with equality only for the
tetrahedral partition. Thus a non-tetrahedral stationary tuple at or
above the tetrahedral value necessarily has unequal radii. Proving a
positive direction in (RH) for every such unequal-radius tuple is a
stronger, targeted sufficient test for the obligation. Neither this
sign nor the equal-radius stability obligation has been proved.
This bridge concerns first-Hermite mass and balanced quantization;
it supplies no uniform full-noise comparison.

### Exact tetrahedral spectrum

Write $z_i^*=c s_i$, with $s_i$ the even sign vectors as below, and let
$\theta=\arctan\sqrt2$, $t=\sin(2\theta)/(2\theta)$.
The conditional mean on a tetrahedral facet is $\eta u$, where $u$
is its shared-sign coordinate axis and
$\eta=\sqrt{\pi/3}/\theta$. This follows from the polar wedge
coordinates: the mean radius is $\sqrt{\pi/2}$ and the mean axial
angle factor is $\sin\theta/\theta$. Consequently
$D_i=(c\eta/2)(s_i s_i^{\mathsf T}-I)$.

For a symmetric perturbation $V=ZS$, the balancing correction is
$\eta^2/2$ times $\|V\|_F^2$ on off-diagonal $S$ and zero on
traceless diagonal $S$. Combining this with the unbalanced facet
calculation in the local proof gives the five exact eigenvalues

$$
\lambda_{\rm off}=-1+t-\frac{\pi}{3\theta^2}
       =-1.653996781499\ldots\quad(\text{multiplicity }3),
$$
$$
\lambda_{\rm diag}=-\frac12-\frac32t
       =-1.240180551423\ldots\quad(\text{multiplicity }2).
\tag{TS}
$$

Their sum is $-4-\pi/\theta^2=-7.442351447343\ldots$.
In (ST), each $\tau_{ij}=2$ and the correction trace is
$\pi/(2\theta^2)$. This checks the constants and shows why the
balancing correction materially changes the averaged Hessian. The
trace criterion is not a positive quantity throughout a neighborhood
of the tetrahedron; any use away from it must exploit stationarity
and the value threshold.

### Evaluating the interface moments

The face quantities require only one bivariate Gaussian probability,
not a numerical derivative of the balance solve. On a face write
$x=hn+y$, $y\perp n$. Let its two winning inequalities in that
plane be $u\cdot y\ge a$, $v\cdot y\ge b$, for tangent unit vectors
$u,v$ with $\rho=u\cdot v\in(-1,1)$ and $s=\sqrt{1-\rho^2}$.
Put

$$
p=\Phi_2(-a,-b;\rho),\quad
f=\phi(a)\Phi((\rho a-b)/s),\quad
g=\phi(b)\Phi((\rho b-a)/s),\quad
e=\phi(a)\phi((b-\rho a)/s).
$$

The unnormalized tangent moments are

$$
M=fu+gv,
$$
$$
N=p(I-nn^{\mathsf T})+afuu^{\mathsf T}+bgvv^{\mathsf T}
 +\frac e s\bigl(uv^{\mathsf T}+vu^{\mathsf T}
                    -\rho(uu^{\mathsf T}+vv^{\mathsf T})\bigr).
\tag{FM}
$$

Gaussian integration by parts on the two bounding lines proves these
formulas. The line at $u\cdot y=a$ has Gaussian area $f$ and
first moment $af u+e(v-\rho u)/s$; the other line gives the symmetric
expression. The ambient area, first moment and second moment are
$\phi(h)p$, $\phi(h)(hnp+M)$ and
$\phi(h)(h^2nn^{\mathsf T}p+h(nM^{\mathsf T}+Mn^{\mathsf T})+N)$.

The [shape diagnostic](check_balanced_shape_hessian.py) implements (FM),
the full balanced Hessian and (SF). At the tetrahedron it reproduces
(TS) and (ST). At the fixed asymmetric tuple it compares all nine
columns of the analytic Hessian with centered finite differences of
the rebalanced gradient: maximum entry discrepancy $7.06\cdot10^{-10}$
at step $2\cdot10^{-5}$. Omitting the correction changes that Hessian
by Frobenius norm $1.1746\ldots$. The
[receipt](balanced_shape_hessian_verification.json) records Python 3.12.3,
NumPy 1.26.4 and the 80-node rule. These are floating-point controls;
the formulas have the analytic proofs above, and no global stationary
classification or interval sign certificate is supplied.
The same diagnostic now checks (RM), compares (RH) with compression
of the full nine-dimensional Hessian at the tetrahedron and the fixed
asymmetric tuple, and recovers all three tetrahedral radius eigenvalues.
The maximum compression discrepancy was $8.89\cdot10^{-16}$.
The asymmetric tuple is balanced but not stationary, so its calculation
uses the full $L^+$ and supplies no evidence for the missing stationary
positive-sign assertion.

## Lifting a planar stationary configuration

Let $z_1,\ldots,z_4\in\mathbb R^2$ be distinct, centered sites with no
three collinear. Choose their balanced affine cells and assume the
centroidal equations $m_i=z_i$. Edges with no positive-length interface
have $a_{ij}=w_{ij}=0$. Gaussian integration by parts still gives
$LZ=Z$, where $Z$ has rows $z_i$. For any $v\in\mathbf1^\perp$, lift
the sites to $z_i(t)=(z_i,tv_i)$ in $\mathbb R^3$ and rebalance. Then

$$
\left.\frac{d^2}{dt^2}J(z(t))\right|_{t=0}
=2v^\mathsf T(L-I)v.
\tag{PL}
$$

Here is a justification at the rank-deficient tuple, where a full-rank
orthant chart is unavailable. With the planar slopes fixed, put
$h(b)=\mathbb E_X\max_i(z_i\cdot X+b_i)$. Its bias Hessian is $L(b)$.
Indeed each equality line contains a winning interval, ray, or empty
set; differentiating a cell mass gives $a_{ij}/d_{ij}$. There are no
triple ties along a line because no three sites are collinear. The
interval endpoints vary continuously, including when a facet shrinks
to a point. Hence the face integrals, and the Hessian, are continuous.
They are globally bounded: each Gaussian line area is at most
$1/\sqrt{2\pi}$ and all $d_{ij}>0$ are fixed. The adjacency graph is
connected, since a generic segment between cell interiors crosses
facets and avoids the finitely many multiple-tie points. All cells have
positive area by balance. Consequently $L(b)$ is positive definite
on $\mathbf1^\perp$ at the balanced biases.

For an independent standard normal $Y$, the lifted unreduced objective
is $\mathbb E_Y h(b+tvY)-\frac14\sum_i b_i$. Its mixed $t,b$
derivative at zero vanishes because $\mathbb E Y=0$, and its second
$t$ derivative is $v^\mathsf TLv$ because $\mathbb E Y^2=1$.
The bounded bias Hessian justifies differentiation under the expectation.
The balance equations and the implicit function theorem give $b'(0)=0$.
The resulting biases are minimizing biases by convexity. Differentiating
the site norm proves (PL).

**Convex quadrilateral exclusion.** Under these hypotheses, if all four
sites are vertices of their convex hull, the remaining eigenvalue
$\lambda$ of $L$ on $\mathbf1^\perp$ satisfies $\lambda\ge2$.
In particular there is a unit lift with $J''(0)\ge2$, so this tuple
cannot be a local maximum in dimension three.

For the proof, let $v$ be a unit vector satisfying
$v\perp\mathbf1$ and $Z^\mathsf Tv=0$. The affine dependence of a
strict convex quadrilateral has two positive and two negative
coefficients, with each sign pair on opposite vertices. This follows
by writing the intersection of the diagonals as a strict convex
combination of each pair of endpoints. Symmetry of $L$, together with
$LZ=Z$ and $L\mathbf1=0$, gives

$$
L=P+(\lambda-1)vv^\mathsf T,\qquad
P=I-\mathbf1\mathbf1^\mathsf T/4.
$$

At least one of the two same-sign pairs has zero interface area. To see
this directly, set $s_i(x)=z_i\cdot x+b_i$ and $k=\sum_i v_i b_i$,
so $\sum_i v_i s_i(x)=k$ for every $x$. On a winning tie of the
positive pair this sum is nonnegative, and it is strictly positive
except at an all-four tie. Thus a positive-length interface for that
pair requires $k>0$. An interface for the negative pair similarly
requires $k<0$. An all-four tie is at most one point because the slopes
span the plane. Both same-sign pairs therefore cannot have positive
interface length.

For an absent same-sign pair $i,j$, the off-diagonal entry of $L$ yields

$$
0=-\tfrac14+(\lambda-1)v_i v_j,\qquad
\lambda-1=\frac1{4v_i v_j}\ge1.
$$

For the last inequality, write $s=v_i+v_j=-(v_k+v_l)$. The unit norm
gives $1\ge s^2$, and $v_i v_j\le s^2/4\le1/4$.
Equation (PL) now proves the strict improvement. Since $v$ is outside
the column space of $Z$, the lifted tuples have rank three for $t\ne0$.
Their balanced assignments have $W_1\ge J(z(t))>J(z)=W_1(A)$ for
sufficiently small nonzero $t$, by completing squares.

The square is an exact control. Its centroid sites are
$(\pm a,\pm a)$ with $a=1/(2\sqrt{2\pi})$. Adjacent sites have
distance $1/\sqrt{2\pi}$ and interface area $1/(2\sqrt{2\pi})$.
Thus $L$ is one half the four-cycle Laplacian, with eigenvalues
$0,1,1,2$. The lift has $J''(0)=2$. The spectral bound above is sharp.

For an internal floating-point control, the following command uses the
[variation diagnostic](check_balanced_variations.py) and the existing
facet integrator without changing the integrator. On Python 3.12.3 and
NumPy 1.26.4 the three printed quotients were respectively
$0.9784197060$, $0.9935889388$, and $0.9978691748$, approaching the
analytic value $(J(t)-J(0))/t^2\to1$. The existing full-rank facet
finite-difference controls also passed. At the fixed asymmetric tuple
retained in the diagnostic, the balancing correction was $0.1439216576$;
the Schur-complement Hessian and a direct finite difference of the
rebalanced gradient differed by $2.47\cdot10^{-10}$ at step $2\cdot10^{-5}$.
The rule has 80 Gauss--Legendre nodes. These are diagnostics; the proof
is the analytic argument above.

```sh
python3 research/plurality-stability/check_balanced_variations.py
```

This exclusion assumes centroidal stationarity and quadrilateral
position. It does not bound $W_1$ for every planar affine assignment.
The next argument handles all non-extreme first-moment configurations,
including an interior centroid, without using the lifting formula.

## A first moment inside the triangle of the other three

Let $A_1,\ldots,A_4$ be any measurable balanced partition of a standard
Gaussian space, in any dimension. Suppose, after relabeling,
$m_4\in\operatorname{conv}\{m_1,m_2,m_3\}$. The convex hull is closed:
an edge point, a repeated first moment, or a degenerate triangle is allowed.
No affine-assignment or stationarity hypothesis is imposed. Set

$$
q=\Phi^{-1}(3/4),\qquad c=\phi(q),\qquad C=\frac9{8\pi},\qquad
B=\frac{\sqrt{C^2/9+16Cc^2}-C/3}{2}
 =\frac3{16\pi}\left(\sqrt{1+64e^{-q^2}}-1\right).
$$

**Non-extreme first-moment bound.** Under these hypotheses,

$$
W_1(A)\le B=0.3252948978196\ldots<\frac13<W_1(T).
\tag{NE}
$$

The idea is to combine two estimates which constrain the fourth moment
in opposite ways. Put $W=W_1(A)$ and $t=\|m_4\|^2$. Every cell of mass
$1/4$ has $\|m_i\|\le c$. Indeed, for a unit vector $u$ and
$Y=u\cdot X$, the pointwise inequality
$(Y-q)(1_{A_i}-1_{\{Y\ge q\}})\le0$, together with the equal masses,
gives $u\cdot m_i\le\mathbb E[Y1_{\{Y\ge q\}}]=\phi(q)$.
Choose $u$ in the direction of $m_i$, with the zero case immediate.
Consequently

$$
W\le3c^2+t. \tag{NE1}
$$

On the other hand the fourth linear score is pointwise dominated by
the largest of the first three scores. Thus

$$
W=\sum_i\int_{A_i}m_i\cdot x\,d\gamma_n(x)
 \le\mathbb E\max_i m_i\cdot X
 =\mathbb E\max_{1\le i\le3}m_i\cdot X.
$$

For any three vectors $u_1,u_2,u_3$, with
$\bar u=(u_1+u_2+u_3)/3$, the exact three-score identity and its
Cauchy--Schwarz consequence are

$$
\mathbb E\max_{1\le i\le3}u_i\cdot X
=\frac{\sum_{i<j}\|u_i-u_j\|}{2\sqrt{2\pi}},\qquad
\left(\mathbb E\max_{1\le i\le3}u_i\cdot X\right)^2
\le C\sum_{i=1}^3\|u_i-\bar u\|^2. \tag{3G}
$$

For completeness, for three real numbers their range is half the sum
of the three pairwise absolute differences. The centered Gaussian
score vector is invariant in law under simultaneous negation, so its
expected maximum is half its expected range. The absolute Gaussian
first moment gives the identity. Squaring, applying Cauchy--Schwarz to
the three lengths, and using
$\sum_{i<j}\|u_i-u_j\|^2=3\sum_i\|u_i-\bar u\|^2$ proves the
inequality. This also covers coincident or collinear vectors.
The constant and Gaussian-maximum interpretation agree with
[Khot and Naor, Corollary 3.4 and Lemma 3.7](https://cs.nyu.edu/~khot/papers/kernel.pdf).
The primary statements and their proofs were inspected. Here (3G) is
derived directly; no four-cell propeller theorem or novelty claim is used.

Since $\sum_i m_i=0$, the first three moments have mean $-m_4/3$ and
their centered squared norm sum is $W-4t/3$. Therefore

$$
W^2\le C(W-4t/3). \tag{NE2}
$$

Substitute $t\ge W-3c^2$ from (NE1) into (NE2). This gives
$W^2+(C/3)W-4Cc^2\le0$. Its nonnegative root is $B$, proving the
stated value bound. This substitution remains valid if $W-3c^2<0$;
no division by $W$, $t$, or a first-moment norm is needed.

The strict comparison with $T$ has an elementary rational check, so it
does not depend on the displayed decimal. Using $25/8<\pi<22/7$ and
$e^{-s^2/2}\le1-s^2/2+s^4/8$ gives

$$
\Phi(2/3)-\tfrac12
 <\tfrac25\left(\tfrac23-\frac{(2/3)^3}{6}+\frac{(2/3)^5}{40}\right)
 =\frac{1508}{6075}<\tfrac14.
$$

Hence $q>2/3$. The exponential series gives
$e^{4/9}>1+4/9+(4/9)^2/2=125/81$, so
$c^2<81/(250\pi)$. The positive root $B$ consequently satisfies, with
$u=\pi B$, $u^2+3u/8<729/500$. The polynomial is increasing for
$u\ge0$, and

$$
(25/24)^2+\tfrac38(25/24)-\tfrac{729}{500}
=\frac{637}{36000}>0.
$$

It follows that $B<25/(24\pi)<1/3$. Finally the tetrahedral formula
and $\arcsin(1/3)>1/3$ imply

$$
W_1(T)>\frac{400}{363\pi}>
\frac7{22}\frac{400}{363}=\frac13+\frac{23}{1331}.
$$

Thus (NE) excludes the requested interior-centroid case by value. It
also excludes every configuration with a first moment on an edge of
the other three, and every repeated first moment. For a nonstationary
affine assignment the hypothesis concerns its actual first moments;
the locations of its input score sites alone do not imply it.
The [variation diagnostic](check_balanced_variations.py) checks the
displayed rational margins using exact fractions and prints $B$ using
the standard-library Gaussian quantile. The decimal is illustrative;
the proof of the strict comparisons is the analytic argument above.

## An explicit stationary exclusion near the rank-deficient boundary

**Quantitative stationary boundary theorem.** Let $Z$ be a centered,
rank-three tuple in $\mathbb R^3$ whose exactly balanced affine cells
have first moments $m_i=z_i$. Put $C=Z^{\mathsf T}Z$ and
$W=\operatorname{tr}C=J(Z)$. If $W\ge W_1(T)$, then

$$
\frac{\lambda_{\min}(C)}{\operatorname{tr}C}>10^{-11}. \tag{QB}
$$

Thus the requested improving-perturbation assertion holds on the region
with ratio at most $\eta=10^{-11}$ because there is no high-value
stationary tuple there. This is an exclusion by the first-order
equations, stronger than a positive-Hessian test in this region.
It asserts neither a planar value bound nor a value bound for every
near-planar nonstationary assignment. The proof uses Gaussian strip
and radial tail integrals analytically, with no quadrature.

**1. High value forces a uniform convex-hull margin.** For the first
moments of any balanced partition, in any Gaussian dimension, set
$W=\sum_i\|m_i\|^2$. Then

$$
W\ge\frac7{20}
\quad\Longrightarrow\quad
\operatorname{dist}(m_i,\operatorname{conv}\{m_j:j\ne i\})
 >\frac1{25}\quad\text{for every }i. \tag{RE}
$$

Here is a quantitative extension of (NE), including its constants.
Take $i=4$, choose a closest convex combination $p$ of the other
three moments, and put $\delta=\|m_4-p\|$, $t=\|m_4\|^2$.
Pointwise,

$$
\max_{1\le i\le4}m_i\cdot X
 \le\max_{1\le i\le3}m_i\cdot X+((m_4-p)\cdot X)_+.
$$

The last summand has expectation $\delta/\sqrt{2\pi}$. The assignment
bound, (3G), centering and the cell moment bound $\|m_i\|\le c$
therefore give

$$
W\le\sqrt{\frac9{8\pi}(W-4t/3)}+\frac\delta{\sqrt{2\pi}},
\qquad t\ge W-3c^2.
$$

The already proved bounds $c^2<81/(250\pi)$ and $\pi>25/8$ give
$c^2<324/3125$ and $9/(8\pi)<9/25$. If $W\ge7/20$ and
$\delta\le1/25$, then $W-\delta/\sqrt{2\pi}>33/100$, whereas

$$
\left(W-\frac\delta{\sqrt{2\pi}}\right)^2
 \le\frac9{8\pi}(4c^2-W/3)
 <\frac9{25}\left(\frac{1296}{3125}-\frac7{60}\right)
 <\left(\frac{33}{100}\right)^2.
$$

This contradiction proves (RE). In particular, a centroidal tuple at
the tetrahedral value level satisfies (RE), since the earlier rational
bound $W_1(T)>1/3+23/1331$ exceeds $7/20$. Also
$W\le4c^2<1296/3125<1/2$, and every site norm is less than $1/2$.

**2. Flattening gives a uniformly nondegenerate quadrilateral.**
Suppose, contrary to (QB), that $\lambda/W\le10^{-11}$, where
$\lambda=\lambda_{\min}(C)>0$. Then

$$
\sqrt\lambda<\frac1{400000}, \tag{BS}
$$

because $W<1/2$ and $(400000)^{-2}> (2\cdot10^{11})^{-1}$.
Choose a unit eigenvector $u$ for $\lambda$, set
$v=Zu/\sqrt\lambda$ and let $y_i=z_i-\sqrt\lambda v_i u$.
Then $\|v\|=1$, $\mathbf1^{\mathsf T}v=0$,
$Z^{\mathsf T}v=\sqrt\lambda u$, and $Y^{\mathsf T}v=0$.
The matrix $Y$ has rank two. For any convex coefficients on the other
three sites, the change in the corresponding difference vector is
at most $2\sqrt\lambda$. Thus (RE) and (BS) imply

$$
\operatorname{dist}(y_i,\operatorname{conv}\{y_j:j\ne i\})
 >\frac1{50}. \tag{PQ}
$$

All four projected sites are vertices of a strict convex quadrilateral.
Its affine dependence $v$ has two positive and two negative entries;
none is zero. Its diameter is less than one because orthogonal
projection does not increase distances and $\|z_i\|<1/2$.

In fact

$$
\min_i|v_i|>\frac1{200}. \tag{VC}
$$

To check this, let $s$ be the common sum of the positive entries and
of the absolute negative entries. Since $1=\|v\|^2\le2s^2$,
the larger entry $b$ in either sign pair exceeds $1/4$. Denote the
smaller entry by $a$, with its site $y_i$, and the larger-entry site
by $y_j$. The affine dependence gives
$a y_i+b y_j=s p$ for a convex combination $p$ of the opposite pair.
Consequently $y_j-p=(a/b)(p-y_i)$, so its distance from that opposite
segment is less than $a/b$. By (PQ) it exceeds $1/50$, whence
$a>b/50>1/200$. This works for both sign pairs.

**3. Every site triangle has a uniform altitude.** From (RE), all
pairwise distances exceed $1/25$. In any three-site triangle the
altitude to a longest side lands in that side's segment: both angles
at the longest side's endpoints are acute. The opposite vertex is
at distance greater than $1/25$ from this segment, by (RE). This
altitude is the smallest of the three altitudes, so all three exceed
$1/25$. Collinear triangles are already excluded by (RE). Hence for
any three distinct labels $i,j,k$, if $n$ is the unit normal to the
equality plane of $i,j$, then

$$
\|(I-nn^{\mathsf T})(z_i-z_k)\|>\frac1{25}. \tag{TA}
$$

**4. One stationary interface is too small.** Let $b$ be the exact
balancing biases and write $s_i(x)=z_i\cdot x+b_i$. Reverse $v$ and
$u$ together if necessary so that $k=v\cdot b\le0$, and let $i,j$
be the positive sign pair of $v$. On their winning interface,
$s_i=s_j\ge s_l,s_m$, where $l,m$ are the negative pair. Therefore

$$
|v_l|(s_i-s_l)+|v_m|(s_i-s_m)
 =\sum_r v_r s_r(x)=\sqrt\lambda\,u\cdot x+k.
$$

On the part with $\|x\|\le4$, both gaps on the left are nonnegative,
and (VC) bounds either one by $800\sqrt\lambda$. In particular this
part of the face lies in the strip

$$
0\le (z_i-z_l)\cdot x+b_i-b_l\le800\sqrt\lambda
$$

inside the equality plane $H_{ij}$. By (TA), the tangent gradient of
the middle expression has length greater than $1/25$. On any plane
$H=\{x:n\cdot x=h\}$, a strip of coordinate width $\ell$ and tangent
gradient length $a$ has Gaussian surface integral at most
$\ell/(2\pi a)$. Indeed write $x=hn+y$ with a standard two-dimensional
tangent Gaussian; its one-dimensional marginal density is at most
$(2\pi)^{-1/2}$, and the normal density has the same bound.
The inner part of this face thus has area less than
$20000\sqrt\lambda/(2\pi)<1/120$, by (BS) and $\pi>3$.

Uniformly over the plane's offset, its outside-ball area satisfies

$$
\int_{H\cap\{\|x\|>4\}}\phi_3(x)\,d\mathcal H^2(x)
 \le\frac{e^{-8}}{\sqrt{2\pi}}<\frac1{1000}. \tag{PT}
$$

For $|h|\le4$, the two-dimensional radial tail is
$\exp(-(16-h^2)/2)$, which multiplied by $\phi(h)$ gives the
displayed bound. For $|h|>4$, the entire plane has area $\phi(h)$
bounded by the same expression. Finally
$e^8>\sum_{r=0}^{5}8^r/r!=8551/15>500$ and
$(2\pi)^{-1/2}<1/2$ prove the rational comparison in (PT).

It follows that the interface area obeys

$$
a_{ij}<\frac1{120}+\frac1{1000}=\frac7{750}<\frac1{100}.
$$

Full-rank centroidal stationarity, however, requires
$a_{ij}=\|z_i-z_j\|/4>1/100$, by the six interface identities and
(RE). This is the contradiction proving (QB).

The [exact constant audit](check_stationary_boundary_margin.py) checks
the rational comparisons, including the strip/tail contradiction and
the admissible $\eta$, with Python's standard-library fractions. It
does not certify the geometric argument or replace its analytic proof.
An internal Python 3.12.3 run passed with robust value-squared margin
$2001/1250000$ and face-area contradiction margin $1/1500$; three
deliberately invalid constant choices were rejected. Reproduce it from
the repository root with
`python3 research/plurality-stability/check_stationary_boundary_margin.py`.
No certified balance solve or numerical Gaussian integration is used
in this boundary theorem. The high-value stability gap on the remaining
compact shape region is still open.

## Boundary exclusion and the remaining global obligation

The rank-one case also has a sharper bound. If the first moments are
collinear, write $m_i=\alpha_i u$ with $u$
unit and sort the $\alpha_i$. The balanced linear assignment to these
scores consists of the four consecutive quartile intervals of
$u\cdot X$. For example, expand a sorted score as its smallest value
plus the three nonnegative successive differences. Each difference
multiplies the integral over a union of the highest-score cells of
fixed mass; that integral is maximized by an upper Gaussian tail.
This proves the asserted quartile assignment, also when scores tie.
Writing $q=\Phi^{-1}(3/4)$, its scalar moments are

$$
c=(-\phi(q),\ \phi(q)-\phi(0),\ \phi(0)-\phi(q),\ \phi(q)).
$$

The assignment bound and Cauchy--Schwarz give
$W_1(A)\le\sum_i\alpha_i c_i\le\sqrt{W_1(A)}\|c\|$, whence

$$
W_1(A)\le2\bigl(\phi(q)^2+(\phi(0)-\phi(q))^2\bigr)<\frac1\pi<W_1(T).
\tag{R1}
$$

The zero-moment case is immediate. Strictness follows because both
$\phi(q)$ and $\phi(0)-\phi(q)$ are positive and sum to $\phi(0)$.

**No rank-deficient global maximum.** Every balanced four-cell
first-Hermite global maximizer in dimension $n\ge3$ has first-moment
rank three. To prove this, suppose the four moments span at most a
plane. Four planar points either have a point in the closed convex
hull of the other three or are the four distinct vertices of a strict
convex quadrilateral. This follows by considering the number of
vertices of their convex hull; repeated and collinear cases are
included. In the first case (NE) puts the value below the feasible
tetrahedron. In the second case, a global maximizing partition must
be a balanced affine assignment to its own moments: if $F(m)>W_1$,
then $J(m)>W_1$ would contradict global maximality. Equality in the
pointwise dual bound forces the affine winner almost everywhere.
The distinct planar slopes have no positive-measure ties, so this
partition depends only on the moment plane. It is centroidal, and
the convex-quadrilateral lift in an unused Gaussian coordinate
strictly improves its value. This again contradicts maximality.

The remaining sufficient statement can be restricted to the tetrahedral
value level:

**High-value full-rank stability gap (open):** every centered rank-three
centroidal balanced affine tuple with $J\ge W_1(T)$, other than the
tetrahedral centroid tuple up to orthogonal maps and label permutations,
has an improving site perturbation for $J$. A positive value of (BV)
in one direction is sufficient. This statement remains open.

The immediate targeted route is the equal-radius stability obligation
above. A positive direction in the three-dimensional radius restriction
(RH) whenever the radii differ, at these high-value stationary tuples,
would connect the shape Hessian directly to the established equal-norm
theorem. It would avoid classifying lower-valued stationary tuples.

**Full-rank stability gap (open):** the same improving-perturbation
statement for every non-tetrahedral full-rank stationary tuple, without
restricting its value. This is a stronger sufficient statement; its
restriction gives the high-value gap. No classification of all the
remaining stationary tuples has been proved here.

Either stability statement implies the global first-Hermite inequality
and its tetrahedral equality characterization. The compact affine reduction
supplies a global maximum at its natural scale. Completing squares
makes it centroidal, and the boundary exclusion makes it full rank.
An improving perturbation would contradict maximality, so the remaining
statement forces the tetrahedron. Equality for any original partition
forces a maximizing assignment to tetrahedral centroid sites, whose
linear winner is unique almost everywhere. The quantization identity
transfers this conclusion to balanced four-symbol squared distortion
for every $n\ge3$. No positive-correlation noise-stability comparison
follows from this reduction.

The [balanced reflected-swap results](RELATIVE_HERMITE_TAIL_OBSTRUCTION.md)
separate the full-noise task from this first-Hermite route. Remote tail
swaps falsify uniform relative higher-Hermite control by the first-Hermite
deficit, so that bridge is retired. Nevertheless every proper reflected
swap of a tetrahedral cell pair strictly loses full noise stability for
every $0<\rho<1$, by an exact Gaussian kernel comparison. That theorem
does not reduce arbitrary balanced competitors to reflected swaps.
The high-value gap above can be attacked with the five shape directions;
a uniform full-noise comparison still requires its own argument.

The **planar value gap remains open** as a separate comparison: every
balanced four-cell affine assignment in $\mathbb R^2$, including ties
subdivided into equal masses, should have $W_1<W_1(T)$. It is no longer
a required premise of the route above. A planar configuration violating
this value bound would have to have its first moments in strict convex
quadrilateral position, by (NE). Its instability in dimension three
does not by itself rule out a large planar value.

Conversely, the global first-Hermite inequality together with the proved
boundary exclusion would settle this strict planar comparison. Extend a
planar balanced partition to dimension three by an independent Gaussian
coordinate. Its first moments and $W_1$ are unchanged and have rank at most
two. A value above $W_1(T)$ would violate the global inequality. Equality
would make the extended partition a rank-deficient global maximizer,
contrary to the boundary exclusion. Thus the planar comparison is a
consequence of the global bound, even though it is not a required premise
of the remaining stability route.

For a computer-assisted route, the explicit local theorem below supplies
a polar Frobenius neighborhood of radius $10^{-5}$. A complete interval proof
still needs rigorous Gaussian integral enclosures and an exhaustive cover
of the remaining compact normalized Gram domain.
For upper bounds on the value, **a certified bias solve is unnecessary**:
every fixed mean-zero bias $\beta$ gives $F(z)\le G(z,\beta)$, even
if the resulting cells are unbalanced. On a box of centered sites of
unit Frobenius norm, a rigorous uniform bound
$G(z,\beta)\le\sqrt{W_1(T)}$ suffices. Different boxes may use
different fixed biases. Certification of stationary tuples, in contrast,
must enforce the balance equations before using (BV) or the six
interface identities. Compactness alone supplies no numerical certificate.
The analytic rank-deficient exclusion is sufficient when searching for
global maxima. In addition, (QB) now confines every high-value full-rank
stationary candidate to $\lambda_{\min}(C)/\operatorname{tr}C>10^{-11}$;
(RE) excludes coincident sites with an explicit separation. Together
with the local exclusion (NR), this leaves a compact domain of normalized
shapes on which an exhaustive interval cover can be attempted. The
Gaussian integrals for balance, stationarity and Hessian rejection on
that domain still require certified enclosures. No exhaustive cover or
practical complexity bound has been established.

The balancing biases also have an explicit compact bound after unit
Frobenius normalization. In the gauge $\sum_i b_i=0$, every minimizing
bias satisfies $\|b\|_\infty\le3$: the dual coercivity estimate gives
$G(Z,b)\ge\max_i b_i\ge\|b\|_\infty/3$, whereas
$G(Z,b)\le G(Z,0)\le\mathbb E\|ZX\|\le\|Z\|_F=1$.
This bounds the possible bias solutions; it does not certify the
balance equations or the required Gaussian integrals on a shape box.

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

## An explicit tetrahedral neighborhood

The local exclusion can be quantified without solving any biases. Keep
$z_i^*=c s_i$ at its centroidal scale and $W_*=W_1(T)$, and put
$h(x)=\max_i z_i^*\cdot x$. For $P=I+S$ symmetric with
$\|S\|_F\le10^{-5}$, the unbalanced upper objective satisfies

$$
\mathcal U((Pz_i^*)_i)
 \le W_*-\frac1{100}\|S\|_F^2. \tag{LR}
$$

Thus every four-cell partition, balanced or otherwise, whose centroid
Gram matrix has this polar representative lies strictly below the
tetrahedron unless $S=0$. No claim about the entire compact shape
domain follows from this radius.

Here is a bound for the Taylor remainder. For $P_t=I+tS$, set
$H(t)=\mathbb E h(P_tX)$ and change variables in its Gaussian integral.
For fixed $y$, write $x=P_t^{-1}y$, $B_t=SP_t^{-1}$ and let $l_t(y)$
be the logarithm of its density. Since $S$ commutes with $P_t$,

$$
l'=-\operatorname{tr}B_t+x^{\mathsf T}B_tx,\quad
l''=\operatorname{tr}B_t^2-3x^{\mathsf T}B_t^2x,\quad
l'''=-2\operatorname{tr}B_t^3+12x^{\mathsf T}B_t^3x.
$$

Writing $R=\|X\|$, bounding each trace by three times its operator
norm, and using the third derivative of an exponential gives

$$
|H'''(t)|\le
 a\frac{1+r}{(1-r)^3}\|S\|_F^3\,
 \mathbb E\bigl[R(60+75R^2+18R^4+R^6)\bigr],
\quad a=\max_i\|z_i^*\|,\quad r=\|S\|_F.
$$

For three standard Gaussian coordinates,
$\mathbb ER^3=4\mathbb ER$, $\mathbb ER^5=24\mathbb ER$ and
$\mathbb ER^7=192\mathbb ER$, so the expectation is
$984\mathbb ER$. The tetrahedral formula and
$\arcsin(1/3)<\pi/6$ give $W_*<4/(3\pi)<1$.
We have $a=\sqrt{W_*}/2<1/2$,
$\mathbb ER<2$, and $(1+r)/(1-r)^3<2$ for $r\le1/10$.
Thus $|H'''(t)|\le1968\|S\|_F^3$. The change of variables
also justifies differentiation: on $r\le1/10$ all density derivatives
are bounded by a fixed polynomial times a decaying Gaussian.

The local facet calculation gives
$D^2H(0)\le(5/6)\sum_i\|S z_i^*\|^2$.
Indeed $\theta>\pi/4>\sqrt2/2$ implies
$t=\sqrt2/(3\theta)<2/3$, so both $\mu_u/2$ and
$3\mu_v/4$ are less than $5/6$. Since $W_*>1/3$,
the second derivative of $\mathcal U$ is at most
$-W_*\|S\|_F^2/9\le-\|S\|_F^2/27$.
Its first derivative vanishes, and its third derivative is bounded
by $3936\|S\|_F^3$. Taylor's theorem yields

$$
\mathcal U((Pz_i^*)_i)
 \le W_*-\frac1{54}\|S\|_F^2+656\|S\|_F^3
 \le W_*-\frac1{100}\|S\|_F^2,
$$

where the last strict arithmetic margin at $r=10^{-5}$ is
$1/54-656/100000-1/100=661/337500>0$.
The shape diagnostic checks this rational arithmetic; the radius is
proved by the uniform derivative bound, not by sampled perturbations.

**Normalized shape neighborhood.** For a centered rank-three tuple of
unit Frobenius norm, let $C=Z^{\mathsf T}Z$. If

$$
\delta(Z):=\|\sqrt{3C}-I\|_F\le\frac1{300000}, \tag{NR}
$$

then $F(Z)^2<W_*$ unless its shape is regular tetrahedral. This bound
does not require stationarity. To prove it, rotate the tuple to
$Z=Z_0P$, where $Z_0=z^*/\sqrt{W_*}$ and $P$ is positive definite.
Then $\|P\|_F=\sqrt3$ and $\|P-I\|_F=\delta(Z)=\rho$.
Bessel's inequality gives $W_1(A)\le3/4$ for every balanced partition,
so its support function $F$ is $\sqrt3/2$-Lipschitz in the site norm.
Consequently $|F(Z)-\sqrt{W_*}|\le\rho/2$.
At the optimizing scale the tuple is $F(Z)Z$, whose polar factor
relative to $z^*$ is $P_*=F(Z)P/\sqrt{W_*}$. Since $W_*>1/3$,

$$
\|P_*-I\|_F
 \le\frac{\rho}{2\sqrt{W_*}}\sqrt3+\rho<3\rho\le10^{-5}.
$$

Apply (LR) and $J\le\mathcal U$. The only exception $P_*=I$
forces $P=I$ by the unit normalization. In particular any remaining
non-tetrahedral high-value stationary tuple must satisfy
$\|\sqrt{3C/W}-I\|_F>1/300000$. This supplies an explicit excluded
region near isotropy, while leaving the sign of (ST) or (AT) outside
that region open.

The checked transport and propeller results cited above likewise do not
supply the sign in (ST) or (AT). These formulas and the local radius do
not assert that the global conjecture has been resolved.

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
