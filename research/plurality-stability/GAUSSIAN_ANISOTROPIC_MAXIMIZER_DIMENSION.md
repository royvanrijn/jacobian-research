# Active coordinates for anisotropic Gaussian noise stability

Let \(R\) be a symmetric matrix on \(\mathbb R^d\) with
\(0<R<I\), and let \((X,Y)\) have standard Gaussian marginals
and cross-covariance \(R\). For a partition \(f=(f_i)_{i=1}^k\)
with prescribed positive masses, set
\[
 S_R(f)=\sum_i\langle f_i,T_R f_i\rangle_{\gamma_d},\qquad
 T_Rg(x)=\mathbb E g\bigl(Rx+(I-R^2)^{1/2}Z\bigr).
\]
An actual maximizer is cylindrical over the span of its cell
centroids. This span is invariant under \(R\). If its dimension
is \(m\), the eigenvalues on it are the largest \(m\) eigenvalues
of \(R\), counted with multiplicity. In an eigenspace tied at the
cutoff any invariant subspace of the required dimension has the
same eigenvalue. In particular the optimal value is attained on
the strongest \(\min(d,k-1)\) coordinates. The active rank can be
lower than this upper bound.

The proof adapts the
[scalar bulk translation argument](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#the-bulk-translation-identity)
to matrix correlations. The extra input is stationarity under
ordinary rotations, which makes the matrix of covariance
derivatives commute with \(R\).

## Covariance derivatives and rotation stationarity

For bounded indicator tuples define the symmetric quadratic form
\[
 M_f(v,w)=\sum_i\langle D_vf_i,T_R D_wf_i\rangle_{\gamma_d}.
                                                               \tag{AM1}
\]
The derivatives here are distributions and the expression is the
inner product after smoothing both arguments by \(T_{R^{1/2}}\).
It is finite and positive semidefinite. For a symmetric matrix
\(B\), the covariance derivative at fixed tuple \(f\) is
\[
        dS_R(f)[B]=\operatorname{tr}(B M_f).              \tag{AM2}
\]
Indeed, identify the Hermite chaos of degree \(n\) with symmetric
\(n\)-tensors, using its Gaussian \(L^2\) inner product. The
operator on this chaos is \(R^{\otimes n}\), and \(D_v\) is
\(\sqrt n\) times contraction by \(v\). These identities follow
in an orthonormal eigenbasis of \(R\) from the one-dimensional
Hermite formulas. Differentiating the tensor product replaces
one of its factors by \(B\). Symmetry makes the \(n\) resulting
terms equal, and contraction in that factor gives (AM2).
The operator norm of the degree-\(n\) derivative is bounded by
\(n b^{n-1}\|B\|\) as long as \(\|R\|\le b<1\). The analogous
second bound is \(n(n-1)b^{n-2}\|B\|^2\) for \(n\ge2\).
They justify termwise differentiation and a uniform quadratic
remainder for bounded \(\sum_i\|f_i\|_2^2\).

At a maximizer, rotating the partition by \(e^{tK}\) for any
skew-symmetric \(K\) preserves all masses. Equivalently, keep
the partition fixed and conjugate the covariance by this
rotation. The first derivative of the objective is zero. By
(AM2) it is, up to the sign of the conjugation convention,
\(\operatorname{tr}(K(RM_f-M_fR))\). The commutator is
skew-symmetric and orthogonal in the trace pairing to every
skew-symmetric matrix, so it vanishes.

## Balanced translations with a matrix correlation

The compact, strictly positive operator \(T_R\) gives an ordinary
maximizer by relaxing to simplex-valued assignments with the
prescribed masses. Its chaos operators tend to zero in norm, so
the relaxed objective is weakly continuous. A fractional tuple
would admit a midpoint split by opposite balanced transfers
between two labels, contradicting strict convexity of the
quadratic objective. As in the scalar argument, the first
variation gives prices \(c_i\) such that each cell almost
everywhere wins the adjusted scores \(s_i=T_R f_i-c_i\).

Let \(v\) be orthogonal to all cell centroids, and translate
by \(f_i^t(x)=f_i(x-tv)\). The first mass changes vanish and the
mass errors are \(O(t^2)\) by differentiation of the Gaussian
density. The
[mass repair for bulk translations](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#optimization-and-volume-repair)
applies here. In detail, the essential cell closures connect
the finite label graph. Along an oriented spanning-tree edge
one can transfer mass from a bounded positive-measure patch
where the two adjusted scores differ by at most any prescribed
\(\epsilon>0\), by continuity at a common point of the closures.
The same patches translated by \(tv\) permit exact repairs of
total mass \(O(t^2)\), with a coefficient independent of
\(\epsilon\). The translated scores converge uniformly on each
fixed compact set because the conditional Gaussian convolution
is smooth. Nonnegativity of the quadratic correction therefore
bounds the loss in the Lagrangian by
\(C\epsilon t^2+o_\epsilon(t^2)\). Since the repaired partition
has the prescribed masses, this shows that the second derivative
of \(\mathcal L(f^t)=S_R(f^t)-2\sum_i c_i\int f_i^t\,d\gamma_d\)
at zero is nonpositive, by letting \(\epsilon\downarrow0\).

The derivative can be taken for arbitrary measurable cells.
Factor \(T_R=T_{R^{1/2}}^2\) and use the translated conditional
Gaussian formula
\[
 T_{R^{1/2}}f_i^t(x)=(T_{R^{1/2}}f_i)(x-tR^{-1/2}v).
\]
The smoothed function has bounded derivatives of each fixed
order, so the translated expression is twice differentiable
in Gaussian \(L^2\). In the resulting pairings the first and
second derivatives of the raw translate are \(-D_v f_i\) and
\(D_v^2 f_i\) in distributions. At the optimum
\(\sum_i s_i D_v f_i=0\) by differentiating the identity
\(\max_i s_i=\sum_i f_i s_i\). At differentiability points of
the maximum, each active score has the same gradient as the
maximum. Gaussian integration by parts now gives
\[
\begin{split}
 \tfrac12\frac{d^2}{dt^2}\mathcal L(f^t)\bigg|_{0}
  &= M_f(v,v)+\sum_i\langle D_v^2 f_i,s_i\rangle_{\gamma_d} \\
  &= M_f(v,v)-\sum_i\langle D_v f_i,D_v T_R f_i\rangle_{\gamma_d} \\
  &= v^{\mathsf T}M_f(I-R)v .                              \tag{AM3}
\end{split}
\]
The middle equality uses the Gaussian adjoint \(D_v^*=v\cdot x-D_v\),
whose first term pairs to zero by the distributional winner
identity. Polynomial test factors can be included by using
compact cutoffs and passing to the limit, since the scores and
their derivatives are bounded and the Gaussian has integrable
polynomial tails. The last equality is the derivative commutation
\(D_v T_R=T_R D_{Rv}\), obtained by differentiating the conditional
mean, or by Hermite lowering.

The product \(M_f(I-R)\) is positive semidefinite with the same
kernel as \(M_f\), since the two factors commute and \(I-R\)
is positive definite. Thus (AM3) and the mass repair force
\(v\in\ker M_f\). Positivity of every Hermite multiplier of
\(T_R\) implies \(D_v f_i=0\) for every \(i\), first on Hermite
coefficients and then in distributions by smoothing with an
ordinary scalar Ornstein–Uhlenbeck operator and letting its
correlation tend to one. Conversely, every invariant direction
has zero centroid projection. The active subspace is therefore
exactly the centroid span, of dimension at most \(k-1\), and
is invariant under \(R\) because \(\ker M_f\) is invariant.

If an omitted eigenvalue \(r_b\) is strictly larger than an
active eigenvalue \(r_a\), rotate those two unit eigenvectors
into one another by a small angle \(t\). On the active space
the resulting covariance changes only in the first eigenvector's
diagonal entry, from \(r_a\) to \(r_a+(r_b-r_a)\sin^2t\).
The derivative in that entry is positive by (AM1)–(AM2), since
the eigenvector is active. This rotation strictly improves the
objective, proving the eigenvalue selection assertion.

## Different partitions and nonsymmetric cross-covariance

Suppose two partitions with the same prescribed masses may be
chosen independently, and the cross-covariance \(C\) has all
singular values in \((0,1)\). Write its singular value
decomposition as \(C=U R V^{\mathsf T}\) and rotate the two
Gaussian endpoints separately by \(U^{\mathsf T},V^{\mathsf T}\).
The transformed partitions, say \(f,g\), have bilinear
agreement \(\sum_i\langle f_i,T_Rg_i\rangle\). Cauchy–Schwarz
in the smoothed space bounds this by
\(\sqrt{S_R(f)S_R(g)}\), and equality at the maximum requires
both tuples to maximize \(S_R\) and their smoothed tuples to
coincide. Injectivity gives \(f=g\) up to null sets.
Consequently the bilinear optimum is obtained by a common
maximizer in these coordinates, using at most \(k-1\) coordinates
from the strongest singular values. The two original partitions
use the corresponding left and right singular spaces.

## The two-axis directional derivative at an isotropic tie

Let \(s_*(R)\) be the balanced three-cell maximum, and let
\(\mathcal M_r\) be the set of its planar maximizers at \(R=rI\),
where \(0<r<1\). For \(a\ge b\) and sufficiently small \(t>0\),
put \(R_t=\operatorname{diag}(r+ta,r+tb)\). With \(M_f\) evaluated
at \(rI\) and eigenvalues \(\mu_1(f)\ge\mu_2(f)\), one has
\[
 \lim_{t\downarrow0}\frac{s_*(R_t)-s_*(rI)}{t}
  =\max_{f\in\mathcal M_r}\bigg[
       \frac{a+b}{2}\,S'_r(f)
       +\frac{a-b}{2}\bigl(\mu_1(f)-\mu_2(f)\bigr)\bigg].   \tag{AM4}
\]
The trace of \(M_f\) is \(S'_r(f)\), and rotations align its
eigendirections with the ordered perturbation. For the lower
bound use each fixed rotated scalar maximizer. For the upper
bound take maximizers at \(R_t\). The uniform quadratic remainder
from (AM2) makes them near maximizers at \(rI\). Weak compactness
and compactness of \(T_{rI}\) give a scalar maximizing limit.
That limit is ordinary by strict positivity, and equality of
the Gaussian \(L^2\) norms then upgrades convergence to strong.
The covariance derivative is continuous in that norm by the
uniform derivative operator bound. This proves the reverse
bound and attainment of the maximum in (AM4).

The formula applies to higher ambient dimensions with these
same two strongest correlations by the reduction above. It
allows rank-one as well as rank-two perturbed maximizers.
At an isotropic tie where the simplex is among the maximizers,
its matrix \(M_f\) is scalar by cyclic rotation symmetry. A
traceless two-axis perturbation \((a,-a)\), \(a>0\), therefore
detects whether any tied shape has unequal covariance-gradient
eigenvalues. This derivative concerns the anisotropic optimal
value at a single scalar correlation.

## Logarithmic covariance comparison

Fix any prescribed positive cell masses and denote their optimum
by \(s_*(R)\). It is a symmetric convex function of the
logarithms of the covariance
eigenvalues on \((-\infty,0)^d\). In particular, if one
log-eigenvalue vector belongs to the convex hull of the
permutations of another, its optimum is at most that of the
other. Spreading a scalar covariance while fixing the
determinant gives a strict comparison. More precisely, let
\(R=\operatorname{diag}(re^{x_1},\ldots,re^{x_d})\) satisfy
\(0<R<I\) and \(\sum_j x_j=0\). If \(f\) maximizes at \(rI\) and
\(W_1(f)=\sum_i|\int y f_i(y)\,d\gamma_d(y)|^2\), then
\[
  s_*(R)-s_*(rI)\ \ge\
     rW_1(f)\left(\frac1d\sum_{j=1}^d e^{x_j}-1\right).
                                                               \tag{AM5}
\]
With at least two cells, \(W_1(f)>0\) by active dimension
reduction, so the right side is positive unless all \(x_j=0\).

For a partition \(f\), use orthonormal multivariate Hermite
coefficients and put \(a_\alpha=\sum_i|\widehat f_i(\alpha)|^2\).
At log-eigenvalues \(z_j\) its objective is
\(\sum_{\alpha} a_\alpha e^{\alpha\cdot z}\). Every coefficient is
nonnegative, so termwise convexity gives convexity of each
objective and hence of their supremum. The supremum is
symmetric because permuting coordinates preserves the
Gaussian and the prescribed masses. This proves the convex
hull comparison. For (AM5), evaluate the fixed scalar
maximizer on each permutation of the \(x_j\) and average.
The exponential factor \(e^{\alpha\cdot x}\) averages to at
least one by Jensen and the zero sum of the \(x_j\). For every
multi-index of degree one it averages to \(d^{-1}\sum_j e^{x_j}\).
The optimum at \(R\) bounds the average, proving (AM5).

In the two-axis determinant-preserving direction
\(R_t=\operatorname{diag}(re^t,re^{-t})\), (AM5) gives the
unconditional gain \(rW_1(f)(\cosh t-1)\). Thus when the
first-order test in (AM4) vanishes, determinant-preserving
anisotropy still increases the optimum at least quadratically.

## Two-axis comparison at fixed arithmetic trace

In the plane, let \(R_t=\operatorname{diag}(r+t,r-t)\) satisfy
\(0<R_t<I\). For a bounded partition \(f\), denote rotational
averaging of \(R_t\) by integrating over \(Q\in SO(2)\) with
normalized Haar measure. Then
\[
   \int_{SO(2)}S_{Q R_t Q^{\mathsf T}}(f)\,dQ
     \ \ge\ S_{rI}(f)+
          \frac{t^2}{4}\sum_i
             \bigl\|T_{\sqrt r}\Delta f_i\bigr\|_2^2.             \tag{AM6}
\]
The coefficient of \(t^2\) is positive for every nonconstant
partition. In particular, for at least two cells with positive
masses and any fixed sum of the two eigenvalues, the isotropic
covariance uniquely minimizes the optimal value.

On the complexified planar tensor space use
\(e_\pm=(e_1\pm i e_2)/\sqrt2\). A real rotation is diagonal in
this basis with opposite angular weights, and \(R_t\) has
diagonal entries \(r,r\) and off-diagonal entries \(t,t\).
In degree \(n\), the normalized symmetrized tensor with \(p\) plus
factors and \(q=n-p\) minus factors has distinct rotation weight
for each \(p\). Rotation averaging therefore keeps exactly the
diagonal entries. To stay in the same weight, a term in the
tensor product flips \(j\) factors in each direction. Its
diagonal coefficient is
\[
    \sum_{j=0}^{\min(p,q)}
         \binom pj\binom qj r^{p+q-2j}t^{2j}.
\]
The \(j=0,1\) terms give \(r^{p+q}+pq\,r^{p+q-2}t^2\) where
the second term is zero if \(pq=0\). Mixed Hermite lowering by
one plus and one minus factor has coefficient \(\sqrt{pq}\),
and \(\Delta=2D_{e_+}D_{e_-}\). This gives (AM6) after summation
over the chaos coefficients. If the quadratic coefficient
vanishes, the mixed coefficients of each indicator vanish.
Then each scalar Gaussian smoothing of that indicator is
a bounded harmonic function and hence constant by Liouville's
theorem. Convergence of smoothing to the indicator in \(L^2\)
makes the indicator constant almost everywhere.

Apply (AM6) to a partition maximizing at \(rI\). At each
rotation of \(R_t\) the optimal value is the same as at \(R_t\),
so the averaged objective is a lower bound for \(s_*(R_t)\).
This fixed-trace comparison does not require the first-order
anisotropy in (AM4) to be nonzero.

## Comparison of two arbitrary spreads at fixed trace

Fix two orthonormal coordinate directions and a positive strict
contraction on their orthogonal complement. Keep that complementary
covariance unchanged and put
\[
 R(x,y)\big|_{\mathbb R^2}=
    \begin{pmatrix}r+x&y\\y&r-x\end{pmatrix},\qquad
 x^2+y^2<\min(r,1-r)^2.
\]
Write \(\Delta_P=D_1^2+D_2^2\) for the spatial Laplacian in
the two chosen directions. For each bounded indicator tuple \(f\),
its stability is subharmonic in \((x,y)\). More precisely,
\[
 (\partial_x^2+\partial_y^2)S_{R(x,y)}(f)
        =\sum_i\|T_{R(x,y)^{1/2}}\Delta_P f_i\|_2^2 .        \tag{AM7}
\]
As in (AM1), interpret the derivatives of the indicators by Gaussian
smoothing. To see (AM7), put \(B=\operatorname{diag}(1,-1)\) and
\(C=\begin{pmatrix}0&1\\1&0\end{pmatrix}\) on the chosen plane
and zero elsewhere. Two derivatives in the degree-\(n\) chaos replace
two distinct factors of \(R^{\otimes n}\), with common factor
\(n(n-1)\). On the two affected slots, the compression of
\(B\otimes B+C\otimes C\) to symmetric two-tensors is the
rank-one operator \(\omega\omega^*\), where
\(\omega=e_1\otimes e_1+e_2\otimes e_2\). Indeed it sends each
\(e_j\otimes e_j\), \(j=1,2\), to \(\omega\), and vanishes on
the symmetrized mixed vector and on tensor slots involving the complement.
Contraction by \(\omega\), with coefficient \(\sqrt{n(n-1)}\),
is exactly \(\Delta_P\) on that chaos. The remaining \(n-2\) copies
of \(R\) give the smoothed squared norm. The derivative bounds in
(AM2) justify termwise summation on compact subsets of the disk.

For any allowed radius \(a\), let \(D\) have covariance \((r-a)I_2\)
on the chosen plane and the unchanged covariance on its complement. Set
\[
 c_f=\sum_i\|T_{D^{1/2}}\Delta_P f_i\|_2^2.
\]
On the closed disk of radius \(a\) we have \(R(x,y)\succeq D\).
Tensor powers preserve this order by telescoping. Applying it in
each chaos on the right of (AM7) shows that
\(S_{R(x,y)}(f)-c_f(x^2+y^2)/4\) is subharmonic. Also \(c_f=0\)
exactly when every \(f_i\) is independent of the two directions.
In one direction, if \(c_f=0\), the smooth bounded functions
\(T_{D^{1/2}}f_i\) are harmonic on each planar slice, since
\(\Delta_P T_{D^{1/2}}=(r-a)T_{D^{1/2}}\Delta_P\).
Liouville's theorem makes them constant on such slices, and injectivity
of the multiplier on every Hermite coefficient gives the same
invariance of the indicators. The converse follows directly.

Let \(0\le b<a<\min(r,1-r)\) and choose a maximizer \(f\) for
\(R(b,0)\). The maximum principle on the disk of radius \(a\) gives
\[
 s_*(R(a,0))-s_*(R(b,0))\ge\frac{c_f}{4}(a^2-b^2),          \tag{AM8}
\]
because each boundary covariance is conjugate to \(R(a,0)\) by a
rotation of the two-plane and its stability is at most \(s_*(R(a,0))\).
In particular every increase of the spread strictly raises the
planar optimum for at least two prescribed positive masses. In
any dimension pairwise equalization cannot increase the optimum.
Equivalently the optimum is Schur-convex in the ordinary covariance
eigenvalues, since ordinary majorization is generated by pairwise
equalizations and permutations. Inequality (AM8) is strict whenever
an optimizer at the less spread covariance depends on the two-plane.

There is a derivative version for each individual maximizer. If \(f\)
maximizes at \(R(a,0)\) with \(a>0\), then on the boundary of the
disk \(S_{R(x,y)}(f)\le S_{R(a,0)}(f)\) by rotation and optimality.
The same maximum principle now bounds \(S_{R(b,0)}(f)\) by
\(S_{R(a,0)}(f)-c_f(a^2-b^2)/4\) for every \(b<a\).
Taking \(b\uparrow a\) and using (AM2) yields
\[
 M_f(e_1,e_1)-M_f(e_2,e_2)\ge c_f a/2 .                  \tag{AM9}
\]
Thus the direction with the larger covariance eigenvalue has the
larger covariance derivative unless the maximizer ignores both
directions. These inequalities compare the optimal values at
distinct eigenvalues and the derivative of an optimizer there.

## Anisotropic three-cell winners cannot be conical

Let \(R\) be a planar positive strict contraction with distinct
eigenvalues. There is no partition into three cells of positive
Gaussian masses that is conical about any point \(p\) and almost
everywhere chooses the largest of the adjusted scores
\(T_R f_i-\beta_i\). In particular anisotropic covariance paths
cannot use a conical deformation of the simplex as a stationary
comparison family.

First, each price difference is in \((-1,1)\). Each cell wins on a
set of positive measure, and every conditional Gaussian probability
\(T_R f_i(x)\) lies strictly between zero and one. Use the winner
inequality on each of the two cells in turn. Write \(h_i(x)=f_i(p+x)\)
for the homogeneous cone indicators and \(\Sigma=I-R^2\). At the
observation \(p+tx\) the conditional expectation of the indicator is
\[
 T_R f_i(p+tx)=
    \mathbb E h_i\big(Rx+t^{-1}(Rp-p)+t^{-1}\Sigma^{1/2}Z\big).
\]
At Lebesgue points of \(h_i\) at \(Rx\), this tends to \(h_i(Rx)\)
as \(t\to\infty\), by the Gaussian approximate identity. The strict
price gaps give a unique winning label in the limit, namely the label
at \(Rx\). For almost every direction \(x\), the cone label wins
for almost every radius \(t>0\), by Fubini. Since \(R\) is invertible
and Lebesgue-point exceptions for the cone indicators are scale
invariant, almost every direction is also admissible in the limit.
Consequently \(h_i(x)=h_i(Rx)\) almost everywhere for each label.

The conical partition has only finitely many essential angular
transition points. Indeed, each adjusted pair-score contrast is a
nonzero real analytic function. Otherwise the two indicator functions
would differ by a constant, by injectivity on positive Hermite degrees.
For each contrast pick a ray from \(p\) on which its restriction is
nonzero as an analytic function. Outside a discrete set of radii on
that ray, the contrast therefore cannot vanish identically on the
corresponding circle. Choose one circle good for all pairs and for
the almost-everywhere winner identity by Fubini. On this circle each
analytic contrast has only finitely many zeros, and the winners are
constant on the remaining open arcs. The cone indicators inherit
these angular representatives.

The map on directions induced by \(R\) preserves the finite set of
essential angular transitions, by the invariance just proved and the
fact that it is a smooth bijection. Every non-axis direction has an
infinite orbit, since the ratio of its coordinates is multiplied by
the ratio of the two distinct positive eigenvalues. Hence transitions
occur only on the four eigenrays, and the partition labels the four
quadrants about \(p\). To use three labels, exactly one label occurs
in two quadrants.

Use the eigenbasis with apex \((a,b)\). Write \(P(x_1)\) and \(Q(x_2)\)
for the conditional probabilities that the corresponding sample
coordinates exceed \(a\) and \(b\), respectively. The two coordinates
of the conditional Gaussian are independent. On the line \(x_1=a\),
write \(P=P(a)\in(0,1)\), while \(Q(x_2)\) is strictly increasing.
If the repeated label occupies adjacent quadrants, rotate or reflect
the coordinate description so that it fills the western half-plane.
On the northern ray of its vertical boundary its score contrast with
the northeastern label is \(1-P-PQ(x_2)\), which is nonconstant. If
it occupies opposite quadrants, choose the northeastern and
southwestern quadrants as the two singleton labels. Their scores
are \(PQ\) and \((1-P)(1-Q)\). The repeated label has score
\(P(1-Q)+(1-P)Q\). On the northern ray the northeastern score minus
the repeated one is \((3P-1)Q-P\), whose constancy requires \(P=1/3\).
On the southern ray the repeated score minus the southwestern one is
\((2-3P)Q+2P-1\), whose constancy requires \(P=2/3\). In both cases
these are genuine interfaces, so continuity and the winner rule
require the pair contrasts to be constant at their price differences.
The contradictions prove the claim for arbitrary positive masses.

The same covariance also excludes a weaker form of asymptotic
tripod geometry. Take any regular three-sector tripod about \(p\)
with fixed distinct labels, and fix \(0<\varepsilon<\pi/6\). There
is no adjusted-score winner which agrees almost everywhere with the
tripod beyond some radius from \(p\) on each of its three open
angular sectors trimmed by \(\varepsilon\) at both edges.

To see this, let \(\mathcal R(u)=Ru/|Ru|\) on the unit circle.
If one label fills an open angular arc beyond some radius, it also
fills every compact subarc of its \(\mathcal R\)-preimage beyond a
possibly larger radius. Indeed at observations \(p+tu\) in that
subarc the conditional mean \(Rp+tRu\) goes into the original arc
at a linearly growing distance from its edges and from the bounded
region. The conditional mass of the prescribed label tends uniformly
to one. The strict price gaps force that label to win there. In
particular the open set of directions that eventually have a
neighborhood filled by this label is invariant under taking
\(\mathcal R\)-preimages.

Modulo the angle \(\pi/2\) between consecutive eigenrays, the three
tripod boundary directions are spaced by \(\pi/6\). At least one
of them is therefore at angular distance at least \(\pi/6\) from
both eigenrays of its quadrant. That quadrant contains trimmed arcs
of the two adjacent tripod labels, one near each eigenray. These are
its only two tripod labels because two tripod boundaries cannot lie
in a quadrant. Iterating \(\mathcal R\) from the trimmed arc near the
eigenray of smaller eigenvalue converges to the eigenray of larger
eigenvalue on the other side. An iterate enters the opposite trimmed
arc. Pulling that label back gives a contradiction on a neighborhood
of the starting direction.
