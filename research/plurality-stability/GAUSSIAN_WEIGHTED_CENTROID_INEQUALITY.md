# A sharp weighted centroid inequality at arbitrary masses

Let \(p\) be an interior mass vector for \(k\) Gaussian
cells. Let \(I(p)\) be the minimal Gaussian cluster perimeter,
counting each interface once. Write \(L_p\) for the
Laplacian of the facet areas of the translated regular
simplex with masses \(p\), acting on \(H=\mathbf1^\perp\).
For an ordinary partition \(A\) of these masses in any
finite dimension, put \(m_i=\int_{A_i}x\,d\gamma(x)\) and let
\(M_A\) have rows \(m_i\). Then
\[
              \operatorname{tr}(L_p^{-1} M_A M_A^{\mathsf T})
                    \le I(p).                             \tag{WC1}
\]
In dimensions \(d\ge k-1\) equality holds precisely at the
simplicial partition with those masses, up to orthogonal transformations,
consistent relabeling and null sets. In lower dimension
the inequality is strict. At equal masses \(L_p\) is a
scalar on \(H\), so (WC1) is the ordinary balanced centroid
comparison. For unequal masses the Laplacian weight varies
with \(p\).

## The balanced dual with a perimeter weight

For a centered covariance \(Q\succeq0\) on \(H\), let
\(\mathcal C_p(Q)\) be the optimal assignment support
\[
 \mathcal C_p(Q)=\inf_{\beta\in H}
      \left\{p\cdot\beta+\mathbb E\max_i(Y_i-\beta_i)\right\},
          \qquad Y\sim N(0,Q).
\]
For positive definite \(Q\) on \(H\), the affine assignment
uses all pair facets with positive Gaussian area. Indeed
one can choose any two coordinates as joint strict leaders
among the shifted scores, and the unshifted \(Y\) has
positive density on \(H\). The minimizing prices exist by
coercivity: the expression is at least
\(p\cdot\beta-\min_i\beta_i\), which tends to infinity on
\(H\) since every \(p_i>0\). The positive facet Laplacian
is the price Hessian, giving a unique smooth price in the
centered gauge. The winner assignment has masses exactly
\(p\).

Choose a centered site matrix \(B\) with \(BB^{\mathsf T}=Q\).
Write \(A_{ij}(Q)\) for the Gaussian area between cells
\(i,j\) of this maximizing assignment, and
\(\ell_{ij}(Q)=|B_i-B_j|\) for the corresponding site
distance. The Gaussian divergence theorem gives the
winning moment rows
\(\sum_j A_{ij}(Q)(B_i-B_j)/\ell_{ij}(Q)\). Consequently
\[
 \mathcal C_p(Q)=\sum_{i<j}A_{ij}(Q)\ell_{ij}(Q),\qquad
 D\mathcal C_p(Q)[D]=\frac12\operatorname{tr}(\mathcal L_Q D),
                                                               \tag{WC2}
\]
where \(\mathcal L_Q\) is the Laplacian with edge weights
\(A_{ij}(Q)/\ell_{ij}(Q)\). The derivative follows by
differentiating the score entries using the maximizing
assignment and its balanced masses. In site coordinates
the derivative in \(\dot B\) is
\(\operatorname{tr}(\dot B^{\mathsf T}\mathcal L_Q B)\);
the covariance variation gives (WC2). This is the same
covariance derivative used in the equal-mass comparison.
[Covariance envelope identities and their geometric input](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#external-endpoint-and-application-audit)

Use \(P\) for orthogonal projection onto \(H\), and let
\(g(s)=\mathcal C_p(Q+sP)^2\). For \(s>0\), (WC2),
Cauchy–Schwarz and the Gaussian multi-bubble theorem yield
\[
 g'(s)=2\left(\sum A_{ij}\ell_{ij}\right)
          \left(\sum A_{ij}/\ell_{ij}\right)
        \ \ge\ 2\left(\sum A_{ij}\right)^2
        \ \ge\ 2I(p)^2 .
\]
Every full-rank assignment here is in dimension \(k-1\),
inside the range of the perimeter theorem. At \(Q=P\)
the site distances are all \(\sqrt2\), so
\(\mathcal C_p(P)=\sqrt2 I(p)\) and
\(\mathcal L_P=L_p/\sqrt2\). Homogeneity and
differentiability at \(P\) give
\[
 \lim_{s\to\infty}\bigl(g(s)-2I(p)^2s\bigr)
                 =I(p)\operatorname{tr}(L_p Q).
\]
The expression on the left is nondecreasing. Continuity
at \(s=0\), using a coupling of the Gaussian site vectors
by their matrix square roots, gives
\[
                \mathcal C_p(Q)^2
                     \le I(p)\operatorname{tr}(L_p Q).
                                                               \tag{WC3}
\]
If \(Q\) is not proportional to \(P\), the site distances
of \(Q+sP\) are not all equal, since the squared distances
determine the centered Gram matrix. Cauchy–Schwarz is then
strict for every \(s>0\), because all the facet areas
are positive. Thus (WC3) is strict off the scalar ray.
[Milman and Neeman, Gaussian multi-bubble perimeter and equality](https://annals.math.princeton.edu/2022/195-1/p02)

The weighted norm in (WC1) is the square of the dual
support
\[
 \sup_{\operatorname{tr}(L_pBB^{\mathsf T})\le1}
                  \operatorname{tr}(M_A^{\mathsf T} B).
\]
For each \(B\) the assignment support for its Gaussian
scores bounds the inner product by
\(\mathcal C_p(BB^{\mathsf T})\). Equation (WC3) proves
(WC1). Equality at the norming nonzero \(B\) forces
\(BB^{\mathsf T}\) proportional to \(P\) and the ordinary
winner assignment at those sites. This gives the stated
rigidity. Conversely at translated simplex sites with
Gram \(P\) the moment matrix equals \(L_p B/\sqrt2\),
so the left side of (WC1) equals
\(\operatorname{tr}(L_p)/2=I(p)\).

## Posterior gradients in the perimeter metric

On the mass simplex, the dual at \(Q=P\) is
\(\mathcal C_p(P)=\sqrt2 I(p)\). Differentiating with
respect to \(p\) gives the centered optimal price \(\beta\).
The price Hessian is \(L_p/\sqrt2\), so
\(-D^2 I(p)=L_p^{-1}\) on \(H\).
For \(p_\rho(x)=(T_\rho1_{A_i}(x))_i\), differentiating
the conditional Gaussian gives
\[
 \nabla p_{\rho,i}(x)=\frac{\rho}{\sqrt{1-\rho^2}}\,
   \int_{(A_i-\rho x)/\sqrt{1-\rho^2}} z\,d\gamma(z).
\]
Applying (WC1) at the posterior masses gives the sharp
pointwise comparison
\[
 \operatorname{tr}
 \bigl[(D p_\rho)^{\mathsf T}(-D^2 I(p_\rho))D p_\rho\bigr]
        \ \le\ \frac{\rho^2}{1-\rho^2}\,I(p_\rho).            \tag{WC4}
\]
For a translated regular simplex this is equality at
every observation point. For a nonsimplicial partition
with positive masses it is strict at every observation
point, because each conditional partition is a similarity
of the original one. The exact metric in (WC4) is tied
to the Gaussian perimeter profile. The ordinary
noise-stability derivative instead uses the Euclidean
quadratic form on label probabilities. A comparison
between those two objectives is still needed.

## An exact posterior perimeter deficit

Let \(p_t=p_{e^{-t}}\) and let \(\bar p=\int p_t\,d\gamma\)
be the initial balanced or unbalanced mass vector. Define
\[
 J(t)=\int I(p_t(x))\,d\gamma(x),\qquad
 D_t(x)=I(p_t(x))-\operatorname{tr}(L_{p_t(x)}^{-1}
                       M_t(x)M_t(x)^{\mathsf T}),
\]
where \(M_t(x)\) is the conditional centroid matrix in (WC4)
before multiplication by \(e^{-t}/\sqrt{1-e^{-2t}}\).
Then
\[
 \frac{J(t)}{\sqrt{1-e^{-2t}}}-I(\bar p)
     =\int_t^\infty
        \frac{e^{-2s}}{(1-e^{-2s})^{3/2}}
                     \int D_s(x)\,d\gamma(x)\,ds .          \tag{WC5}
\]
In particular the translated simplex uniquely minimizes
the mean posterior perimeter at every positive time.

To justify differentiation, fix a compact interval of
positive times. Each \(p_{t,i}(x)\) is bounded below by
\(c e^{-C(1+|x|)^2}\) by restricting the conditional
source integral to a fixed compact part of \(A_i\)
of positive measure. For the centered site prices
\(\beta(p)\) at \(Q=P\), the event that \(i\) wins implies
\(Z_i-Z_j\ge\beta_i-\beta_j\). The variance is two, so a
positive price difference is at most
\(2\sqrt{\log(1/p_i)}\). Thus \(\nabla I(p_t)=\beta(p_t)/\sqrt2\)
grows at most linearly in \(|x|\) on this time interval.
The first and second spatial derivatives of \(p_t\) are
uniformly bounded there by the Gaussian kernel derivative
formulas. Since \(\partial_t p_t=(\Delta-x\cdot\nabla)p_t\),
these bounds justify differentiating \(J\) under the
integral. Gaussian integration by parts gives
\[
 J'(t)=\int\operatorname{tr}
  [(D p_t)^{\mathsf T}(-D^2I(p_t))D p_t]\,d\gamma .
\]
The right side is integrable by (WC4), since \(I\) is
bounded by the number of facets times the maximum
Gaussian area of an affine hyperplane. Cutoffs on
large balls remove the boundary term by the linear growth
bound. Substitution of (WC4) in its exact centroid form
shows that the derivative of \(J(t)/\sqrt{1-e^{-2t}}\)
is \(-e^{-2t}\int D_t\,d\gamma/(1-e^{-2t})^{3/2}\).
As \(t\to\infty\), \(p_t\to\bar p\) pointwise, and boundedness
and continuity of \(I\) at \(\bar p\) give the limit in
(WC5). Nonnegativity and strictness follow from (WC1).

## A sharp posterior entropy and centroid tradeoff

Let \(A\) be an ordinary balanced Gaussian partition, \(T\) the
regular simplex, and \(n=k-1\). Write \(\pi=\mathbf1/k\),
\(I_*=I(\pi)\), and \(P\) for orthogonal projection onto the
label plane \(H\). For \(0<s<1\) put \(p(x)=T_s1_A(x)\),
\(J_s=\mathbb E I(p(X))\), and \(\overline L=\mathbb E L_{p(X)}\).
Let \(W_1(A)=\|M_A\|_F^2\) and
\(\Delta=W_1(T)-W_1(A)\). Then

\[
 s^2 W_1(A)+\frac12\|\overline L\|_F^2\le W_1(T),
 \qquad W_1(T)=\frac{2I_*^2}{n}.
 \tag{WC6}
\]

Equality holds exactly for regular simplex partitions, with the
usual cylindrical extension and symmetries. The inequality is
valid in every finite ambient dimension, with strictness when
the dimension is below \(n\). The all-cell statement imports
the balanced first-moment theorem with the source-audit scope
given in the repository. [Balanced simplex first-moment
comparison](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#external-endpoint-and-application-audit).

For the proof, add an independent standard Gaussian \(Z\) on \(H\).
Assign \((x,z)\) the label of the translated model simplex of masses
\(p(x)\) evaluated at \(z\). This is an ordinary balanced partition
with conditional masses \(p(x)\), using the continuous model price
map in the interior mass simplex. Its first-moment matrix has an
\(x\) block \(sM_A\), by self-adjointness and \(T_s x=sx\).
The \(z\) block is \(\overline L/\sqrt2\), by the fixed-fiber
Gaussian flux identity. Apply the balanced first-moment comparison
to the lifted partition. For model masses \(\pi\), the facet
Laplacian is \((2I_*/n)P\) by symmetry and trace, which gives the
stated value of \(W_1(T)\).

If \(A\) is a regular simplex, the lift is the same winning rule
on \(sO x+\sqrt{1-s^2}z\), with \(O\) a coisometry onto \(H\).
Conversely, equality forces the lift to be a centered regular
simplex. Each of its pair facets in a fixed \(z\) fiber contains
an open portion of the model pair hyperplane. Thus the label
components of the simplex site differences are positive multiples
of \(e_i-e_j\). Consistency around label triangles gives a common
multiple, and the equal site distances give a scalar row Gram
matrix for the \(x\) parts of the sites. The lifted price map
is therefore linear and, unless constant, a scalar multiple
of a coisometry. Constant prices would make \(p\) constant and
contradict injectivity of \(T_s\). The nonconstant case makes
\(p=T_r1_{T'}\) for some centered simplex \(T'\) and \(0<r<1\).
The semigroup law and injectivity force \(r=s\) and \(A=T'\):
if \(r\ne s\), one of the two ordinary partitions would be
a strictly smoothed indicator tuple. This is the same geometric
identification as for the entropy equality case.
[Rigidity of the simplicial lift and model price
map](GAUSSIAN_SCORE_SUPPORT_FLOW.md#sharp-comparison-along-the-indicator-flow).

Since \(\operatorname{tr}\overline L=2J_s\), (WC6) gives

\[
 \frac{2}{n}\big[J_s^2-(1-s^2)I_*^2\big]
 +\frac12\left\|\overline L-\frac{2J_s}{n}P\right\|_F^2
 \le s^2\Delta .
 \tag{WC7}
\]

In particular the integrated posterior metric deficit in (WC5)
has a dimension-independent upper bound from the ordinary first
moments. With \(\tau=-\log s\), one obtains

\[
 0\le
 \int_\tau^\infty
       \frac{e^{-2t}\mathbb E D_t}{(1-e^{-2t})^{3/2}}\,dt
   =\frac{J_s}{\sqrt{1-s^2}}-I_*
 \le \frac{n s^2\Delta}{4I_*(1-s^2)} .
 \tag{WC8}
\]

Indeed (WC5) makes the middle expression nonnegative. Drop the
matrix term in (WC7), divide by \(1-s^2\), and factor the
difference of squares, whose other factor is at least \(2I_*\).
The bound controls an integrated perimeter-profile deficit.
The ordinary noise-stability deficit involves the Euclidean
quadratic variation of the posterior.

## The linear response to correlation in the label coordinates

Use the same lift and give the independently sampled Gaussian
coordinate blocks correlations \(\alpha\) in \(x\) and \(\beta\)
in \(z\), with no cross correlations between blocks. At
\(\beta=0\) the agreement of the lift is \(S_{s^2\alpha}(A)\).
Its derivative in \(\beta\) there is

\[
 \frac12\mathbb E\,\operatorname{tr}(\overline L(U)^2),
 \qquad
 \overline L(U)=\mathbb E[L_{p(\sqrt\alpha\,U+\sqrt{1-\alpha}\,V)}\mid U],
 \tag{WC9}
\]

where \(U,V\) are independent standard Gaussians. Indeed only
the first \(z\)-chaos contributes, with the model centroid
matrix \(L_{p(x)}/\sqrt2\). Conditioning the two correlated
\(x\) observations on their common \(U\) gives the formula.
For the original lift that assigns \(s x+\sqrt{1-s^2}\,z\)
according to \(A\), the first \(z\)-chaos coefficient instead
is \((1-s^2)S'_{s^2\alpha}(A)\), since the whole agreement then
equals \(S_{s^2\alpha+(1-s^2)\beta}(A)\). In this original
lift the auxiliary Gaussian has the ambient dimension of
\(A\); the model lift uses \(n\) label coordinates.

There is a sharp conditional comparison in the profile metric.
Let \(0\le\alpha<1\), put \(q(U)=T_{s\sqrt\alpha}1_A(U)\)
and \(\theta^2=s^2(1-\alpha)/(1-s^2\alpha)\), and write \(M(U)\)
for the centroid matrix of \(A\) after conditioning on the
shared \(U\) and standardizing its Gaussian source. Then

\[
 \frac12\operatorname{tr}\big(L_{q(U)}^{-1}\overline L(U)^2\big)
 \ \ge\ \frac{(\mathbb E[I(p(x))\mid U])^2}{I(q(U))}
 \ \ge\ (1-\theta^2)I(q(U))
 \ \ge\ (1-\theta^2)\operatorname{tr}
                    \big(L_{q(U)}^{-1}M(U)M(U)^{\mathsf T}\big).
 \tag{WC10}
\]

The first inequality is matrix Cauchy–Schwarz between
\(L_{q(U)}^{-1/2}\overline L(U)\) and \(L_{q(U)}^{1/2}\), using
the trace identities. The second is (PF4) for the standardized
conditional source with posterior multiplier \(\theta\). The
last is (WC1). The original lift has conditional first
\(z\)-chaos matrix \(\sqrt{1-\theta^2}M(U)\), by Gaussian
linear regression, so the two ends of (WC10) measure the
corresponding model and original conditional responses.
The regular simplex attains equality, since replacing its
conditional partitions by translated model simplices gives
its original lift after identifying the active \(z\) coordinates.
[Conditional sharp entropy
comparison](GAUSSIAN_SCORE_SUPPORT_FLOW.md#sharp-comparison-along-the-indicator-flow).

At \(\alpha=0\), \(q(U)=\pi\) and the metric is scalar, so
(WC10) gives ordinary derivative domination, strict for
nonsimplex \(A\). More explicitly (WC5), the Laplacian trace
identity, and Jensen's inequality give the lower bound
\((1-s^2)W_1(T)\) on (WC9) for every \(\alpha\). Hence the
model derivative strictly exceeds the original one whenever
\(s^2\alpha\le 2/3\) and
\[
 2 s^2\alpha\big(1-1/k-W_1(A)\big)<\Delta .
\]
Here \(m(s^2\alpha)^{m-1}\le2s^2\alpha\) for every \(m\ge2\)
bounds the higher Hermite terms of \(S'_{s^2\alpha}(A)\).
For general \(\alpha\), (WC10) has the observation-dependent
weight \(L_{q(U)}^{-1}\), whereas the agreement derivative
in (WC9) uses the Euclidean metric. The higher powers of
\(\beta\) also enter the pair objective.

## Iteration improves first moments

Write \(\mathcal L_sA\) for the ordinary balanced partition
constructed in the proof of (WC6). Its moment deficit contracts:
\[
       0\le W_1(T)-W_1(\mathcal L_sA)
               \le s^2\bigl[W_1(T)-W_1(A)\bigr].          \tag{WC11}
\]
Indeed \(\operatorname{tr}\overline L=2J_s\) and (WC5) give
\(\|\overline L\|_F^2/2\ge 2J_s^2/n\ge (1-s^2)W_1(T)\).
Use the two centroid blocks of the lift. Both bounds in (WC11)
are strict for a nonsimplex partition, by the equality
characterizations in (WC5) and (WC6).

In particular, fix any \(s\in(0,1)\) and iterate the lift, adding
fresh label coordinates on every application. The first-moment
deficits decrease at least geometrically by (WC11). The
dimension-independent first-moment stability theorem makes
the summed Gaussian label distances to the regular simplex
orbit at most a constant times \(s^m\sqrt{\Delta}\) after \(m\)
steps. One may use the first lift when the original ambient
dimension is below \(n\).
[Dimension-independent Gaussian label distance](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#gaussian-margin-and-partition-distance)

This gives a reduction for the three-cell value theorem. For
any fixed \(s\in(0,1)\), it suffices that every actual maximizer
at the tested positive correlation admit a nondecreasing
iteration step \(A\mapsto\mathcal L_sA\) for ordinary noise
stability. In fact the lifted partition would also be a
maximizer. Iteration would then give maximizing partitions
arbitrarily close to a regular simplex in Gaussian label
distance. Uniform local simplex rigidity across ambient
dimensions would force some iterate to be a simplex, which
would force the original partition to be a simplex by (WC6).
The monotonicity of full noise stability needed for this
reduction remains open.
[Uniform local simplex isolation](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md#uniform-repair-and-comparison)

## Uniqueness of the sharp local posterior metric

Fix an interior mass vector \(p\) and use the label
hyperplane \(H=\mathbf1^\perp\). Among symmetric matrices
\(Q\) on \(H\), the only ones for which the model simplicial
cluster is stationary for the quadratic centroid objective
\(\operatorname{tr}(M^{\mathsf T}QM)\) are multiples of
\(L_p^{-1}\). Here rows of \(M\) are the cell first moments
and \(L_p\) is the Gaussian facet-area Laplacian of the
model cluster. Global maximization requires the multiple to
be nonnegative.

In active simplex coordinates the model centroid matrix is
\(M=L_p/\sqrt2\), by Gaussian divergence across its facets.
The centroid objective has first response given by the affine
scores with slopes \((QM)_i\). If the difference of two such
scores were nonconstant on their shared facet, exchanging
equal small Gaussian masses from opposite sides at two
points of different score would strictly improve the objective
to first order. Thus stationarity requires
\[
                L_p Q(e_i-e_j)\ \in\ \mathbb R(e_i-e_j)
                                    \qquad (i\ne j).
\]
The simplex has every pair facet. The relation
\(e_i-e_j+e_j-e_\ell=e_i-e_\ell\) and linear independence
of the first two differences force all these eigenvalues
to coincide. The edge directions span \(H\), so
\(Q=cL_p^{-1}\). If \(c<0\), the model value is negative,
whereas a symmetric partition into Gaussian radial bands of
the prescribed masses has all centroids zero. The weighted
centroid theorem supplies the converse for \(c\ge0\).

Suppose \(k\ge3\) and \(G\) is a \(C^3\) function on the
interior probability simplex whose Hessian gives a maximizing
centroid metric for the model at every mass vector. Then
\[
                  G(p)=-cI(p)+\ell(p),\qquad c\ge0,
\]
where \(\ell\) is affine. Indeed \(D^2G=-\theta(p)D^2 I\)
by the preceding uniqueness and the model profile Hessian.
Equality of mixed third derivatives gives
\((\partial_a\theta)I_{bj}=(\partial_b\theta)I_{aj}\) in affine
coordinates on \(H\). Multiplying in the last coordinate
by the inverse Hessian, and choosing two distinct coordinate
directions, gives \(\nabla\theta(p)=0\).

This constrains a pointwise matched Bellman approach to noise
stability. Suppose a smooth function \(G(t,p)\) gives an
upper bound on the squared-posterior objective by a local
generator inequality against every admissible conditional
centroid matrix at \(p\). Suppose its generator inequality
is an equality along the regular simplex posterior process.
For each interior time the simplex posterior has full
support in the interior mass simplex, by the continuous
bijective Gaussian price map. Therefore its model centroid
matrix must maximize the Hessian-weighted quadratic objective
at every interior \(p\), and each \(G(t,\cdot)\) must belong
to the span of \(I\) and affine functions. A continuous
terminal value \(\sum_i p_i^2\) is impossible in this class.
Indeed the span is finite dimensional and pointwise closed,
and \(I\) is not quadratic. A constant Hessian for \(I\)
would give constant \(L_p\) by the profile Hessian identity
and hence constant \(I(p)=\operatorname{tr}(L_p)/2\), a
contradiction to strict concavity. This obstruction concerns
a smooth pointwise calibration using all conditional
centroid matrices, rather than a comparison that keeps
the geometric compatibility of different posterior states.
