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
