# Stationary and locally maximizing balanced first moments

Write \(W_1(A)=\sum_{i=1}^k|m_i|^2\), where
\(m_i=\int_{A_i}x\,d\gamma_d\) and every cell has mass \(1/k\).
Assume \(d\ge k-1\). Two classifications follow from the sharp
first-moment and Gaussian perimeter inequalities.

- Every balanced centroidal affine winner with moment rank \(k-1\)
  is the centered regular simplex. Centroidal means that the
  sites of the affine winner rule are exactly the cell moments.
- Every balanced measurable local maximizer of \(W_1\) in Gaussian
  label distance is the centered regular simplex. In particular,
  the rank-deficient affine winners that survive first variation
  have a strictly improving variation using an unused coordinate.

The inputs are Zhang's frozen Theorem 1 and the Gaussian multi-bubble
theorem, with the balanced model perimeter normalization
\(S_*^2=(k-1)F_k/2\) in equation (9) of the frozen source.
[Retained first-moment theorem and flux formulas](sources/zhang-equal-gaussian-cells-ebfd905.md)
[Milman–Neeman, Theorems 1.1–1.2](https://arxiv.org/html/1805.10961v3)

## Facet Laplacian and full-rank rigidity

For a centroidal affine winner with distinct sites let
\(\ell_{ij}=|m_i-m_j|\), let \(A_{ij}\) be the Gaussian surface
measure of the shared facet (zero if there is no facet), and set
\[
 L=\sum_{i<j}\frac{A_{ij}}{\ell_{ij}}d_{ij}d_{ij}^{\mathsf T},
 \qquad d_{ij}=e_i-e_j .
\]
With rows \(m_i\) in \(M\), Gaussian divergence on the polyhedra,
whose inward facet normals are \((m_i-m_j)/\ell_{ij}\), gives
\(LM=M\). One may truncate by large balls and send the radius to
infinity since the spherical Gaussian boundary contribution
tends to zero. In particular
\[
 W_1(A)=\sum_{i<j}A_{ij}\ell_{ij},\qquad
 S(A)^2\le\left(\sum_{i<j}A_{ij}/\ell_{ij}\right)W_1(A)
           =\tfrac12\operatorname{tr}(L)W_1(A),                \tag{MS1}
\]
where \(S(A)=\sum_{i<j}A_{ij}\) counts each interface once.
If \(M\) has rank \(k-1\), its columns span \(\mathbf1^\perp\),
so \(L=P=I-\mathbf1\mathbf1^{\mathsf T}/k\). The endpoint bounds
\(W_1(A)\le F_k\) and \(S(A)\ge S_*\) make (MS1) an equality.
Therefore \(W_1(A)=F_k\), and the equality case of Theorem 1
identifies the simplex. In particular no lower-valued full-rank
stationary centroidal tuple is possible.

## Instability of a rank-deficient centroidal affine winner

Suppose the sites are distinct and their span has rank less
than \(k-1\). Take an unused Gaussian coordinate \(z\) orthogonal
to their span. On a facet oriented out of cell \(i\) take a
normal speed \(c_{ij}z\). It has zero first-order volume change
on each facet. The first derivative of the moment tuple is
\(\sum_{i<j}A_{ij}c_{ij}d_{ij}\) in the \(z\) coordinate,
since \(\int z^2\,d\gamma_1=1\). The pure second-order term
from the linear adjusted score along the moving facet is
\(-\ell_{ij}\int c_{ij}^2 z^2\,d\nu\), as its outward
normal derivative is \(-\ell_{ij}\). Thus half the constrained
second derivative is
\[
 \left|\sum_{i<j}A_{ij}c_{ij}d_{ij}\right|^2
       -\sum_{i<j}A_{ij}\ell_{ij}c_{ij}^2 .                   \tag{MS2}
\]
Indeed expand \(W_1(A_t)\) around \(M\). The linear part in
the changed indicators equals the integral of the adjusted
affine scores, since the masses are fixed. The change along
a moved two-label facet is the integral of its linear
normal score difference over the transfer strip. This gives
half the local term for the linear part, and the factor of two
from expanding \(W_1\) gives the term in (MS2).

To justify the speeds, first truncate them by smooth cutoffs
on each active facet away from its junctions and by even
cutoffs in the unused Gaussian coordinates. Oddness in \(z\)
keeps the linear mass flux zero. Three distinct affine scores
cannot tie along a full facet: their sites would be collinear
and the middle score would be a convex combination of the
other two everywhere, giving its cell zero mass.
For the compact regular speeds
use additional small displacements on fixed regular facets
of a spanning tree of the cell adjacency graph to make the
masses exact by the implicit function theorem. This graph is
connected since line segments between cell interiors can avoid
the lower-dimensional strata. The correction starts at second
order and contributes zero through the linear adjusted scores
on its facets. Gaussian integrability lets these truncated
coefficients in (MS2) tend to those of the full speeds. A
strictly positive limiting coefficient therefore gives an
actual exact-mass improving deformation at arbitrarily small
Gaussian label distance.

In (MS2) take \(c_{ij}=d_{ij}^{\mathsf T}v/\ell_{ij}\). The
value becomes \(v^{\mathsf T}(L^2-L)v\). If no improving
variation exists, positivity of \(L\) forces \(0\le L\le P\),
so (MS1) again gives \(W_1(A)=F_k\), which contradicts rank
deficiency. Equivalently every such rank-deficient centroidal
winner has an eigenvalue of \(L\) above one.

## From measurable local maximality to affine scores

At a local maximum, two cell moments cannot coincide. If
\(m_i=m_j\), swap two sufficiently small sets of equal mass
in those cells with different vector integrals. Such sets
can be chosen in disjoint bounded balls with separated
projections, using distinct density points in the two cells
and nonatomicity. The linear term of the \(W_1\) change is
zero and the squared change of the moment tuple is positive.

The label rule at a local maximum maximizes the linear
assignment with scores \(m_i\cdot x\) under its exact masses.
Otherwise let a better balanced assignment \(B\) have strict
linear gain. On each ordered transfer set from \(A\) to \(B\),
select the same small fraction of its measure from the upper
quantile of the integrable score gain. Nonatomicity handles
ties at the quantile. The gain on each selected set is at least
the selected fraction of the full gain on that transfer.
Matching the same fraction of every directed transfer preserves
all masses. The positive linear gain and nonnegative quadratic
term contradict local maximality as the fraction tends to zero.

The usual price minimization
\(\inf_\lambda [k^{-1}\sum_i\lambda_i+
\mathbb E\max_i(m_i\cdot X-\lambda_i)]\) attains its minimum:
normalize \(\min_i\lambda_i=0\) and lower-bound the expression
by \(\max_i\lambda_i/k\) using a zero-priced mean-zero score.
Distinct sites make ties null. The derivative in the prices
gives the prescribed masses, so its winner achieves the
optimal linear value. Equality in the pointwise score bound
identifies the labels almost everywhere. This reduces a
measurable local maximum to the centroidal affine class above.

These classifications concern \(W_1\). The full-rank argument
discharges the four-cell centroidal stationarity obligation
as a special case. The finite positive-correlation score
equation involves all Hermite levels.
