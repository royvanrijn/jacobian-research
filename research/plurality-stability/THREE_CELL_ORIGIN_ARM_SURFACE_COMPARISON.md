# A Gaussian surface comparison for three arms from the origin

Suppose a balanced planar Gaussian noise-stability maximizer at
\(0<r<1\) has an unbounded embedded rectifiable arc of each pair
type issuing from the Gaussian origin. Then
\[
                         S'_r(A)\ge S'_r(T),              \tag{OA1}
\]
with strict inequality unless \(A\) is the centered simplex
\(T\) up to an orthogonal transformation, relabeling and null
sets. The derivatives keep the partitions fixed. The three arcs
may curve or backtrack and additional interfaces are allowed.

In particular this class cannot contain a nonsimplex tie at a
correlation approached from the right by correlations where the
simplex is optimal. If there is a strict failure of simplex
optimality, a nonsimplex maximizer at the last such transition
must lack at least one of the three origin-based arms.
On each compact interior correlation band the same exclusion
holds for sufficiently long arcs of the three types whose
starting points lie in a sufficiently small fixed disk.

## An unoriented surface bound from stability

Orient the pairs cyclically, so their label vectors \(d_i\)
satisfy \(|d_i|^2=2\) and \(d_i\cdot d_j=-1\) for distinct
pairs. Write \(\nu_i\) for the Gaussian surface measure of pair
\(i\). Define the geometric surface functional
\[
 \mathcal J_r(A)=2\sum_i\iint M_r\,d\nu_i\,d\nu_i
       +\tfrac12\sum_{i\ne j}\iint M_r\,d\nu_i\,d\nu_j .  \tag{OA2}
\]
Here and below the off-diagonal sum is ordered. Choose three
planar unit vectors \(n_i\) at \(120^\circ\). Assign the speed
\(n_i\cdot v\) on the whole pair \(i\), for one vector \(v\) of
an orthonormal basis of the plane at a time. Add the two
first extra-coordinate Hessian inequalities (G1). Their
diagonal term is \(\sum_i\int a_i\,d\nu_i=rS'_r(A)\) by (G3).
The two interaction terms add to \(r\mathcal J_r(A)\), since
\(d_i\cdot d_j\) is multiplied by \(n_i\cdot n_j\). Thus
\(S'_r(A)\ge\mathcal J_r(A)\). At the simplex these speeds
are translations when the \(n_i\) are the reference normals,
so (G2) gives equality \(\mathcal J_r(T)=S'_r(T)\).
[Translation, stability, and derivative identities (G1)–(G3)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

## A pointwise arclength comparison

The kernel of (OA2) against Euclidean arclength, apart from
a common positive factor, is
\[
 H_r(x,y)=\exp\left[-\frac{|x|^2+|y|^2-2r x\cdot y}
                                 {2(1-r^2)}\right].      \tag{OA3}
\]
This follows by multiplying the Mehler kernel by the two
Gaussian surface densities. Parametrize the three arcs from
the origin by arclength, and let \(x_i\) and \(y_i\) be the
positions on arm \(i\) at lengths \(s>0\) and \(t>0\).
An unbounded arc has infinite length, so every such parameter
is available. Put
\[
 d_i=s^2-|x_i|^2,\qquad e_i=t^2-|y_i|^2,\qquad
 \delta_i=st-x_i\cdot y_i.
\]
These quantities are nonnegative, and
\(2\delta_i\le d_i+e_i\) because
\(|x_i-y_i|\le |s-t|\). Relative to straight reference arms
at \(120^\circ\), the own-arm exponent at \((x_i,y_i)\)
therefore gains at least
\((1-r)(d_i+e_i)/(2(1-r^2))\). Each own-arm kernel is
at least its reference value.

For the cross-arm terms, let \(X=\sum_i x_i\), \(Y=\sum_i y_i\).
Assume \(t\ge s\) by symmetry and put \(a=t/s\). Since
\(2\delta_i\le d_i+e_i\),
\[
 |y_i-a x_i|^2=2a\delta_i-a^2d_i-e_i\le(a-1)e_i .
\]
With \(Z=Y-aX\), Cauchy–Schwarz and completion of the square give
\[
 X\cdot Y\ge-|Z|^2/(4a)\ge-\tfrac34\sum_i e_i .         \tag{OA4}
\]
The sum of all six ordered cross-arm exponent gains relative
to the tripod is
\[
 \frac{2\sum_i(d_i+e_i)+2r(X\cdot Y+\sum_i\delta_i)}
      {2(1-r^2)}\ge0 .                                  \tag{OA5}
\]
The equality on the left follows by summing (OA3), since the
reference cross inner product is \(-st/2\). The sign follows
from (OA4) and nonnegativity of the deficits. All six reference
cross exponents are equal, so Jensen's inequality for the
exponential bounds the cross-kernel sum below by its tripod
value at the same \((s,t)\).

Integration yields \(\mathcal J_r(A)\ge\mathcal J_r(T)\).
All extra interfaces only add nonnegative terms. Equality
forces \(d_i(s)=0\) for all \(s\) by the own-arm bound and
continuity. An arclength-parametrized arc with \(|x_i(s)|=s\)
is a straight ray by the equality case in the length bound.
For three straight rays, (OA5) then has numerator
\(2rst|\sum_i x_i(1)|^2\). Its vanishing forces the three
unit ray vectors to have pairwise angle \(120^\circ\).
There can be no additional interface of positive length.
The actual strict winner cells of a maximizer have constant
labels on the resulting three complementary sectors, giving
the centered simplex. This proves (OA1).

## A uniform strict term from the translation gap

Let \(g_v\) be the pair-constant speed \(n_i\cdot v\) used
above, and let \(P\) be the weighted orthogonal projection onto
the actual translation speeds. The extra-coordinate gap (G17)
sharpens the first comparison to
\[
 S'_r(A)-\mathcal J_r(A)\ge\frac{1-r}{r}\mathcal E_r(A),
 \qquad \mathcal E_r(A)=\sum_{v\in\mathcal B}\|g_v-Pg_v\|_a^2,
                                                               \tag{OA6}
\]
where \(\mathcal B\) is an orthonormal basis and
\(\|f\|_a^2=\int a f^2\,d\nu\). On each compact interior band
\(I\), this residual has a positive lower bound over the
nonsimplex balanced global maximizers, if there are any.
[The first extra-coordinate gap (G17)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)

Suppose the residual tends to zero along maximizers
\(A^{(m)}\) at correlations in \(I\). Pass to a strong planar
maximizer limit and a limiting correlation. The projected
translation speeds have joint values \(B_m N^{(m)}\) on the
oriented interface for a real matrix \(B_m\). The matrices
are bounded. The Gram matrix of translations is positive
definite and continuous by (G24), and the unprojected speeds
are bounded by (G3). Take \(B_m\to B\). On each compact disk
(G20) bounds the reciprocal gradient integrals uniformly,
so Cauchy–Schwarz gives
\[
 \sum_i\int_{\Sigma_i^{(m)}\cap D}
         |n_i-B N^{(m)}|\,d\mathcal H^1\longrightarrow0.  \tag{OA7}
\]
To pass from \(B_m\) to \(B\), use the uniform local surface
bound from (G20) and the upper bound on score gradients.
Gaussian surface and arclength are equivalent on the disk.
[Continuity of the translation Gram matrix (G24)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-uniform-gap-from-conical-geometry)
[Reciprocal gradient bound (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

At least two pair types retain positive local arclength in
a subsequence. Indeed take local weak limits of all three
arclength measures, using their uniform local bounds. If at
most one pair survived, a label not incident to it would have
zero distributional derivative in the limit by strong label
convergence and the interface formula for its derivative.
Its mass is \(1/3\), a contradiction. Thus for two independent
reference normals \(n_i\) there are limits of unit normals
\(N_i\) with \(BN_i=n_i\) by (OA7), and \(B\) is invertible.

Orient the pairs cyclically and choose vertex vectors \(b_i\)
with \(b_j-b_i=B^{-1}n_{(i,j)}\). They exist since the three
edge vectors sum to zero. The distributional Jacobian of
\(\sum_i b_i1_{A^{(m)}_i}\) is the sum of the tensors
\((b_j-b_i)\otimes N^{(m)}\) on the interfaces against
arclength. By (OA7) its difference from the positive symmetric
normal tensors \(N^{(m)}\otimes N^{(m)}\) tends to zero in local
total variation. The limiting field has symmetric positive
semidefinite distributional Jacobian. This argument needs no
persistence of individual arcs at dormant score ties.

Convolving the limit field with smooth mollifiers gives
curl-free fields with convex potentials and uniform Lipschitz
bounds. Normalize the potentials at the origin and take a
locally uniform limit. Its gradient equals the three vectors
\(b_i\) on their respective cells. Supporting planes with a
given gradient have the same intercept at every differentiable
point of that gradient by the two supporting-plane inequalities.
The potential thus equals the maximum of three affine functions
and its cells are convex up to null sets. The convex-maximizer
rigidity theorem gives a simplex. Uniform local growth on \(I\)
precludes nonsimplex maximizing sequences converging to that
orbit. This proves the uniform strict term.
[Convex-maximizer rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-and-conical-maximizers)
[Uniform local growth](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md)

## Arcs close to the origin and endpoint implications

Suppose the three arms start at points \(p_i\) with
\(\delta=\max_i|p_i|\le1\). Parametrize each arm relative to
its starting point, so the norm and length bounds hold for the
translated positions. In (OA3), replacing \((x,y)\) by
\((p_i+x,p_j+y)\) increases the quadratic expression in the
negative exponent by at most \(C_I(\delta^2+\delta(s+t))\).
The resulting common exponential lower bound can be used
before (OA5). Integration gives
\[
               \mathcal J_r(A)\ge S'_r(T)-C_I\delta .    \tag{OA8}
\]
Here the constant is uniform on \(I\) since multiplying the
reference integrand by \(1+s+t\) leaves a uniformly integrable
Gaussian. Equations (OA6)–(OA8) give a strictly positive uniform
lower bound on \(S'_r(A)-S'_r(T)\) for nonsimplex maximizers
with such arms whenever \(\delta\) is below a sufficiently
small fixed positive threshold depending on \(I\).

If \(A\) ties \(T\) at a correlation \(r\) approached from
the right by simplex-optimal correlations, the fixed-shape
difference \(S_s(T)-S_s(A)\) is zero at \(s=r\) and nonnegative
along that right sequence. Thus \(S'_r(A)\le S'_r(T)\).
This excludes the configurations above. If the failure set
is nonempty, its last endpoint is interior by the exact
near-one interval. A sequence of improving planar maximizers
approaching that endpoint has a nonsimplex tied limit by
compactness and local growth. The comparison applies there.
It also applies at an interior maximum of the strict stability
advantage. At such a correlation a fixed maximizing partition
has equal correlation derivative to \(T\), since their difference
attains its scalar maximum there.
[Exact near-one comparison and endpoint compactness](THREE_CELL_NEAR_ONE_EXACT.md#two-endpoints-of-any-strict-failure-interval)

Three finite arms of sufficient length give the same exclusion.
Indeed, if their lengths are at least \(L\), apply the comparison
only on \([0,L]^2\) and discard all other terms of (OA2). The
omitted reference tails are \(O_I(e^{-c_I L^2})\) for some
\(c_I>0\), by the positive definite Gaussian exponent of (OA3),
uniform on \(I\). Hence the lower bound in (OA8) loses at most
\(C_Ie^{-c_I L^2}\). Equation (OA6) absorbs both errors when the
starting points are sufficiently central and the arcs are long.
This conclusion places no restriction on the network beyond the
specified three arcs.
