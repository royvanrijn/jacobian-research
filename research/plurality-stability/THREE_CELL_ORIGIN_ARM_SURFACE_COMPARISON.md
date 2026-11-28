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
For a common starting point the displacement loss is quadratic.
A weighted average condition on the arm positions also gives
the exclusion for some junctions at arbitrary distance from the
origin.

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

## Quadratic loss at a common displaced junction

If all three arms start at the same point \(p\), the loss in
(OA8) is quadratic in \(|p|\). On a compact interior correlation
band \(I\), there is a constant \(C_I\) such that
\[
                \mathcal J_r(A)\ge S'_r(T)-C_I|p|^2.       \tag{OA9}
\]
For arms of length at least \(L\ge2\), the bound loses only an
additional \(C_Ie^{-c_I L^2}\). The geometric estimate itself
uses neither balance nor the score equations.

Write \(x_i(s)=p+u_i(s)\) for the arclength parametrizations,
\(S(s)=\sum_i u_i(s)\), and \(D(s)=\sum_i(s^2-|u_i(s)|^2)\).
Let \(J_d(s,t)\) and \(J_o(s,t)\) be the same-arm and cross-arm
reference kernels of the centered tripod in the normalization
of \(\mathcal J_r\), including the common positive Gaussian
factor omitted in (OA3). Put
\[
 w(s)=\int_0^\infty [2J_d(s,t)+J_o(s,t)]\,dt.
\]
For \(t\ge s>0\), set \(a=t/s\). The arclength estimate in the
proof of (OA4) gives
\[
        |S(t)-aS(s)|^2\le3(a-1)D(t).
\]
Write \(X=S(s)\), \(Y=S(t)\), and \(Z=Y-aX\). Both identities
\(X\cdot Y=a|X|^2+X\cdot Z\) and
\(X\cdot Y=|Y|^2/a-Y\cdot Z/a\) give square-completion bounds.
Keeping one quarter of the first square in each and averaging gives
\[
 X\cdot Y\ge \tfrac18(a|X|^2+a^{-1}|Y|^2)-|Z|^2/(3a)
          \ge \tfrac18(a|X|^2+a^{-1}|Y|^2)-D(t).           \tag{OA10}
\]
In (OA5) the unshifted cross-arm exponent sum is therefore
at least
\[
 \frac{D(s)+D(t)}{1+r}
       +\frac{r}{8(1-r^2)}(a|S(s)|^2+a^{-1}|S(t)|^2).
\]
The sum of the three unshifted same-arm exponent gains is at
least \((D(s)+D(t))/(2(1+r))\).

Adding \(p\) to both kernel arguments changes the exponent by
\(-[|p|^2+p\cdot(u_i(s)+u_j(t))]/(1+r)\). Apply
\(e^z-1\ge z\) to each total exponent change relative to its
reference kernel, multiply by the weights in (OA2), and
integrate. The linear contribution of the shift is exactly
\[
 -\frac{3|p|^2}{1+r}\int_0^\infty w(s)\,ds
 -\frac{2}{1+r}p\cdot\int_0^\infty w(s)S(s)\,ds .
\]
The unshifted contributions dominate a fixed positive constant
times
\[
 G=\int_0^\infty w(s)D(s)\,ds+\int_1^2|S(s)|^2\,ds .
\]
For the radial term, the same-arm and cross-arm lower bounds
integrate to \((1+r)^{-1}\int w D\) after applying the two weights.
For the second term, restrict the positive square terms to
\(1\le s\le t\le2\), where \(J_o\) and \(r/(1-r^2)\) are
uniformly bounded below on \(I\). The two squares together
control all of \(\int_1^2|S(s)|^2\,ds\).

To control the remaining linear term, fix \(u\in[1,2]\). For
\(s\ge u\), the arclength estimate gives
\(|S(s)-(s/u)S(u)|^2\le3(s/u-1)D(s)\).
For \(s\le u\), multiply the same estimate with the roles
reversed by \((s/u)^2\) to bound this square by \(3D(u)\).
In both ranges,
\[
 |S(s)|\le C\big[(1+s)|S(u)|+\sqrt{D(u)}
                              +\sqrt{(1+s)D(s)}\big].
\]
Average over \(u\in[1,2]\), then integrate against \(w\).
Cauchy–Schwarz and the uniformly bounded Gaussian moments
of \(w\), together with its positive uniform lower bound on
\([1,2]\), give \(\int w|S|\le C_I\sqrt G\).
The resulting bound on the energy difference is
\(c_I G-C_I|p|\sqrt G-C_I|p|^2\), which proves (OA9)
by completing the square. Extra interfaces contribute only
nonnegative terms. For finite arms restrict all integrals to
\([0,L]\), including the marginal defining \(w\). Its lower
bound on \([1,2]\) and upper moment bounds are uniform for
\(L\ge2\). The omitted centered reference energy has the
Gaussian tail bound used in (OA8).

Consequently (OA6) excludes a right-contact nonsimplex
maximizer whenever the squared displacement of the common
junction and the finite-arm error are together smaller than
the strict first extra-coordinate margin after multiplication
by the constants above.

The radial gain and the shift terms combine to give a condition
independent of the size of the displacement. Substituting
\(u_i=x_i-p\) gives the bound
\[
\begin{split}
 \mathcal J_r(A)-S'_r(T)\ \ge\
 &\frac{1}{1+r}\int_0^\infty w(s)
                   \sum_i\bigl(s^2-|x_i(s)|^2\bigr)\,ds\\
 &+\frac{r}{16(1-r^2)}\iint_{(0,\infty)^2} J_o(s,t)
     \left(\frac{t}{s}|S(s)|^2+\frac{s}{t}|S(t)|^2\right)ds\,dt .
\end{split}                                                     \tag{OA11}
\]
Thus (OA6) excludes a nonsimplex right tie whenever the arms'
weighted second radial moment about the Gaussian origin is at
most that of the centered reference. For a simpler sufficient
condition, define the weighted average position
\[
 m_{\rm arm}
     =\frac{\int_0^\infty w(s)\sum_i x_i(s)\,ds}
                  {3\int_0^\infty w(s)\,ds}.
\]
The first term in (OA11) is the sum of
\((1+r)^{-1}\int wD\) and
\((1+r)^{-1}(3\int w)(|p|^2-2p\cdot m_{\rm arm})\).
In particular it is nonnegative if
\(p\cdot m_{\rm arm}\le |p|^2/2\), which includes a zero weighted
average. For finite arms restrict the integrals in (OA11) to
the available interval and subtract \(C_Ie^{-c_I L^2}\).
Deriving the radial moment condition from the maximizer
equations remains open.

## Optimal pair directions and strictness at a maximizer

For an actual balanced planar maximizer let
\(m_{ij}=\iint M_r\,d\nu_i\,d\nu_j\), where the indices run
over the three pair types in cyclic order. The \(m_{ij}\)
are positive by complete adjacency. For each distinct
\(i,j,k\) set
\[
 c_i=\sqrt{m_{ij}m_{ik}/m_{jk}},\qquad
 d=\left(2\max_i c_i-\sum_i c_i\right)_+ .
\]
The optimal choice of unit vectors constant on each pair
in the test (OA2) has energy
\[
 \mathcal J_r^{\rm opt}(A)
        =2\sum_i m_{ii}+\sum_i c_i^2-d^2.                 \tag{OA12}
\]
On each compact interior correlation band \(I\), there
is \(\epsilon_I>0\) such that every nonsimplex maximizer
satisfies
\[
                 S'_r(A)\ge \mathcal J_r^{\rm opt}(A)
                                      +\epsilon_I.       \tag{OA13}
\]
In particular, a nonsimplex right tie requires
\(\mathcal J_r^{\rm opt}(A)\le S'_r(T)-\epsilon_I\).

For any choice of unit vectors \(u_i\) in a plane,
the energy after summing the two tests in (G1) equals
\[
 2\sum_i m_{ii}-2\sum_{i<j}m_{ij}u_i\cdot u_j
       =2\sum_i m_{ii}+\sum_i c_i^2-\left|\sum_i c_i u_i\right|^2 .
\]
The minimum length of the resultant is \(d\). If all three
triangle inequalities hold, take the directions of a
closed triangle with side lengths \(c_i\). Otherwise the
largest vector points opposite the other two, attaining
the lower bound from the triangle inequality. This proves
(OA12). The gain over (OA2) is
\(\frac12\sum_{i<j}(c_i-c_j)^2-d^2\).
For \(d=0\) it is positive unless all \(c_i\) are equal.
For \(c_1>c_2+c_3\) the gain is
\(c_1(c_2+c_3)-3c_2c_3>0\), and the other cases follow by
permutation.

We first verify that all pair arclength measures have
uniformly positive mass on one fixed disk for maximizers
in \(I\). If this failed, use increasing disks and a
subsequence so that one pair arclength tends locally to
zero. Strong compactness gives a planar maximizing limit.
The uniform polynomial bound on interface length gives
uniform Gaussian tails for all the Mehler surface
pairings. In Euclidean arclength their kernel is a
positive Gaussian on the product plane, as in (OA3).
In particular the pairings involving the vanishing pair
tend to zero, also after multiplication by bounded normal
components.
[Polynomial interface length on compact bands](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#polynomial-growth-of-the-network)

Orient the two remaining pairs outward from their common
label and take weak limits locally of their Euclidean
arclength measures marked by the unit normals. For each
translation speed \(s_v=v\cdot N\),
\(Q^{(1)}(|s_v|)\ge Q^{(1)}(s_v)=0\). The local terms
cancel. The label inner products of the remaining pairs
are positive. Thus in the product of the limiting marked
measures the integral of the positive kernel times
\(|v\cdot N||v\cdot N'|-(v\cdot N)(v\cdot N')\) is zero.
This uses weak convergence on compact sets and the uniform
Gaussian product tails. For a countable dense family of
\(v\), this says the surviving normals all have the same
direction. Indeed two distinct unit normals admit strict
opposite signs against some vector of that family. At
least one measure survives, since otherwise every label
would have zero distributional gradient in the limit.
Passing the distributional gradient identities to the
limit now makes every indicator invariant in the
direction perpendicular to the surviving normal. This
contradicts full active dimension. The claimed pair mass
bound follows. The Gaussian kernel has a positive
minimum on the fixed disk product, so every \(m_{ij}\)
is bounded below. Its uniform upper bound follows from
the same Euclidean arclength growth estimate.
[Full active dimension and complete adjacency](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#complete-adjacency)

To prove (OA13), suppose a sequence of nonsimplex maximizers
on \(I\) has gap tending to zero. Choose optimal \(u_i\)
at every term and pass to limits in the correlations,
partitions, \(m_{ij}\), and \(u_i\). The limiting vectors
are optimal for the limiting positive weights, by
continuity of the finite-dimensional optimization. Write
\(A\) for the limiting partition and \(u_i\) for these
limiting directions. The first extra-coordinate gap
(G17) implies that each component of the chosen direction
fields approaches its weighted translation projection
in \(L^2(a\,d\nu)\). Collect the two projection vectors
in a real matrix \(B_n\), so the projected direction
field is \(B_n N^{(n)}\). The translation Gram matrix
stays positive definite in the limit by (G24) and full
active dimension, and the direction fields have uniformly
bounded weighted norm by (G3). Thus \(B_n\) is bounded.
Take a subsequence with \(B_n\to B\). The reciprocal
gradient bound (G20), as in (OA7), gives
\[
 \int_{\Sigma_i^{(n)}\cap D}
       |u_i-BN_i^{(n)}|\,d\mathcal H^1 \longrightarrow0
 \quad\text{for each pair and compact disk }D.             \tag{OA14}
\]
The limiting partition cannot be a simplex, by uniform
local growth about the simplex orbit.

If the \(u_i\) are noncollinear, their optimal resultant
vanishes and \(\sum_i c_i u_i=0\) with positive limiting
\(c_i\). Every \(u_i\) lies in the range of \(B\) by (OA14)
and the uniform pair mass bound. Hence \(B\) is invertible.
Choose label sites \(b_i\) whose cyclic differences on
the pairs \((i,j)\) are \(b_j-b_i=c_e B^{-1}u_e\).
The vanishing resultant supplies consistency around the
cycle. In the distributional Jacobian of
\(\sum_i b_i 1_{A_i^{(n)}}\), (OA14) allows replacement
of each jump tensor by \(c_e N_e^{(n)}\otimes N_e^{(n)}\)
with error tending locally to zero in total variation.
As in the proof of (OA6), smoothing the limiting field
therefore gives convex potentials. The resulting limit
cells are convex up to null sets. Convex maximizer
rigidity forces a simplex, a contradiction.

In the collinear case the \(u_i\) have mixed signs on
their common unit direction \(u\). All equal signs
would be strictly improved by the equilateral choice for
positive \(m_{ij}\). Reorient each pair so the selected
direction is \(u\) on all three. A mixed sign pattern
in the cyclic orientation gives a total order on the
labels, denoted \(L,M,H\). Relabel the pair normals so
they point from the lower to the higher label. Put
\(v=B^{\mathsf T}u\). By (OA14),
\(v\cdot N_e^{(n)}\to1\) in local arclength \(L^1\)
on every directed pair. Let \(\lambda_{LM},\lambda_{LH},
\lambda_{MH}\) be the locally weak limits of the three
unoriented arclength measures. In distributions,
\[
 D_v1_{A_L}=-(\lambda_{LM}+\lambda_{LH}),\qquad
 D_v1_{A_H}=\lambda_{LH}+\lambda_{MH}.                     \tag{OA15}
\]
In particular \(v\ne0\) and almost every \(v\)-parallel
fiber has labels in increasing order. The fiber statement
follows by slicing the distributional inequalities in
(OA15), since a one-dimensional indicator with derivative
of one sign is monotone up to null sets.

On a regular \(L,M\) arc away from genuine triple
junctions there is a neighborhood where the \(H\) label
is absent up to a null set. This follows from the
definition of a genuine triple point, including when
the unused score equals the winning values along the arc.
The second identity in (OA15) makes both measures
involving \(H\) vanish there. The approximating \(L,M\)
normals in that neighborhood converge to the level-set
normal. Indeed the pair scores converge in \(C^1\) on
compacts, after passing to a limit of the bounded
centered prices, and the pair gradient is nonzero on
the arc by (G20). Passing the full distributional label
gradients to the limit consequently identifies
\(\lambda_{LM}\) with the arclength of the limiting arc.
The first identity in (OA15) then gives
\(v\cdot N_{LM}=1\) on the arc. The same argument treats
\(M,H\) using the absence of \(L\). On any regular arc
where the third score is strictly lower, the score
convergence directly removes the other two pair measures
locally. Thus \(v\cdot N_e=1\) on all such arcs.

Any arc with this normal property is straight, since its
continuous unit normal takes values in the intersection
of the unit circle with \(v\cdot N=1\), a set of at most
two points. At a noncollinear genuine triple point all
three outgoing arcs have a strictly losing third label
on their interiors near the vertex. Their unoriented
directions cannot consist of only two lines. Three rays
on two lines would contain an opposite pair and give a
sector of opening at least \(\pi\), whereas the
noncollinear linearization has all openings below
\(\pi\). A collinear genuine triple point would have a
cusp between two analytic arcs with the same tangent.
Both arcs have a strictly losing third label on their
interiors in the cusp description. Their straightness
would make them coincide, leaving no cusp. Hence there
is no genuine triple point.
[Regular arcs and the two types of genuine triple](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Complete adjacency now gives a straight \(L,M\) line
and a straight \(M,H\) line. The regular jump arcs
cannot stop at any finite point without a genuine
triple junction, so both lines are entire. They cannot
intersect and therefore are parallel. Each is transverse
to \(v\), and on almost every \(v\)-parallel fiber they
give the \(L\)-to-\(M\) and \(M\)-to-\(H\) transitions.
Fiber monotonicity implies that the three cells are
parallel slabs, with the outer two unbounded in the
normal coordinate. This contradicts full active
dimension and completes the proof of (OA13).
[Continuation of regular jump edges](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#the-local-graph-and-its-three-faces)

## Balance alone permits violation of the arm moment condition

Fix \(r\), so the reference marginal \(w\) in (OA11) is fixed,
and let \(W=\int_0^\infty w(s)\,ds\) and
\(M=\int_0^\infty s w(s)\,ds\).
For an arbitrary vertex \(p\), the direction of \(X-p\) for
a standard Gaussian \(X\) has a strictly positive continuous
density on the circle. Choose three consecutive angular intervals
of equal mass for this density. They form an ordinary balanced
fan with connected cells and three straight arms
\(x_i(s)=p+s v_i\). Its moment margin in (OA11) is
\[
          -3 W|p|^2-2M p\cdot\sum_i v_i
                  \le -3W|p|(|p|-2M/W).
\]
So the arm moment condition fails for sufficiently distant
vertices, even with connected cells and no other interfaces.

If extra interfaces are allowed, it can fail for arbitrarily
small nonzero \(p\) as well. Start with the three regular rays
translated by \(p\), so \(\sum_i v_i=0\). The mass errors tend
to zero as \(p\to0\). In each cell choose in advance two disjoint
interior disks of positive mass with positive distance from all
the reference rays, so the disks stay inside that cell for small
\(p\). Transfer smaller concentric disks from surplus cells
to deficit labels. A surplus can be allocated across at most
two deficits, and each required disk mass is attained by
continuity as its radius varies. For small \(p\) all errors
are below the fixed disk capacities. These changes preserve
the three selected infinite rays as actual pair interfaces,
and produce exact balance. The selected arms have
\(\sum_i(x_i(s)-p)=0\), so their angular term in (OA11)
vanishes, and their radial term before the prefactor is
\(-3W|p|^2<0\). This is an obstruction to deducing the
arm moment condition from balance alone.
