# The Gaussian adaptive surface limit for any number of cells

At high correlation the scaled squared adaptive surface energy
converges, in the variational sense, to \(\sqrt{2/\pi}\) times
Gaussian cluster perimeter. This holds in every fixed finite
dimension and for every finite number of cells, even with scalar
directions. The sharp lower bound uses Gaussian isoperimetry in
each conditional Gaussian law. Two label rows of the adaptive
tensor retain that lower bound at each limiting interface,
including when intermediate labels occur in the approximants.
[The conditional Gaussian surface formula](GAUSSIAN_CONDITIONAL_SURFACE_CALIBRATION.md)
[A separate planar proof using Jordan paths](THREE_CELL_ADAPTIVE_SURFACE_HIGH_CORRELATION_LIMIT.md)

Let \(A=(A_1,\ldots,A_k)\) be an ordinary partition of
\(\mathbb R^d\) with \(d\ge1\) and \(k\ge2\). Write
\(P_\gamma(A)\) for total Gaussian interface measure,
counting each pair once. Put
\(\Sigma_{ij}=\partial^*A_i\cap\partial^*A_j\),
\(\nu_{ij}=\phi_d\mathcal H^{d-1}|_{\Sigma_{ij}}\), and
\(d_{ij}=e_i-e_j\). For a positive direction rank \(D\) or \(D=\infty\),
define
\[
 F_{D,r}(A)
 =\sup_{|u_{ij}|\le1}
       \left\|T_{\sqrt r}\sum_{i<j}d_{ij}\otimes u_{ij}\nu_{ij}\right\|_2 .
\]
The measurable directions take values in \(\mathbb R^D\) if \(D\)
is finite and in \(\ell^2\) otherwise. The norm is Gaussian
\(L^2\) of the label tensor, and the noise operator on a measure
uses the relative Mehler kernel. Discard partitions outside the
finite Gaussian perimeter class by assigning them infinite energy.
Within that class the displayed norm is allowed to be infinite.
Use Gaussian label distance \(\sum_i\gamma_d(A_i\triangle B_i)\).

For every \(r_n\to1\) from below and any sequence of direction ranks
\(D_n\), the functionals \(\sqrt{1-r_n^2}\,F_{D_n,r_n}^2\) have
precompact bounded-energy sequences and \(\Gamma\)-converge in
Gaussian label distance to \(\sqrt{2/\pi}\,P_\gamma\). The same
statement holds at fixed Gaussian cell masses. In particular
recovery sequences preserve the prescribed masses exactly.

## The conditional isoperimetric lower bound

Put \(s=\sqrt r\), \(\tau=\sqrt{1-r}\), and \(\sigma=\sqrt{1-r^2}\).
For one cell \(E\) with Gaussian perimeter measure
\(\nu_E=\phi_d|D1_E|\), write \(f=T_s1_E\), \(h=T_s\nu_E\).
At each \(z\), let \(E_z=(E-sz)/\tau\). Then
\[
 \gamma_d(E_z)=f(z),\qquad P_\gamma(E_z)=\tau h(z),\qquad
                          |\nabla f(z)|\le s h(z).
                                                               \tag{AH1}
\]
The surface formula follows by changing variables in the
Euclidean perimeter measure. The factor \(\tau\) is the
product of \(\tau^d\) from the Gaussian kernel and
\(\tau^{-(d-1)}\) from the perimeter change of variables.
The gradient formula follows by differentiation of the Gaussian
convolution and integration by parts against \(D1_E\).
For each fixed \(z\) the conditional density divided by \(\phi_d\)
is bounded in the source variable, so finite Gaussian perimeter
justifies these identities by truncation.

Let \(I(t)=\phi_1(\Phi^{-1}(t))\) on \((0,1)\), with value zero
at the endpoints, and \(\Psi(t)=\int_0^t I(v)\,dv\).
The scalar Gaussian isoperimetric inequality gives
\(\tau h\ge I(f)\). Combining it with (AH1) yields the
pointwise inequality
\[
              \sigma h^2\ge
              \frac{\sigma}{s\tau}\,|\nabla\Psi(f)|.
                                                               \tag{AH2}
\]
Here \(\sigma/(s\tau)=\sqrt{(1+r)/r}\to\sqrt 2\), and
\(\Psi(1)=\int_{\mathbb R}\phi_1(t)^2\,dt=1/(2\sqrt\pi)\).
Thus the lower bound is normalized to give
\(P_\gamma(E)/\sqrt{2\pi}\) for a single cell.
[Borell, Gaussian isoperimetry, Theorem 3.1](https://dornsife.usc.edu/sergey-lototsky/wp-content/uploads/sites/211/2023/06/Borell-Subspace-original.pdf)

## Compactness and localization at the limiting tags

At finite energy, write \(h_{ij}=T_s\nu_{ij}\) and
\(h_i=\sum_{j\ne i}h_{ij}\), the smoothed perimeter of cell \(i\).
Test with independent uniform scalar signs across all
\(q=\binom{k}{2}\) tags. Since \(|d_{ij}|^2=2\), the expected
squared norm is \(2\sum_{i<j}\|h_{ij}\|_2^2\). In particular
\[
 F_{D,r}(A)^2\ge 2\sum_{i<j}\|h_{ij}\|_2^2
       \ge \frac{2}{q}\left\|\sum_{i<j}h_{ij}\right\|_2^2 .
                                                               \tag{AH3}
\]
Each tag lies in \(L^2\) by testing that tag alone.
Consequently \(\sigma\|h_i\|_2^2\) is bounded for every cell
on bounded-energy sequences and \(P_\gamma(A)=O(\sigma^{-1/2})\).

Write \(f_{n,i}=T_{s_n}1_{A^{(n)}_i}\). By (AH2), the
functions \(\Psi(f_{n,i})\) have bounded local \(BV\) norm.
Also \(\int I(f_{n,i})^2\,d\gamma_d\le
\tau_n^2\|h_{n,i}\|_2^2=O(\tau_n)\). Since \(I\) is
continuous and strictly positive on every compact subset
of \((0,1)\), this implies
\(\int f_{n,i}(1-f_{n,i})\,d\gamma_d\to0\).
For every indicator \(1_E\) and \(f=T_s1_E\), the Hermite
expansion gives
\[
 \|1_E-f\|_2^2\le\|1_E\|_2^2-\|f\|_2^2
                           =\int f(1-f)\,d\gamma_d .
                                                               \tag{AH4}
\]
Indeed \((1-s^m)^2\le1-s^{2m}\) for every degree \(m\).
The bounded continuous inverse of \(\Psi\) on its range
converts local \(BV\) compactness of \(\Psi(f_{n,i})\)
to local \(L^1\) compactness of \(f_{n,i}\). Equation (AH4)
then gives local label compactness of the original
partitions. The Gaussian tail yields global label compactness.
At each limit \(A\), \(\Psi(f_{n,i})\) converges to
\(\Psi(1)1_{A_i}\). Lower semicontinuity of variation in (AH2)
gives finite Gaussian perimeter of each limit cell.

To prove the sharp cluster bound, take pairwise disjoint
compact sets \(K_{ij}\subset\Sigma_{ij}(A)\) missing
arbitrarily little limiting tagged perimeter. This is
possible by mutual singularity and inner regularity of
the tagged measures. Enclose the compact sets in disjoint
bounded open source neighborhoods \(U_{ij}\), with
open observation neighborhoods \(K_{ij}\subset V_{ij}\)
and \(\overline{V_{ij}}\subset U_{ij}\). Omit empty
compact sets. On the source region \(U_{ij}\), choose
a total order of the labels with \(i\) first and \(j\)
last. Set scalar directions to orient every label
vector from the earlier to the later label in this
order. The \(i\)-row of the smoothed tensor is the
sum of the smoothed incident tags restricted to
\(U_{ij}\), and the \(j\)-row is the negative of its
incident sum. Independent uniform signs across the
different source regions cancel the cross terms in
the expected squared norm. Therefore \(F_{D_n,r_n}^2\)
bounds the sum over these regions of the squared
norms of the two indicated smoothed cell perimeters.

For a cell \(A_i^{(n)}\), let \(h_{n,i}^U\) denote the
smoothing of its Gaussian perimeter measure restricted
to \(U=U_{ij}\). If \(z\in V=V_{ij}\), the distance from
\(s_nz\) to \(U^{\rm c}\) is bounded below for \(n\) large.
For such observations and omitted source points, the
relative kernel obeys
\[
 M_{s_n}(z,x)=\tau_n^{-d}
   \exp\!\left(-\frac{|x-s_nz|^2}{2\tau_n^2}+\frac{|x|^2}{2}\right)
                 \le C\tau_n^{-d}e^{-c/\tau_n^2}.
\]
The estimate is uniform in \(x\notin U\) by
\(|x|^2\le2|x-s_nz|^2+2|z|^2\). Equation (AH3) bounds
the omitted source mass by \(C\sigma_n^{-1/2}\).
Thus on \(V\), replacing \(h_{n,i}^U\) by \(h_{n,i}\)
changes its squared norm after multiplication by
\(\sigma_n\) by \(o(1)\), by Cauchy–Schwarz and (AH3).
The same holds for \(j\).

Apply (AH2) and weighted local \(BV\) lower semicontinuity
on each \(V_{ij}\). Both limiting cell perimeters charge
\(K_{ij}\) in the amount \(\nu_{ij}(K_{ij})\). Summing
the two selected bounds and then exhausting the tags gives
\[
 \liminf_n \sigma_n F_{D_n,r_n}(A^{(n)})^2
                       \ge \sqrt{\frac2\pi}\,P_\gamma(A).
                                                               \tag{AH5}
\]
For a liminf that is finite, the argument applies to a
subsequence with bounded scaled energy. Infinite-perimeter
limits are excluded on such subsequences by (AH2)–(AH4).

## Recovery with prescribed volumes

First consider a fixed polyhedral partition \(B\) with a
finite number of facets, which may be unbounded. Divide
them into facets with disjoint relative interiors, each
with a fixed tag. The joint Gaussian Lebesgue density is
\[
 j_r(x,y)=(2\pi)^{-d}\sigma^{-d}
       \exp\!\left(-\frac{|x|^2+|y|^2}{2(1+r)}
                   -\frac{r|x-y|^2}{2\sigma^2}\right).
                                                               \tag{AH6}
\]
For a fixed source point \(x\) and target facet \(L\),
\(\sigma\int_L j_r(x,y)\,d\mathcal H^{d-1}(y)\)
is bounded by \(C e^{-|x|^2/4}\) for \(r\ge1/2\) by
integrating the displacement Gaussian on the target
supporting hyperplane. It tends to zero when \(x\) is
at positive distance from \(L\), by taking half the
displacement exponent for that distance. At almost
every point \(x\) of \(L\), the change of variables
\(y=x+\sigma v\) in its hyperplane gives the limit
\(\phi_d(x)/\sqrt{2\pi}\). Distinct facets have
intersections of zero area, including when they lie
in the same hyperplane. Dominated convergence on
each source hyperplane therefore gives this coefficient
on the diagonal facets and zero on all cross terms.
These estimates also show finite energy at every interior
correlation for a finite polyhedral partition.

The incidence coefficients have absolute value at
most two, and every admissible direction inner product
has absolute value at most one. The sum of the diagonal
facet bounds is \(2P_\gamma(B)/\sqrt{2\pi}\) after scaling.
Together with (AH5) on the constant sequence, this proves
\[
             \sigma F_{D,r}(B)^2
                        \longrightarrow \sqrt{\frac2\pi}\,P_\gamma(B)
                                                               \tag{AH7}
\]
uniformly in the direction rank. The Mehler semigroup
identifies the quadratic form in (AH6) with the original
squared norm, so the estimate applies to all measurable
directions at once.

For a general finite-Gaussian-perimeter partition, fill
the exterior of a large cube by one label of positive
mass. At good radii this adds at most \(2d\phi_1(R)\)
perimeter on the boundary of the cube \((-R,R)^d\).
The truncated partitions converge in label measure
and perimeter by the Gaussian tail and lower
semicontinuity. Each has finite Euclidean perimeter.
Encode its labels by the orthonormal vectors \(e_i\),
so that the total variation of its derivative equals
\(\sqrt2\) times the Euclidean interface measure.
The polyhedral density theorem produces approximating
partitions with derivatives close in total variation
to the derivative after a diffeomorphism tending to
the identity in \(W^{1,\infty}\). Since \(\phi_d\) is
bounded and uniformly continuous, the area formula
for that diffeomorphism gives convergence of the weighted
perimeters. The approximants also converge in label
measure by local convergence and the Gaussian tail.
Truncate each approximant at a sufficiently large good
cube and fill the exterior by one label. The same cube
wall estimate and perimeter tail bound make the change
arbitrarily small. The polyhedral mesh of the approximation
is locally finite, so the truncations have finitely many
facets.
[Braides, Conti and Garroni, Theorem 2.2 and Lemma 2.6](https://arxiv.org/pdf/2105.14530)

For \(d\ge2\), repair the cell masses of this approximating
sequence by transfers in small disjoint cubes. Choose a
distinct density point of each donor label in the target
partition for each possible recipient, using only
positive-mass labels throughout. In a cube centered at
the chosen point, the Gaussian measure of the donor
part varies continuously with its radius. For every
fixed small positive radius it is eventually bounded
away from zero, by label convergence. Thus a cube of
radius tending to zero supplies each required transfer
exactly. Use at most one transfer for each ordered
donor–recipient pair, so the cubes can be disjoint.
Relabeling the donor part inside a cube does not increase
the total perimeter in its interior. The new interfaces
on the cube walls have perimeter at most their Gaussian
area, which tends to zero. The corrected approximants
are still polyhedral.

For \(d=1\), truncation already gives finitely many
jumps. The adjacency graph of the positive-mass labels
in a sufficiently large fixed interval is connected,
as one sees by traversing between points of all labels.
Choose fixed jumps there whose label pairs form a
spanning tree. Varying their positions gives volume
derivatives \(\phi_1(x)(e_i-e_j)\) up to orientation
and hence an invertible mass derivative on the centered
mass plane. The inverse function theorem corrects the
truncation error by displacements tending to zero,
with a vanishing change of Gaussian perimeter.

We have therefore obtained fixed polyhedral partitions
with the exact desired masses, converging to \(A\) in
label measure and in Gaussian perimeter. Apply (AH7)
to them and take the usual diagonal recovery sequence,
using its uniformity in rank. For an infinite-perimeter
target the recovery upper bound is automatic. This proves
the stated \(\Gamma\)-convergence.

## Consequences for optimal values and moving maximizers

For prescribed positive masses \(p\) in \(d\ge k-1\), write
\(P_*(p)\) for the minimum Gaussian cluster perimeter.
The Gaussian multi-bubble theorem gives the unique
simplicial orbit at those masses. Compactness, (AH5)
and recovery give
\[
 \lim_{r\uparrow1}\sqrt{1-r^2}\inf_{A:\,\gamma_d(A_i)=p_i}
                         F_{D,r}(A)^2
              =\sqrt{\frac2\pi}\,P_*(p).
                                                               \tag{AH8}
\]
Every sequence attaining this leading coefficient
converges in Gaussian label distance to that orbit.
For equal masses the orbit consists of centered
regular simplex partitions and their cylindrical extensions.
[Milman and Neeman, Theorems 1.1–1.2](https://arxiv.org/html/1805.10961v3)

For an actual noise-stability maximizer at positive
correlation with prescribed positive masses in \(d\ge k-1\),
the identity
\(F_{\infty,r}(A)^2=S'_r(A)\) also holds for arbitrary \(k\).
Indeed its cylinder in one extra coordinate is a
maximizer by dimension independence. Test a smooth compactly
supported ambient field in the original coordinates
times an odd cutoff of the extra coordinate. The resulting
normal speed is \(g(x)\) times the cutoff, so all first
volume variations vanish. The constrained second variation
and the limit to the degree-one Gaussian eigenfunction in
that coordinate bound the scalar interaction by
\(r^{-1}\int a g^2\,d\nu\), with \(a_{ij}=|\nabla
(T_r1_{A_i}-T_r1_{A_j})|\). The linear map from ambient
fields to the first volume variations is onto the centered
label plane. Indeed, coefficients annihilating this map
give a linear combination of the indicators with zero
distributional gradient, hence a constant. Positivity of
each mass makes all the coefficients equal. The implicit
function theorem therefore supplies the constraint correction.
On a fixed compact source set the interaction kernel is
bounded. To obtain any bounded measurable speed there,
approximate the ambient field \(gN_{ij}\) on the tagged
supports by bounded smooth fields in surface \(L^2\).
The surface measures of distinct tags are mutually
singular, so this prescribes each tag separately.
Test one tag at a time and exhaust the source to see that
each smoothed surface measure is square-integrable. This uses the
finite Gaussian perimeter of maximizing partitions
and the boundedness of \(a\) at fixed \(r\). For fixed
positive masses the perimeter proof uses only that no
adjusted pair score is identically zero. Choose an anchor
where the score is nonzero and apply the ray zero-count
argument at the fixed correlation.
The absolute kernel integral is then finite, so the
scalar bound extends to all bounded measurable speeds.
Sum the scalar bound over direction coordinates.
The physical normal field attains equality by the
translation identity. The Gaussian integration by
parts identity \(r S'_r(A)=\sum_{i<j}\int a_{ij}\,d\nu_{ij}\)
identifies the resulting value. In this calculation the
same regular-interface and second-variation inputs apply
for any number of labels, using their two-set versions
when \(k=2\). No full rank or nonvanishing lower bound
on the score gradients is needed.
[Heilman and Tarter, Theorem 1.9 and Lemmas 2.4, 4.8 and 5.4](https://arxiv.org/pdf/2011.05583v2)
[The finite Gaussian perimeter argument](FINITE_REFLECTION_TEST_FOR_MAXIMIZERS.md#uniform-gaussian-perimeter-and-a-hermite-tail)
[The translation and derivative identities in three cells](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Let \(V_r(p)\) be the optimal stability at the fixed
masses. Its known small-noise asymptotic is
\(V_r(p)=1-c(p)\sqrt{1-r^2}+o(\sqrt{1-r^2})\),
where \(c(p)=\sqrt{2/\pi}P_*(p)\). For every fixed
partition, \(S_r(A)\) is convex in \(r\in(0,1)\)
by its Hermite expansion with nonnegative coefficients.
For an actual optimizer \(A_r\), its derivative lies
between the backward and forward secants of \(V\)
through \(r\), because \(S_r(A_r)=V_r\) and
\(S_t(A_r)\le V_t\) at each other correlation \(t\).
Use \(t=r\pm\eta(1-r)\), take \(r\uparrow1\), and
then send \(\eta\downarrow0\). Both secants after scaling
by \(\sigma\) converge in this two-step limit to
\(c(p)\). Thus every moving family of actual optimizers
satisfies \(\sigma F_{\infty,r}(A_r)^2\to c(p)\).
[The small-noise optimal value asymptotic (SP2)](GAUSSIAN_ALL_CELL_SMALL_NOISE_PRESSURES.md#perimeter-convergence-at-moving-masses)

There is also a local saturation consequence. For an
asymptotically minimizing sequence, write
\[
 \begin{aligned}
 \sigma h_i^2-\frac{\sigma}{s\tau}|\nabla\Psi(f_i)|
 &=\frac{\sigma}{\tau}h_i(\tau h_i-I(f_i))\\
 &\quad+\frac{\sigma}{s\tau}I(f_i)(s h_i-|\nabla f_i|).
 \end{aligned}                                                   \tag{AH9}
\]
Both terms are nonnegative. Along a convergent sequence
attaining (AH5) with equality in the limit, their
integrals tend to zero for the two incident labels on
every compact subset of an open region where the
limiting perimeter has a single tag. To see this,
choose the source and observation neighborhoods
in the proof of (AH5) so that the neighborhood for
that tag also contains the chosen compact set.
This is possible because compact sets carrying
other tag perimeters can be taken outside the open
region. The same proof bounds the limsup of the
sum of the local nonnegative defects by
\(\sqrt{2/\pi}\) times the limiting perimeter omitted
by the chosen compact tag sets. Send the omitted
perimeter to zero. The two terms thus quantify local
posterior isoperimetric loss and cancellation of the
smoothed normals, even when the raw tagged measures
are too large to converge.

## Residual pair energies on single-tag neighborhoods

The localization also retains information about the
other labels. For a source neighborhood \(U=U_{ij}\)
of the limiting tag \(i,j\), orient the scalar fields
on tags incident to \(i,j\) as above. On each tag
\(a,b\) not incident to either \(i\) or \(j\), use an
independent uniform scalar sign. Write \(h_{ab}^U
=T_s(\nu_{ab}|_U)\), ignoring index order.
The \(i,j\) rows have squared norms \(\|h_i^U\|_2^2\)
and \(\|h_j^U\|_2^2\). In each other row \(a\),
the deterministic term is \(-h_{ia}^U+h_{aj}^U\).
The random terms have mean zero, and each internal
tag contributes its squared norm in both of its
rows. Therefore the expected squared norm also includes
\[
 \sum_{a\notin\{i,j\}}\|h_{ia}^U-h_{aj}^U\|_2^2
        +2\sum_{\substack{a<b\\a,b\notin\{i,j\}}}\|h_{ab}^U\|_2^2.
\]
As before, multiply the whole field in each disjoint
source region by an independent uniform sign before
summing the fields, so the expected regional energies
add. This test still uses only one direction coordinate,
uniformly in the available ranks.

Suppose a sequence converging in label distance to \(A\)
attains (AH5) with equality in the limit. Let \(C\) be
a compact subset of an open region whose limiting
perimeter has only the tag \(i,j\). Then, with full
smoothed tagged functions,
\[
 \sigma_n\int_C\left[
 \sum_{a\notin\{i,j\}}|h_{n,ia}-h_{n,aj}|^2
       +2\sum_{\substack{a<b\\a,b\notin\{i,j\}}}h_{n,ab}^2
                 \right]d\gamma_d
                         \ \longrightarrow\ 0 .
                                                               \tag{AH10}
\]
To verify this, require \(C\subset V_{ij}\) in the
compact-tag localization. The omitted source terms on
\(V_{ij}\) have exponentially small norms by the same
kernel bound used for (AH5), and (AH3) controls the total
squared norms by \(O(\sigma_n^{-1})\). This permits
replacement of all the restricted smoothed measures by
the full ones in the residual on \(C\). Keep that
nonnegative residual in the regional lower bound, and
apply (AH2) on the observation neighborhoods to the two
selected cell rows. The limsup of the scaled residual
is at most \(\sqrt{2/\pi}\) times the perimeter omitted
by the limiting compact tag sets. That omission is
arbitrarily small, proving (AH10). In particular
\(\sigma_n\int_C|h_{n,i}-h_{n,j}|^2d\gamma_d\to0\),
because their difference is the sum of the first
differences in (AH10). The same proof applies to
arbitrary rank sequences \(D_n\ge1\).

## Two-phase conditional profiles at transitions of actual maximizers

Let \(A^{(n)}\) be actual noise-stability maximizers at fixed
positive masses and \(r_n\to1\), aligned so that they converge
to a simplicial limit. Fix \(C\) compact in an open region
where the limiting perimeter has only tag \(i,j\), and fix
\(0<\epsilon<1/2\). For an observation center \(z\in C\),
write \(E_{n,a}(z)=(A^{(n)}_a-s_nz)/\tau_n\).
Let \(\mathcal H_{ij}\) be the partitions into complementary
half-spaces labeled \(i,j\), with all other cells empty.
For every \(\delta>0\),
\[
 \int_{\substack{z\in C,\ \epsilon\le f_{n,i}(z)\le1-\epsilon\\
        \operatorname{dist}_\gamma(E_n(z),\mathcal H_{ij})\ge\delta}}
                  \sigma_n h_{n,i}(z)^2\,d\gamma_d(z)
                               \ \longrightarrow\ 0 .
                                                               \tag{AH11}
\]
Here \(\operatorname{dist}_\gamma\) is the infimum of Gaussian
label distance to that class. This classifies the conditional
partitions at transition-weighted centers on limiting regular
facets. The conclusion concerns an energy-weighted family
of centers. The separate microscopic cut theorem treats
arbitrary bounded centers and in particular full-score junctions.
[Microscopic cut minimizers at full-score junctions (JD6)](GAUSSIAN_NEAR_ONE_JUNCTION_DEGREE_AND_COMPONENTS.md#microscopic-minimizing-partitions-at-the-ties)

To prove (AH11), the first term of (AH9) gives
\[
     1-\frac{I(f_{n,i})}{\tau_n h_{n,i}}\ \longrightarrow\ 0
\]
in measure weighted by \(\sigma_n h_{n,i}^2\,d\gamma_d\)
on \(C\). The same statement holds for \(j\) with its own
weight. Equation (AH10) implies \(h_{n,j}/h_{n,i}\to1\)
in measure under the \(i\)-weight. On the set where this
ratio is at least \(1/2\), the \(i\)-weight is at most four
times the \(j\)-weight, so both isoperimetric losses tend
to zero in measure under the \(i\)-weight there. Ratios
are used only where the denominator is positive.

If the integral in (AH11) stays bounded below on a
subsequence, choose centers \(z_n\) in the indicated set
where both relative losses tend to zero and
\(h_{n,j}/h_{n,i}\to1\). Along a subsequence,
\(f_{n,i}(z_n)\to t\in[\epsilon,1-\epsilon]\).
Then \(\tau_n h_{n,i}(z_n)\to I(t)\), and
\(I(f_{n,j}(z_n))\to I(t)>0\). The conditional masses of
both \(i,j\) stay away from zero. By (AH1), the Gaussian
perimeter of \(E_{n,i}(z_n)\) tends to \(I(t)\). Gaussian
perimeter compactness and equality in scalar isoperimetry
make the limit of this cell a nontrivial half-space \(H\).
Extract a further subsequence so that the other conditional
indicators converge weakly star locally to functions \(u_a\).
Their limits are nonnegative and sum to one, with
\(u_i=1_H\). The \(j\)-limit is nonzero, since local
weak convergence and the Gaussian tail preserve its
limiting conditional mass.
[Milman and Neeman, scalar equality case of Theorem 1.2](https://arxiv.org/html/1805.10961v3)

Let \(\beta_{n,a}\) be the centered score prices.
At \(x=s_nz_n+\tau_n w\), the adjusted score of label
\(a\) is exactly
\[
 \int 1_{E_{n,a}(z_n)}\!
       \left(r_n w-\tau_n s_n z_n+\sqrt{1+r_n}\,y\right)
                         \,d\gamma_d(y)-\beta_{n,a}.
                                                               \tag{AH12}
\]
The prices tend to zero by (SP1). Since \(z_n\) is bounded,
the kernels in (AH12) converge to Gaussian kernels of
mean \(w\) and covariance \(2I_d\). On every compact
observation set the scores converge uniformly to these
Gaussian convolutions of \(u_a\). Indeed local weak-star
convergence tests each fixed truncated kernel, the
Gaussian tails make the truncation error uniform, and
the bounded kernel derivatives give equicontinuity.
The winner rule passes to the resulting fractional limit.
For labels \(a,b\), the integral on any bounded set of
the \(a\)-indicator times the positive part of the score
of \(b\) minus that of \(a\) is zero. Uniform score
convergence and weak-star convergence give the same
identity with \(u_a\) and the limiting scores.
[The vanishing price estimate and score convention (SP1)](GAUSSIAN_ALL_CELL_SMALL_NOISE_PRESSURES.md)

At a point of \(\partial H\) the \(i\)-score is \(1/2\).
Every nonzero other label has a strictly positive score
there by positivity of the Gaussian kernel, while their
scores sum to \(1/2\). If two such other labels were
nonzero, every other score would be strictly below \(1/2\).
Continuity and the finite label set would make \(i\) the
unique winner in a ball crossing \(\partial H\), contrary
to \(u_i=1_H\) and the limiting winner rule. Thus \(j\)
is the only other nonzero label, and \(u_j=1_{H^{\rm c}}\).
Weak convergence to these indicators gives strong local
label convergence, hence strong Gaussian label convergence.
This contradicts the chosen positive distance from
\(\mathcal H_{ij}\) and proves (AH11).

The transition weights see the limiting facet. For a
nonnegative continuous cutoff \(\chi\) supported in
\([\epsilon,1-\epsilon]\), multiply (AH2) by
\(\chi(f_{n,i})\). The right side becomes
\(\frac{\sigma_n}{s_n\tau_n}|\nabla\Psi_\chi(f_{n,i})|\)
for \(\Psi_\chi(t)=\int_0^t \chi(v)I(v)\,dv\).
Thus for an open observation set \(O\) contained in the
single-tag region,
\[
 \liminf_n\int_O \chi(f_{n,i})\sigma_n h_{n,i}^2\,d\gamma_d
      \ \ge\ \sqrt2\,\nu_{ij}(O)\int_0^1\chi(v)I(v)\,dv .
                                                               \tag{AH13}
\]
This is weighted \(BV\) lower semicontinuity applied to
\(\Psi_\chi(f_{n,i})\), whose limit is
\(\Psi_\chi(1)1_{A_i}\).

The limiting comparisons and weighted defect convergence
leave the finite-correlation strict adaptive comparison open.
