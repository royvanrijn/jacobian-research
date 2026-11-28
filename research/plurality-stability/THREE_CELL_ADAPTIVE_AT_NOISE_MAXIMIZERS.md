# Adaptive surface energy near a Gaussian noise maximizer

At an actual planar three-cell noise-stability maximizer,
every unrestricted adaptive optimizing direction field is
the physical normal field followed by a single isometric
embedding. The gap in the first extra-coordinate block
quantifies this assertion. Consequently along a compact
smooth shape flow, the squared adaptive energy and the
correlation derivative of noise stability agree to first
order, with an error of order \(t^2\).

In particular, if an actual maximizer is also a local
balanced minimum of the adaptive functional, the mixed
shape–correlation pairing must vanish in all balanced
directions. If the three radial second moments agree,
this forces the simplex.

## Normal frames at all actual maximizers

Fix \(0<r<1\), with the masses positive, and let \(A\) be
an actual planar Gaussian maximizer. Use its tagged
interfaces \(\Sigma_e\), Gaussian surface measures \(\nu_e\),
outward normals \(N_e\), and positive normal score
multiplier
\(a_e=|\nabla T_r(1_{A_i}-1_{A_j})|\) for \(e=(i,j)\).
Write \(K\) for the scalar tagged Mehler kernel including
the factor \(d_e\cdot d_f\), where \(d_e=e_i-e_j\).
For an admissible Hilbert-valued direction field \(u\),
denote the squared smoothed energy by \(\mathcal E(u)\).
Put \(F=F_\infty(A)\). Normal saturation gives
\[
       F^2=S'_r(A)=r^{-1}\int a\,d\nu,\qquad
       K N=(a/r)N .
\]
The kernel inequalities are in the fixed weighted space
\(L^2(a\,d\nu)\), summed over direction coordinates when
needed. The first extra-coordinate block has kernel
\(V=\{v\cdot N:v\in\mathbb R^2\}\) and gap at least
\(1-r\) on its orthogonal complement.
[Normal saturation](THREE_CELL_ADAPTIVE_NORMAL_SATURATION.md#normal-saturation-and-the-scalar-interface-form)
[The first extra-coordinate gap, G17](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)

The essential support of the normal directions contains
at least three distinct unoriented lines. To prove this,
suppose it were contained in two lines. By continuity on
each regular analytic jump arc, the normal there takes
one fixed direction, so every such arc is straight.
At a noncollinear genuine triple, the three linearized
winner wedges have positive angles strictly less than
\(\pi\), giving three unoriented normal lines. This
excludes such a triple under the assumption. At a
collinear genuine triple the two cusp walls are tangent
along one side of the common tangent line and separate
from one another there. They cannot both be straight.
This excludes the collinear case as well.

In the absence of genuine triples, a straight actual
open arc extends along its entire affine line. A finite
endpoint would lie in the closures of both incident
strict winners, with the third absent in a neighborhood.
The pair gradient stays nonzero there by (G20). Its
analytic zero curve then continues the straight line,
and the two labels on its sides show that it is still
an actual arc. Any two nonparallel complete jump lines
would intersect. For the same pair this contradicts
nonvanishing of the pair gradient. For different pairs
the intersection would lie in the closures of three
strict winners. All jump lines would consequently have
to be parallel. Their common tangential derivative of
each cell indicator vanishes distributionally, making
the partition depend on one coordinate, contrary to
full active dimension. This proves the normal-direction
claim.
[Local analytic interfaces and both triple types, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)
[Two active directions and full moment rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

There is a constant depending on \(A,r\) such that for
every admissible field in a direction Hilbert space of
dimension at least two,
\[
 \inf_{L^*L=I_2}\|u-LN\|_{L^2(a\nu)}
                  \le C_{A,r}\sqrt{F^2-\mathcal E(u)}.      \tag{DSM1}
\]
Here \(L\) ranges over linear isometric embeddings of
the active plane into the direction space. To prove this,
write the coordinatewise projection onto \(V\) as \(DN\)
for a linear map \(D\) from the active plane. The weighted
normal moment matrix is positive definite, so \(D\) is
bounded uniformly in the operator norm for \(|u|\le1\).
For \(\delta=F^2-\mathcal E(u)\), the gap gives
\[
 r\delta\ \ge\
  \int a(1-|u|^2)\,d\nu
    +(1-r)\|u-DN\|_{L^2(a\nu)}^2 .                          \tag{DSM2}
\]
With \(C=D^*D\), Cauchy–Schwarz and (DSM2) bound
\[
         \int a\,|N^\mathsf T(C-I_2)N|\,d\nu
                                           \le C_{A,r}\sqrt\delta .
\]
Indeed use the absolute difference between \(|DN|^2\)
and \(|u|^2\), then add the nonnegative term \(1-|u|^2\).
The functional \(H\mapsto\int a|N^\mathsf T HN|\,d\nu\)
is a norm on symmetric two-by-two matrices: a nonzero
homogeneous quadratic form in the plane can vanish on
at most two lines. Its equivalence to the matrix norm
gives \(\|C-I_2\|\le C_{A,r}\sqrt\delta\). For small
\(\delta\), \(L=D C^{-1/2}\) is an isometry close to
\(D\), proving (DSM1). For larger \(\delta\) compare
with any fixed isometry and use the common weighted
norm bound. In particular all attaining fields are
exactly the isometric normal frames.

## Second-order agreement along shape paths

Let \(A_t=\Phi_t(A)\), where \(\Phi_t\) is a smooth
family of diffeomorphisms with \(\Phi_0\) the identity,
supported in one compact set relative to the identity.
Transport the interface measures with their Gaussian
tangential Jacobians and pull back all directions to
the original tagged interfaces. Write \(\mathcal E_t(u)\)
for their squared energy after this transport with
direction values held fixed.

On the product of the original tagged supports the
kernel of \(\mathcal E_t-\mathcal E_0\) is bounded by
\(C|t|\) uniformly. The difference vanishes when both
points are outside the transport support. If one point
is in that compact set, the Mehler kernel and its
source derivatives there are bounded uniformly over the
other point, since their polynomial factors are
dominated by Gaussian decay in the other point.
The Gaussian tangential Jacobians have uniformly
bounded derivatives on the compact set. The same
argument gives uniform second-order expansion of the
kernel. The global reciprocal bound \(\int a^{-1}d\nu<\infty\)
then makes these kernel bounds bounds on the bilinear
forms in \(L^2(a\,d\nu)\), by Cauchy–Schwarz from this
space into \(L^1(d\nu)\).
[Global reciprocal bound, IB1](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)

For any \(u\) compare with the isometric frame \(LN\) in
(DSM1). The energy of \(LN\) after transport is
\(\mathcal E_t(N)\), by invariance under the direction
isometry. Thus
\[
 \mathcal E_t(u)-\mathcal E_t(N)
 \le -\delta+C |t|\sqrt\delta \le C' t^2 .
\]
The first inequality uses the bound on the difference
form and the uniformly bounded weighted norms of the
two fields. Taking the supremum shows
\(F_\infty(A_t)^2=\mathcal E_t(N)+O(t^2)\).

At the deformed partition the ordinary normal field
pulled back to the old measures is
\(N_t=(D\Phi_t)^{-\mathsf T}N/|(D\Phi_t)^{-\mathsf T}N|\).
It satisfies \(N_t-N=O(t)\) in the weighted space.
The first-order change in \(\mathcal E_0\) on replacing
\(N\) by \(N_t\) vanishes, because \(K N=(a/r)N\)
and the derivative of a unit field is orthogonal to
that field. Combining this with the same form bound
gives \(\mathcal E_t(N_t)=\mathcal E_t(N)+O(t^2)\).
Gaussian integration by parts gives the normal
energy identity also at the moved partition,
\(\mathcal E_t(N_t)=S'_r(A_t)\). For this identity,
differentiate the spectral stability sum, or
differentiate the Gaussian kernel and integrate by
parts twice against the cell indicators. The signed
interface measure of a label is the Gaussian
distributional gradient of its indicator, up to the
outward sign. Smoothing these measures by \(T_{\sqrt r}\)
and squaring gives the identity. The measures retain
finite energy under the compact flow because its
changed part has compact support.
[Surface derivative identity for ordinary partitions](THREE_CELL_ADAPTIVE_NORMAL_SATURATION.md#normal-saturation-and-the-scalar-interface-form)

In particular
\[
         0\le F_\infty(A_t)^2-S'_r(A_t)
                                  \le C_\Phi t^2 .          \tag{DSM3}
\]
The lower bound uses the deformed normals as a
competitor. The same upper bound holds for the
original two-direction surface functional, which
is at most the unrestricted one and also contains
the deformed normal field.

## Balanced stationarity and the dilation obstruction

For a compact smooth flow let \(h=X\cdot N\) be its
initial normal speed. Differentiating the ordinary
shape variation with respect to \(r\) gives
\[
        D_X S'_r(A)=2\int g_e h_e\,d\nu_e,\qquad
         g_e=\partial_r T_r(1_{A_i}-1_{A_j}) .
                                                               \tag{DSM4}
\]
The interchange is justified by the smooth conditional
kernel and compact transport, and volume corrections
make no contribution to first derivatives for an
initially balanced speed. Formula (DSM3) equates
this first variation to that of the squared adaptive
energy, for both ranks.

If \(A\) is also a local balanced minimum of either
adaptive functional at \(r\), then (DSM4) vanishes
for every balanced compact regular speed. Conversely
this vanishing is equivalent to shape stationarity
of the adaptive functionals at an actual maximizer.
The mixed dilation resolvent applies: write \(b=x\cdot N\),
\(d_i=\int_{A_i}(2-|x|^2)\,d\gamma_2\), and let \(w\) be
the translation speed with mass flux \(d\). The
resolvent shows that a nonsimplex simultaneous
stationary partition must have \(d\ne0\), and its
translation-orthogonal part of \(b\) must obey the
coupling bound (MR3). In particular if the conditional
radial second moments are equal, then \(d=0\) by balance.
The resolvent forces \(b=0\), so the cells are centered
cones. The conical classification of actual maximizers
makes this partition the regular simplex.

Equivalently, any nonsimplex actual maximizer with
equal conditional radial second moments has a balanced
compact smooth deformation strictly decreasing the
adaptive energy at the correlation of maximality. To
see the descent directly, use \(b\) as speed in (G19),
where radial balance gives zero flux, and approximate
by compact regular speeds with two fixed flux corrections.
The bound is strictly negative for \(b\ne0\), and
(DSM3)–(DSM4) give the decreasing deformation. This
criterion concerns actual maximizers with the stated
moment condition. A nonsimplex contact could have
unequal radial moments, or fail to minimize the
auxiliary surface functional.
[Mixed dilation resolvent and its flux correction (MR1–MR3)](THREE_CELL_MIXED_DILATION_RESOLVENT.md#resolvent-on-the-mass-plane)
[Dilation sign and the first-block gap, G19](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)
[Conical rigidity for actual maximizers](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)
