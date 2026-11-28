# Reflection symmetry of one maximizing Gaussian cell forces corners on its mirror

Let \(A\) be an actual planar Gaussian three-cell noise-stability
maximizer at \(0<r<1\), with three positive prescribed masses.
Suppose one cell is invariant up to null sets under reflection
through a line through the Gaussian origin. Its essential closed
cell cannot have a one-sided \(C^1\) boundary at any point of the
mirror line. In particular no cell of an actual maximizer can be
a disk, even if the other two labels have no reflection symmetry.

If the invariant cell has a bounded convex representative, the
two endpoints of its intersection with the mirror are genuine
noncollinear triple junctions. At both points the convex cell has
a corner of aperture strictly less than \(\pi\). These aperture
bounds are uniform for bounded junction points on compact
positive-correlation bands and compact interior mass sets.

## Sorting the two labels on reflected orbits

Write \(i\) for the invariant label, \(j,k\) for the other labels,
\(R\) for reflection in the line \(L\), and \(H\) for one of the
two open half-planes. Let \(v\) be the unit normal to \(L\)
pointing into \(H\). Write \(f_l=1_{A_l}\),
\(U_l=T_rf_l\) and \(V_l=U_l-\beta_l\) for the unadjusted
scores and the actual adjusted scores, with constant volume
prices \(\beta_l\). The strict-winner sets \(O_l\) agree with
the cells almost everywhere.
[Fixed-mass score assignment](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality)

On a reflected pair \(\{x,Rx\}\), either both points are in
\(i\), both have the same other label, or there is one
occurrence each of \(j,k\). Only the last pairs may be
reordered. Reordering preserves every cell mass since a
reflection preserves Gaussian measure, and there is still
one occurrence of each of \(j,k\) on every reordered pair.

Write \(g=f_j-f_k\). Its even part stays fixed under this
reordering. The odd part on \(H\) is
\(s(x)=f_j(x)-f_j(Rx)\), since \(f_j+f_k=1-f_i\)
is even. Sorting puts \(j\) in \(H\) whenever a pair is
mismatched, replacing \(s\) by \(|s|\). The pair energy is
one half the sum of the energies of \(f_j+f_k\) and \(g\),
and the Mehler operator preserves reflection parity. Thus
the change of total stability under this sorting is
\[
 \iint_{H^2} \bigl[M_r(x,y)-M_r(x,Ry)\bigr]
          \bigl(|s(x)s(y)|-s(x)s(y)\bigr)\,d\gamma_2(x)d\gamma_2(y).
                                                               \tag{SM1}
\]
The first bracket is strictly positive for \(x,y\in H\) since
the radial factors are unchanged and
\(x\cdot y-x\cdot Ry=2(x\cdot v)(y\cdot v)>0\).
Consequently at a maximizer \(s\) has one weak sign. It is
not zero almost everywhere, since otherwise all cells are
reflection-invariant and their Gaussian first moments lie
on \(L\), in contradiction to full planar moment rank.
By exchanging \(j,k\), take \(s\ge0\) and nonzero.
[Reflecting and sorting a globally reflected label pair, RG1](REFLECTION_SYMMETRY_OF_MAXIMIZERS.md#sorting-a-reflected-pair)
[Full active moment rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

Conditional Gaussian convolution yields
\[
 U_j(x)>U_j(Rx),\qquad U_k(x)<U_k(Rx) \qquad (x\in H),
                                                               \tag{SM2}
\]
by integrating \(s\) against the positive reflected kernel
difference. The sign for \(k\) also follows from evenness of
\(U_j+U_k=1-U_i\). On the mirror itself, differentiation
of the conditional Gaussian density gives, with
\(\sigma^2=1-r^2\) and \(\phi_\sigma\) the planar Gaussian
density of covariance \(\sigma^2 I\),
\[
 \partial_v U_j(p)
    =\frac{r}{\sigma^2}\int_H(y\cdot v)s(y)
                          \phi_\sigma(y-rp)\,dy>0,\qquad
 \partial_v U_k(p)=-\partial_v U_j(p),\qquad p\in L .
                                                               \tag{SM3}
\]
Indeed \(rp\) lies on \(L\), so its conditional source density
at \(y\) equals the density at \(Ry\). The volume prices
are constant in space and do not affect these derivatives.

## Boundary labels and the mirror obstruction

Use \(C_i=\overline{O_i}\) as the essential closed cell;
equivalently it is the support of \(1_{A_i}\,dx\).
It is invariant under \(R\) even if a strict representative
has isolated dormant ties deleted. By continuity label \(i\)
is a maximal adjusted-score label throughout \(C_i\).

Take \(x\in C_i\cap H\). By reflection \(Rx\in C_i\), so
\(V_k(Rx)\le V_i(Rx)=V_i(x)\). Equation (SM2) implies
\(V_k(x)<V_i(x)\). If \(x\in\partial C_i\), label \(i\)
cannot win strictly, since then \(x\) would lie in the
interior of \(O_i\). Hence on \(\partial C_i\cap H\) one
has \(V_i=V_j>V_k\). On the opposite half-plane
the same argument gives \(V_i=V_k>V_j\).
In particular every regular open boundary arc of this
cell off the mirror has its adjacent label determined by
its side of the mirror.

Suppose \(p\in L\cap\partial C_i\) has a neighborhood in
which \(C_i\) is the closure of one side of an embedded
\(C^1\) curve. By reflection invariance its inward normal
at \(p\) is fixed by \(R\) and lies in \(L\).
The tangent is therefore parallel to \(v\), and the
curve crosses the mirror transversely. The boundary
equalities just obtained hold on the approaching arcs
on the respective sides. In particular they hold on
the regular interface subarcs given by local analytic
stratification; in this neighborhood no portion of
the actual one-sided boundary consists only of dormant
interior ties. Differentiating the boundary equalities
along the two approaching portions of the \(C^1\) curve
shows
\(\partial_v(V_i-V_j)(p)=\partial_v(V_i-V_k)(p)=0\).
This contradicts (SM3). The argument does not require
the other two cells to be exchanged globally by \(R\).
[Local analytic jump geometry](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

## Convex mirror corners and the disk exclusion

If \(C_i\) is bounded and convex, its nonempty interior
meets the mirror by convexity and symmetry, so \(C_i\cap L\)
is a nondegenerate closed segment. At each endpoint \(p\)
the boundary approaches from both open half-planes.
The side-specific boundary labels imply that both \(O_j\)
and \(O_k\) approach \(p\). Indeed at a boundary point
off the mirror the third score is strictly below the
two maximal scores, and the open complement of \(C_i\)
near that point contains strict wins of the unique
remaining adjacent label. Thus \(p\) is a genuine
triple junction.

Use coordinates normal to and along \(L\). At \(p\)
the horizontal components of the gradients of
\(V_j-V_i\) and \(V_k-V_i\) are \(d\) and \(-d\)
respectively, with \(d=\partial_v U_j(p)>0\).
If the three score gradients were collinear, the
gradient of \(V_i\) would lie between the other two
because of these horizontal components. Its tangent
winner cone would then have empty interior. This is
impossible for a convex cell containing an interior
ball: the convex hull of \(p\) and that ball gives
an open cone of directions entering \(C_i\).
Hence the triple is noncollinear. The two incident
boundary rays of \(i\) reflect to each other. Their
normals are the two contrast gradients, whose horizontal
components have equal nonzero magnitude. It follows that
the tangential components are equal. The gradients
therefore have the form \((d,a)\) and \((-d,a)\) with
\(a\ne0\); the aperture of cell \(i\) is
\[
                              2\arctan(|a|/d)<\pi .
                                                               \tag{SM4}
\]

On compact interior parameter sets this upper bound
on the aperture is uniform at mirror junctions in
any fixed disk. To see this, write \(m_l=\int x f_l\,d\gamma_2\).
The \(v\)-components of the three moments are
\(0,m_j\cdot v,-m_j\cdot v\), and
\(m_j\cdot v=\int_H (y\cdot v)s(y)\,d\gamma_2(y)\).
Uniform full moment rank therefore bounds this
nonnegative integral below. Its part outside a
sufficiently large fixed disk is uniformly small
by Gaussian tail integrability and \(0\le s\le1\).
In (SM3) the relative Mehler kernel on the remaining
fixed source disk has a uniform positive lower bound
for the bounded observations \(p\) and the correlations
in the band. Thus \(d\) stays bounded below. Gaussian
score derivative bounds keep \(|a|\) bounded above,
and (SM4) gives the uniform aperture loss.
[Uniform full moment rank on compact positive parameter sets](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md#uniform-full-moment-rank-for-three-cells)

Every disk is reflection-invariant along a line
through its center and the Gaussian origin and has
a smooth boundary crossing this line. The mirror
obstruction thus excludes any disk cell. Moreover
the minimum Gaussian symmetric difference from a
disk, over all cells of actual maximizers on the stated
compact parameter sets and over all disks, is positive.
Otherwise take a sequence converging strongly to an
actual maximizer and disks giving vanishing symmetric
difference for a fixed label. The limiting disk masses
are bounded away from zero and one. If centers and
radii stay bounded, a subsequence converges in Gaussian
label distance to a disk of positive radius, already
excluded. Otherwise along a subsequence both the center
norm \(c_n\) and the radius \(R_n\) tend to infinity;
one cannot escape while the other stays bounded,
because disk mass would tend to zero or one. The
offsets \(c_n-R_n\) stay bounded by the same mass
condition and the triangle inequality. Take their
limit \(\tau\) and the limit of the unit center
directions \(u\) along a subsequence. The disks then
converge in Gaussian measure to the half-plane
\(\{x:x\cdot u>\tau\}\), by expanding their squared
distance inequalities and using the null boundary
of that half-plane. This would give an individual
cylindrical cell of a maximizer, also excluded.
[Strong maximizing compactness](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)
[Every maximizing cell uses both active directions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions)
