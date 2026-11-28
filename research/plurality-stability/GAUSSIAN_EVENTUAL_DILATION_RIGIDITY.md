# Discrete dilation rigidity for Gaussian winner sets

Let \(A\) be a balanced three-cell Gaussian noise-stability maximizer
in its active plane at \(0<\rho<1\). Fix \(p\in\mathbb R^2\) and
\(0<\lambda<1\). Suppose
\[
        1_A(p+\lambda x)=1_A(p+x)
        \quad\text{for almost every } |x|>R_0/\lambda       \tag{ED1}
\]
for some \(R_0>0\), with label indicators taken as a vector. Then
\(A\) agrees outside a compact set with a finite conical partition
about \(p\). If all three labels occur in the exterior, the cone is
centered at the Gaussian origin. Its sectors have either three equal
angles or six equal angles, with the labels in cyclic order. Every
bounded strict-winner component is contained in the compact core.
If all three labels occur in the exterior and at least one strict
winner set is connected, \(A\) is the centered regular simplex.
The proof combines the stationary score equation with Euclidean
first variation on the exterior.
[Outward threshold recurrence](GAUSSIAN_OUTWARD_THRESHOLD_LIMIT.md)
[Gaussian maximizer interface structure](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

## The exterior model acquires the second scale

Extend the exterior labels to a globally \(\lambda\)-invariant
partition \(B\) in coordinates relative to \(p\). For example, for
each \(x\ne0\) take the label eventually assigned to
\(p+\lambda^{-n}x\). Assumption (ED1) gives an eventually constant
sequence for almost every \(x\), and \(B(x)=A(p+x)\) for almost
every \(|x|>R_0\).

Use the adjusted scores \(V_i=T_\rho1_{A_i}-\beta_i\), with
\(\sigma^2=1-\rho^2\) and \(\eta=(\rho-1)p\). No two prices
differ by as much as one. Indeed, each label wins on a set of
positive measure and all Gaussian conditional cell probabilities
belong to \((0,1)\). At observation \(p+\lambda^{-n}x\) the
conditional source equals
\[
               p+\lambda^{-n}[\rho x+\lambda^n(\eta+\sigma Z)] ,
\]
where \(Z\) is standard planar Gaussian. Up to a probability tending
to zero, this lies outside the core. There its label agrees with
\(B(\rho x+\lambda^n(\eta+\sigma Z))\) by discrete invariance.
At almost every \(x\) the Gaussian average converges to the label
vector of \(B(\rho x)\), by the approximate identity at its
Lebesgue points. A unit-vector probability has a unique adjusted
winner at its own label since the price differences have magnitude
less than one. The observation label is \(B(x)\) for all large
\(n\), so \(B(x)=B(\rho x)\) almost everywhere. If
\(\log\lambda/\log\rho\) is irrational, the two scales already give
conicality, by local \(L^1\) continuity of the dilation action and
the density of the subgroup they generate. The proof below also
covers commensurable scales.

## The exterior network is stationary for Euclidean length

Fix a regular exterior interface point \(p+x\). Its scaled copies
\(x_n=p+\lambda^{-n}x\) have the same adjacent pair. A fixed physical
offset \(u\) at \(x_n\) has conditional source
\[
 p+\lambda^{-n}[\rho x+\lambda^n(\eta+\rho u+\sigma Z)] .
\]
Outside the core, the source labels are those of
\(B(x+\delta_n(\eta+\rho u+\sigma Z))\) with
\(\delta_n=\lambda^n/\rho\). This uses the two exact scales of \(B\).
The core error is \(O(e^{-c/\lambda^{2n}})\) for fixed \(x\ne0,u\).

At \(u=0\), use tangent and normal coordinates at \(x\), with
normal \(N\) pointing into label \(i\). The regular wall is a graph
\(t=g(s)\) with \(g(0)=g'(0)=0\), and \(i\) wins on its upper
side and \(j\) below. Put \(W=\eta+\sigma Z\) and write its independent
tangential and normal components as \(W_T,W_N\). The pair contrast
at \(x_n\) is, up to an exponentially small error,
\[
   2\mathbb E\,\overline\Phi
      \left(\frac{g(\delta_n W_T)/\delta_n-\eta\cdot N}{\sigma}\right)-1 .
\]
The graph can be extended with \(|g(s)|\le C s^2\) when estimating
this expectation. The change from that extension and any labels
beyond the local graph both lie in the Gaussian tail. The contrast
is the same price difference for all \(n\). Its limit is
\(2\Phi(\eta\cdot N/\sigma)-1\). Subtracting that limit and dividing
by \(\delta_n\) gives
\[
        -\frac{\phi_1(\eta\cdot N/\sigma)}{\sigma}\,
                     g''(0)\,\mathbb E W_T^2 .
\]
This follows by dominated convergence using Taylor's formula for
\(g\) and the bounded derivative of \(\Phi\), with domination by a
constant times \(W_T^2\). The coefficient of \(g''(0)\) is nonzero.
Hence every regular exterior arc is straight.

A genuine exterior triple cannot have collinear distinct score
gradients. In that case its local graph has distinct tangent
arcs by the collinear-junction description, whereas straight arcs
through the same point with the same tangent lie on a single line.
At any genuine exterior triple, the full translated windows about
\(x_n\) therefore converge to its three-sector tangent fan. Indeed
the regular noncollinear score linearization gives the observation
labels at \(x_n+u\) from \(B(x+\lambda^n u)\). The conditional source
formula above gives the score limit from the same tangent fan at
conditional mean \(\eta+\rho u\) and variance \(\sigma^2 I\).
Passing the winner inequalities to this limit shows that the tangent
fan satisfies the adjusted winner rule for these shifted conditional
Gaussians, with the original prices.

Each incident pair of the fan has sector angles strictly less than
\(\pi\). Its shifted pair score contrast is constant on its complete
supporting line by analyticity. Along its forward ray it tends to
the binary halfspace contrast \(2\Phi(\eta\cdot N/\sigma)-1\).
Along the reverse ray it tends to zero, since that ray is interior to
the third sector. Hence \(\eta\cdot N=0\) for each incident pair.
The normals span the plane, so \(p=0\). All pair prices at this
vertex are equal. The fan at \(u=0\) has its conditional cell
probabilities equal to its angular fractions, so these fractions
are \(1/3\). Thus every exterior triple has three angles of
\(120^\circ\). When \(p\ne0\) there are no exterior triple vertices.
[Shifted pair-line comparison for flat partitions](https://arxiv.org/pdf/1403.0885#page=11)

The locally finite exterior network \(\Gamma\) is stationary for
Euclidean arclength against every field supported outside the core.
Integrating derivatives along the straight edges gives only endpoint
terms, and the three outgoing unit tangents sum to zero at each
vertex. Local finiteness on compact sets justifies summing this
identity. The small-variance step above is the regular-interface
calculation for Gaussian threshold dynamics, with the shifted mean
retained at each scale.
[Gaussian threshold dynamics and its local interface limit](https://websites.umich.edu/~esedoglu/Papers_Preprints/esedoglu_otto.pdf)

## An expanding-period variation forces radial tangents

Let \(r=|y-p|\), \(e_r=(y-p)/r\) and \(\tau\) denote the unit
tangent along \(\Gamma\). Choose a smooth nonincreasing \(\psi(r)\)
equal to one on \(r\le R\), with \(R>R_0\), and zero beyond some
larger radius. Set \(\chi(r)=\psi(\lambda r)-\psi(r)\). This
nonnegative function is supported outside the core. The tangential
divergence of the compact radial field \(\chi(r)e_r\) is
\[
     \frac{\chi(r)}r(1-(\tau\cdot e_r)^2)
                     +\chi'(r)(\tau\cdot e_r)^2 .
\]
The integral of the second term is zero. In
\(\chi'(r)=\lambda\psi'(\lambda r)-\psi'(r)\), the two integrals
identify by centered dilation on the support of \(\psi'\) and its
outward copy. The arclength factor is exactly \(\lambda\), and the
squared tangential projection is unchanged. Stationarity and
nonnegativity now give radial tangents wherever \(\chi>0\).
Varying the cutoffs covers the exterior. Consequently there are no
exterior vertices and every active edge is part of a ray from \(p\).
The network meets an exterior circle in finitely many points by
local analytic stratification and radial transversality. Since
there are no exterior vertices, the labels are constant radially
between these rays. This proves the claimed finite exterior cone.
Inside each open exterior sector, the sector's label wins almost
everywhere. Its non-strict winner points are isolated for each
competing label by (FC4). Removing this locally finite set leaves
the open exterior sector connected and unbounded. A bounded
strict-winner component therefore cannot meet the exterior.
[Isolation of dormant interior ties, (FC4)](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)

## Centering and the exterior labels

Suppose that all three labels occur in this cone. An incident pair
contrast is constant on the full supporting line of its exterior
ray by analyticity. Tangential Gaussian deconvolution then gives a
constant normal slice contrast with the shifted Gaussian mean
\((\rho-1)p\cdot N\), as in (G15). For sufficiently large absolute
tangential coordinate, the entire normal slice agrees with the
exterior cone. Its finite step formula is analytic on each sign
of the tangential coordinate, and hence gives the same constant
on each whole side. If the normal mean is nonzero, the
linear-independence argument in (G15) leaves only rays parallel or
perpendicular to the line. A pure reverse ray gives a contrast of
zero or \(\pm1\), incompatible with the shifted halfspace forward
contrast. A different active pair on the reverse ray gives two
independent constant score contrasts on the line. This would force
the same forward and reverse score vectors, despite their different
two-label supports. The reverse ray must have the same active pair
and orientation. The exterior cone would then use just these two
labels, contradicting the hypothesis. Each active line therefore
passes through the Gaussian origin; with three exterior labels
there are nonparallel lines, so \(p=0\).
[Shifted normal slices and (G15)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

For the centered cone, the normal slice equation and Laplace
uniqueness in (G15) give reflection in each active line swapping
the two adjacent labels and fixing the third, now for the exterior
cone. Adjacent widths are consequently equal. As in the centered
conical argument, the labels repeat in the order \(1,2,3\) up to
permutation, since a run \(i,j,i\) would propagate as a binary
pattern by the reflections. If this three-label cycle repeated at
least three times then each pair would occur on at least two
nonparallel lines through zero. Its actual adjusted score contrast
vanishes on both lines, so its gradient at zero is zero. All three
scores would then tie at zero with identical gradients. The
reciprocal-gradient condition forbids regular walls accumulating
there, and the resulting dormant interior tie is excluded at zero.
Thus the three-label exterior has either three or six sectors.
[Nonzero gradients at actual walls and (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)
[Central dormant-tie exclusion, (FC5)](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)

With six sectors, each label has two opposite open exterior sectors.
A full pair-score equality line involving that label separates them
and cannot meet its strict-winner set. Each strict winner therefore
has at least two unbounded components. With three sectors, the
exterior gives all three full regular simplex pair lines. The
exact-line connectedness obstruction (EL1) then forces the simplex
if any one strict-winner set is connected.
[Exact score lines and connectedness, (EL1)](THREE_CELL_EXACT_LINE_GEOMETRY.md#each-strict-winner-cell-is-disconnected-away-from-the-simplex)

This argument leaves compact changes of the exterior cone, including
binary and constant exteriors, and does not obtain eventual discrete
scaling from ordinary maximality.

## Central components and rotation in the six-sector case

Suppose the exterior cone has six sectors. To locate the components,
number the six consecutive
open sectors so that their exterior labels are
\((1,2,3,1,2,3)\). All three adjusted scores tie at zero by the full
pair-line equalities. Two maximal scores at zero have distinct
gradients by the central separation result. The gradient contrasts
are normal to their respective pair lines. Their cyclic sum is zero,
so they have a common nonzero proportionality factor to those of a
regular tripod with the same pair lines. The central strict-winner
germ is the corresponding tripod, with one of the two orientations
according to the sign of this factor.
[Central maximal-tie separation](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#uniform-separation-at-every-maximal-tie)
[Exact germ and pair factorization](GAUSSIAN_FINITE_SIMPLEX_GERM_PROPAGATION.md#consecutive-circular-arcs-force-the-planar-simplex)

With the displayed exterior labeling, the central labels in the six
sectors are \((3,3,2,2,1,1)\) or \((2,1,1,3,3,2)\).
Indeed the active central rays alternate, and the label on the sector
between two consecutive active central rays is the common label of
their two pair types. The two full equality lines involving label \(i\)
miss its strict-winner set. Their four complement chambers include
two opposite \(60^\circ\) chambers where \(i\) occurs in the exterior
and two opposite \(120^\circ\) chambers with no exterior \(i\).
The central \(i\)-germ is in one of these \(120^\circ\) chambers,
so its component stays there and is bounded.

There is a second bounded component for each label, in its opposite
\(120^\circ\) chamber. On an active central ray of pair \(i,j\)
write \(h=V_i-V_k=V_j-V_k\), where \(k\) is the remaining label.
This function is positive on the ray near zero and near infinity.
The normal derivative of \(V_i-V_j\) has opposite signs at those
two ends, as seen from the central and exterior label lists.
If \(h\) were nonnegative on the whole ray, its isolated zeros
and the isolated zeros of the normal derivative could be
avoided by nearby points of that ray. The points so obtained would
be actual regular \(i,j\) walls accumulating at a zero of the
normal derivative. This contradicts (G20), since the tangential
derivative along the full equality line is zero. Thus \(h<0\)
on an interval of the active central ray. That interval lies in
the \(120^\circ\) chamber opposite the central chamber of label
\(k\). Its strict-winner component belongs to that chamber and is
bounded, since \(k\) has no exterior there. The three active
central rays have distinct remaining labels. Every label therefore
has one bounded component in each of its \(120^\circ\) chambers
and one unbounded component in each of its \(60^\circ\) chambers.

The six-sector case also requires at least ten genuine triple
junctions. Each of the three active central rays has at least two
distinct sign-changing zeros of \(h\), as it starts and ends
positive and is negative somewhere. On each of the other three
rays the analogous \(h\) starts negative and ends positive, giving
at least one sign change. At any such sign change, nearby points
where \(h>0\) and the pair's normal derivative is nonzero lie on
actual regular \(i,j\) walls. Strict wins of both \(i,j\) accumulate
there. Strict wins of \(k\) accumulate from where \(h<0\).
Thus each sign change is a genuine triple on a distinct nonzero
position of the ray. Including the origin gives the count.
These ten vertices are trivalent. In the collinear normal form
put \(t=V_{\rm high}-V_{\rm low}\) and write the middle score
relative to the low score as \(q(s)+\ell(s,t)t\) with
\(0<\ell<1\). On the high-low equality the spectator contrast
is \(-q(s)\). On either other pair equality its sign is the sign
of \(q(s)\), by solving \(t=q/(1-\ell)\) or
\(t=-q/\ell\), respectively. A sign change on the exact pair
line therefore makes \(q\) change sign, giving the one-sided
collinear cusp with three arcs. Noncollinear vertices also
have three arcs.
[Nonvanishing pair gradients at accumulating walls, (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

In fact there must be at least sixteen genuine triple junctions.
For a counterclockwise loop transverse to the regular interfaces,
give each cyclic transition \(i\to i+1\) weight \(+1\) and each
reverse transition weight \(-1\). Its index is the sum of these
weights divided by three. The index on a large circle equals the
sum of the indices on small counterclockwise circles around all
genuine triple vertices. To see this, remove the small disks
from the large disk. Orient each regular arc using its cyclic
label transition. Its two boundary intersections have opposite
contributions to the oriented boundary of the resulting domain,
and a closed component contributes none. The boundary index
of the domain is zero. At a trivalent junction the index on a
small circle is \(+1\) or \(-1\), according to the cyclic label
order. At a four-arc collinear junction the middle label occurs
twice, so the index is zero. The large-circle index is \(+2\)
and the origin index is \(-1\).

At each outward ray let \(\tau\) point outward, \(n=J\tau\)
with \(J\) a counterclockwise \(90^\circ\) rotation, and write
its pair as \(i,j=i+1\) in the exterior cyclic order. Set
\(b=\partial_n(V_i-V_j)\). Then \(b<0\) on the far exterior
ray and \(b>0\) near zero when the ray is active in the central
germ. At a sign change of \(h\), the active \(i,j\) arc at the
junction points toward the side with \(h>0\). Its local index
is \(+1\) if that side is toward decreasing radius and \(b>0\),
or toward increasing radius and \(b<0\). In these two cases a
counterclockwise small circle crosses the \(i,j\) arc from \(i\)
to \(j\). In the other cases its index is \(-1\).

An intermediate positive stretch of \(h\) contributes opposite
indices at its two ends, because the sign of \(b\) cannot
change there. Isolated non-sign-changing zeros of \(h\) inside
such a stretch do not permit a sign change of \(b\), by the
same accumulating-wall argument using (G20). Every active
central ray has a first positive stretch ending at an index
\(+1\) junction and a last positive stretch beginning at an
index \(+1\) junction. Every inactive central ray has only the
last such endpoint. Any other genuine triple on a pair line
has no sign change of \(h\). It is collinear, since in the
noncollinear case the derivative of \(h\) along the exact
pair line is nonzero. Its collinear normal form then has an
even-order cusp and index zero. Thus the sum of the indices at
the nonzero vertices on the six rays is \(+9\). The remaining
sum is \(2-(-1+9)=-6\), so there are at least six negatively
oriented trivalent vertices off the pair lines. Together with
the origin and the nine ray vertices, this proves the count.

Each open \(60^\circ\) sector contains at least one of the
off-line negative vertices. Truncate the sector at a small
central radius and a large exterior radius. On any compact
segment of a bounding \(i,j\) line, a sufficiently thin
neighborhood contains no off-line \(i,j\) interface. Otherwise
these interfaces would accumulate at a point of the segment.
If \(b\ne0\) there, the local pair equality is uniquely the
line by the implicit function theorem. If \(b=0\), the
accumulation contradicts (G20). The two radial sides of the
truncation can therefore be pushed a little into the sector,
transverse to the interfaces and avoiding all vertices, with
no crossing of the respective line pair.

Let \(i\) be the exterior label in the sector. Reading the
boundary counterclockwise, the left side runs outward and
misses \(i-1\leftrightarrow i\). The right side runs inward
and misses \(i\leftrightarrow i+1\). The allowed transitions on
each side form a path on the three labels. Their signed sum
depends only on the endpoints, because traversing an edge
and traversing it in reverse cancel. If the central label is
\(i-1\), the two sides contribute \(-2\) and \(-1\), respectively.
If it is \(i+1\), they contribute \(-1\) and \(-2\). The inner
and outer arcs make no transitions. The boundary index of
this truncated sector is consequently \(-1\), which forces a
negative trivalent vertex in the open sector.

The central and exterior parts of the rotational zero mode have
strictly negative mutual interaction before the intervening
interfaces are included. Let \(J\) rotate vectors counterclockwise
by \(90^\circ\). Orient each pair from \(i\) to \(i+1\), with labels
read cyclically in the exterior, and put
\(\zeta(x)=Jx\cdot N(x)\). For a small radius \(r_1>0\) and a
large radius \(r_2>r_1\), let \(I\) be the interface inside
\(B_{r_1}\), \(E\) the interface outside \(B_{r_2}\), and \(M\)
the remaining interface. Choose \(r_1,r_2\) so that \(I\) consists
of the three straight germ rays and \(E\) of the six exterior rays.
The cyclic label order of the germ is opposite to the exterior
order. Thus \(\zeta(x)=-|x|\) on \(I\) and
\(\zeta(x)=|x|\) on \(E\).

Write \(d(x)=e_i-e_{i+1}\) for the oriented pair's label vector,
\(d\nu=\phi_2\,ds\) for Gaussian arclength, and
\(K(x,y)=(d(x)\cdot d(y))M_\rho(x,y)\) for the signed interaction
kernel. Let \(C_{XY}=\int_{X\times Y}K(x,y)\zeta(x)\zeta(y)
\,d\nu(x)d\nu(y)\) for two of the pieces. Then
\[
 C_{IE}= -\frac{12}{(2\pi)^2\sigma^2}
    \int_0^{r_1}\!\int_{r_2}^{\infty}
        uv\,e^{-(u^2+v^2)/(2\sigma^2)}
   \left[\cosh\!\left(\frac{\rho uv}{\sigma^2}\right)
      -\cosh\!\left(\frac{\rho uv}{2\sigma^2}\right)\right]dv\,du
   <0.                                                        \tag{ED2}
\]
For a fixed central ray at radius \(u\), the two exterior rays of
its pair type are aligned and antialigned, with label dot product
two. The remaining four make angles of \(60^\circ\) or \(120^\circ\)
with the central ray and have label dot product \(-1\). In the
joint Gaussian density at radii \(u,v\) the resulting angular sum
is \(4(\cosh(t)-\cosh(t/2))\), where \(t=\rho uv/\sigma^2\).
Summing the three central rays and using the opposite signs of
\(\zeta\) gives (ED2).

Each of the three cut speeds \(1_X\zeta\) is balanced. For smooth
radial cutoffs this is Gaussian divergence of the field
\(x\mapsto\chi(|x|)Jx\), whose Gaussian divergence vanishes
identically. Approximation gives the three speeds, since their
Gaussian arclength integrals with polynomial weights are finite.
The rotation equation \(K\zeta=a\zeta\) then gives the constrained
form
\[
 Q_0\!\left(\sum_{X\in\{I,M,E\}}\alpha_X1_X\zeta\right)
   =\sum_{X<Y} C_{XY}(\alpha_X-\alpha_Y)^2 .
\]
Put \(L=-C_{IE}>0\), \(X=C_{IM}\) and \(Y=C_{ME}\).
Nonnegativity on the coefficients with \(\alpha_M=0\) requires
\[
        X>L,\qquad Y>L,\qquad XY\ge L(X+Y),
                 \qquad X+Y\ge4L .                           \tag{ED3}
\]
Indeed the diagonal entries are \(X-L,Y-L\) and the
off-diagonal entry is \(L\). The determinant bound and
\(XY\le (X+Y)^2/4\) give the two last inequalities.
The intervening network must also have sufficient radial
interface length in a weighted sense. The rotation equation on
\(M\) gives
\[
 \int_M a\zeta^2\,d\nu
   =C_{IM}+C_{ME}+C_{MM}\ \ge\ 4L+C_{MM}\ \ge\ 4L.
\]
The nonnegativity of \(C_{MM}\) follows from the positive
definiteness of the scalar Mehler kernel applied coordinate by
coordinate to the vector-valued interface measure. Since
\(|\zeta(x)|=|x||\tau(x)\cdot e_r(x)|\) for the unit interface
tangent \(\tau\), this is a lower bound on the squared radial
tangent component integrated against \(a|x|^2d\nu\).
The central and exterior cuts alone give no upper bound on it.
[Rotation Jacobi equation and the cut form](THREE_CELL_THETA_ROTATION_CUT_MATRIX.md#the-exact-matrix-and-its-mass-constraints)

## Tangential simplex contacts have positive curvature

Suppose \(A\) has the six-sector exterior and ties the centered
regular simplex \(T\) at \(r\), with simplex-optimal correlations
accumulating at \(r\). For the fixed partition deficit
\(D(s)=S_s(T)-S_s(A)\), suppose \(D'(r)=0\). Then \(D''(r)>0\).

The mixed shape test (G10) implies that \(D''(r)=0\) would make
the mixed score derivative \(g_{ij}=\partial_r
T_r(1_{A_i}-1_{A_j})\) constant on the whole regular \(i,j\)
interface. Test with any compact regular speed of zero Gaussian
arclength mean on that pair, including speeds made from two bumps
on different components. Such a speed preserves both independent
mass constraints to first order. On any exterior ray of this pair,
the mixed derivative tends to zero. In fact reflection in its line
swaps the two cone labels and fixes the observation point, so for
the exterior reference cone the score contrast is zero there for
every correlation. The actual cell indicators differ from the
exterior cone on a bounded set, and the conditional Gaussian
formula for their correlation derivative gives Gaussian decay
on the ray. Hence the constant for every pair vanishes.

With \(g=0\) on the entire actual interface, the dilation test
(G19) forces \(x\cdot N=0\) almost everywhere there. This means
\(x\cdot D1_{A_i}=0\) in distributions, so the actual partition
is centered conical. Its classification gives the simplex,
contrary to the six-sector exterior. Taylor expansion of \(D\) at
the approaching simplex-optimal correlations gives \(D''(r)\ge0\),
which completes the strict curvature conclusion.
[Mixed shape and correlation test, (G10)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)
[Dilation test, (G19)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)

## A global single-scale symmetry suffices for conicality

For any finite number of analytic winner scores whose pair contrasts
are nonidentically zero, exact global invariance under a nonunit
positive dilation about \(p\) implies invariance under all positive
dilations about \(p\), modulo null sets. For each contrast take its
first nonzero homogeneous Taylor term at \(p\), allowing order
zero. Away from the union of their homogeneous zero sets, all
pairwise signs at \(p+\lambda^n x\) eventually agree with the signs
of those terms at \(x\). These determine a unique winner and are
constant along each positive ray. Exact discrete invariance gives
the same label almost everywhere at \(p+x\).

For an actual balanced three-cell maximizer, the arbitrary-center
conical classification now gives the centered regular simplex.
There is also a uniform exclusion of approximate global symmetry
over a compact positive-correlation band, for centers and contraction
ratios in compact subsets of \(\mathbb R^2\) and \((0,1)\),
respectively. In fact the nonsimplex maximizers on the band form
a strongly compact class in the active plane by the uniform simplex
isolation. Gaussian label distance from the transform under a
specified dilation is jointly continuous in the partition and
these parameters. For continuity, restrict the affine changes of
variables to a fixed ball, where the transported Gaussian densities
are uniformly comparable, and obtain a uniform tail bound from
compactness of the parameters. Approximation by continuous compactly
supported functions proves continuity in the affine maps. If the
disagreement had infimum zero, a limit would exhibit exact global
discrete invariance for a nonsimplex maximizer, which has just been
excluded.
[Arbitrary-center conical classification](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)
[Uniform local simplex isolation](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md#uniform-repair-and-comparison)
