# Balanced planar cones with four or more cells are cylindrically unstable

Fix \(k\ge4\) and \(0<\rho<1\). An ordinary balanced
planar adjusted-score winner conical about a point
cannot satisfy the balanced second-variation inequality
after extension by one independent Gaussian coordinate.
Consequently no actual balanced \(k\)-cell maximizer in
ambient dimension at least three can have a conical
active plane. The instability is a second angular mode
of the cycle of interface rays.

## Reduction to a regular sector fan

The positive-mass winner scores are real analytic and
no adjusted pair contrast vanishes identically, by
Gaussian injectivity. Around a center of conicality,
choose for each contrast a radial direction on which
it is not identically zero. For almost every radius
all pair contrasts restrict to nonzero analytic
functions on the circle and therefore have only
finitely many zeros there. Conicality then gives
finitely many sectors. The full supporting lines
of their boundary rays define a finite polyhedral
partition about the same center, allowing unions of
its chambers as cells.

The adjusted pair scores agree on each shared open
facet. The single-correlation centering lemma for
balanced polyhedral cones therefore makes the
partition conical about the Gaussian origin. The
planar conical classification with an extra stable
coordinate then makes it the regular \(k\)-sector
fan. That classification excludes adjacent pair
contrasts of Taylor degree greater than one by
reciprocal-gradient integrability at the vertex.
It then uses the exact pair lines and their
reflections to make the sectors equally wide.
[Centering balanced polyhedral cones](AFFINE_SIMPLEX_FAN_RIGIDITY.md#centering-balanced-polyhedral-cones-from-one-correlation)
[Planar conical classification under cylindrical stability](AFFINE_SIMPLEX_FAN_RIGIDITY.md#planar-origin-conical-classification-with-cylindrical-stability)

## The signed cycle form

Write \(\theta=2\pi/k\). Label the cells cyclically and
orient ray \(j\) from cell \(j\) into cell \(j+1\),
with indices modulo \(k\). Its normal \(N_j\) rotates
by \(\theta\) from one ray to the next. The incidence
vectors \(d_j=e_j-e_{j+1}\) satisfy
\[
 d_i\cdot d_j=
 \begin{cases}
  2,&i=j,\\
 -1,&i,j\text{ are neighboring rays in the cycle},\\
  0,&\text{otherwise}.
 \end{cases}
\]
Let \(a_j\) be the magnitude of the normal score
gradient on ray \(j\), \(d\nu_j=\phi_2\,d\mathcal H^1\),
and \(M_\rho\) the Gaussian-relative Mehler kernel.
By rotation and cyclic relabeling the integral
\(c=\int a_j\,d\nu_j\) is independent of \(j\).
Write \(m_0,m_1\) for the integrals of \(M_\rho\)
over an equal-ray product and a neighboring-ray
product, respectively, with their surface measures.
Both are strictly positive and finite. Indeed
\(M_\rho(x,y)\phi_2(x)\phi_2(y)\) is bounded by a
constant times \(e^{-c_\rho(|x|^2+|y|^2)}\) with
\(c_\rho>0\). Each \(a_j\) is bounded by Gaussian
convolution derivative estimates.

For speeds \(f_j\) that are constant along each ray,
the first extra-coordinate form is
\[
 Q^{(1)}(f)
   =c\sum_j f_j^2
     -2\rho\left(m_0\sum_j f_j^2
                      -m_1\sum_j f_j f_{j+1}\right).       \tag{CS1}
\]
To justify this test, first cut off the speeds near
the vertex and at large radius. Multiply by a
compactly supported odd function of the independent
coordinate. This realizes compact smooth normal
speeds on regular wall interiors with zero linear
volume flux. Exact volume correction can use
spanning-tree bumps on the rays. Let the odd
functions tend in Gaussian \(L^2\) to the linear
coordinate and then remove the radial cutoffs by
the integrability bounds. Cylindrical stability
therefore requires \(Q^{(1)}(f)\ge0\) in (CS1).
[Heilman–Tarter, Lemma 4.8](https://arxiv.org/html/2011.05583#S4)
[The first extra-coordinate form, G1](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

The normal translation identity \(\rho K s_v=a s_v\),
where \(s_v=v\cdot N\), follows by Gaussian integration
by parts and stationarity on the rays. In particular
its first extra-coordinate form vanishes. Choose the
translation whose normal components are
\(f_j=\cos(j\theta)\). Since
\(\sum_j f_jf_{j+1}=\cos\theta\sum_j f_j^2\), (CS1)
gives \(c=2\rho[m_0-m_1\cos\theta]\).
[Heilman–Tarter, Lemma 5.4](https://arxiv.org/html/2011.05583#S5)
[The normal translation identity, G2](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Now take \(f_j=\cos(2j\theta)\). Its local term still
uses \(c\), and its cyclic neighbor pairing is
\(\cos(2\theta)\sum_jf_j^2\). Substitution gives
\[
 Q^{(1)}(f)
   =2\rho m_1[\cos(2\theta)-\cos\theta]\sum_jf_j^2<0 .
                                                               \tag{CS2}
\]
Indeed \(0<\theta\le\pi/2\) and cosine strictly decreases
between \(\theta\) and \(2\theta\). This contradicts
cylindrical stability. For \(k=3\), the second mode
has the same cyclic pairing as the first, so this
calculation is consistent with the two mixed
rotation zero modes of the tripod.

If a global maximizer in dimension at least three
is cylindrical over an active plane, adding a
single independent coordinate to that planar
partition must give a maximizer in dimension
three. An improvement there would extend to the
original ambient space. Thus its planar partition
has the cylindrical stability used above.
The argument excludes conical active-plane shapes.
It leaves general nonconical rank-two maximizers and
higher-dimensional rank-deficient cones open.
[Heilman–Tarter, Theorem 1.9](https://arxiv.org/html/2011.05583#S1)

## Homogeneous planar cones for the Euclidean Gaussian cut

Let \(G_\varepsilon\) be the centered Euclidean Gaussian
density of covariance \(\varepsilon^2I\). A partition
that minimizes the \(G_1\) cut against compact changes,
with free volumes, may be conical over a planar factor.
If it is, that planar cone has at most three positive
phases. Three phases force one \(120^\circ\) tripod.
Two phases force a straight line. This also applies
to planar cones extended over any number of independent
Euclidean coordinates.

A compact-replacement minimizer assigns almost every
point to a largest convolution score. Otherwise a
transfer on a sufficiently small bounded set of
strictly positive score gap improves the cut by a
linear term that dominates its quadratic interaction
error. The pair contrasts \(G_1*(1_{A_i}-1_{A_j})\)
are real analytic and are nonzero for distinct
nonempty cells by Gaussian convolution injectivity.
Restricting to a generic circle in the active plane
therefore gives finitely many genuine switching rays.
[Euclidean score condition for homogeneous minimizers](EUCLIDEAN_THREE_CELL_JUNCTION.md#homogeneous-minimizing-junctions-and-stationary-compact-defects)

By translation, take the cone's vertex to be the origin.
It is invariant under simultaneous dilation of all ambient
coordinates. For a compact change \(B\) of \(A\), write
\[
 \mathcal J_\varepsilon(B,A)=\iint G_\varepsilon(x-y)
     \left[1_{\{B(x)\ne B(y)\}}-1_{\{A(x)\ne A(y)\}}\right]\,dx\,dy ,
\]
where \(A(x)\) and \(B(x)\) denote the respective labels.
Changing variables in the cut of the dilated competitor
gives \(\mathcal J_1(\varepsilon^{-1}B,A)
=\varepsilon^{-d}\mathcal J_\varepsilon(B,A)\).
Thus \(\mathcal J_\varepsilon(B,A)\ge0\) for every
\(\varepsilon>0\). For a piecewise polyhedral \(B\),
\(\varepsilon^{-1}\mathcal J_\varepsilon(B,A)\) tends
to \(\sqrt{2/\pi}\) times the relative Euclidean
interface area. To see this, first restrict both
endpoints to a box containing the compact change with
positive clearance. Contributions with one endpoint
outside this box have a Gaussian tail. Inside the box,
the normalized cut for a fixed displacement
\(\varepsilon z\) tends to the sum over facets of
their areas times \(|z\cdot N|\), since a crossing of
a flat facet lies in a slab of width
\(\varepsilon|z\cdot N|\) and edges contribute
vanishing relative measure. The finitely many
supporting hyperplanes bound the normalized cut by
a constant times \(|z|\). Integration in \(z\) gives
the stated limit. A compact negative relative
interface area would therefore contradict cut minimality.
[Small-kernel perimeter comparison in the Euclidean problem](EUCLIDEAN_THREE_CELL_JUNCTION.md#homogeneous-minimizing-junctions-and-stationary-compact-defects)

Consider two consecutive rays at angle \(\alpha<\pi\)
with unit endpoints \(z_1,z_2\). Move their common
endpoint to \(p\), at distance \(t>0\) from the origin
along the interior bisector. Relabel the triangle
with vertices \(0,z_1,p\) by the label beyond the
first ray, and relabel the triangle with vertices
\(0,p,z_2\) by the label beyond the second ray.
The old segments \(0z_1,0z_2\) disappear, and the
new segments are \(pz_1,pz_2\), together with \(0p\)
if the two outer labels differ. All other interfaces
are unchanged. The derivative of the perimeter
change at \(t=0\) is
\[
 \mathbf 1_{\{\text{outer labels differ}\}}-2\cos(\alpha/2).
                                                               \tag{CS3}
\]
With at least four rays some consecutive gap is at
most \(\pi/2\), making (CS3) negative. With three rays,
every gap must be at least \(2\pi/3\), so all three
are equal. With two rays the outer labels in this
test coincide. If the rays are not opposite, the
smaller gap is less than \(\pi\), again giving a
negative derivative.

For a cone over a planar factor with \(m\ge1\) extra
coordinates, use the planar change over the transverse
box \([-L,L]^m\). Its planar interface improvement
is multiplied by \((2L)^m\). The new walls on the
transverse box boundary have area at most the
planar changed area times \(2m(2L)^{m-1}\).
For large \(L\) the net interface change is negative.
Thus the same exclusion holds in every ambient
dimension. In particular a microscopic full tie with
at least four positive phases cannot have this conical
planar structure.
[Microscopic unit-Gaussian minimizers at full ties, JD6](GAUSSIAN_NEAR_ONE_JUNCTION_DEGREE_AND_COMPONENTS.md#microscopic-minimizing-partitions-at-the-ties)
