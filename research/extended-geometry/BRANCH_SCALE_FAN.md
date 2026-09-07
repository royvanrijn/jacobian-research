# Branch-scale fans and the missing H1/H2 graph

## Result of the first three-cluster experiment

The degree-five obstruction is the first radial chart of a general
branch-scale construction, but the hyperplanes

\[
 \mu_i v(x_i)=\mu_j v(x_j)
\]

do **not** by themselves determine the full normalized branch-stabilization
graph.

The exact degree-six profile `(2,2,2)` gives eight conclusions.

1. The moving critical-value scales are monomials of the predicted degrees.
   Their pairwise initial ratios produce the weighted braid fan.  In the
   equal-weight degree-six chart this is the unimodular \(A_2\) fan with six
   maximal cones.
2. On the locus where three leading critical values coincide, a second-scale
   cross-ratio varies while all radial braid-fan data remains fixed.  The full
   stable-target graph therefore strictly refines the radial toric graph.
3. In degree six that refinement is now explicit.  The radial target is the
   blowup of \(\mathbf P^2\) at its three coordinate points; the complete
   stable target is the further blowup at \([1:1:1]\), hence
   \(\overline M_{0,5}\).  Its pullback has four reduced triple-Maxwell
   branches over the root-scale chart.
4. All thirteen radial source scale types have explicit degree-six component
   decompositions, matching node-index partitions, and saturated Kummer
   monoids.
5. The nonradial target arrangement adds no hidden local center.  Each
   pairwise Maxwell divisor pulls back to two Kummer branches with source-node
   partition \((1,1,2,2)\); the three residual radial--Maxwell crossings are
   transverse, and the coordinate intersections are already radial equality
   faces.
6. The central Hurwitz problem has two ambient algebraic cover classes, but
   the labelled source-root cross-ratio separates them.  Together with ACV
   deformation theory, this proves that the full admissible-cover graph and
   the explicit branch graph have the same labelled coarse normalization in
   this chart.
7. Normalization choices and inertia are different.  Full-chain
   centralizer matching shows that connector flips permute the normalized
   radial lifts, so the labelled radial inertia is trivial even when the
   nodal special fiber has connector automorphisms.  Pairwise
   and triple Maxwell nodes have respectively two and four normalization
   branches, but each normalized branch retains only one diagonal
   \(\mu_2\).
8. The resulting algebraic candidate is canonical: take second roots of
   the three pairwise Maxwell boundary divisors and the exceptional
   triple-Maxwell divisor.  This root construction is \(S_3\)-equivariant;
   the separate pair-root quotient adds
   \((S_2)^3\rtimes S_3\).  Smooth tame stack reconstruction identifies it
   with the labelled ACV graph on this chart.  The global \(S_6\)-action is
   canonical on graph closure; an explicit recursive presentation on all
   other root-partition types remains.

Thus the correct candidate for the H1/H2 object is recursive: a weighted
braid subdivision at the first scale, followed by nested resonance
modifications on equal-scale strata.  Equivalently, it should be the
logarithmically saturated pullback of the tropical stable-target or tropical
admissible-cover fan, not a single monomial blowup and not merely the first
weighted braid fan.

The general local normalization/inertia theorem and the corrected formal H2
quotient are isolated in
[Labelled node saturation and the corrected H2
quotient](LABELLED_NODE_SATURATION.md).  Full root-label gluing is canonical
on the global labelled normalized graph closure.  The stable-target graph is
now identified with the normalized wonderful pullback, and
[source-vertex rigidity](SOURCE_VERTEX_RIGIDITY.md) reconstructs its
component maps from compatible target-flag fibers.  The
[finite-normalization theorem](SOURCE_GRAPH_FINITE_NORMALIZATION.md)
identifies the resulting coarse polynomial source graph with the wonderful
target graph, so corrected H2 is unconditional.  The
[general radial source atlas](GENERAL_RADIAL_SOURCE_ATLAS.md) completes the
source construction on every ordered first-scale stratum and for arbitrary
cluster multiplicities.  The
[polynomial monodromy-forest theorem](POLYNOMIAL_MONODROMY_FORESTS.md)
also fixes every simple-branch resonance source dual graph and node
partition.  The
[monodromy-centralizer theorem](MONODROMY_INERTIA_CHARACTERS.md) closes all
radial and generic simple-resonance inertia.  The
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) supplies all
higher-codimension flag-divisor coordinates, contraction maps, and
simultaneous characters; its selected-factor/log-étale comparison supplies
the global coverage and descent required for `H1-STACK`.

The exact checker is
[`scripts/verify_branch_scale_fan.py`](../scripts/verify_branch_scale_fan.py),
and its machine-readable fan and valuation data is
[`artifacts/generated-results/branch_scale_fan_degree6.json`](../artifacts/generated-results/branch_scale_fan_degree6.json).
The complete target-arrangement checker is
[`scripts/verify_degree_six_branch_target_graph.py`](../scripts/verify_degree_six_branch_target_graph.py),
with artifact
[`artifacts/generated-results/branch_target_graph_degree6.json`](../artifacts/generated-results/branch_target_graph_degree6.json).

## 1. The local scale lemma

Let \(c_1,\ldots,c_r\) be distinct cluster centers and suppose the roots in
cluster \(i\) have the one-scale form

\[
 c_i+x_i\alpha_{i1},\ldots,c_i+x_i\alpha_{i\mu_i}.
\]

Write

\[
 Q_i(X)=\prod_{a=1}^{\mu_i}(X-\alpha_{ia}).
\]

Spectator roots and all other clusters contribute a unit at \(c_i\).  For a
polynomial \(F\) having these roots,

\[
 F(c_i+x_iX)
   =x_i^{\mu_i} U_i Q_i(X)+
     o(x_i^{\mu_i}),                              \tag{1.1}
\]

where

\[
 U_i=
 \prod_{j\ne i}(c_i-c_j)^{\mu_j}
 \prod_{\text{spectators }s}(c_i-s)
\]

is nonzero.  Differentiating with respect to \(W\) gives

\[
 F'(c_i+x_iX)
   =x_i^{\mu_i-1}U_iQ_i'(X)+
     o(x_i^{\mu_i-1}).                            \tag{1.2}
\]

If \(\xi\) is a simple root of \(Q_i'\) and \(Q_i(\xi)\ne0\), Hensel
lifting gives a critical point

\[
 r_{i,\xi}=c_i+x_i\xi+o(x_i)
\]

and a moving critical value

\[
 \lambda_{i,\xi}
   =x_i^{\mu_i}U_iQ_i(\xi)(1+o(1)).               \tag{1.3}
\]

An affine target normalization multiplies all differences of critical values
by the same factor.  It therefore cancels from ratios.  For any two
nonzero local critical values,

\[
 \frac{\lambda_{i,\xi}}{\lambda_{j,\eta}}
   =u_{ij,\xi\eta}
      \frac{x_i^{\mu_i}}{x_j^{\mu_j}},
 \qquad u_{ij,\xi\eta}\in\mathcal O^\times.        \tag{1.4}
\]

This proves the proposed exponent rule on a one-scale, noncaustic cluster
chart.  If \(Q_i'\) has multiple roots, if \(Q_i(\xi)=0\), or if the internal
shape \((\alpha_{ia})\) also degenerates, (1.3) is only the first layer and
must itself be refined.

The multiplicity in (1.4) is the **number of roots in the cluster**, not the
ramification index \(\mu_i-1\).  A cluster of \(\mu_i\) roots has
\(\mu_i-1\) local critical points, but its nonzero critical values have
radial order \(\mu_i v(x_i)\).

## 2. Recovery of the degree-five blowup

For the degree-five chart in
the [corrected Hurwitz--LL comparison](HURWITZ_LL_COMPACTIFICATION.md),
the local cluster polynomials are

\[
 Q_0(X)=X^2(X-1),\qquad Q_1(X)=X(X-1).
\]

Their moving critical points are \(2/3\) and \(1/2\), with nonzero critical
values \(-4/27\) and \(-1/4\).  After the common target normalization,

\[
 \lambda_0=-\frac4{27}\frac{x^3}{y}(1+o(1)),
 \qquad
 \lambda_1=-\frac14y(1+o(1)).
\]

Hence

\[
 \frac{\lambda_1}{\lambda_0}
   =\frac{27}{16}\frac{y^2}{x^3}(1+o(1)).          \tag{2.1}
\]

The wall is \(3v(x)=2v(y)\).  The primitive new ray is \((2,3)\), exactly
the fan subdivision of the normalized blowup of \((x^3,y^2)\).  In general,
the normalized graph of

\[
 [x_i^{\mu_i}:x_j^{\mu_j}]
\]

has wall \(\mu_i v(x_i)=\mu_jv(x_j)\); in its two-coordinate plane the
primitive ray is

\[
 \frac1{\gcd(\mu_i,\mu_j)}(\mu_j,\mu_i).           \tag{2.2}
\]

This explains the degree-five modification without treating it as an
isolated blowup.

## 3. The degree-six three-cluster chart

Take

\[
 F(W)=
 W(W-x)(W-1)(W-1-y)(W-3)(W-3-z).                  \tag{3.1}
\]

At \(x=y=z=0\), the roots form three pair clusters centered at \(0,1,3\).
The limiting local models are

\[
\begin{array}{c|c|c}
\text{center}&x_i^{-2}F(c_i+x_iX)&
 \lambda_i/x_i^2\\ \hline
0&9X(X-1)&-9/4\\
1&4X(X-1)&-1\\
3&36X(X-1)&-9.
\end{array}                                       \tag{3.2}
\]

Thus the three moving critical values are

\[
 \lambda_x=-\frac94x^2+o(x^2),\qquad
 \lambda_y=-y^2+o(y^2),\qquad
 \lambda_z=-9z^2+o(z^2).                          \tag{3.3}
\]

The special derivative factors as

\[
 F_0'(W)=
 2W(W-1)(W-3)(3W^2-8W+3).
\]

The three linear factors give the moving critical values in (3.3).  The
quadratic factor gives two critical points whose critical values remain
nonzero, so (3.3) is the complete moving list.

The valuation vectors in the coordinates \((x,y,z)\) are

\[
\begin{array}{c|c}
\text{function}&\text{valuation vector}\\ \hline
\lambda_x&(2,0,0)\\
\lambda_y&(0,2,0)\\
\lambda_z&(0,0,2)\\
\lambda_x/\lambda_y&(2,-2,0)\\
\lambda_x/\lambda_z&(2,0,-2)\\
\lambda_y/\lambda_z&(0,2,-2).
\end{array}                                       \tag{3.4}
\]

As an independent check, the script eliminates \(W\) from
\(F'(W)=0\) and \(L-F(W)=0\).  For each of the six traits obtained by
permuting

\[
 (v(x),v(y),v(z))=(1,2,3),
\]

the lower Newton polygon of the degree-five critical-value eliminant gives

\[
 \{0,0,2,4,6\}.                                   \tag{3.5}
\]

The two zeroes are the stationary critical values and the other three are
exactly \(2v(x),2v(y),2v(z)\).

## 4. The radial normalized graph

Consider the graph of all three initial ratio maps

\[
 [x^2:y^2],\qquad [x^2:z^2],\qquad [y^2:z^2].
\]

The fan of its normalization is the common refinement by

\[
 2v_x=2v_y,\qquad 2v_x=2v_z,\qquad 2v_y=2v_z.      \tag{4.1}
\]

It has six maximal cones, one for each total ordering of
\(2v_x,2v_y,2v_z\).  For example,

\[
 0\le v_x\le v_y\le v_z
\]

is generated by

\[
 (1,1,1),\qquad(0,1,1),\qquad(0,0,1).              \tag{4.2}
\]

In general the cone for
\(v_{i_0}\le v_{i_1}\le v_{i_2}\) is generated by the full indicator,
the indicator of \(\{i_1,i_2\}\), and the indicator of \(\{i_2\}\).
All six determinant checks are \(1\).  The rays in the positive octant are

\[
\begin{split}
 &(1,0,0),(0,1,0),(0,0,1),\\
 &(1,1,0),(1,0,1),(0,1,1),(1,1,1).
\end{split}                                       \tag{4.3}
\]

Modulo the common diagonal ray this is the ordinary \(A_2\) braid fan.  For
unequal multiplicities it is pulled back by

\[
 (v_1,\ldots,v_r)\longmapsto
 (\mu_1v_1,\ldots,\mu_rv_r),                       \tag{4.4}
\]

and the lattice must be saturated.  This is the precise sense in which the
first correction is a weighted braid-fan modification.

On a nonresonant trait,

\[
 v(\lambda_i-\lambda_j)
   =\min\{\mu_i v(x_i),\mu_jv(x_j)\}                \tag{4.5}
\]

when the two weights differ, and the same formula holds on an equality wall
when the initial residues are distinct.  Consequently the six chambers give
the six radial caterpillar types of the stabilized target containing
\(0,\lambda_x,\lambda_y,\lambda_z,\infty\).

This is exactly the level at which existing tropical admissible-cover theory
supports the proposal: Cavalieri--Markwig--Ranganathan construct the
tropicalization of admissible covers and prove compatibility of the
classical and tropical tautological, including branch, maps
([arXiv:1401.4626](https://arxiv.org/abs/1401.4626)).  Tropical moduli of
stable rational curves supplies the target fan of metric trees
([arXiv:0708.2268](https://arxiv.org/abs/0708.2268)).  These results support a
skeleton-level comparison; they do not identify the repository-specific
normalized graph or its completed local rings.

## 5. Why the braid fan is not the whole graph

The initial Maxwell differences in (3.3) are

\[
\begin{aligned}
 \operatorname{in}(\lambda_x-\lambda_y)
   &=-\frac14(3x-2y)(3x+2y),\\
 \operatorname{in}(\lambda_x-\lambda_z)
   &=-\frac94(x-2z)(x+2z),\\
 \operatorname{in}(\lambda_y-\lambda_z)
   &=-(y-3z)(y+3z).
\end{aligned}                                      \tag{5.1}
\]

These are binomial resonance loci inside the equal-weight strata.  They are
not torus-invariant boundary strata of the radial fan.

The decisive test is the triple resonance

\[
\begin{aligned}
 x&=t+At^2,\\
 y&=\frac32t+Bt^2,\\
 z&=\frac12t+Ct^2.
\end{aligned}                                      \tag{5.2}
\]

All three critical values have the same first term

\[
 \lambda_x=\lambda_y=\lambda_z
   =-\frac94t^2+O(t^3).                            \tag{5.3}
\]

Exact implicit expansion of the three critical points gives

\[
\begin{aligned}
 \lambda_x-\lambda_y
   &=-\frac34(6A-4B+1)t^3+O(t^4),\\
 \lambda_z-\lambda_y
   &= \frac32(2B-6C+1)t^3+O(t^4).
\end{aligned}                                      \tag{5.4}
\]

After the common first-scale target component is extracted, the next target
bubble has cross-ratio

\[
 R=
 \frac{\lambda_x-\lambda_y}{\lambda_z-\lambda_y}
 =
 -\frac{6A-4B+1}{2(2B-6C+1)}.                     \tag{5.5}
\]

For \((A,B,C)=(0,0,0)\), this is \(-1/2\); for \((1,0,0)\), it is
\(-7/2\).  The scale valuations, the leading scale ratios, and hence the
entire radial braid-fan point are identical in these two cases, while the
stabilized second-scale four-pointed target is different.

Therefore:

\[
 \boxed{\Gamma_{\mathrm{branch}}
   \text{ strictly refines }
   \Gamma_{\mathrm{radial}}
   \text{ at the triple-resonance locus}.}
\]

This does not contradict the degree-five calculation.  With only two moving
scales there is no three-point residue configuration carrying a further
cross-ratio.  The normalized blowup of \((x^3,y^2)\) is the complete first
radial correction in that chart.

## 6. The complete degree-six stable-target graph

The resonance refinement in Section 5 is a standard Kapranov blowup, not an
unspecified additional modification.

Fix the target markings \(0,\infty\).  Modulo the remaining target scaling,
three distinct nonzero branch values define

\[
 (0,\infty,\lambda_x,\lambda_y,\lambda_z)
 \longmapsto
 [A:B:C]=[\lambda_x:\lambda_y:\lambda_z]\in\mathbf P^2. \tag{6.1}
\]

Consequently

\[
 M_{0,5}\cong
 \mathbf P^2\setminus
 \{ABC(A-B)(A-C)(B-C)=0\}.                         \tag{6.2}
\]

The six-line arrangement in (6.2) has four triple points

\[
 p_A=[1:0:0],\quad p_B=[0:1:0],\quad
 p_C=[0:0:1],\quad p_\Delta=[1:1:1],               \tag{6.3}
\]

and three remaining double points

\[
 [0:1:1],\qquad[1:0:1],\qquad[1:1:0].             \tag{6.4}
\]

Kapranov's construction specializes here to

\[
 \overline M_{0,5}
 \cong
 \operatorname{Bl}_{p_A,p_B,p_C,p_\Delta}\mathbf P^2.
                                                               \tag{6.5}
\]

Kapranov's Chow-quotient construction identifies
\(\overline M_{0,n}\) with an iterated blowup of projective space
([arXiv:alg-geom/9210002](https://arxiv.org/abs/alg-geom/9210002)).
In the present case (6.5) can also be checked directly.  Each of the six
lines passes through exactly two points of (6.3), so its strict transform has
self-intersection \(1-2=-1\).  Together with the four exceptional curves,
these give ten \((-1)\)-curves.  Their dual incidence graph is the Petersen
graph, and

\[
 K^2=K_{\mathbf P^2}^2-4=9-4=5,                   \tag{6.6}
\]

the degree-five del Pezzo model of \(\overline M_{0,5}\).

The first three blowups in (6.5) are toric:

\[
 L_3=\operatorname{Bl}_{p_A,p_B,p_C}\mathbf P^2.   \tag{6.7}
\]

This is the permutohedral surface with the six \(A_2\) braid-fan chambers.
Thus \(L_3\) is exactly the radial target compactification of Section 4.
The only target blowup missing from the radial construction is

\[
 \overline M_{0,5}=\operatorname{Bl}_{p_\Delta}L_3.
                                                               \tag{6.8}
\]

The point \(p_\Delta\) records simultaneous equality of the three first-scale
branch residues, precisely the phenomenon detected in (5.2)--(5.5).

### Pullback to the root-scale chart

On the leading exceptional projective plane of the source chart, the branch
map is

\[
 [x:y:z]\longmapsto
 \left[-\frac94x^2:-y^2:-9z^2\right].             \tag{6.9}
\]

The inverse image of \(p_\Delta\) consists of four reduced points:

\[
 [x:y:z]
 =
 \left[1:\epsilon\frac32:\delta\frac12\right],
 \qquad \epsilon,\delta\in\{\pm1\}.                \tag{6.10}
\]

Indeed, its tangent-cone equations are

\[
 4y^2-9x^2=0,\qquad x^2-4z^2=0.                   \tag{6.11}
\]

On the chart \(x=1\), the Jacobian of (6.11) with respect to \(y,z\) is
nonzero at all four points.  Hence the four directions are reduced.

This leading calculation is exact on the critical-point etale chart.  Let
\(r_i\) be the critical point specializing to the center \(c_i\).  Substituting

\[
 r_i=c_i+x_iR_i
\]

into \(F'(r_i)=0\), dividing by \(x_i\), and applying the implicit-function
theorem at \(R_i=1/2\) gives a regular \(R_i\).  Substitution into \(F\)
then gives

\[
 \lambda_i=x_i^2u_i,\qquad u_i\in\mathcal O^\times,              \tag{6.12}
\]

with special units \(-9/4,-1,-9\).  After a finite etale extension adjoining
square roots of the units, (6.12) becomes a diagonal Kummer map.  The four
directions (6.10) therefore lift uniquely to four smooth formal branches of
the exact triple-Maxwell locus

\[
 I_\Delta=
 (\lambda_x-\lambda_y,\lambda_z-\lambda_y).        \tag{6.13}
\]

Let \(\widehat U\) denote this labelled critical-point chart.  The complete
degree-six stable-target graph over \(\widehat U\) is consequently the
normalization of the following explicit sequence:

1. resolve
   \([\lambda_x:\lambda_y:\lambda_z]\) by the normalized blowup of
   \((\lambda_x,\lambda_y,\lambda_z)\);
2. blow up the three strict transforms over the coordinate points
   \(p_A,p_B,p_C\), equivalently the three radial nested-scale axes;
3. blow up the four now-disjoint strict transforms of the formal branches
   of (6.13).

Steps 1 and 2 give the weighted braid-fan graph.  Step 3 is the exact
pullback of (6.8).  Thus the target-stabilization part of the requested
degree-six normalized graph is no longer conjectural.

### The equal-scale admissible source

Part of the source comparison can also be constructed exactly.  On the
equal-scale trait \(x=y=z=t\), the main target component carries the central
degree-six map

\[
 F_0(W)=W^2(W-1)^2(W-3)^2.                         \tag{6.14}
\]

Its ramification contributions are

\[
 5\quad\text{at infinity},\qquad
 1+1+1\quad\text{at }0,1,3,\qquad
 1+1\quad\text{at the two stationary critical points}.
\]

Their sum is \(10=2\cdot6-2\), as required by Riemann--Hurwitz.

The first target bubble receives three source tails with maps

\[
 9X(X-1),\qquad4X(X-1),\qquad36X(X-1).             \tag{6.15}
\]

They have degrees \(2+2+2=6\), branch values
\(-9/4,-1,-9\), and one point of index two over the attaching target node
at infinity.  Thus (6.14)--(6.15) form an admissible special fiber of total
degree six over both target components.

Let \(\tau\) smooth the target node and let \(s_x,s_y,s_z\) smooth its three
source nodes.  Admissibility gives

\[
 \tau=s_x^2=s_y^2=s_z^2.                            \tag{6.16}
\]

The completed node-deformation ring is therefore

\[
 R_{\rm node}=
 k[[\tau,s_x,s_y,s_z]]/
 (s_x^2-\tau,s_y^2-\tau,s_z^2-\tau).               \tag{6.17}
\]

Eliminating \(\tau\) gives the radical decomposition

\[
 (s_y^2-s_x^2,s_z^2-s_x^2)
 =
 \bigcap_{\epsilon,\delta=\pm1}
 (s_y-\epsilon s_x,s_z-\delta s_x).                \tag{6.18}
\]

Consequently the normalization of (6.17) is the product of four formal
discs

\[
 \tau=q^2,\qquad
 (s_x,s_y,s_z)=(q,\epsilon q,\delta q).             \tag{6.19}
\]

These are exactly the four Kummer directions (6.10).  Thus the normalized
source-node deformation and the pullback of the diagonal target center agree
on the equal-scale chart.

The first diagonal refinement agrees as well.  Centering each tail at its
critical point gives the exact equations

\[
\begin{aligned}
 S-\lambda_x&=9U_x^2,\\
 S-\lambda_y&=4U_y^2,\\
 S-\lambda_z&=36U_z^2.
\end{aligned}                                      \tag{6.20}
\]

When the three branch values collide and the target sprouts the
\(p_\Delta\)-bubble, (6.20) supplies three source nodes of index two above
the new target node, with the same saturated four-branch pattern (6.17)--(6.19).

This proves the admissible source comparison on the generic equal-scale face
and its first triple-Maxwell refinement.

### The complete radial admissible atlas

The source combinatorics and saturated node monoids can in fact be checked
on every radial cone and equality face.  Write an ordered scale partition as

\[
 B_0\mid B_1\mid\cdots\mid B_{\ell-1},              \tag{6.21}
\]

where the clusters in \(B_a\) become visible on the \(a\)-th target bubble.
On that bubble, the contribution of cluster \(i\), whose level is \(h_i\),
has one of three forms:

\[
\begin{array}{c|c|c}
 h_i>a & S=V_i^2 & \text{one degree-two connector},\\
 h_i=a & S=u_iX_i(X_i-1) & \text{one degree-two natural tail},\\
 h_i<a & S=V_{i,+},\ S=V_{i,-} & \text{two degree-one strands}.
\end{array}                                         \tag{6.22}
\]

Every row contributes total degree two, so every target bubble has source
degree six.  The degree-two components have ramification contribution two
and each degree-one strand has contribution zero, verifying
Riemann--Hurwitz component by component.

Suppose \(m\) clusters have activated before the node following a target
bubble.  The source-node indices there are

\[
 \underbrace{1,1,\ldots,1,1}_{2m\ {\rm entries}},
 \underbrace{2,\ldots,2}_{3-m\ {\rm entries}},       \tag{6.23}
\]

whose sum is always six.  If \(k=3-m>0\), saturation of the common smoothing
relation has \(2^{k-1}\) normalization branches; if \(k=0\), it has one.
For example, the strict chamber \(x\mid y\mid z\) has successive index
partitions

\[
 (2,2,2),\qquad(1,1,2,2),\qquad(1,1,1,1,2),         \tag{6.24}
\]

and Kummer branch counts \(4,2,1\).  The equal-scale face has the single
partition \((2,2,2)\) and four branches, recovering (6.17)--(6.19).

There are exactly thirteen ordered set partitions of three labelled
clusters: one one-level face, six two-level faces, and six strict chambers.
The verifier enumerates all thirteen, builds (6.22), checks degree and
Riemann--Hurwitz on every component, checks (6.23) at every node, and
computes every Kummer saturation count.  Thus the radial admissible source
atlas is complete at the combinatorial and characteristic-monoid level.

This does **not** yet construct one algebraic family over the entire radial
modification.  The deformation-theoretic discussion below separates the
ambient overlap gluing, which is standard, from the still-missing extension
that selects the polynomial admissible-cover component.  Label descent and
the global comparison remain separate.

### The complete Maxwell source atlas

The nonradial boundary can also be exhausted in this degree-six chart.  The
three pairwise Maxwell divisors are

\[
 A=B,\qquad A=C,\qquad B=C.                        \tag{6.25}
\]

After an etale extension absorbing the units in (6.12), write
\(\lambda_i=q_i^2\).  The pullback of a pairwise divisor is

\[
 q_i^2-q_j^2=(q_i-q_j)(q_i+q_j)=0,                 \tag{6.26}
\]

so it has exactly two reduced Kummer branches.

When \(\lambda_i=\lambda_j\), the new target bubble receives the two exact
quadratic tails

\[
 S-\lambda_i=c_iU_i^2,\qquad
 S-\lambda_j=c_jU_j^2,                             \tag{6.27}
\]

and two degree-one strands from the remaining cluster.  Its total source
degree is \(2+2+1+1=6\), and the source nodes over its attaching node have
indices

\[
 (1,1,2,2).                                        \tag{6.28}
\]

The completed node ring is, up to the two index-one variables that equal
\(\tau\),

\[
 k[[\tau,s_i,s_j]]/(s_i^2-\tau,s_j^2-\tau).
                                                               \tag{6.29}
\]

Its normalization has the two branches \(s_j=\pm s_i\), exactly matching
(6.26).  For the triple collision, the same construction gives
\((2,2,2)\) and the four sign branches of (6.17)--(6.19).

It remains to check where these strata meet the radial boundary.  The three
unblown double points of the target arrangement are

\[
 [0:1:1],\qquad[1:0:1],\qquad[1:1:0].
\]

Their source preimages are respectively

\[
\begin{aligned}
 &[0:\pm3:1],\\
 &[\pm2:0:1],\\
 &[1:\pm3/2:0].
\end{aligned}
\]

The Jacobian determinants of the radial and Maxwell equations at these six
points are \(\pm6,\pm9,\pm3\), so every intersection is reduced and
transverse.  At a coordinate triple point, the radial blowup records the
ratio of the two vanishing values.  The Maxwell equality has two reduced
directions on that exceptional curve:

\[
 [x:y]=[1:\pm3/2],\quad
 [x:z]=[1:\pm1/2],\quad
 [y:z]=[1:\pm3].                                   \tag{6.30}
\]

These are precisely the two-cluster equality faces already included in the
thirteen-type radial atlas.  At \(p_\Delta\), the additional target blowup
and its four source branches are exactly (6.8)--(6.20).  Therefore every
stratum and every intersection of the six-line target arrangement now has
the expected degree-six source components, node indices, and normalized
Kummer branches.  No further local blowup center is forced in this chart.

What is still missing is not another combinatorial stratum.  The next
question is whether standard admissible-cover deformation theory glues these
formal models and, separately, whether it selects the closure of the
polynomial family.

### What admissible-cover deformation theory now supplies

The ambient gluing problem is actually standard once an algebraic admissible
cover point has been selected.  If one target node has source nodes of
indices \(e_1,\ldots,e_r\), the Harris--Mumford completed ring has node
factor

\[
 k[[\tau,s_1,\ldots,s_r]]/
 (s_1^{e_1}-\tau,\ldots,s_r^{e_r}-\tau).           \tag{6.31}
\]

Its normalization has

\[
 M=\frac{\prod_i e_i}{\operatorname{lcm}(e_1,\ldots,e_r)}
                                                               \tag{6.32}
\]

branches.  Thus (6.32) gives \(M=2\) for \((1,1,2,2)\) and \(M=4\)
for \((2,2,2)\), exactly the pairwise and triple Maxwell counts above.

Abramovich--Corti--Vistoli identify the twisted-cover stack with the
normalization of the Harris--Mumford admissible-cover space
([arXiv:math/0106211](https://arxiv.org/abs/math/0106211)).
Cavalieri--Markwig--Ranganathan make the consequence used here explicit:
the normalized admissible-cover stack is smooth with toroidal boundary, its
branch morphism is toroidal, and the normalization multiplicity at a target
edge is (6.32)
([arXiv:1401.4626](https://arxiv.org/abs/1401.4626), Sections 4.2 and 5.2).
Therefore the rings (6.17), (6.29), and the radial node rings are not merely
candidate equations; after choosing the relevant algebraic cover point they
are the exact toroidal completed-ring factors, and their overlaps glue in
the ambient normalized admissible-cover stack.

This does **not** finish the repository comparison.  The same tropical
admissible-cover type can correspond to several irreducible algebraic
strata, and several normalization branches can map to the same tropical
cone.  This is precisely the distinction emphasized in Section 4.2.4 of
the cited tropicalization paper.  Consequently equality of fans,
component types, and node rings does not select the closure of the polynomial
family.

The degree-six local gap is now sharply reduced to one morphism problem:
construct the extension

\[
 \widetilde\Gamma_{(2,2,2)}
 \longrightarrow
 \overline{\mathcal H}^{\,\mathrm{ACV}}_{0\to0,6}              \tag{6.33}
\]

from the normalized branch graph, extending the smooth polynomial-cover
map, and prove that its image is the intended irreducible component.  Once
(6.33) exists, the ambient ACV deformation theorem supplies the overlap
gluing and completed local rings automatically.

### Central Hurwitz selection and the coarse comparison

There is one remaining possible source of nonuniqueness in (6.33): the
algebraic cover on the central target component.  Its ramification profiles
are

\[
 (6),\qquad(2,2,2),\qquad
 (2,1,1,1,1),\qquad(2,1,1,1,1).                  \tag{6.34}
\]

An exact \(S_6\) branch-cycle enumeration gives twelve factorizations after
fixing the six-cycle at infinity, and two simultaneous-conjugacy classes.
The two classes remain distinct after the Hurwitz move that forgets the
ordering of the simple branch points.  Thus the ambient target and
ramification data alone do **not** select the polynomial cover.

The source-root marking does.  A polynomial with profiles \((6)\) and
\((2,2,2)\) is the square of a cubic.  Normalize its three labelled roots to
\(0,1,a\):

\[
 H_a(W)=\bigl(W(W-1)(W-a)\bigr)^2.                 \tag{6.35}
\]

If \(u,v\) are the critical values of the cubic, then the two simple branch
values of \(H_a\) are \(u^2,v^2\).  The invariant

\[
 J(a)=\frac{(u^2+v^2)^2}{u^2v^2}
 =
 \frac{4(8a^6-24a^5+21a^4-2a^3+21a^2-24a+8)^2}
 {729a^4(a-1)^4}                                  \tag{6.36}
\]

forgets their ordering and common target scale.  The numerator of
\(J(a)-J(3)\) factors as

\[
\begin{aligned}
 64&(a-3)(a+2)(2a-3)(2a+1)(3a-2)(3a-1)\\
 &\cdot
 (9a^6-27a^5+79a^4-113a^3+79a^2-27a+9).          \tag{6.37}
\end{aligned}
\]

The six linear roots are exactly the \(S_3\)-relabeling orbit

\[
 3,\quad\frac13,\quad-2,\quad-\frac12,\quad
 \frac32,\quad\frac23.
\]

The final sextic is irreducible over \(\mathbf Q\) and accounts for the
second Hurwitz class.  Moreover \(a=3\) is a simple root of (6.37) and is not
a root of the sextic.  Hence the labelled source-root cross-ratio selects
the polynomial central cover etale-locally and with multiplicity one.

This closes the labelled **coarse** comparison in the degree-six chart.
Indeed, let \(G_6\) be the coarse normalization of the full polynomial
admissible-cover graph closure and let \(B_6\) be the explicit normalized
branch graph constructed above.  Forgetting the source cover gives a proper
birational map

\[
 \pi:G_6\longrightarrow B_6.                       \tag{6.38}
\]

The Kummer normalizations separate every node-gluing choice.  Every
noncentral vertex cover has degree one or two and is unique up to
isomorphism with the prescribed markings.  Equations (6.34)--(6.37) select
one reduced central cover.  Finally ACV deformation theory says that, after
the vertex covers are fixed, deformations are controlled by the target and
the node equations.  Thus \(\pi\) is quasi-finite.  It is therefore finite;
since it is birational and \(B_6\) is normal, it is an isomorphism.

Consequently (6.33) exists on the labelled coarse degree-six chart.  The
next two subsections distinguish
normalization choices from actual inertia and construct the canonical local
root-stack candidate.  The later recursive resonance theorem glues all
root-partition charts and proves the global comparison.

### Normalization choices are not automatically inertia

For one target node with source indices \(e_1,\ldots,e_r\), put
\(L=\operatorname{lcm}(e_1,\ldots,e_r)\).  Choosing roots of unity in the
normalization gives \(\prod e_i\) parametrizations, but simultaneous
\(\mu_L\)-scaling reparametrizes the same normalized component.  Hence

\[
 \#\{\text{normalized branches}\}
   =\frac{\prod_i e_i}{L}.                          \tag{6.39}
\]

The diagonal \(\mu_L\) in this count is not, by itself, inertia.  Actual
inertia is the subgroup realized by automorphisms of the cover that preserve
the chosen normalized lift and every source-root label.

This test changes the radial answer.  A connector \(S=V^2\) has a deck
flip, and a labelled tail on a different source component does not
componentwise remove it.  One must match the actions on both halves of every
source node along the complete target chain.  The deepest target screen
contains the labelled zero fiber, while the outer vertex is rigidified by
the stationary polynomial monodromy.  The full-centralizer calculation
shows that every nontrivial connector action changes at least one
simultaneous sign-ratio branch.  Thus

\[
 I_{\rm radial}^{\rm labelled}=1                   \tag{6.40}
\]

on all thirteen radial types.

At a Maxwell node the situation is different.  Its quadratic square tails
have no zero-root labels, so all their deck flips are permitted.  For \(m\)
colliding tails, \((\mu_2)^m\) acts on the \(2^{m-1}\) normalized
sign-ratio branches and the stabilizer of one branch is the diagonal
\(\mu_2\).  Therefore

\[
\begin{array}{c|c|c|c}
 \text{stratum}&\text{square-tail sign group}&
 \text{normalized branches}&\text{labelled inertia}\\ \hline
 \text{pairwise Maxwell}&(\mu_2)^2&2&\mu_2\\
 \text{triple Maxwell}&(\mu_2)^3&4&\mu_2\\
 x\mid y\mid z&(\mu_2)^3&8&1.
\end{array}                                         \tag{6.41}
\]

This also corrects the earlier interpretation of the tropical weights.
Formula (6.39) is the \(W3\) normalization multiplicity of
Cavalieri--Markwig--Ranganathan.  One must compute genuine cover
automorphisms only after selecting a normalized lift; the raw source-node
sign group cannot be inserted as \(W1\).

### The Maxwell-boundary root stack

Let \(D_{xy},D_{xz},D_{yz}\subset B_6\) be the strict transforms of the
three pairwise Maxwell divisors, and let \(D_{xyz}\) be the exceptional
triple-Maxwell divisor.  Define

\[
 \mathcal B_6^{\rm Max}
 =
 \sqrt[2]{(B_6,D_{xy})}
 \mathop{\times}_{B_6}
 \sqrt[2]{(B_6,D_{xz})}
 \mathop{\times}_{B_6}
 \sqrt[2]{(B_6,D_{yz})}
 \mathop{\times}_{B_6}
 \sqrt[2]{(B_6,D_{xyz})}.                          \tag{6.42}
\]

This is an algebraic Deligne--Mumford stack, not merely a complex of
groups.  At the generic point of any one Maxwell divisor its completed
chart is

\[
 [\operatorname{Spec}k[[s]]/\mu_2],
 \qquad \tau=s^2,\qquad s\longmapsto-s,             \tag{6.43}
\]

which is exactly the diagonal square-tail automorphism in (6.41).  The
three pairwise strict transforms are mutually disjoint after blowing up the
triple point; each meets \(D_{xyz}\) once.  Hence a pair--triple crossing
has inertia \((\mu_2)^2\), one factor for each target node.
This is the standard divisor root construction of Cadman
([arXiv:math/0312349](https://arxiv.org/abs/math/0312349)).

No radial root is added in (6.42).  There the branch value already has the
form \(\lambda_i=x_i^2u_i\), so the root-scale graph carries the saturated
source smoothing parameter.  Moreover (6.40) shows that source labels kill
the would-be diagonal stabilizer.  Along a Maxwell divisor, by contrast,
the target smoothing is first order in a branch-value difference and
(6.43) supplies the missing square root.

The \(S_3\)-action permutes \(D_{xy},D_{xz},D_{yz}\) and fixes
\(D_{xyz}\), so (6.42) descends equivariantly.  The four radial orbit types
still have sizes \(1,3,3,6\), but all have trivial labelled node inertia.
The separate quotient forgetting the order inside each source-root pair
has stabilizer

\[
 (S_2)^3\rtimes S_3,                                \tag{6.44}
\]

and should not be conflated with the Maxwell root factors.

There is no extra codimension-two inertia.  At a pair--triple crossing, five
local square-tail/connector flips satisfy three independent compatibility
equations, leaving \((\mu_2)^2\).  At a radial--Maxwell crossing the
full-chain radial condition kills the radial sign and the two Maxwell signs
agree, leaving one \(\mu_2\).  These exhaust the boundary crossings of the
surface.

The local algebraization ambiguity is therefore gone.  The selected ACV
component is smooth and tame in characteristic zero, has trivial generic
stabilizer, and has coarse space \(B_6\) by (6.38).  Its ramification divisor
and all crossing stabilizers are exactly those of (6.42).  The standard
bottom-up characterization of smooth tame Deligne--Mumford stacks
([arXiv:1503.05478](https://arxiv.org/abs/1503.05478)) therefore gives

\[
 \mathcal G^{\rm ACV}_{6,(2,2,2)}
 \cong \mathcal B_6^{\rm Max}.                     \tag{6.45}
\]

Thus the labelled coarse and stack comparisons are complete on the
`(2,2,2)` chart.  The full \(S_6\)-action is canonical on the global graph;
the recursive resonance atlas constructs and contracts the source flag
divisors on every other collision type.

## 7. Stable-target and coarse-source theorems

The recursive target construction is now a theorem.

> **Wonderful-pullback target theorem.**  
> On a simultaneous one-scale root-cluster chart, the tropical branch map
> begins with the integral-linear weights
> \(w_i=\mu_i v(x_i)\).  The first toroidal modification is the saturated
> pullback of the braid fan by \(v\mapsto w\).  On every equality stratum,
> take the initial critical-value residues.  Whenever a subset of those
> residues collides, repeat the construction with their differences.  The
> normalization of the full stable-target graph is the resulting recursive
> nested-set modification: equivalently, the normalized wonderful pullback
> of the labelled branch-diagonal building set.

The word “recursive” is essential.  Valuations of the individual branch
values determine the first scale; valuations and residues of their
differences determine subsequent scales.  Combinatorially this replaces the
radial caterpillar subfan by the relevant pullback of the full tropical
\(\overline M_{0,r+2}\) tree fan.  Algebraically it replaces a purely
monomial toric modification by a wonderful/nested-set sequence whose later
centers include strict transforms of binomial resonance loci.
The complete construction and its permutation-equivariant building-set
certificate are in
[The stable-target graph as a wonderful
pullback](BRANCH_GRAPH_WONDERFUL_PULLBACK.md).

The coarse source comparison is also a theorem: the normalized polynomial
closure in the finite admissible-cover pullback is finite birational over
the normal wonderful graph, hence has that graph as its coarse space.  By
[source-vertex rigidity](SOURCE_VERTEX_RIGIDITY.md), it can be stated
without an open-ended Hurwitz enumeration:

> **Explicit polynomial flag-atlas theorem.**  
> Over the wonderful stable-target graph, the labelled polynomial root
> degeneration has a unique dominant source tree with weighted inverse-image
> divisors for every target flag, compatible with nested-set edge
> contractions.  The two-fiber formula then reconstructs the unique
> source-vertex maps.  After labelled-node saturation, the resulting
> normalized tame stack is the polynomial ACV graph.

The [general radial source theorem](GENERAL_RADIAL_SOURCE_ATLAS.md) proves
this theorem on the complete weighted-braid face poset for arbitrary
cluster multiplicities.  The monodromy-forest theorem determines the source
combinatorics on the additional simple-branch residue-resonance centers.
Monodromy centralizers compute their generic inertia, and the recursive
resonance theorem supplies their higher-codimension algebraic flag-divisor
placement, descent, and simultaneous characters.

## 8. Consequences for H1-COARSE, H1-STACK, and H2

The labelled coarse and stack target/source comparisons are complete in the
degree-six `(2,2,2)` chart.  In general `H1-COARSE`, `H1-STACK`, H2, and
coarse H3 are complete.  The two H1 claims remain separate because the
all-degree stack atlas has substantially greater proof obligations than the
finite-normalization comparison.

The logarithmic comparison object is built as follows.

1. Start on labelled root-cluster atlases.
2. Form the normalized wonderful pullback of all branch-diagonal ideals;
   its first layer is the saturated weighted braid subdivision and its later
   layers resolve equal residues and internal caustic/Maxwell strata.
3. Inside the finite pullback of the fully marked admissible-cover stack,
   normalize the closure of the generic polynomial cover.  Its coarse space
   is the graph from step 2 by finite normalization.
4. Use the general connector/local-tail/identity rule on radial strata;
   use monodromy subforests for resonance source trees and node indices;
   construct their actual flag-divisor positions, then reconstruct every
   source-vertex map from two fibers and one third-fiber equation.
5. Add the root indices forced by source-node smoothing after subtracting
   the roots already present in the saturated root-scale graph; in the
   sextic chart this is precisely (6.42).
6. Use the canonical \(S_n\)-action on the global labelled graph closure;
   no separately chosen label-gluing cocycle is needed.
7. Apply the formal subgroup quotient and normalization/Stein contraction
   to the finite selected-root incidence used by H2.  This step does not
   wait for the explicit atlas in steps 4--5.

The degree-six result rules out a proof that identifies the full normalized
graph using only the walls \(\mu_i v(x_i)=\mu_jv(x_j)\).  It supports those
walls as the universal first layer, identifies the first missing centers,
and constructs their complete stable-target blowup.  The radial admissible
source atlas and all its saturated node monoids are explicit on all thirteen
ordered scale types.  The pairwise and triple Maxwell source-node rings, and
every intersection with the radial boundary, are also explicit.  There every
cluster tail begins with

\[
 T=u_iX_i(X_i-1)
\]

and has expansion factor \(2\) at its attaching node.  Accordingly the
target and source smoothing parameters should satisfy
\(\tau_{\rm target}=\tau_{\rm source}^2\), and equations
(6.16)--(6.20) and (6.26)--(6.29) verify the resulting normalizations.
The labelled coarse extension (6.38) selects the polynomial
admissible-cover component.  Equations (6.39)--(6.41) separate
normalization multiplicity from inertia, and (6.42) constructs the
Maxwell-boundary root stack; (6.45) identifies it with the labelled ACV
graph.  The local \(S_3\) and \((S_2)^3\rtimes S_3\) descent data is
explicit.  The general subgroup-quotient theorem identifies the pair-block
construction with the corrected H2 marked/unmarked quotient over the global
labelled graph.  The target fan and label descent are therefore complete in
general.  Two-fiber rigidity also removes any residual vertexwise
Hurwitz-class choice once the source flag divisors are known, and the
general radial atlas constructs those divisors and maps on all ordered
first-scale types.  Monodromy subforests also construct the source dual
graphs and node indices on nested simple-branch resonance types.  The
recursive resonance atlas completes the algebraic flag-divisor construction,
contraction compatibility, and simultaneous inertia characters at every
nested resonance intersection.

## 9. Newton-boundary interface

The output is already in a form usable by the plane Newton-boundary
programme:

- the branch-scale exponent rows are
  \((2,0,0),(0,2,0),(0,0,2)\);
- the ratio rows are their pairwise differences;
- the first resonance polynomials are the three binomials in (5.1);
- the triple-resonance second layer is the pair of linear forms
  \(6A-4B+1\) and \(2B-6C+1\).

This suggests a compiler interface based on successive initial forms:
monomial Newton data chooses a radial cone, then binomial initial ideals
choose nested resonance strata.  The degree-six artifact records exactly
those inputs without claiming that the resulting plane-boundary
compactification has yet been constructed.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_h1_h2_comparison_obstruction.py
.venv/bin/python scripts/verify_branch_scale_fan.py
.venv/bin/python scripts/verify_degree_six_branch_target_graph.py
.venv/bin/python scripts/verify_degree_six_admissible_equal_scale.py
.venv/bin/python scripts/verify_degree_six_admissible_radial_atlas.py
.venv/bin/python scripts/verify_degree_six_admissible_maxwell_atlas.py
.venv/bin/python scripts/verify_degree_six_central_hurwitz_selection.py
.venv/bin/python scripts/verify_degree_six_stack_inertia.py
.venv/bin/python scripts/verify_degree_six_stacky_fan_descent.py
```

The first command checks the degree-five constants and the second checks:

- all three degree-six local critical-value scales;
- completeness of the moving critical-value list;
- six eliminant Newton polygons;
- all rays, cones, and unimodularity determinants of the radial fan;
- the three Maxwell initial forms; and
- the varying second-scale cross-ratio at fixed radial data.

The third command independently checks:

- all intersections of the six-line target arrangement;
- the four-point Kapranov blowup and its ten boundary curves;
- the Petersen boundary-incidence graph;
- the four reduced preimages of the diagonal point; and
- the exact factorization \(\lambda_i=x_i^2u_i\) on the critical chart.

The fourth command checks the central and tail Riemann--Hurwitz data, the
degree-six admissible special fiber, the exact simple-ramification squares at
the diagonal bubble, and the four-branch normalization of the source-node
deformation ring.

The fifth command enumerates all thirteen radial ordered scale types and
checks the degree-six source over every target bubble, Riemann--Hurwitz on
every source component, all node-index partitions, and all Kummer
normalization counts.

The sixth command checks the three pairwise and one triple Maxwell
source-node rings, their two- and four-branch normalizations, all six
radial--Maxwell intersection points, and the exceptional directions over the
three coordinate target points.

The seventh command enumerates the central \(S_6\) branch cycles, detects the
two ambient Hurwitz classes, factors the exact square-cubic branch invariant,
and proves that the labelled source-root cross-ratio \(a=3\) selects the
polynomial class as a reduced local branch.

The eighth command separates normalization branches from genuine
label-preserving cover inertia.  It proves that all thirteen
equal-multiplicity degree-six radial inertias are trivial and that pairwise
and triple Maxwell branches each retain one diagonal \(\mu_2\).  Unequal
multiplicities obey the separate full-chain product formula in the recursive
resonance atlas.

The ninth command constructs the four-divisor Maxwell root complex, checks
all pair--triple face inclusions and the \(S_3\)-action, confirms trivial
inertia on the four radial quotient orbits, verifies the pair--triple and
radial--Maxwell crossing groups, and keeps the
\((S_2)^3\rtimes S_3\) pair-root quotient separate.
