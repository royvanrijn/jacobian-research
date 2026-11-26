# Meridional generation through covers of one link component

Let `L` be a tame link in the three-sphere. Use `b(L)` for bridge
number and `μ(L)` for meridional rank. The number of generators on an
individual component can obstruct a proposed small meridional generating
set even when that component and all the other components are unknotted.
For example, if `K_r` is the connected sum of `r ≥ 1` nontrivial
knots, its two-component untwisted Bing double satisfies

```text
b(B(K_r)) ≥ μ(B(K_r)) ≥ ceil((4r + 6)/3).
```

The lower bound uses successive double branched covers. Each proper
nonempty sublink of `B(K_r)` is an unknot, so the separate meridional
ranks of the proper sublinks give only the bound two.

## Component counts and plain spheres

For indexed components `K_i` of `L`, let `m_i` denote the number of
meridians of `K_i` in a generating tuple. For every nonempty index set
`I`, write `L_I` for the corresponding sublink. Meridional Dehn filling
on the deleted components induces an epimorphism from the link group
to the group of `L_I`. Under it a meridian of a deleted component
becomes the identity and each retained meridian remains a meridian of
the same component. Therefore

```text
Σ_{i∈I} m_i ≥ μ(L_I).
```

All these inequalities hold for the same tuple. In particular, the
minimum total of any nonnegative integer counts satisfying the
inequalities from known proper-sublink lower bounds is a lower bound
for `μ(L)`. A tuple with one meridian per component forces every
proper sublink to have meridional rank equal to its number of
components.

Blair, Kjuchukova and Pfaff's Theorem 1.1 says that the plain sphere
number, minimized over diagrams of a link in the three-sphere, is
`b(L)`. In Definition 2.2 each move uses a simple transverse plane
loop avoiding crossing neighborhoods. The loop meets the new strand
exactly once and all other intersections are on previously colored
strands. A counterexample to the meridional rank conjecture would thus
require a generating tuple of size `k` and a proof that every diagram
needs more than `k` plain-sphere seeds. A bridge lower bound greater
than `k` supplies the latter proof. A computation on one diagram
only constrains that diagram. [Published Theorem 1.1 and Definition 2.2](representations/sources/bkp-plain-sphere-published-2026.pdf)

## A bound after branching over one component

Write `L = U ∪ K_1 ∪ ⋯ ∪ K_s`, where `U` is unknotted, and suppose a
meridional generating tuple has `u` entries on `U` and `m_i` entries on
`K_i`. Fix `d ≥ 2` and take the cyclic `d`-sheeted cover branched along
`U` alone. Its ambient manifold is again the three-sphere. For a sublink
`P` of the preimage of `L \ U`, let `e_i` be the sum of the covering
degrees of the chosen components over `K_i`. One chosen component over
`K_i` has degree `d / gcd(d, lk(U,K_i))`, since the preferred longitude
of `K_i` acts on the sheets by translation by `lk(U,K_i)` modulo `d`
and its meridian acts trivially.

The group of `P` has a generating tuple with at most

```text
t = (d − 1)(u − 1)             additional arbitrary generators,
a_i = e_i m_i                 meridians over K_i,
```

where each chosen component over `K_i` of degree `e` receives at most
`e m_i` meridians. In particular, if `u = 1` and `P` consists of
one degree-one lift of `K_i`, then `μ(P) ≤ m_i`.

**Proof.** Present the link group on the proposed tuple. Its
presentation complex has one vertex and one edge per generator.
Take the cover of this complex associated to the homomorphism to
`Z/d` given by linking with `U`. Its fundamental group is that of
the corresponding exterior cover. Each positive meridian on `U` has
one cycle of `d` edges in the lifted graph. An inverse gives the
same undirected cycle. The loop around this cycle is a lifted
meridional filling slope for the branched completion over `U`.
Conjugation of the original meridian changes the lifted base path
and can permute the sheets, leaving this slope unchanged up to
conjugacy and inversion. Thus one can impose the filling relation
and eliminate one edge from each such cycle. The cycles for
different generators use disjoint sets of edges. For a single
`U` meridian, its remaining `d − 1` edges form a spanning tree.
Collapsing this tree leaves `(d − 1)(u − 1)` other edges from
meridians on `U`.

Every generator on `K_i` lifts to a loop at each sheet. With its
base path lifted, the loop is a meridian of the component of the
preimage reached by that path. Each component of degree `e` receives
exactly `e` of these loops per original generator. Collapsing the
tree gives based meridians from these loops. Filling the preimage
components outside `P` kills their meridians and gives the claimed
generating tuple. Remaining lifted relators or filling relators
can only reduce the required counts. This is the partial-filling
version of the graph argument for [arbitrary permutation branched covers](representations/GENERAL_BRANCHED_COVER_RANK_BOUNDS.tex).

**A second cover.** Suppose a link group in the three-sphere has a
presentation on `a` specified meridians and `t` other generators.
Take its connected double branched cover over the whole link.
Each of the `a` meridians has a two-edge cycle whose filling relation
eliminates one edge. The other generators supply at most `2t` edges.
The graph stays connected because each removed edge had an alternative
path within its own cycle. After collapsing one tree edge, the
fundamental group of the completed double branched cover therefore
has ordinary rank at most

```text
a + 2t − 1.
```

Its first homology over any field satisfies the same upper bound.
This step treats the additional generators as arbitrary.

**Bridge positions.** Suppose a bridge sphere has one upper bridge
on `U` and `b_i` upper bridges on `K_i`. Each side of the sphere
is a ball with one boundary-parallel branch arc. Its cyclic branched
cover is a ball, as seen by modeling the arc as the axis of rotation
in a disk times an interval. The original bridge disks for the
other arcs can be chosen disjoint from the branch arc, so those
disks lift to disjoint bridge disks. A chosen component of degree
`e` over `K_i` uses `e b_i` lifted upper bridges. Deleting the
other lifts gives

```text
b(P) ≤ Σ_i e_i b_i.
```

This compares the algebraic and geometric counts for the same
cover when the branch component uses one meridional generator or
one bridge, respectively.

For a split union of `U` with a knot `K`, a ball containing `K`
lifts to `d` disjoint balls. Choosing all lifts gives `d` split
copies of `K`, each of degree one. The sublink inequalities applied
to the copies give `μ(P) ≥ d μ(K)`, matching the first-cover bound
for a tuple with `u = 1` and `m_1 = μ(K)`.

## Two unknotted components with growing meridional rank

For a knot `K`, let `B(K)` denote its untwisted Bing double and
`K^r` the knot with its string orientation reversed. Both components
of `B(K)` are unknotted. In the double cover branched along either
component alone, each of the two lifts of the other component has
degree one and is isotopic, up to orientation, to `J = K # K^r`.
This is Livingston and Van Cott's Covering Theorem A in the two-leaf
case. Their proof identifies a degree-one lift by placing the
unbranched component in a solid torus that lifts to two copies, and
the deck involution exchanges the two copies. Either depth-one leaf
can be marked. [Livingston and Van Cott, Theorem 4 and Figure 4](https://arxiv.org/pdf/1009.3164v1)

Let `g = rank π_1(Σ₂(J))`, where `Σ₂` denotes the double branched
cover. Suppose some generating tuple for `B(K)` has `u` meridians on
one component and `v` on the other. Branch over the first component
and retain just one lift of the second. The resulting knot group
has `a ≤ v` specified meridians and `t ≤ u − 1` additional generators.
Apply the second-cover bound, then exchange the roles of the components
to obtain

```text
g ≤ min(v + 2u − 3, u + 2v − 3).
```

Thus `u + v ≥ ceil(2(g + 3)/3)`. Also, any tuple with `u = 1`
must have `v ≥ μ(J)` by the first-cover bound alone.

The preimage of a connected-sum sphere meeting a knot in two
points is a sphere.
Capping its two sides by the double branched cover of a ball with
one boundary-parallel arc gives the covers of the two summands.
Consequently `Σ₂(J) = Σ₂(K) # Σ₂(K^r)`. Reversing knot orientation
does not change the double branched cover. Van Kampen gives the
free product of two copies of `π_1(Σ₂(K))`, so Grushko's theorem
on rank additivity for free products gives `g = 2R`, where
`R = rank π_1(Σ₂(K))`. It follows that

```text
μ(B(K)) ≥ ceil((4R + 6)/3).
```

Grushko obtains rank additivity by putting a generating tuple
of a free product into the factors by a change of free generating
basis. [Grushko, §2, p. 170](https://www.mathnet.ru/eng/sm6023)

For every nontrivial `K`, Grove and Wilking's Theorem C ensures
`R ≥ 1`, which yields `μ(B(K)) ≥ 4`. For a connected sum of
`r` nontrivial knots, the same branched-cover connected-sum
decomposition and rank additivity give `R ≥ r`.
[Grove and Wilking, Theorem C](https://arxiv.org/pdf/1304.4827v1)

A homological version uses `h = dim_F H_1(Σ₂(K); F)`, for any
field `F`. Since `R ≥ h`, it gives

```text
μ(B(K)) ≥ ceil((4h + 6)/3).
```

For a trefoil, the meridional presentation
`⟨a,b | aba = bab⟩` becomes `⟨a,b | a² = b² = (ab)³ = 1⟩`
after imposing the two meridian squares. The even subgroup, which
is the fundamental group of its double branched cover, is generated
by `ab` and has order three, as seen by mapping `a,b` to the two
adjacent transpositions in `S_3`. Hence the connected sum of `r`
trefoils has `h = r` over `F_3`, giving the claimed unbounded lower
bound.

These are lower bounds for unrestricted meridional generating tuples
and hence for bridge number. A separation requires both an explicit
unrestricted generating tuple and a strictly larger bridge lower
bound on the same link. For a proposed three-meridian generating
tuple of a two-component link of unknots, the only possible counts
are `(1, 2)` and `(2, 1)`. To rule out either count, branch over the
singly represented component and exhibit a degree-one lift of the
other component whose meridional rank is at least three. The
Bing-double argument excludes both counts for every nontrivial
companion.
