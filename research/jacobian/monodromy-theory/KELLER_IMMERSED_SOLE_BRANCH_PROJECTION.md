# Projection restrictions for an immersed sole Keller branch

Let \(F\) be a complex plane Keller map of generic sheet degree \(N>1\),
with full finite normalization \(\pi:W\to\mathbb A^2_\mathbb C\) in its
source field. Suppose that its reduced ramified target curve \(B\) is
irreducible and that every affine analytic branch of \(B\) is smooth.
For a polynomial target coordinate \(p\), let \(w(p)\) be the degree of
its restriction to \(B^\nu\simeq\mathbb A^1\).

**Theorem KellerImmersedSoleBranchCompositeWidth.** The width \(w(p)\)
is composite for every polynomial target coordinate \(p\).

By
[KellerOrevkovInvariantObstruction](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#local-degree-and-the-finite-branch-quotient),
the Orevkov invariant of \(B\) is nonabelian. Choose a companion \(q\)
forming a polynomial target coordinate system. Write \(p(t),q(t)\)
for the restrictions on the normalization. Each has degree at least
two. Finiteness of the projection excludes a constant restriction,
and a degree-one restriction makes \(B\) a polynomial graph, contrary
to Chau's exclusion of a nonproperness component isomorphic to the
affine line.

Put \(n=\deg p(t)\). When the leading degree of \(q(t)\) is a multiple
of \(n\), replace \(q\) by \(q-cp^j\) for the scalar \(c\) and positive
integer \(j\) canceling that term. The replacement is a polynomial
target automorphism. It strictly decreases the degree of the
restriction of \(q\), and its new degree is still at least two by the
same graph argument. After these cancellations the companion has degree
\(m\) not divisible by \(n\).

If \(n\) were prime, \(n,m\) would be coprime. Swap the two coordinates
if needed so that \(x(t),y(t)\) have degrees \(M=\max(n,m)\) and
\(L=\min(n,m)\), respectively. A local parameter \(z\) at infinity
on the normalization can be chosen with \(x=z^{-M}\). Then
\(y=c' z^{-L}+\cdots\), with \(c'\ne0\), so the first exponent \(L/M\)
of \(y\) in \(x\) already has the full reduced denominator \(M\).
The knot at infinity has just one characteristic pair, \((L,M)\),
in the conventions of
[Neumann–Norbury, Section 4, equation (2)](sources/neumann-norbury-orevkov-invariant.pdf).
Their terminal number is \(d_1=L+M>1\) by equation (7).
Since this is a rational immersed component with one point at
infinity, their Lemma 5.4 with
\(d_1=-K\cdot\overline B-1\) shows that its affine complement group
is abelian, as in their proof of Theorem 2. In particular its
Orevkov invariant is abelian. A polynomial target automorphism
carries affine local complement groups to affine local complement
groups, so this contradicts the original obstruction.

## Genus of slices with immersed index-three ramification

**Theorem KellerImmersedTripleRamificationGenus.** Suppose in addition
that \(W\) is smooth and that every finite ramified source germ of
\(\pi\) is a cyclic power map of index three along its immersed
image branch. All remaining source germs are unramified. The
compactified generic \(p\)-slice has genus at least two in every
polynomial target coordinate.

The
[smooth-normalization polar argument](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#compactification-over-the-p-line)
gives a projective smooth-surface model \(f:S\to\mathbb A^1_p\) over
the finite \(p\)-line. Its fibres are reduced and connected.
Its polar boundary consists of disjoint sections, and \(W\) is the
complement of these sections. The smooth immersed index-three local
power map gives the fibre germ \(t^m+v^3=0\) at a critical point of
multiplicity \(m\ge2\) of \(p\) on the branch normalization, where
\(m\) is the local degree of that projection. Everywhere else the
finite fibre is smooth, including the points at the poles. This
description concerns each finite source germ separately and applies
when two image branches meet. The genus-zero case is already excluded
by the non-normal-crossings germ in the smooth-normalization argument.

Assume that the generic slice genus is one. Flatness gives arithmetic
genus one on every finite fibre. For a fibre \(C\) with \(k\)
components and normalization genera \(g_i\), let \(\delta_x\)
be the delta invariant of a singular point \(x\) and \(b_x\) the
number of its local analytic branches. The normalization sequence gives
\[
 1=\sum_{i=1}^k g_i + \sum_x\delta_x-k+1 .
\]
Connectivity gives \(\sum_x(b_x-1)\ge k-1\), since a singular point
joins at most \(b_x\) components and can reduce their number of
connected groups by at most \(b_x-1\). Consequently
\(\sum_x(\delta_x-b_x+1)\le1\).

For \(t^m+v^3=0\), this summand is one when \(m=2\) or \(m=3\),
and is at least three otherwise. In fact, if \(3\nmid m\) the one
branch has coprime characteristic orders \(3,m\), giving
\(\delta_x=m-1\). If \(m=3s\) the three smooth branches meet pairwise
with multiplicity \(s\), giving \(\delta_x=3s=m\). Thus a singular
fibre has exactly one singularity. A cusp has one analytic branch
and cannot connect distinct components, so its fibre is integral.
With an ordinary triple point, the genus formula gives \(k\ge3\)
and connectivity gives \(k\le3\). The fibre consists of three
smooth rational curves meeting pairwise transversely at the point.

On a smooth or integral cuspidal fibre the dualizing bundle is
trivial by duality and degree zero. In the reducible case, a
component \(C_i\) meets the other two components with total
intersection number two. Since the whole fibre has intersection zero
with \(C_i\), its self-intersection is \(-2\). Adjunction gives
degree zero for the dualizing bundle on \(C_i\simeq\mathbb P^1\).
Duality gives a nonzero section of the fibre dualizing bundle,
since the arithmetic genus is one. On each component such a
section is identically zero or nowhere zero. The restrictions have
the same value at a point shared by the components, so this nonzero
section is nowhere zero on the whole fibre.

The relative dualizing line bundle is therefore trivial on all finite
fibres. Cohomology and base change and the evaluation map pull it
back from a line bundle on \(\mathbb A^1_p\). This line bundle is
trivial. Since the base has trivial canonical bundle, restriction
to \(W\) gives \([K_W]=0\), contrary to the
[ramification class calculation](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#canonical-class-on-the-finite-affine-normalization).

## Six-sheet consequence

**Corollary SixSheetImmersedAlternatingCompositeWidth.** For a six-sheet
plane Keller map with monodromy \(A_6\) and two unramified omitted
primes, the sole ramified curve has composite width at least six in
every polynomial target coordinate, and every compactified generic
slice has genus at least two. If a coordinate has width six, its
slice has \((g_p,r_p)=(2,4)\) or \((3,2)\).

The
[six-sheet cusp classification](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#consequences-in-generic-sheet-degree-six)
gives smooth \(W\), and zero excess on the sole index-three prime
gives the cyclic local power maps at every point of that prime.
All other finite source germs are étale by purity. The
composite-width theorem excludes width five, and the
[prior polar bound](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#six-sheet-consequences)
gives width at least five. For width six the directional formula
\(2w(p)=6+2g_p+r_p-2\) gives \(2g_p+r_p=8\). The genus bound
and \(r_p\ge1\) give the two listed pairs.

**Corollary SixSheetImmersedTranspositionFoldDegree.** In the
six-sheet sole-transposition profile, if the ramified branch is
immersed, its width is at least eight in every polynomial target
coordinate and its total degree is at least nine in every target
frame. If this profile also has an exceptional \((2,1,5)\) fold, both
longitudinal degrees are at least \(48\) and at least one is at least
\(54\) in the original normalized fold chart.

The
[six-sheet smooth-normalization bound](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#six-sheet-consequences)
makes the directional width at least seven. Composite width
raises it to eight. Canceling the highest parameter power by a
nonzero linear combination of the two coordinates in a fixed
frame strictly reduces the width below the original total degree,
which must therefore be at least nine. In the exceptional-fold
case, the coordinatewise bound in that original frame gives
\(\deg_a P\ge6\,w(P)\) and \(\deg_a Q\ge6\,w(Q)\), where the
widths here are the restriction degrees of the respective target
coordinates. The immersed hypothesis holds, in particular, when
the sole transposition and three unramified omitted primes
saturate the ramification budget.

The [degree-pair refinement](KELLER_IMMERSED_FOLD_DEGREE_PAIR.md)
raises this immersed-slice bound to 60 and describes the equality case.
