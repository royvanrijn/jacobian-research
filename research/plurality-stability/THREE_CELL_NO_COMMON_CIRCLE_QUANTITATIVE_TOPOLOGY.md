# Quantitative topology without a common winner circle

Let \(A\) be an actual planar Gaussian noise-stability maximizer with
three positive prescribed masses at \(0<\rho<1\). Assume no centered
circle meets all three open strict-winner sets. Then some label has
a uniformly positive amount of Gaussian mass outside every one of
its connected components. Some label also misses a centered annulus
and has positive Gaussian mass on both sides of it. The two masses
and the annular width are bounded below on compact interior
parameter sets.

If all three labels are unbounded, the number of bounded connected
winner components inside \(B_R\) grows at least logarithmically.
The proof uses the positive radial overlap between labels in the
bounded region and the contraction of occupied phases at infinity.

## Fragmentation and an annulus that one label misses

Write \(O_i\) for the strict-winner sets and
\(\mathcal R_i=\{s>0:O_i\cap\{|x|=s\}\ne\varnothing\}\)
for their open radial images. They give the ordinary cells
up to Gaussian null sets. Let \(p_i(s)\) denote the normalized
angular mass of label \(i\) at radius \(s\), and write
\(d\mu(s)=s e^{-s^2/2}ds\). Set
\[
 \delta_{ij}=\int \min(p_i,p_j)\,d\mu,\qquad
 \delta=\min_{i<j}\delta_{ij}.
\]
Radial overlap of actual maximizers gives \(\delta>0\),
uniformly bounded below on compact positive-correlation bands
and compact sets of positive mass triples.
[Radial overlap and its uniform lower bound](THREE_CELL_RADIAL_FIBER_EXCLUSIVITY.md#a-stability-bound-from-the-radial-overlap)

Let \(C_i\) be a largest-Gaussian-mass connected component of
\(O_i\). Every component is open and has positive mass, so there
are at most countably many components and a largest one exists.
Put \(\eta_i=\gamma_2(A_i\setminus C_i)\).
The radial overlap of \(C_i,C_j\) is at least
\(\delta_{ij}-\eta_i-\eta_j\), since deleting a slice part
decreases the minimum by at most its deleted angular mass.
If all three overlaps of these components were positive,
their interval radial images would have pairwise nonempty
intersection. Intervals with pairwise nonempty intersection
have a common point, which would be a common winner radius.
Consequently some pair satisfies
\[
                  \eta_i+\eta_j\ge\delta_{ij}.
                                                               \tag{QC1}
\]
In particular some label has Gaussian mass at least
\(\delta/2\) outside every one of its connected components.

To find a radial separator, set \(\varepsilon=\delta/8\).
For each label take an interval \(I_i=[a_i,b_i]\) of radial
quantiles with Gaussian mass \(\varepsilon\) of that label on
each of the two exterior sides. These quantiles exist and are
finite and positive: the radial measure \(p_i\,d\mu\) has
no atoms and has total mass at least \(\delta\).
The interval endpoints satisfy \(a_i<b_i\) as the total mass
exceeds \(2\varepsilon\). The overlap measure
\(\min(p_i,p_j)d\mu\) loses at most \(4\varepsilon\) of its
mass upon restricting to \(I_i\cap I_j\). Thus every pair
intersection has \(\mu\)-mass at least \(\delta/2\).
The interval \(J=\bigcap_i I_i\) has endpoints realized by two
of the intervals \(I_i\), or is itself one of them contained
in the other two. In either case it is a pair intersection,
so \(\mu(J)\ge\delta/2\).

On \(\operatorname{int}J\), the three relatively closed sets
\(\{s:s\notin\mathcal R_i\}\) cover the interval. If each had
empty interior, their three open dense complements would
have nonempty intersection. Hence some label misses an open
subinterval of \(\operatorname{int}J\), and also misses any
smaller closed subinterval \([a,b]\) of positive length.
Since \([a,b]\subset I_i\) for every label, the missing
label \(i\) satisfies
\[
 O_i\cap\{a\le |x|\le b\}=\varnothing,\qquad
 \gamma_2(A_i\cap B_a)\ge\varepsilon,\qquad
 \gamma_2(A_i\setminus B_b)\ge\varepsilon .
                                                               \tag{QC2}
\]
This gives two separated groups of the same label without a
hypothesis on its junctions. It does not require either group
to be connected.

## A uniform annular width on compact parameter sets

Fix a compact positive-correlation band and a compact set of
positive mass triples. Choose a common lower bound
\(\delta_0>0\) for the pair overlaps and use
\(\varepsilon=\delta_0/8\) in the quantiles.
There is a width \(w>0\), depending only on these parameter
sets, such that (QC2) holds with \(b-a\ge w\) and both masses
at least \(\varepsilon/2\).

Suppose otherwise, and take a sequence of maximizers for which
all such annular separators have widths tending to zero.
Pass to limiting parameters and a strongly converging
sequence of ordinary maximizing partitions in the plane.
Center the prices by \(\sum_i\beta_i=0\). They are bounded
because each label wins on a positive-mass set where
\(0\le U_i\le1\). After taking a subsequence the prices converge,
and the adjusted scores converge locally uniformly to those
of the limiting partition. Its limiting prices give the
winner rule: a strict losing comparison persists on a
neighborhood and prevents positive limiting mass of the
losing label there. A common winner circle in the limit
would have three strict score gaps at three points of that
circle, persisting in the sequence. Thus the limit also has
no common winner circle.
[Strong maximizer compactness](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)
[Local regularity of maximizing strict winners](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Construct the core interval \(J\) for the limit. If a radius
in \(\operatorname{int}J\) meets two strict winners, both
strict inequalities persist at their respective fixed
angles on a smaller closed radial interval inside
\(\operatorname{int}J\). Shrink that interval so their
gaps have a positive minimum. The same two labels then
strictly win somewhere at each of its radii throughout
the sequence for all sufficiently large indices.
Absence of a common circle forces the remaining label
to miss that entire interval.

If no radius in \(\operatorname{int}J\) meets two strict
winners, choose one radius there with a strict winner.
Such a radius exists because the score tie sets have
Lebesgue measure zero. The same label strictly wins at
one fixed angle on a smaller open radial interval.
There is no other strict label at any radius in that
interval, so the corresponding open annulus belongs
almost everywhere to the same cell. A maximal tie there
is a dormant interior score tie. For each pair it is
isolated by (FC4), since the contrast has positive
definite Hessian at its zero. There are only finitely
many such ties on any smaller compact annulus:
the maximal-tie set there is relatively closed, and
an accumulation point would itself be a nonisolated
dormant interior tie. Take a closed radial subinterval
of positive length avoiding all their radii. The
same label strictly wins at every point of this
closed annulus, with a positive minimum gap by
compactness. This again persists in the sequence,
so either other label is absent.
[Dormant interior ties, FC4](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)

In both cases the closed annulus is fixed and lies
inside the core of every limiting label. Each label
therefore has at least \(\varepsilon\) limiting Gaussian
mass on each side. Strong convergence gives at least
\(\varepsilon/2\) on each side for every label in the
sequence at all sufficiently large indices. The fixed
positive annular width contradicts the choice of the
sequence.

## A logarithmic count when every label is unbounded

Suppose in addition that all three \(O_i\) are unbounded.
Set \(L=R_{\mathrm{source}}/(1-\rho)\), where
\(R_{\mathrm{source}}\) is the source radius in (PP2).
An unbounded label meets the radial interval
\[
                         [\rho a-L,a+L]
                                                               \tag{QC3}
\]
for every \(a>L/\rho\). Indeed start from one of its
strict winners \(p\) with \(|p|>a\), and take the first
integer \(n\ge1\) with \(\rho^n|p|\le a\), so that
\(\rho^n|p|>\rho a\). The iterated source estimate
provides a strict winner of the same label within \(L\)
of \(\rho^n p\). In particular all three labels meet
every block \((S,\rho^{-4}S)\) once \(S\ge S_0\).
To see this, use (QC3) with \(a=\rho^{-5/2}S\);
the additive errors fit inside the block for large \(S\).
The threshold \(S_0\) can be chosen uniformly on the
compact parameter sets above.
[Iterated contraction of occupied phases, PP2](THREE_CELL_PRESSURE_CONTRACTION_PROPAGATION.md#contractive-recentering-and-isolated-ends)

Take two consecutive blocks
\((S,\rho^{-4}S)\) and \((\rho^{-4}S,\rho^{-8}S)\)
and select one connected strict-winner component per
label meeting the first block. If two selected
components extended to radii beyond \(\rho^{-8}S\),
their interval radial images would cover the entire
second block. The remaining label meets the second
block, contradicting absence of a common circle.
Therefore at least two selected components, with
different labels, have outer radii in
\((S,\rho^{-8}S]\) and are bounded.

Applying this to the disjoint successive bands with
endpoints \(S_0\rho^{-8n}\) gives distinct components
on distinct bands because each has a single outer radius.
For \(R\ge S_0\) it proves the lower bounds
\[
 \begin{split}
 \#\{\text{bounded strict-winner components contained in }B_R\}
       &\ge 2\left\lfloor\frac{\log(R/S_0)}
                                 {8\log(1/\rho)}\right\rfloor,\\
 \#\{\text{such components from either of any two fixed labels}\}
       &\ge \left\lfloor\frac{\log(R/S_0)}
                                 {8\log(1/\rho)}\right\rfloor .
 \end{split}                                                   \tag{QC4}
\]
An open bounded component cannot contain a point at its
outer radius, so containment in the open \(B_R\) holds
also at an endpoint of a band.

In this situation at most one label has an unbounded
connected component. If two labels had such components,
their radial images would each contain a full radial tail.
The third label would be absent throughout that tail,
contrary to its unboundedness. The logarithmic count
therefore describes a necessary form of fragmentation in
the three-unbounded-label case. It establishes neither
an upper bound on the number of bounded components nor an
exclusion of this case.
