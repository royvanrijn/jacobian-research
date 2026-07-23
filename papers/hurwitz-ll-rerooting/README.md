# Hurwitz--LL compactification of rerooting

This paper contains the admissible-cover closure, rerooting quotient stack,
marked-zero-fiber LL degree, formal collision comparison, and the boundary
Picard and caustic/Maxwell calculations.

The missing global comparison is resolved negatively as stated by the
[global comparison obstruction](global-comparison-obstruction.tex).  In
degree five, the stable zero-root compactification forgets the target
cross-ratio with leading term `x^3/y^2`; the admissible-cover model resolves
it by the normalized blowup of `(x^3,y^2)`.  Thus there is no global morphism
from the former compactification to the latter.  The correct package is their
normalized graph correspondence followed by the normalization/Stein
contraction to the finite incidence.  The standalone
[DVR marking audit](dvr-marking-audit.tex) remains the canonical source for
the abstract finite-cover valuative theorem and quotient collision models.

The next local chart is computed in the
[branch-scale fan note](../../extended-geometry/BRANCH_SCALE_FAN.md).  The
general node-normalization/inertia formula and the corrected formal H2
quotient over the labelled graph are isolated in the
[labelled-node saturation note](../../extended-geometry/LABELLED_NODE_SATURATION.md).
[The stable-target graph as a wonderful
pullback](../../extended-geometry/BRANCH_GRAPH_WONDERFUL_PULLBACK.md)
constructs the target graph for arbitrary labelled branch collisions.  The
[source-vertex rigidity
theorem](../../extended-geometry/SOURCE_VERTEX_RIGIDITY.md) proves that,
once the source tree and weighted inverse-image divisors of target flags are
known, two fibers and one third-fiber equation uniquely reconstruct each
component map.  The
[general radial source
atlas](../../extended-geometry/GENERAL_RADIAL_SOURCE_ATLAS.md) constructs
those components and their saturated node partitions on every ordered
first-scale stratum for arbitrary multiplicities.  The
[polynomial monodromy-forest
theorem](../../extended-geometry/POLYNOMIAL_MONODROMY_FORESTS.md) determines
all nested simple-branch resonance source trees and node indices.  The
[source-graph finite-normalization
theorem](../../extended-geometry/SOURCE_GRAPH_FINITE_NORMALIZATION.md)
shows that the wonderful target graph is already the complete coarse
polynomial admissible-cover graph and makes corrected H2 and coarse H3
unconditional.  The
[monodromy-centralizer
theorem](../../extended-geometry/MONODROMY_INERTIA_CHARACTERS.md) proves
the full-chain radial inertia formula: equal cluster multiplicities give
trivial inertia, whereas unequal multiplicities can retain a finite
subgroup.  It also computes generic simple-resonance inertia from anchors
and node partitions.  The
[recursive resonance atlas
theorem](../../extended-geometry/RECURSIVE_RESONANCE_ATLAS.md) constructs
all source flag divisors by normalized initial forms in framed residue
screens, proves contraction descent, and extracts the full simultaneous
centralizer characters.  This completes corrected H1.  The
degree-six `(2,2,2)` boundary has the predicted six-cone radial braid fan, but
a triple leading-value resonance retains a second-scale cross-ratio.  Thus
the radial fan is only the first toric layer; the full logarithmic graph needs
the recursive nested resonance refinement described there.  In this chart
the target refinement is completely explicit: it is the pullback of the
four-point blowup of `P^2` giving `Mbar_0,5`, with four reduced
triple-Maxwell branches above the diagonal blowup center.  The source-cover
and saturated-node comparison is exact on the generic equal-scale face and
its first diagonal refinement.  The admissible component types and saturated
node monoids are also verified on all thirteen radial scale types; exact
source-node rings are verified on all pairwise and triple Maxwell strata and
their radial intersections.  Standard admissible-cover deformation theory
then supplies ambient toroidal gluing.  The central Hurwitz problem has two
ambient classes, but the labelled source-root cross-ratio selects the
polynomial class as a reduced branch; this identifies the labelled coarse
graphs in the degree-six chart.  The remaining stack discrepancy is also
explicit: labelled radial inertia is trivial, while each pairwise or triple
Maxwell target node retains one diagonal \(\mu_2\).  The matching candidate
is the iterated second-root stack along the four Maxwell boundary divisors;
its face maps and codimension-two inertia are explicit, and smooth
tame-stack reconstruction identifies it with the labelled ACV graph on this
chart.  The three-pair wreath-product quotient remains separate.  Gluing
under the full root-label action is canonical on the global labelled graph,
and the corrected H2 marked/unmarked quotient is formal there.  Corrected
H1, H2, and coarse H3 are therefore complete; the comparison with the
unmodified root-stable quotient remains false.

ACV and Deopurkar provide the compactification technology cited in the paper.
They do not externally review the repository's specialized LL-degree or
boundary-class calculations; no such external review is recorded.

Build from this directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```
