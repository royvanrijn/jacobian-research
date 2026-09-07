# Latent Mordell--Weil lattice calibration

## Outcome

The concise method and target closeout is
[`LATENT_LATTICE_REVERSE_ENGINEERING_REPORT.md`](LATENT_LATTICE_REVERSE_ENGINEERING_REPORT.md).
After a second independently frozen failure on E29 and ICARM 398--400, this
research line is parked pending new record fibres.  The control technique is
retained; neither target failure is promoted to a subgroup nonexistence
theorem.

### Frozen E29 and ICARM 398--400 run (2026-09-01)

The unchanged secondary tag is
`LATENT-LATTICE-E29-398-400-FROZEN-2026-09-01-v1`, with manifest hash
`8795cdd203ba1c698e0f0534c14a45a91c37a6f4d405e795a9f2f295f86bfcba`.
E29 failed closed because dimensions 19 and 20 were absent from at least one
persistence cutoff.  ICARM 398 hit the unchanged 600-second PARI resource
bound at the first height candidate.  Independent invocations selected 12 on
ICARM 399 and 16 on ICARM 400.  Hence the required four-of-four recurrence
failed and all component and later stages remained gated off.  The consolidated
artifact hash is
`e89f6156a22f7a5b762e6c525db7b95b4d625ac1a96d712ca8f4053aa05f62cb`.

### Frozen wgxli run (2026-09-01)

The user subsequently authorized one cautious target run after an immutable
method freeze.  The content-addressed tag is
`LATENT-LATTICE-WGXLI-FROZEN-2026-09-01-v1`; its manifest hash is
`ef6f8b7be7a14095efa7529fb795d237e06465ba1cec023dcb4845287609c9f4`.
The manifest fixes dimensions `10..20`, adaptive cloud bounds, all seeds and
beam widths, the three-cutoff persistence score, finite-code ensembles, the
component sampler, equal structural-channel weights, and fivefold hold-out
rules.  It explicitly forbids target-informed tuning under the tag.

The independent dimension stage used every displayed public generator and
returned:

| ICARM fibre | ambient rank | bound | rays | frozen dimension result |
|---|---:|---:|---:|---:|
| 351 | 25 | 38 | 1,806 | 10 |
| 356 | 29 | 58 | 2,655 | fail: dimensions 19/20 missing at one or more cutoffs |
| 376 | 22 | 30 | 2,054 | 10 |
| 377 | 23 | 34 | 3,300 | 13 |
| 385 | 29 | 62 | 2,099 | fail: dimensions 18/19/20 missing at every cutoff |

Thus no dimension recurs in at least four fibres.  The frozen result is
`FAIL_FROZEN_DIMENSION_RECURRENCE`, with artifact hash
`c5ca3a8746ce2d889ad0f9ff0f1eb9015f320ae0f820147da966ff0b45b5732e`.
By the precommitted protocol, relation-component matching, primitive-Hermite
matching, finite-index matching, held-out component prediction, and equation
interpolation were not run.  In particular this run returns no candidate
abstract wgxli lattice and no embeddings.

This is a bounded negative result for the tagged proposal/persistence method,
not a theorem that the five curves have no common primitive subgroup of
dimension `10..20`.  The two failed fibres reflect proposal omissions inside
the fixed beam.  The values 10 and 13 are statistical selector outputs, not
proved generic ranks.

The bounded Phase-0 graph-walk and dimension controls pass.  Exact graph-walk consensus
recovers the five fixed-dimension control spaces, and the cross-dimension
persistence audit described below changes the Fermigier estimate from 13 to
12 at both independent height bounds.  A separate center-free component
calibration also recovers the exact primitive rank-16 R17 component in the
rank-25 control and completes it to rank 17.  The two later target applications
described above produce no candidate generic lattice or target embedding.

The earlier failed selectors remain bounded negative results, not
nonexistence results for a common wgxli lattice.  The new control pass also
does not satisfy target success level C from the research brief: no exhaustive
search of all primitive target sublattices of dimensions 10 through 20 has
been performed.

The original `v1` selector artifact is retained as a historical first pass.
It predates two correctness fixes: non-unit identities `a +/- b = m*c` are no
longer collapsed to unit ternary relations, and proposal indices now refer to
the caller's height-ordered record list rather than canonical vertex order.
The active replay artifact is `v2`; `v1` must not be refreshed in place.

## Exact withheld truth

The calibration-truth artifact recovers coordinates by a 120--150 digit
height-dual solve and then accepts them only after exact elliptic group-law
replay.  It contains:

- primitive rank-17 embeddings of the published R17 sections in the public
  rank-at-least 25, 26, 27, and 28 subgroups;
- the exact rank-12 Fermigier--Mestre embedding in ICARM 245;
- the primitive closure of that rank-12 rational space in the displayed
  rank-20 group; and
- the Smith factors of every embedding.

The curve-245 generic subgroup has index `2^11` in its primitive closure.
This distinction matters: a selector of primitive rational spaces cannot by
itself recover the actual specialized generic subgroup.

## Blind calibration result

Every cloud below is a complete enumeration of primitive unoriented vectors
through the displayed height bound in the full public independent subgroup.
Relations are the complete unoriented hypergraph of visible
`a +/- b = c` triples.  The color-refinement digest is coordinate-free but is
only an isomorphism invariant, not a complete canonical-labelling theorem.

| control | ambient rank | bound | lines | selected/truth intersection |
|---|---:|---:|---:|---:|
| R17 rank-at-least 25 | 25 | 40 | 2,155 | 14/17 |
| R17 rank-at-least 26 | 26 | 43 | 1,921 | 17/17 |
| R17 rank-at-least 27 | 27 | 52 | 2,313 | 17/17 |
| R17 rank-at-least 28 | 28 | 60 | 2,423 | 17/17 |

Thus three of four positive controls are recovered exactly by the active v2
selector.

On ICARM 245, the corrected v2 maximum-integrality-likelihood scan selects
dimension 13, so it fails the requirement to recover approximately dimension
12.  Its selected space intersects the withheld Fermigier space in dimension
6.  A second rank-12 search grown from 3,000 exact additive hyperedges has
intersection dimension 5.  For comparison, the withheld true space contains
144 retained lines and 112 integral lines; recognizing its statistics after
supplying truth does not constitute blind recovery.

Exact finite `E(F_p)/2E(F_p)` and `E(F_p)/3E(F_p)` codes and exact component
codes at the declared multiplicative places `2,5,13,19,37` were added to the
curve-245 complex.  They do not repair the failed subspace proposal stage.

## Post-v1 high-recall audit

The corrected library also contains a bounded enclosure/core-extension
proposal channel.  This materially improves recall but does not yet supply a
valid joint selector:

- all four R17 controls contain the exact rank-17 truth in both the
  arithmetic-priority and relation-only 3,000-seed ledgers; in the
  arithmetic channel it is the final (lowest-scoring) proposal in every
  fibre, showing that the old score direction was wrong for generic sections;
- ICARM 245 improves from a 6/12 selected intersection to a direct 11/12
  proposal; the exact primitive Fermigier space occurs once in a 5,385-entry
  rank-12 extension ledger and at rank 65 after a two-cutoff enclosure-
  intersection/arithmetic/relation score;
- ICARM 282 reaches 11/12 at direct blind rank 1 and contains the exact
  Fermigier space once at refined blind rank 34 of 4,904; and
- the `u=28917/20` sibling has 68 retained truth rays spanning rank 12, but
  the direct 3,000-seed channels reach only 8/12.  Truth-containing pairwise
  enclosures first occur at rank 18, not rank 15.  Exhaustive coordinate-
  subset height, relation-count, normalized-Gram, and LLL two-generator-shell
  controls do not select the truth (best ranks 491 before LLL and 781 after
  LLL among the filtered coordinate proposals).

These are bounded diagnostic observations, not a new successful calibration
artifact.  They establish high recall on R17 and substantial improvement on
the Fermigier controls, while also showing that single-fibre support,
integrality, truncated theta data, and elementary reduced-Gram signatures are
still inadequate selectors.  Therefore the target gate remains closed.

## Finite-aware calibration

The reusable finite layer now separates two roles that the earlier pass had
conflated:

- a **source-local proposal key**, which uses actual quotient classes only to
  branch from rare finite fingerprints inside one fibre; and
- a **source-free candidate signature**, which retains candidate image ranks,
  unoriented class multiplicities, cyclic element orders, and induced
  unit/scaled relation types.  Its digest forgets public point labels,
  quotient bases, component orientations, reduction-prime names, and the
  ordering of equal-type blocks.

For each control fibre, the finite calibration uses the first three
one-dimensional quotient blocks for each of `ell=2,3` as development data and
the next three as an untouched validation set.  All reduction primes are at
most 251.  Finite-priority proposal generation contains the exact R17 space
in all four positive controls, at blind ranks 1792, 1666, 1227, and 1067.  It
reaches maximum truth intersections 11/12 on ICARM 245, 11/12 on ICARM 282,
and 8/12 at `u=28917/20`.  It therefore passes R17 proposal recall but fails
the required Fermigier calibration.

Finite profile matching is not a selector by itself.  In the explicit
leave-rank-25-out diagnostic, all rank-25 proposals are compared with the
rank-26--28 R17 development profiles.  The true rank-25 R17 space improves
from source rank 1792 to finite-profile rank 188, but a false candidate is
still selected.  Disjoint held-out blocks do not turn this into recovery.

One necessary nuisance separation emerged from the negative controls.  On
the `u=28917/20` fibre all six sampled mod-2 quotient maps vanish on the known
rank-12 subgroup, although they do not vanish on the full displayed rank-20
subgroup.  This is exact fibre-specific divisibility/saturation information,
not an abstract height-lattice invariant.  The artifact consequently reports
both strict profiles and profiles conditioned on active quotient blocks; it
does not match raw finite classes or silently discard inactive blocks.

The finite-aware artifact remains `FAIL_FINITE_PROPOSAL_RECALL`.  It proves
the exact finite calculations within the declared ensembles and a bounded
failure of this generator/selector.  It neither proves that finite codes are
useless in a later joint method nor that a common generic lattice is absent.

## Cross-bound finite-aware shape calibration

The next control-only replay separates proposal recall from selection.  On
ICARM 245 it independently constructs rank-15 enclosure ledgers from the
complete height-28 and height-29 clouds, retains the leading 200 enclosures
at each bound, and tests all 40,000 pairs.  Two finite-field annihilator keys
are rejection filters only.  Every survivor is regrouped by its exact
rational annihilator, and every retained candidate is primitively saturated.

The fixed score

```text
arithmetic LLR + 0.1 * induced ternary relations
               + 2 * exact cross-bound occurrence count
```

gives 3,799 surviving pairs and 2,939 distinct exact rank-12 spaces.  There
are no two-prime collisions in this run.  The exact primitive Fermigier
subgroup occurs twice, has 144 retained rays and 535 induced ternary
relations, and ranks 65th.  Thus the earlier blind v1 dimension scan selects
12 rather than a fake forced rank-17 core, and the exact truth is moved into a
bounded top-128 ledger.  Rank 65 is materially better proposal recall, but it
is not blind recovery.

For the four R17 controls, finite-seeded ledgers contain the exact truth at
source ranks 1792, 1666, 1227, and 1067.  A scale-free cloud-height shortlist
followed by the intrinsic Hermite statistic recovers rank 25 exactly when the
rank-26--28 truth lattices are the training controls.  The symmetric
leave-one-out experiment, however, selects truth in only one of four fibres:
the other truth cloud ranks are 1266, 482, and 200, outside the top-64 Hermite
stage.  An exact additive complex on an intrinsic complete shortest shell was
also tested; all four true R17 specializations have different shell digests
at the declared 128-vector minimum, so that invariant is too brittle under
specialization.

Accordingly `latent_lattice_shape_calibration_v1.json` has status
`PASS_PROPOSAL_CALIBRATION_SELECTOR_FAIL`.  Exact R17 recall, a held-out R17
recovery, blind dimension 12, and top-128 Fermigier recall pass.  Symmetric
joint selection does not.  The wgxli gate remains closed.

### Post-artifact recurrence and relation-consensus audits

A third cutoff does not repair the selector.  The height-27 ledger has 3,281
rank-15 enclosures, with truth-containing enclosures only at blind ranks 490,
1456, and 2012.  Under the declared adaptive bounds 512 x 200, 512 x 200,
and 200 x 200, the exact rank-12 truth ranks 184 on the 27--28 pair and 65 on
the 28--29 pair, and is absent on 27--29.  Among all pre-closure exact keys it
ranks 165 under the lexicographic multi-bound recurrence score; 138 false
cores recur in all three pair ledgers.  Strict cutoff recurrence is therefore
rejected within these bounds.

The Fermigier height comparison also exposed an essential saturation
distinction.  On ICARM 245 the published generic rank-12 subgroup has Smith
factors `2,2,...,2,1` and index `2^11` in its primitive closure; on ICARM 282
the Smith factors are `12,2,...,2`.  A primitive intersection candidate must
be compared with the closure Gram, not the stored generic-subgroup Gram.  The
corrected primitive-shell profile still fails (the ICARM 245 truth ranks 209
among the first 256 two-bound candidates), so it is retained only as a
diagnostic.

The exact coefficient-relation benchmark is substantially stronger.
After expressing every retained truth ray in the corresponding generic basis,
157 rays occur in at least three controls and span rank 17.  In every
leave-one-fibre-out split, the training two-of-three core has 113--137 rays
visible on the held-out fibre, again spanning rank 17.  The replayed
`latent_lattice_relation_consensus_v1.json` therefore has status
`PASS_CONTROL_EXACT_RELATION_SIGNALS`.  Rational-ray normalization now clears
denominators before canonicalization and records those denominators as
finite-index metadata.  This is necessary on ICARM 245, where all retained
primitive-closure rays have denominator 2 relative to the displayed generic
subgroup, and on ICARM 282, where denominators 2, 4, 6, and 12 occur.  The
aligned ICARM 282 / `u=28917/20` sibling pair has 30 common rational rays
spanning rank 12.  ICARM 245 is a different Fermigier--Mestre family and is
not asserted to share their labelled rays.

The relation benchmark remains deliberately supervised: published embeddings
align the coefficient systems.  The exact hypergraph validator independently
lifts the active held-out cores to primitive rectangular matrices of shapes
`17 x 25`, `17 x 26`, `17 x 27`, and `17 x 28`.  Its Smith factors are all
one, and exact full-cloud replay finds respectively 238, 266, 291, and 304
training-core rays, each spanning rank 17.  Thus cutoff loss does not destroy
the signal once an injection is supplied.

The current blind unequal-cloud proposal generator still fails.  It combines
relation-star seeds, robust scale-free height-angle pruning, a bounded beam,
codimension-one metric reseeding, exact rectangular lifting, and exact global
replay.  In the frozen rank-25 box it tests 256 center pairs, retains a beam of
500 for 80 steps, expands 39,714 states, and makes 500 exact-lift attempts.
It reaches rank 17 but no candidate passes the 100-ray global-replay gate; the
maximum is 49 versus the supervised truth value 238.  Its artifact status is
`FAIL_BLIND_R17_RECOVERY_GATE_CLOSED`.

This failure is informative.  Exact relations plus height angles can produce
primitive false rank-17 subspaces; a sampled false core had only
9-dimensional intersection with R17.  Disjoint finite quotient profiles
separate that false core (`0.1530` development distance) from truth
(`0.04583`), but finite scoring cannot select a truth candidate absent from
the proposal ledger.  The next selector must match whole relation components
jointly or score exact partial-subspace replay before full rank.  No wider
edgewise beam is authorized by this calibration.

## Exact proper-subspace replay calibration

The partial-subspace alternative is now implemented without completing a
proper component to a square ambient basis.  Given a partial ray injection,
the code takes the primitive closure of its source span, rewrites every mapped
ray in intrinsic integer coordinates, lifts the induced relation injection,
and replays every supplied source ray lying in that rational subspace.  The
rectangular target matrix is Smith-tested.  The same matrix can be restricted
to good-reduction quotient and bad-component codes through the existing
source-free finite signature.

On the supervised rank-25 R17 control, the deterministic truth path after 287
relation-frontier steps has 103 rays and rank 16.  Exact saturation and lift
find two global-sign choices.  The source hyperplane contains 362 supplied
rays; 194 replay into the held-out target cloud, span rank 16, and support 318
target ternary relations.  The target Smith factors are all one.  Its finite
signature uses the disjoint three-development plus three-held-out quotient
blocks and retains 279 target-cloud rays in the target rational subspace.

This strong validator does not yet repair proposal selection.  In the frozen
oracle-center search, all 128 retained seed lineages receive a descendant
slot, and at most 400 proper subspaces of ranks at least 10 are audited before
beam truncation.  The run expands 44,697 states, but its best audited partial
map replays only 29 rays and 22 relations; its best full-rank lift replays 47
rays, and no embedding crosses the 100-ray gate.  Therefore
`latent_lattice_partial_replay_v1.json` has status
`PASS_EXACT_PARTIAL_REPLAY_SELECTOR_FAIL`.

What this changes is the failure diagnosis: a large correct rank-16 component
is exactly distinguishable once proposed, but one- or two-vertex edgewise
continuations lose it even when the correct center pair is supplied.  The next
bounded generator must propose whole components or maintain richer
within-seed branch diversity, with finite codes used to prune those branches.
Simply widening the old global beam is not justified.  No wgxli record was
loaded.

## Finite-aware center-star proposal

The first whole-component proposal unit assigns a center together with several
incident ternary relations.  Each matched arm adds both noncentral rays at
once; candidate arms are height-ordered, but every accepted partial map is
checked for equal total ranks over `F_2` and `F_3`.  This is a necessary exact
condition for the primitive integral embeddings sought here.  The beam is
stratified by matched-edge count and rational rank so cutoff-heavy components
are not discarded solely for having fewer visible arms.

On the supervised rank-25 control, the complete visible truth star at source
vertex 116 and target vertex 15 has 10 incident relations, 21 rays, and rank
11.  Its exact primitive lift replays 32 of the 40 source-hyperplane rays and
26 ternary relations; the associated six-block finite signature retains 35
target rays.  Thus a correct star is exactly certifiable once supplied.

The bounded proposal ledger processes 20 star layers with a beam of 512,
expands 9,282 states, rejects 4,152 branches by mod-2/mod-3 rank, and makes 500
exact partial-replay audits.  It returns 476 distinct target subspaces, but the
largest postselection intersection with the rank-11 truth star is only 9.
Checking every newly introduced pair and triple over both finite fields takes
151 seconds in the same box, rejects no additional branch, and leaves the same
9/11 ceiling.  It is retained as an explicit expensive option, not the
default.

Accordingly `latent_lattice_star_component_v1.json` has status
`FAIL_STAR_COMPONENT_RECALL_GATE_CLOSED`.  The next generator needs a
multi-center component or a joint cross-fibre invariant capable of ranking
partial stars before local beam loss.  Neither a wider edgewise beam nor the
current single-center star is calibrated for target use.

## Primitive-shape intersection consensus

The next invariant compares whole candidate subspaces rather than individual
rays or stars.  A truth-free replay ledger stores every finite-seeded R17
proposal, its exact primitive embedding matrix, its induced finite relation
fingerprint, and the 80-digit ambient height Gram.  All 5,887 proposal
matrices are Smith-audited; every one is primitive.  The primitive closure is
nevertheless recomputed before restricting the height form, so the same code
also handles nonsaturated future proposals correctly.

For each fibre, the numerical prefilter retains the 64 largest scale-free
Hermite invariants.  Every one of the `64 choose 2 = 2016` pairs is then
ranked exactly over `Q`.  A candidate receives the equal-weight sum of:

```text
within-pool normalized Hermite extremality
  + fraction of pooled candidates meeting it in dimension at least k-1.
```

This fixed score selects the exact withheld R17 subgroup in all four positive
controls.  The selected source indices are 1792, 1666, 1227, and 1067, and
the pinned artifact includes the four exact primitive `17 x r` embedding
matrices.  It also expresses each selected basis in the published section
basis, Smith-certifies the four changes as unimodular, and transports the
exact pinned R17 Gram into every selected basis.  Thus the positive control
recovers the subgroup and its abstract lattice, rather than merely a
17-dimensional rational span.  The result is materially stronger than
nearest-neighbour matching
of the compact relation fingerprints, primitive Hermite extremality alone,
theta profiles, or intrinsic-shell relation profiles; each of those
intermediate selectors fails and remains recorded as such.

The same invariant does not complete the Fermigier calibration.  The
truth-free ICARM-245 replay fixes the earlier height-28/29, top-200 by
top-200, 40,000-pair bound and serializes all 128 primitive rank-12
candidates.  Thus later selector experiments do not rerun the expensive
elliptic arithmetic.  The earlier dimension scan selects 12, not the forced
rank-17 core.  Intersection consensus moves the exact Fermigier subspace from
zero-based proposal rank 65 to zero-based rank 12 in the fixed top-64 pool,
but selects candidate 11, whose intersection with truth has dimension 9.
This is improved identification and correct dimension calibration, not exact
recovery.

Finite codes remain audit data rather than a positive identity score.  When
the ICARM-245 candidates are compared with the two other Fermigier truth
profiles using source-free active-block matching, the exact truth ranks 113
of 128 on development blocks and 73 of 128 on disjoint held-out blocks
(zero-based ranks).  Folding that adverse signal into the selector would be
post-hoc overfitting, so it is explicitly excluded.

At this intermediate stage the positive artifact had status
`PASS_R17_EXACT_INTERSECTION_CONSENSUS`, while the negative-control artifact
has status `PARTIAL_FERMIGIER_DIMENSION_PASS_SELECTOR_GATE_CLOSED`.  Phase 0
was not yet complete.  The exact graph-walk refinement below supersedes that
selector result; no wgxli curve is loaded in either calculation.

## Exact graph-walk Phase-0 calibration

The calibrated refinement keeps the same primitive Hermite prefilter and
exact codimension-one candidate graph, but replaces degree alone by two exact
graph statistics.  For adjacency matrix `A`, every candidate receives:

- its exact triangle count; and
- its exact length-four walk count, the corresponding entry of `A^4 1`.

Ties are assigned average rank percentiles.  If the two largest Hermite
invariants differ by at least `0.005`, the separated extremum is selected.
Otherwise the score is

```text
1.5 * (triangle percentile + length-four-walk percentile)
    + Hermite percentile.
```

The default pool contains the 64 largest primitive Hermite candidates.  The
selector recovers all four exact R17 subgroups and the exact primitive
rank-12 Fermigier rational space.  The selected R17 bases are transported by
exact unimodular matrices to the published section bases and carry the exact
pinned R17 Gram.  The selected Fermigier candidate is source index 65, the
withheld primitive closure itself, rather than merely a high-overlap
candidate.

This result is stable throughout the declared calibration box:

- pool sizes `64,80,96,112,128`;
- Hermite-gap thresholds `0.004,0.005,0.006,0.007`; and
- graph weights `1.25,1.5,1.75`.

All 300 conditional fixed-dimension control selections in that box are exact.
Recomputing all ambient
height Grams and primitive Hermite signatures at 120 digits also leaves all
five default selections unchanged.  Thus
`latent_lattice_graph_walk_calibration_v1.json` has status
`PASS_PHASE0_GRAPH_WALK_CONTROL_CALIBRATION_WITH_DIMENSION_WINDOW`.

The pass has two important boundaries.  The graph, its rational intersection
edges, triangle counts, walk counts, embedding matrices, and R17 transports
are exact within the pinned proposal ledgers.  Hermite ordering, percentile
weights, gap threshold, and extrapolation to an unseen family remain
heuristic.  Moreover the Fermigier selector recovers the correct primitive
rational space but does not blindly infer which index-`2^11` sublattice is
generated by the published generic sections.  That index remains a withheld
postselection audit.  The corrected blind scan estimates 13 rather than
uniquely 12, so the graph result by itself is explicitly conditional on
retaining the rank-12 ledger from that window.  The cross-dimension audit
below closes this bounded control gap without changing the graph artifact.

## Cross-dimension persistence calibration

The dimension audit replaces the active rank-12 proposal by the truth-free
graph-walk selection and compares its local dimension profile at complete
height-28 and height-29 clouds.  It includes every distinct rank-11 face
obtained by intersecting that core with the 128 replay candidates and every
distinct rank-13 extension obtained by adjoining one retained short ray.  It
then maximizes the threshold-free three-level bottleneck

```text
persistence(k) = min(best integrality LLR at k-1, k, and k+1).
```

Raw integrality enrichment still prefers a false rank-13 extension: its LLR
is `24.20898` at height 28 and `40.25603` at height 29.  Persistence instead
selects rank 12 uniquely at both bounds.  At height 28 its score is `18.36234`
versus at most `7.43459` for every other dimension; at height 29 it is
`34.03847` versus at most `11.36279`.  The selected rational space equals the
withheld Fermigier primitive closure in both postselection audits.  The two
clouds contain 1,928 and 2,618 rays, with respectively 293 and 351 distinct
rank-13 retained-ray extensions.

Accordingly `latent_lattice_dimension_persistence_v1.json` has status
`PASS_FERMIGIER_CROSS_DIMENSION_PERSISTENCE`.  Candidate coordinates, cloud
membership, rational faces and retained-ray extensions, integrality counts,
and additive relation counts are exact in the declared ledger.  The LLR and
its use as a dimension score are statistical.  The rank-13 exhaustion is over
one-retained-ray extensions, not all rational rank-13 superlattices, and the
published generic subgroup's index `2^11` is still not inferred blindly.  No
unseen-family claim follows from this control calibration.

## Center-free relation-component calibration

The next positive-control calculation replaces oracle centers by whole dense
components.  Its generic source proposal is obtained from 2,000 deterministic
rank-16 samples of the 472-ray R17 training complex.  Exact full-complex replay
selects a unique largest hyperplane with 362 rays and 1,521 induced ternary
relations; the next support sizes are 348, 337, and 337.  It has rank 16 over
`Q`, `F_2`, and `F_3`.  No center label is supplied.

On the held-out rank-25 fibre, the first 16 height-prefiltered rank-17 parent
proposals are each split internally into the lowest three quarters of their
rays for component development and the remaining quarter for replay.  Four
dense rank-16 components per parent are retained from 400 deterministic
samples.  Equality of `F_2` and `F_3` ranks, a 17-quantile absolute
height-angle profile, and primitive Hermite distance are rejection filters.
Among the four survivors, selection is lexicographic by replayed held-out rays
per added rank and then total exact replay; local relation count is not in the
selection objective.

The selected component replays 71 held-out rays and 279 rays in the full
cloud.  It is the exact rational image of the blindly selected 362-ray source
hyperplane and has Smith factors `1^16` in the public rank-25 subgroup.  Its
parent is proposal 1792.  Adding the one missing rank gives a primitive
rank-17 embedding with Smith factors `1^17`; its exact basis change to the
published R17 specialization is also unimodular.  Published coordinates and
the truth proposal index are used only in this postselection audit.

The default numerical slacks are `0.0021` for the angle-profile distance and
`0.007` for primitive Hermite distance.  All nine combinations of angle slack
`0.0021,0.0023,0.0025` and Hermite slack `0.0065,0.007,0.008` select the same
exact parent.  Two narrower development experiments are retained in the
research history rather than presented as passes: rank-2/3 per-arm replay had
no discrimination, and Hermite slack `0.005` excluded the exact component.

Thus `latent_lattice_relation_component_calibration_v1.json` has status
`PASS_CENTER_FREE_R17_RANK16_COMPONENT_AND_RANK17_COMPLETION`.  The exact
claims are conditional on the frozen 16-parent, four-component, sampled box.
The parent prefilter, deterministic sampling, height filters, and thresholds
remain heuristic.  The source coefficient complex itself comes from
supervised positive-control alignment, so this is a calibrated matcher, not
an unseen-family recovery.  No wgxli fibre was loaded.

## What is proved and what is heuristic

Exact within the recorded bounds:

- vector enumeration coordinates and primitivity;
- rational-span saturation and intersection dimensions;
- additive relations;
- rational point arithmetic, integrality, and coordinate complexity;
- finite-reduction quotient codes;
- declared multiplicative component codes, including pair-sum replay; and
- all withheld point identities and embedding Smith factors.

Numerical or heuristic:

- canonical heights and all scores derived from them;
- beam and relation-seeded subspace selection;
- color refinement as a proxy for hypergraph matching; and
- any interpretation as a generic family or height lattice.

No equation interpolation, displayed-label sign/permutation search,
unrestricted `GL(17,Z)` search, target-family identification, K3 assumption,
or target Gram reconstruction was performed.

## Reproduction

```sh
PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 elliptic-curves/cas/build_latent_lattice_calibration_truth.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 elliptic-curves/cas/calibrate_latent_lattice_method.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 elliptic-curves/cas/calibrate_finite_aware_latent_lattice.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 elliptic-curves/cas/calibrate_latent_lattice_shape.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 elliptic-curves/cas/calibrate_latent_lattice_relation_consensus.py --check

sage -python \
  elliptic-curves/cas/calibrate_latent_lattice_hypergraph_matcher.sage --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_metric_relation_search.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_partial_replay.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_star_component.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_joint_fingerprints.py --check

PYTHONPATH=elliptic-curves \
  python3 elliptic-curves/cas/calibrate_latent_lattice_joint_shape.py --check

PYTHONPATH=elliptic-curves \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_intersection_consensus.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 elliptic-curves/cas/build_latent_lattice_fermigier_replay.py --check

PYTHONPATH=elliptic-curves \
  python3 elliptic-curves/cas/calibrate_latent_lattice_fermigier_replay.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_graph_walk_consensus.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_dimension_persistence.py --check

PYTHONPATH=elliptic-curves \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_relation_components.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 -m unittest elliptic-curves/tests/test_latent_lattice.py -v
```

The two pinned outputs are
[`latent_lattice_calibration_truth_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_calibration_truth_v1.json)
and
[`latent_lattice_calibration_v2.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_calibration_v2.json).
The finite-aware replay is
[`latent_lattice_finite_calibration_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_finite_calibration_v1.json).
The cross-bound shape replay is
[`latent_lattice_shape_calibration_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_shape_calibration_v1.json).
The exact supervised relation benchmark is
[`latent_lattice_relation_consensus_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_relation_consensus_v1.json).
The exact rectangular validator is
[`latent_lattice_hypergraph_matcher_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_hypergraph_matcher_v1.json).
The bounded blind unequal-cloud failure is
[`latent_lattice_metric_relation_search_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_metric_relation_search_v1.json).
The exact proper-subspace validator and bounded selector failure is
[`latent_lattice_partial_replay_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_partial_replay_v1.json).
The finite-aware center-star proposal failure is
[`latent_lattice_star_component_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_star_component_v1.json).
The truth-free R17 proposal and compact-fingerprint replay is
[`latent_lattice_joint_fingerprint_ledger_v1.json.gz`](../../artifacts/generated-results/elliptic-curves/latent_lattice_joint_fingerprint_ledger_v1.json.gz),
with its failed nearest-neighbour summary in
[`latent_lattice_joint_fingerprints_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_joint_fingerprints_v1.json).
The primitive shape replay and its failed one-statistic summary are
[`latent_lattice_joint_shape_ledger_v1.json.gz`](../../artifacts/generated-results/elliptic-curves/latent_lattice_joint_shape_ledger_v1.json.gz)
and
[`latent_lattice_joint_shape_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_joint_shape_v1.json).
The exact positive-control consensus is
[`latent_lattice_intersection_consensus_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_intersection_consensus_v1.json).
The truth-free Fermigier replay and its partial calibration are
[`latent_lattice_fermigier_replay_v1.json.gz`](../../artifacts/generated-results/elliptic-curves/latent_lattice_fermigier_replay_v1.json.gz)
and
[`latent_lattice_fermigier_consensus_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_fermigier_consensus_v1.json).
The current Phase-0 passing control artifact is
[`latent_lattice_graph_walk_calibration_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_graph_walk_calibration_v1.json).
The bounded cross-dimension Fermigier pass is
[`latent_lattice_dimension_persistence_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_dimension_persistence_v1.json).
The center-free rank-16 component and rank-17 completion calibration is
[`latent_lattice_relation_component_calibration_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_relation_component_calibration_v1.json).
The superseded `v1` bytes remain available for provenance but are not the
active replay target.
