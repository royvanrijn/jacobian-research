# Latent Mordell--Weil lattice calibration

## Outcome

The calibration gate fails.  Consequently the method was **not** applied to
ICARM 351, 356, 376, 377, or 385, and it produces no candidate generic lattice
or target embedding.

This is a bounded negative result about the selector, not a nonexistence result
for a common wgxli lattice.  In particular it does not satisfy success level C
from the research brief: no exhaustive search of all primitive sublattices of
dimensions 10 through 20 has been performed.

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
The superseded `v1` bytes remain available for provenance but are not the
active replay target.
