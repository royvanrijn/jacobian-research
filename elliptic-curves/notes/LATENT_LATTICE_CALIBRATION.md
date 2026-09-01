# Latent Mordell--Weil lattice calibration

## Outcome

The calibration gate fails.  Consequently the method was **not** applied to
ICARM 351, 356, 376, 377, or 385, and it produces no candidate generic lattice
or target embedding.

This is a bounded negative result about the selector, not a nonexistence result
for a common wgxli lattice.  In particular it does not satisfy success level C
from the research brief: no exhaustive search of all primitive sublattices of
dimensions 10 through 20 has been performed.

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
| R17 rank-at-least 25 | 25 | 40 | 2,155 | 13/17 |
| R17 rank-at-least 26 | 26 | 43 | 1,921 | 16/17 |
| R17 rank-at-least 27 | 27 | 52 | 2,313 | 15/17 |
| R17 rank-at-least 28 | 28 | 60 | 2,423 | 17/17 |

Thus only one of four positive controls is recovered exactly.

On ICARM 245, the maximum-integrality-likelihood dimension is 12, so the new
scan avoids the old forced-rank-17 error.  It nevertheless selects a rank-12
space whose intersection with the withheld Fermigier space has dimension only
6.  A second search grown from 3,000 exact additive hyperedges has intersection
dimension 5.  For comparison, the withheld true space contains 144 retained
lines, 112 integral lines, and has integrality likelihood-ratio statistic
18.3623; recognizing that statistic after supplying truth does not constitute
blind recovery.

Exact finite `E(F_p)/2E(F_p)` and `E(F_p)/3E(F_p)` codes and exact component
codes at the declared multiplicative places `2,5,13,19,37` were added to the
curve-245 complex.  They do not repair the failed subspace proposal stage.

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
  python3 -m unittest elliptic-curves/tests/test_latent_lattice.py -v
```

The two pinned outputs are
[`latent_lattice_calibration_truth_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_calibration_truth_v1.json)
and
[`latent_lattice_calibration_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_calibration_v1.json).

