# Six higher-parameter masked controls recover their withheld directions

All six fixed controls recover the withheld generic direction through an
exactly verified integer combination of the original 17 independent points.
All 72 point-search boxes and their rational-map replays completed. None of
the searches directly recovered either sign of the original withheld basis
representative. The [relation certificate](../../artifacts/generated-results/elliptic-curves/higher_r17_masked_relations_v1.json)
contains each returned point and its exact 17-coordinate group word.

This tests a specific gap left by the
[paired score and visibility audit](PAIRED_PRODUCT_SCORE_2026-09-06.md): known
generic directions can remain recoverable on larger-coefficient fibres even
when the original representatives are outside the completed boxes. It does
not prove adequate sensitivity to exceptional directions or explain the
low ranks detected in the paired candidate cohort.

The frozen selector takes the smallest retained index in each family of the
already completed higher-parameter cohort, without filtering its point result:

| Control | Completed boxes | Returned points up to sign | Withheld coefficient in exact recovery word |
| --- | ---: | ---: | ---: |
| `074d9-010` | 12 | 104 | -1 |
| `07ca9-001` | 12 | 101 | -1 |
| `08234-005` | 12 | 89 | -1 |
| `08f72-011` | 12 | 77 | 1 |
| `103b2-159` | 12 | 76 | -1 |
| `11952-057` | 12 | 64 | -1 |

Section zero is withheld uniformly. The preparer separates its coordinates
and the original independence proof from the sixteen-point worker input.
Only the retained principal metric block selects the twelve deepest of 256
fixed sampled parity classes, using the earlier ordinary-control geometry.
These are sampled classes in a rounded specialized metric, not an exhaustive
maximum stratum or a proof of canonical-height optimality. Exact PARI quartic
transforms then supply the search coordinates. Every box has height 125,000
and a ten-second cap; there is one worker, no replacement and no adaptive wave.

Oracles are opened only after all six searches and map replays terminate.
Finite reductions modulo 3, 5, 7 and 11 propose small rational coefficient
words for at most the first 128 returned points per curve. An exact rational
group-law equality verifies each reported recovery. In these six cases the
successful words are integral with nonzero withheld coefficient. Since the
original 17 points are separately certified independent, each witness lies
outside the rational span of the retained sixteen points. No inference of
rational dependence is made merely from a finite congruence.

Preparation took 15.091 seconds, the six searches and their replays 83.948
seconds, and relation construction/replay 1.509/1.047 seconds on the recorded
machine. The endpoint is recovery of known directions under a changed
reference subgroup, not a new rank result. Six successful controls do not
provide a calibrated population success rate, saturation or an upper bound.

`higher_r17_masked_controls.py replay --index i` checks a retained attempt for
`i=0,...,5`. `audit_higher_r17_masked_relations.sage --check` independently
checks the original independence proofs and exact recovery words, without
performing a new point search or repeating the finite word proposals.
