# F7 cubic-cube coefficient obstruction

**GAOF7CUBE1:** no Gao F7 eliminant has a cubed cubic divisor, including
repeated or shared roots. Consequently sum floor(m_r/3)<=2; at most two
distinct roots have multiplicity at least three, and every multiplicity is
at most eight by this obstruction alone. The later
[GAOF7MULT1 proof](../gao-f7-linear-cube/HIGH_ROOT_PROOF.md) sharpens the
maximum to five and proves that it occurs.

This excludes three nonsimple multiplicity lists for a three-point fibre:
(9), (6,3), (3,3,3). GAOF7MULT1 also excludes(7,2). The later
[GAOF7FIVE1 proof](../gao-f7-quintuple/PROOF.md) classifies all118 normalized
fivefold targets and excludes(5,4),(5,2,2). Only(4,3,2),(3,2,2,2) and
global three-point existence remain **OPEN**.

[Proof and full scope](PROOF.md), [761-check receipt](receipt-integer.json),
[streaming exact verifier](verify_rank_integer.py),
[compressed certificate, about2.08GB](cubic-six-rank-dag.json.gz),
[coefficient input](cubic-six-input.json),
[tested dependencies](requirements-replay.txt).

The verifier checks734 rational predecessor nodes after exact integer
denominator clearing, then a full402x402 integer coefficient matrix
containing every actual monomial, including1. Its nonzero determinant
modulo1000000007 proves invertibility over Q. The replay passed in724.77
seconds; it runs no Groebner-basis search or rational reconstruction.
The [direct-rational alternative](verify_rank.py) reached its600-second cap
after last reporting700 nodes and has no complete PASS receipt.

The result has no independent end-to-end review, formal verification or
established novelty priority and does not settle JC2. The
[next proof](../gao-f7-linear-cube/README.md) uses quadratic cancellation
and exact resultant interpolation rather than the superseded four-variable
necessary-equation route.
