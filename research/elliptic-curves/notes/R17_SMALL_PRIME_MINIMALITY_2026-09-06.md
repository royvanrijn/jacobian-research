# Universal small-prime minimality for the six compact R17 families

For every primitive integer parameter pair, each displayed integral short model
is minimal at every prime **5 through 997 except 13**, whenever its fibre is
nonsingular. The [separate 13-adic classification](R17_INTEGRAL_13_PARAMETER_CHARTS_2026-09-06.md)
gives its unique removable scale and proves that reduction remains bad there.
This closes a model-normalization gap; it supplies no new point or rank bound.

The [exact certificate](../../artifacts/generated-results/elliptic-curves/r17_small_prime_minimality_v1.json)
accounts for all 990 family/prime pairs. The integer Sylvester adjugate in the
[resultant proof](R17_SCALING_PRIME_SUPPORT_2026-09-06.md) implies that a removable
scale requires resultant valuation at least four. This excludes 964 pairs.

For the remaining 26 pairs, complete projective residue trees test
A_h = 0 modulo p^4 and B_h = 0 modulo p^6. Affine parameters and the infinity
chart are both included. Every tree becomes empty by depth three; the largest
modulus is 47^3 = 103823. An independent checker enumerates each full residue
ring at every recorded depth, using direct modular polynomial evaluation.
All 26 candidate pairs are excluded; no branch remains UNKNOWN or censored.

For primes at least five, c4 = -48 A_h and c6 = -864 B_h make these coefficient
conditions necessary and sufficient for nonminimality of the short model.
Their universal failure proves the stated local minimality theorem.
Primes 2 and 3, all primes above 997, exact conductors and ranks remain outside
this result. In particular, the large unfactored resultant cofactors still
prevent a complete global scaling-prime classification.

## Replay

```sh
python3 elliptic-curves/cas/report_r17_small_prime_minimality.py --check
python3 elliptic-curves/cas/verify_r17_other_small_prime_scalings.py
```

The aggregate recomputes the exact resultants and complete residue trees and
runs the independent residue checker. It binds the passed independent Sage
resultant verification as well. All local checks pass. All22 isolated supplement stages now pass; see the replay completion below. The finite residue protocol
allows at most depth six, 4096 live residues per chart and 200000 candidates per
level; the computation closed well inside those bounds without point searches.

## Portable replay completion

All22 stages pass in the [standalone supplement replay](../../artifacts/generated-results/elliptic-curves/compact192_followup_supplement_portable_replay_v2.json), covering five exact histories, ten point-cloud proofs, rational geometry, the follow-up summary and five small-prime proof checks. The [7005716-byte archive](../../artifacts/generated-results/elliptic-curves/compact192_followup_supplement_evidence_v2.zip) needs no base archive. Its frozen notes precede this completion update. Version one omitted a command-invoked cloud checker; that failed package and replay are preserved. Version two includes the checker and all recorded cloud-proof sources, with unchanged arithmetic inputs.
