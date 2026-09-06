# Exact resultant support for possible R17 model scalings

The six homogeneous coefficient resultants are computed exactly and independently
verified with Sage. Trial division through997 leaves an unfactored cofactor
of1657through1780 bits in every family. **The complete prime support is UNKNOWN.**
No large-factorization campaign or new parameter/point search is part of this audit.

The [exact audit](../../artifacts/generated-results/elliptic-curves/r17_scaling_prime_support_v1.json)
retains every integer resultant, small-prime valuation and residual cofactor.

## Necessary condition for a removable scale

For degree8 and12 binary forms A_h,B_h, multiplication by them gives the
integer20-by20 Sylvester matrix on forms of degree19. Its adjugate expresses
the resultant times any degree19 monomial as an integer polynomial combination
of A_h and B_h. If a primitive pair has p^4 dividing A_h and p^6 dividing B_h,
one coordinate is a p-adic unit, so **p^4 divides the resultant**.

For p at least5 the c4,c6 invariants make these coefficient conditions
necessary and sufficient for a removable short-model p scale. The resultant
condition is only necessary: it supplies candidate primes, not actual scales.
Primes2and3 need the integral-model translation analysis separately.

| Family | Known primes at least5 with resultant valuation at least4 | Unfactored cofactor bits |
| --- | --- | ---: |
| 103b2 | 5, 7, 11, 13, 17, 19, 47 | 1657 |
| 11952 | 5, 13 | 1779 |
| 074d9 | 5, 7, 11, 13, 17, 19, 23 | 1698 |
| 07ca9 | 5, 7, 13, 17, 19, 29 | 1717 |
| 08234 | 5, 7, 13, 19, 23, 37 | 1727 |
| 08f72 | 5, 13, 17, 29 | 1780 |

The [complete13-adic classification](R17_INTEGRAL_13_PARAMETER_CHARTS_2026-09-06.md)
already settles13. The [complete small-prime audit](R17_SMALL_PRIME_MINIMALITY_2026-09-06.md)
now excludes every other displayed candidate by independent exhaustive residue
checks. Every prime5through997 other than13 is universally minimal.
The unfactored cofactors remain unresolved; small-prime trial division does
not establish that all remaining primes occur with valuation less than4.

Replay `audit_r17_scaling_prime_support.py check` and
`verify_r17_scaling_prime_support.sage`. The first uses exact Bareiss elimination
on the homogeneous Sylvester matrix; the second uses Sage polynomial resultants.
Both pass. All22 isolated supplement stages now pass; see the replay completion below.

## Portable replay completion

All22 stages pass in the [standalone supplement replay](../../artifacts/generated-results/elliptic-curves/compact192_followup_supplement_portable_replay_v2.json), covering five exact histories, ten point-cloud proofs, rational geometry, the follow-up summary and five small-prime proof checks. The [7005716-byte archive](../../artifacts/generated-results/elliptic-curves/compact192_followup_supplement_evidence_v2.zip) needs no base archive. Its frozen notes precede this completion update. Version one omitted a command-invoked cloud checker; that failed package and replay are preserved. Version two includes the checker and all recorded cloud-proof sources, with unchanged arithmetic inputs.
