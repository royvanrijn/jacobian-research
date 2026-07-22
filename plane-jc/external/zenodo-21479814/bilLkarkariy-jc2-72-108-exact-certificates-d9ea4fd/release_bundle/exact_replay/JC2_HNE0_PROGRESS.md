# JC2 (72,108) — hard branch progress report

Date: 2026-07-21

## Scope

This report concerns the remaining hard branch obtained from Proposition 4.3, Case 1, after the exact Laurent reduction and the exhaustive branch/cascade reductions contained in the earlier certificate package.

## Exact regeneration completed

The branch was rebuilt independently from `case1_checkpoint.pkl`.

It reduces exactly to six equations in three variables `(h,u1,u2)` over a degree-35 number field. Their term counts are

    30, 30, 40, 40, 40, 30

and their total degrees are

    7, 7, 8, 8, 8, 7.

A Z/7 grading was verified term by term. After scaling the variables with weights `(5,6,5)`, the coefficient field descends exactly to

    L = Q[w]/(w^5 - w^4 + 3*w^3 + 3*w^2 + 26).

The transformed six-equation system is stored in `hne0_polred.pkl`.

## Finite-field certificate and fixed-support audit

At p=71, Singular constructs and directly checks

    h = T1*E1 + T2*E2 + T3*E3 + T4*E4

inside the degree-5 finite-field extension.

After reduction modulo the syzygy module, the four multipliers contain 385 monomials, split as

    100, 99, 93, 93.

For this fixed multiplier support, the coefficient map has

    402 rows x 385 columns

and rank

    385

over F_(71^5). Hence the p=71 multiplier vector is unique within this fixed support.

This is still not a characteristic-zero certificate.

## CRT experiment: important negative result

Certificates were computed at 800 good primes near 5*10^8. Their product has about 6960 decimal digits.

Naive rational reconstruction returned 1165 apparent rational coordinates and left 760 unresolved. A fresh holdout prime, not used in the CRT, rejected every one of the 1165 apparent reconstructions.

Therefore those values were accidental rational reconstructions of modular residues, not exact coefficients. The CRT approach has not yet produced any verified characteristic-zero multiplier coefficient.

Conditionally, if an exact solution exists with this fixed support, its coefficient heights exceed the current ordinary rational-reconstruction range, and the explicit identity is much larger than expected.

## Exact CAS attempt

An exact `liftstd` run over the degree-5 number field was allowed 10 minutes under a 1.8 GB virtual-memory cap. It did not finish and produced no certificate.

## Honest status

The exact target remains

    h = T1*E1 + T2*E2 + T3*E3 + T4*E4

in `L[h,u1,u2]`.

No characteristic-zero identity has yet been obtained, so the hard branch and the (72,108) case are not yet excluded.

## Recommended next move

Blind CRT on Singular's prime-dependent syzygy normal forms should not be continued as the primary strategy. The next rigorous route is to build a fixed square linear subsystem from the 402x385 coefficient matrix, then use a deterministic p-adic/Dixon solve or exact sparse number-field linear algebra. Every lifted solution must be verified against all 402 rows exactly.
