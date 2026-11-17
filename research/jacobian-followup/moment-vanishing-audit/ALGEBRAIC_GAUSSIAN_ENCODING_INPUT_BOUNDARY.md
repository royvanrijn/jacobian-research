# Rank labels require a minimal branch polynomial

16 November2026. Validation addendum. The universal
[Gaussian realization proof](ALGEBRAIC_FIXED_POINT_LIFTS_AND_UNIVERSAL_GAUSSIAN_IMAGES.md)
states the correct rank [K(t,f):K(t)]. Its branch-normalization identities
work with any polynomial F having nonzero derivative on the selected branch.
The finite-transcription helper `encode` also records the degree of F as
`actual_generic_rank`; that label presupposes that F is minimal. All five
retained models meet this condition, but the helper did not enforce it.
Its sealed evidence is preserved.

The current rational entry point is
[encode_checked](algebraic_gaussian_lift_encoding_guard.py). It verifies
irreducibility over Q(t) before calling the retained constructor and rejects
nonminimal input. This is a scope guard, not a new minimization algorithm.
Over other effective fields use their exact minimal branch polynomial.

The exact rejected control is F=(x-t)(x-t-1), branch f=t. The normalization
identities remain valid and yield g=0, H=t*z^2, so the actual lift is safe
with rank1. The degree2 of the nonminimal F cannot be used as its rank.
This correction does not alter any retained proof, model or assurance flag.

[Guard verification](algebraic_gaussian_lift_encoding_guard_verification.json)
checks all five retained inputs and rejects this nonminimal polynomial,
a repeated polynomial, a constant and an unsupported coefficient parameter.
