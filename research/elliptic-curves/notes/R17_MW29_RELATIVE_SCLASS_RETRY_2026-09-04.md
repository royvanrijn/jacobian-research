# MW29-relative S-class retry: historical stopping boundary

This completed bounded closure experiment on record fibres 356 and 385 has no
claim in [`MATH_STATUS.json`](../../MATH_STATUS.json).  It starts after the 29
known point half-ideals have been quotiented out, so it cannot calibrate an
MW17-only candidate gate or prove a global class-group, Selmer, or rank bound.
The retained [relative Selmer filter](ELKIES_R17_RELATIVE_2SELMER_PIPELINE.md)
states the current `UNKNOWN` boundary.

## Retained result

Three exact relation routes gave zero quotient-rank gain: dense reduced-ideal
collisions (14,341 and 13,213 attempts), bounded short-vector candidates
(88,648 and 94,134), and a 200-residual-ideal targeted tranche per curve.
The final relation graphs had no dependencies while new residual relations
introduced roughly three fresh outside prime-ideal vertices per edge.  The
particular family is therefore underdense; increasing its generic budget is
not supported.  All one-sided PARI class-quotient attempts timed out in
relation generation before yielding a class-group result.  These are bounded
observations, not a class-group upper bound or a proof of absence.

The durable method is to preserve incomplete factorization as exact residual
ideal-HNF vertices and add a relation only after its principal identity and
prime-ideal valuations verify.  A repeated reduced ideal is useful only with
the retained reduction multipliers; projective normalization prevents rational
multiples from becoming tautological rows.  The full derivation, hashes,
commands and local-artifact locations are preserved byte-for-byte in the
[archive](../../archive/elliptic-curves/notes/R17_MW29_RELATIVE_SCLASS_RETRY_2026-09-04.md.txt).

The input and output ledgers are ignored local artifacts and are absent from
this checkout, so no portable replay is available here.  Do not relaunch its
collectors or widen their budgets by default.  A future closure route needs a
certified F2 global quotient generator/upper bound, a ray-class or equivalent
complete provider, or an independently justified relation construction.
