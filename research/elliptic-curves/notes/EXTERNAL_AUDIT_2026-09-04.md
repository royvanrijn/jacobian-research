# Historical external mathematical and machinery audit — 4 September 2026

This audit examined the September 4 working tree under bounded replay limits.
Its inventory, command outcomes, environment failures, and proposed next steps
are historical. The exact report is retained in
[`archive/elliptic-curves/notes/EXTERNAL_AUDIT_2026-09-04.md.txt`](../../archive/elliptic-curves/notes/EXTERNAL_AUDIT_2026-09-04.md.txt)
(SHA-256 `a66605b10d3c66d375c2051cf1b5aa9d12856690bcad19b9a41d54d179ec50a7`),
with its [execution record](../../artifacts/generated-results/elliptic-curves/external_audit_20260904.json)
and [historical replay inputs](../../archive/elliptic-curves/external-audit-2026-09-04/README.md).

[MATH_STATUS.json](../../MATH_STATUS.json) and its canonical sources are the
current mathematical authority. Use [REPRODUCE.md](../REPRODUCE.md) for active
verification guidance; the archived commands do not authorize a replay,
descent, or search campaign.

## Preserved implementation rule

When ranking finite local or fingerprint signatures modulo known
Mordell--Weil images, compute a quotient normal form by clearing every known
echelon pivot. A membership test may stop at a first free coordinate, but that
partial remainder is not a canonical quotient representative: equivalent
cosets can otherwise appear independent. Check the residual rank against
`rank(known + candidates) - rank(known)` on small exhaustive fixtures.

This finite signature is only a selected target. A zero residual says that
target does not distinguish the candidate from the known span; it does not
establish a global square, a Selmer conclusion, or a rank bound. The current
implementation and its regression are indexed in [ALGORITHMS.md](../../knowledge/ALGORITHMS.md).

The audit's pairing-attestation, primality, coordinate, and protocol repairs
remain represented by their current canonical sources and tests. The archive
retains their historical context without keeping a competing status narrative.
