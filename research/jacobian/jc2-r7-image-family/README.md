# The complete gap-21 image family has two rational coefficient components

14 November 2026. PF2R7FAM1 classifies every normalized birational
degree-(6,10), gap-21 pair: after b=1, all coefficients lie on one
polynomial A1 family or one explicit Laurent Gm family. The proof covers
all exceptional branches and replays three literal graph memberships.

PF2R7FAM2 gives the exact finite nonimmersion loci in both charts.
Every immersed r7 image, including candidates with additional ramified
primes, now has one of these two explicit parameter descriptions.

The subsequent [collision classification](../jc2-r7-collisions/README.md)
identifies every affine singularity stratum of the polynomial family L,
including its34 exceptional parameters. The
[Laurent classification](../jc2-r7-laurent-collisions/README.md) completes
K as well:237 exceptional Laurent parameters and268 distinct nonnodal
coefficient points in the normalized union. Global matching remains open.

This is a set-theoretic classification, with an authored written proof
and exact algebra. No formalization or independent external review is
claimed. Global second-image equations, complement relations and Keller
realization remain open.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-image-family/verify.py
```

[Proof and formulas](PROOF.md) · [exact definitions](family.py) ·
[coefficient catalogue](families.json) · [certificates](certificates.json) · [replay](receipt.json) ·
[validation](validation.json) · [current omitted ledger](../jc2-omitted-divisors/README.md).
[Author-side audit](REVIEW.md) lists the completeness and assurance boundaries.
