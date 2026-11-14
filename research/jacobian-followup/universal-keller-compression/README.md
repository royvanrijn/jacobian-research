# Universal Keller maps in one fewer dimension

For every N>=5, an explicit determinant-one map on A^(N-1) realizes every rank-N finite etale algebra over every characteristic-zero field as a complete fibre. It retains degree N, monodromy S_N and absolute stable atomicity. The former construction used A^N.

The [proof](PROOF.md) removes one coefficient using the coprime weights five and six. It supplies explicit target formulas, a short Jacobian and full-fibre calculation, Zariski-local universality for families, and sufficient positive-characteristic and finite-field bounds. Its quotient section removes scaling choices; arbitrary primitive-presentation non-descent persists.

- UACOMP1: compressed universal maps and improved target-dimension upper bound.
- UAFAM1: complete-fibre universality locally on every Q-scheme.
- UAPOS1: characteristic p>N, with finite-field universality when q>4N-13.
- UAQUOT1: root-free scaling quotient and the surviving non-descent obstruction.

[Exact checker](verify.py) · [receipt](receipt.json) · [sources](sources.json) · [mathematical authority](../../MATH_STATUS.json).

From the repository root:

    .venv/bin/python research/jacobian-followup/universal-keller-compression/verify.py

The all-rank statements have written proofs; exact controls are not formalization or independent whole-proof review. Minimum dimension and a different canonical compiler remain open.

