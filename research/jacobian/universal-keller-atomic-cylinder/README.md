# Atomic cylinders and the affine Keller boundary

Atomic nonfinite maps on a geometric torus-quotient cylinder have base degree one. Their canonical marked data therefore restrict to an actual affine fibre, without the earlier metabelian hypothesis. Primitive self-map atomicity also needs no acyclicity when units are constant and the canonical bundle is torsion. On the one fixed fourfold Gm times A3, the exact nonfinite degree spectrum is every integer at least3, with infinitely many stably inequivalent atomic volume-preserving representatives in each degree.

- [Written general proofs](PROOF.md), including the exact normal-affine category boundary and infinitely many stably inequivalent nonconstant dilation families.
- [Literature and limitations](LITERATURE_AND_GAPS.md).
- [Exact controls](verify.py), [receipt](receipt.json), [source pins](sources.json), [integration](integration.json), [focused validation](validation.json).

Replay `python verify.py` with [requirements](requirements.txt). The controls do not formalize Ax's theorem, Zariski Main, geometric monodromy or relative descent. The [trace filter](../universal-keller-canonical-filter/PROOF.md) now supplies canonical non-atomic whole-affine-fibre bundle descent in every rank. Its atomic version, stronger trivial-bundle requirements and publication priority remain unresolved.

The [arithmetic Weil extension](../universal-keller-feedback-rigidity/PROOF.md) gives a separate base-field result: evaluating a primitive seed in A is k-stably self-atomic exactly when the finite etale algebra A is a field. Stable factor length after any constant extension is its field-component count; absolute atomicity fails when the rank exceeds one.
