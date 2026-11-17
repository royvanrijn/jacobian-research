# Exact prime-power degree spectrum

In characteristic zero, the nontrivial degrees of strongly absolutely stably atomic finite etale self-maps of smooth affine geometric torus quotients are exactly the prime powers. Every degree can be attained with a volume form and field versality for a specified finite etale group scheme. This closes the formerly remaining odd binary degree-existence gate.

- [Written proofs](PROOF.md): holonomy-twisted isogenies, the ideal-class equation, its prime-conductor solution, an explicit nonprincipal degree2048 control, and punctured canonical-line volume lifting.
- [Literature and limitations](LITERATURE_AND_GAPS.md).
- [Exact controls](verify.py), [receipt](receipt.json), [sources](sources.json), [integration](integration.json), [focused validation](validation.json).

Replay `python verify.py` with [requirements](requirements.txt). The arbitrary-degree proof uses classical arithmetic, not the bounded control table. Sources have nonconstant units; full symmetric universality, canonical affine-space landing and effective gluing remain open. Independent proof review and publication-priority review remain needed.

The [atomic-category extension](../universal-keller-atomic-cylinder/PROOF.md) shows that allowing nonfinite maps and affine identity cylinders does not enlarge this normal-affine-atomic degree spectrum. Self-map atomicity admits separate nonfinite dilation families and must not be conflated with that stronger category.
