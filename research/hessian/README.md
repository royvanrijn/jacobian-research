# Hessian and Schur reductions

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Restricted Hessian classes, polynomial graphs and formal termination.

**Authored assessment:** The selected quintic, Ricci-flat and trace classes are closed; global coverage and polynomial termination remain open.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| Complete quintic graph exclusion with an explicit nilpotent algebra | Internal result | proved | [Source](../jacobian/affine-complements/FULL_COMPLEX_QUINTIC_GRAPH_EXCLUSION_AND_CUBIC_NILPOTENCE.md) |
| Ricci-flat HC4 inversion and the plane shear lemma | Internal result | proved | [Source](hc4-triangular-graphs/RICCI_FLAT_HC4_AND_ONE_DIMENSIONAL_JACOBIAN_IMAGES.md) |
| Even-moment lattices and the sixth-moment trace criterion | Internal result | open | [Source](hc4-triangular-graphs/EVEN_MOMENT_LATTICES_AND_SIXTH_MOMENT_TRACE_TEST.md) |

**Computed obligations:** 10 unresolved; 0 resolved within their exact scope.


Scope limits:

- The degree-five graph exclusion concerns a fixed v2 parent and graph degree, not all potentials or embeddings.
- Ricci-flat HC4 is a restricted sufficient theorem using the attributed Tang structural alternative. Scalar curvature does not substitute for Ricci flatness. Retained invertible cotangent examples have nonzero Ricci, so Ricci flatness is not necessary for invertibility.
- The exact moment lattice/trace criterion retains LW0, LY0 and KN0. Automatic T3 perfection remains the missing implication to JC2; HC4 and JC2 are linked by the stated Tang/cotangent equivalence.

Tracked obligations (computed states):

- **mathematical / open:** [Higher graph degrees, other embeddings, global coverage and unrestricted HC4 remain open.](README.md) — GAP-HESSIAN-TERMINATION
- **mathematical / open:** [Full HC4 and global coverage beyond the Ricci-flat subclass remain open; Ricci flatness is not necessary for invertibility, and the attributed Tang structural theorem is a separate input.](README.md) — GOAL-HC4
- **mathematical / open:** [Automatic perfection or full codifferent multiplication is not proved; HC4 and the plane Jacobian conjecture remain open.](README.md) — GAP-AUTOMATIC-T3-PERFECTION

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
