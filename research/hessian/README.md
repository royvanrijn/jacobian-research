# Hessian and Schur reductions

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Restricted Hessian classes, polynomial graphs and formal termination.

**Authored assessment:** The selected quintic, Ricci-flat and trace classes are closed; global coverage and polynomial termination remain open.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| Complete quintic graph exclusion with an explicit nilpotent algebra | Internal result | proved | [Source](../jacobian/affine-complements/FULL_COMPLEX_QUINTIC_GRAPH_EXCLUSION_AND_CUBIC_NILPOTENCE.md) |
| Ricci-flat HC4 inversion and the plane shear lemma | Internal result | proved | [Source](hc4-triangular-graphs/RICCI_FLAT_HC4_AND_ONE_DIMENSIONAL_JACOBIAN_IMAGES.md) |
| Even-moment lattices and the sixth-moment trace criterion | Internal result | proved | [Source](hc4-triangular-graphs/EVEN_MOMENT_LATTICES_AND_SIXTH_MOMENT_TRACE_TEST.md) |
| First source derivatives and finite source reconstruction | Internal result | proved | [Source](hc4-triangular-graphs/REPLICA_DEGREE_AND_SOURCE_DERIVATIVE_FIELD_RECONSTRUCTION.md) |
| Reduced branch trace defects and symmetric matrix factorizations | Internal result | proved | [Source](hc4-triangular-graphs/REDUCED_BRANCH_PAIRING_AND_SIXTH_MOMENT_DEFECT.md) |
| Regular source-coordinate completions and the primitive obstruction | Internal result | proved | [Source](hc4-triangular-graphs/REGULAR_SOURCE_COORDINATE_COMPLETIONS_AND_PRIMITIVE_RESIDUES.md) |
| Generic coordinate order and exact closed-collision conductor | Internal result | proved | [Source](hc4-triangular-graphs/GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md) |

**Computed obligations:** 11 unresolved; 0 resolved within their exact scope.


Scope limits:

- The degree-five graph exclusion concerns a fixed v2 parent and graph degree, not all potentials or embeddings.
- Ricci-flat HC4 is a restricted sufficient theorem using the attributed Tang structural alternative. Scalar curvature does not substitute for Ricci flatness. Retained invertible cotangent examples have nonzero Ricci, so Ricci flatness is not necessary for invertibility.
- The exact moment lattice/trace criterion retains LW0, LY0 and KN0. Automatic T3 perfection remains the missing implication to JC2; HC4 and JC2 are linked by the stated Tang/cotangent equivalence.
- Two first source velocities already reconstruct the full function field in every source frame; their integral landing is still open. Reduced branch residue pairings and matrix factorizations are proved by finite-normal-surface algebra independently of KQ0. The analytic moment point hypothesis is now supplied by the audited Flenner-Zaidenberg quotient theorem; the separate differential source-layer input in KU2 remains unverified.
- In generic source directions two first velocities give the finite normalization of the source first-jet image and R=S[P_v,Q_v] exactly. The open target-integrality step remains. Smooth exact-volume completions with a regular actual source coordinate and positive canonical class exist in every n>=2,q>=2 of the explicit control family; their displayed finite plane maps have Jacobian zeros.
- Global directional coordinate-conductor degree2 and4 are fully excluded, including degree drops and tacnodes. Noninvertibility forces even global degree>=6; special fibre degrees2 or4 and all automatic integrality/termination obligations remain open.
- The generic hyperelliptic genus-three ambient branch at global conductor degree6 is excluded by three-period rigidity. This does not exclude intrinsic quartic/cable embeddings or raise the unconditional global lower bound above6.

Tracked obligations (computed states):

- **mathematical / open:** [Higher graph degrees, other embeddings, global coverage and unrestricted HC4 remain open.](README.md) — GAP-HESSIAN-TERMINATION
- **mathematical / open:** [Full HC4 and global coverage beyond the Ricci-flat subclass remain open; Ricci flatness is not necessary for invertibility, and the attributed Tang structural theorem is a separate input.](README.md) — GOAL-HC4
- **mathematical / open:** [Automatic perfection or full codifferent multiplication is not proved; HC4 and the plane Jacobian conjecture remain open.](README.md) — GAP-AUTOMATIC-T3-PERFECTION
- **mathematical / open:** [Global directional conductor degrees2 and4 are fully excluded, so noninvertibility forces even global degree>=6. At degree6 the generic hyperelliptic genus-three ambient branch is excluded, including persistent contacts and special degree drops. The other two source-classification branches, intrinsic quartic and degree-six cable, remain open, as do higher degrees and any unexpanded family-coverage obligations. Individual special fibres may have degree2 or4. Automatic first-velocity integrality, full codifferent multiplication, trace perfection, HC4 and JC2 remain unresolved.](hc4-triangular-graphs/HYPERELLIPTIC_SEXTIC_CONDUCTORS_AND_THREE_PERIOD_RIGIDITY.md) — GAP-HESSIAN-HIGHER-COORDINATE-CONDUCTORS

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
