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
| Generic coordinate order and exact closed-collision conductor | Internal result | proved | [Source](hc4-triangular-graphs/GLOBAL_SEXTIC_CONDUCTOR_EXCLUSION_AND_PERIOD_POLE_REPAIR.md) |
| Given-replica index from a descended differential | Internal result | proved | [Source](hc4-triangular-graphs/DIRECT_REPLICA_INDEX_FROM_A_DESCENDED_DIFFERENTIAL.md) |
| A regular source coordinate gives a global source translation | Internal result | proved | [Source](hc4-triangular-graphs/HAMILTONIAN_COMMON_COORDINATE_AND_ADDITIVE_COMPLETION.md) |
| Canonical full completion and whole unramified boundary coverage | Internal result | proved | [Source](hc4-triangular-graphs/CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md) |
| Valuation and exact-volume data do not force a source coordinate | Internal result | proved | [Source](hc4-triangular-graphs/VALUATION_MATRIX_AND_CONTRACTED_SOURCE_AXES_CONTROL.md) |
| A given regular source coordinate forces smooth boundary and whole root charts | Internal result | proved | [Source](hc4-triangular-graphs/COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md) |
| The actual intermediate source algebra is smooth affine etale | Internal result | proved | [Source](hc4-triangular-graphs/INTERMEDIATE_SOURCE_ALGEBRA_AND_REPLICA_FIBRE_COVERAGE.md) |
| Given-replica intermediate source images have no point gaps | Internal result | proved | [Source](hc4-triangular-graphs/INTERMEDIATE_SOURCE_ALGEBRA_AND_REPLICA_FIBRE_COVERAGE.md) |
| Individual source coordinates and complete algebraic boundary charts | Internal result | proved | [Source](hc4-triangular-graphs/INDIVIDUAL_SOURCE_COORDINATES_AND_UNCONDITIONAL_BOUNDARY_CHARTS.md) |
| Exact volume forces one common polynomial source pencil | Internal result | proved | [Source](hc4-triangular-graphs/EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md) |
| Unramified boundary residue degree and target normalization | Internal result | proved | [Source](hc4-triangular-graphs/UNRAMIFIED_BOUNDARY_RESIDUE_DEGREE_AND_TARGET_NORMALIZATION.md) |
| Codifferent preimages under target derivations | Internal result | proved | [Source](hc4-triangular-graphs/CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md) |

**Computed obligations:** 10 unresolved; 3 resolved within their exact scope.


Scope limits:

- The degree-five graph exclusion concerns a fixed v2 parent and graph degree, not all potentials or embeddings.
- Ricci-flat HC4 is a restricted sufficient theorem using the attributed Tang structural alternative. Scalar curvature does not substitute for Ricci flatness. Retained invertible cotangent examples have nonzero Ricci, so Ricci flatness is not necessary for invertibility.
- The exact moment lattice/trace criterion retains LW0, LY0 and KN0. Automatic T3 perfection remains the missing implication to JC2; HC4 and JC2 are linked by the stated Tang/cotangent equivalence.
- Two first source velocities already reconstruct the full function field in every source frame; their integral landing is still open. Reduced branch residue pairings and matrix factorizations are proved by finite-normal-surface algebra independently of KQ0. The analytic moment point hypothesis is now supplied by the audited Flenner-Zaidenberg quotient theorem; the separate differential source-layer input in KU2 remains unverified.
- In generic source directions two first velocities give the finite normalization of the source first-jet image and R=S[P_v,Q_v] exactly. The open target-integrality step remains. Smooth exact-volume completions with a regular actual source coordinate and positive canonical class exist in every n>=2,q>=2 of the explicit control family; their displayed finite plane maps have Jacobian zeros.
- Global directional coordinate-conductor degrees2,4 and6 are fully excluded, including all genus-three ambient classes and degree drops. Noninvertibility forces even global degree>=8; each individual source line has only the lower bound2.
- The full genus-three exclusion combines the hyperelliptic theorem with intrinsic quartic/cable action gates, filtered generic descent and period-pole repair. Neither a uniform family target automorphism nor a special-fibre lower bound8 is assumed. Higher conductors and automatic trace integrality remain open.
- Direct given-replica differential trace proves index divisibility without arbitrary-LND triangularization. Every actual intermediate source intersection algebra is smooth affine etale; replica stability fills all point gaps and proves disjoint smooth A1 reduced fibres. None makes velocities integral over B.
- Individual regular source coordinates now follow from the inspected Borisov minimality proof and an elementary positive-pullback sign argument. Actual unramified boundaries are smooth algebraic A1 with entire root charts. Exact volume gives a strict pole inequality; all-length multiplicity and cross-intersection bounds force ONE common source pencil, without adopting Xie or arbitrary-LND classification. This concerns the entire unramified completion, not ramified-prime integrality.
- Ordinary codifferent differentiation now characterizes the full finite order: S={f in L:D_Pf,D_Qf in Lambda}, and one generic target D suffices uniformly on source polynomials. In source-finite frames, one generic first velocity in Lambda forces inversion. Its automatic membership remains open; common-coordinate regularity on the unramified completion does not supply it.

Tracked obligations (computed states):

- **mathematical / open:** [Higher graph degrees, other embeddings, global coverage and unrestricted HC4 remain open.](README.md) — GAP-HESSIAN-TERMINATION
- **mathematical / open:** [Full HC4 and global coverage beyond the Ricci-flat subclass remain open; Ricci flatness is not necessary for invertibility, and the attributed Tang structural theorem is a separate input.](README.md) — GOAL-HC4
- **mathematical / open:** [Automatic perfection or full codifferent multiplication is not proved; HC4 and the plane Jacobian conjecture remain open.](README.md) — GAP-AUTOMATIC-T3-PERFECTION
- **mathematical / open:** [Global directional coordinate-conductor degrees2,4 and6 are fully excluded for actual noninvertible complex plane Keller maps in a monic source frame, including every genus-three ambient class and special degree drops. Noninvertibility therefore forces even global degree>=8. Every specialized source line still has only the lower bound2; special-line degrees2,4 or6 are not excluded. Higher conductor degrees, automatic first-velocity integrality over the original full target normalization, full codifferent multiplication, trace perfection, HC4 and JC2 remain unresolved.](hc4-triangular-graphs/GLOBAL_SEXTIC_CONDUCTOR_EXCLUSION_AND_PERIOD_POLE_REPAIR.md) — GAP-HESSIAN-HIGHER-COORDINATE-CONDUCTORS
- **mathematical / open:** [Can actual Keller/source constraints force one generic first velocity into the ordinary codifferent in a source-finite frame, eliminating the exact polar marked traces at original omitted ramified primes?](hc4-triangular-graphs/CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md) — GAP-HESSIAN-AUTOMATIC-ONE-VELOCITY-CODIFFERENT

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
