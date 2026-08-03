# Plane Jacobian constraint program

This directory contains two complementary JC(2) programmes: the audit and
local reproduction of the current plane degree frontier, and a structural
programme based on the canonical finite normalization of an arbitrary
hypothetical counterexample.  It is separate from the repository's
three-dimensional counterexample construction.

| Document | Purpose |
| --- | --- |
| [PROVENANCE.md](PROVENANCE.md) | Exact Zenodo/arXiv versions, files, licenses, and hashes |
| [DEGREE_FRONTIER_125.md](DEGREE_FRONTIER_125.md) | Theorem scope, reduction chain, and historical frontier |
| [FINITE_NORMALIZATION_PROGRAM.md](FINITE_NORMALIZATION_PROGRAM.md) | Unconditional surface finite-flatness theorem; canonical branch/missing-boundary cover; arbitrary-puncture rigidity, bounded Pareto signature atlas, residual-different identity, and log-surface programme |
| [JC2_FINITE_NORMALIZATION_FRONTIER.md](JC2_FINITE_NORMALIZATION_FRONTIER.md) | Cubic cusp countermodel to automatic residue immersion; clean-packet classification; Orevkov Euler-budget closure of the cusp and all geometric degree three |
| [JC2_QUARTIC_PACKET_FRONTIER.md](JC2_QUARTIC_PACKET_FRONTIER.md) | Orevkov's exact quartic jump/two-boundary dichotomy; \(3+1\) cusp and \(2+2\) collision atlas; monodromy exclusion of the lone-cusp packet |
| [JC2_GLOBAL_COX_PACKET_ATTACK.md](JC2_GLOBAL_COX_PACKET_ATTACK.md) | Global quartic Cox and boundary-deletion reduction: raw pole pairs are unbounded and equal numerical matrices can have different residue cancellation; the exact finite Laurent compiler exhausts every lowering triangular polynomial, and the consecutive pole-change inequality proves marked multi-pole peak reduction for reduced alternating Jung words; the remaining gates are the invariant degree-four height and conductor-pairing tests |
| [QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md](QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md) | Conductor-decorated endpoint-semigroup no-finiteness theorem and bounded lattice compiler: the cusp semigroup is \(\langle2,3\rangle\), every connector has signed semigroup \(\{(u,v):u+v\le0\}\), and arbitrary connector count/contact prevents an unconditional finite quartic enumeration |
| [KELLER_PENCIL_AT_INFINITY_EXPERIMENT.md](KELLER_PENCIL_AT_INFINITY_EXPERIMENT.md) | Exact target-shear calibration of the full linear pencil: the same \(3+1\)/\(2+2\) finite packet and the same Jacobian determinant support generic zeta functions \(1\) and \((1+t^2)^{-1}\), proving that packet data do not determine pencil topology |
| [PAIR_72_108_REPRODUCTION.md](PAIR_72_108_REPRODUCTION.md) | Newton-to-coefficient reconstruction, exact ideals, certificates, and division audit |
| [WEIGHTED_WRONSKIAN_FIRST_BLOCK.md](WEIGHTED_WRONSKIAN_FIRST_BLOCK.md) | Hyperelliptic/de Rham interpretation and residual-scaling quotient of the audited first block |
| [SUPERELLIPTIC_DERHAM_ENGINE.md](SUPERELLIPTIC_DERHAM_ENGINE.md) | Reusable character-wise Hermite reduction, Gauss--Manin matrices, exact scalar Picard--Fuchs extraction, and frontier experiment design |
| [NEWTON_BOUNDARY_DICTIONARY.md](NEWTON_BOUNDARY_DICTIONARY.md) | Qualified comparison with boundary/valuation language |
| [BOUNDARY_LATTICE_PREFILTER.md](BOUNDARY_LATTICE_PREFILTER.md) | Chart-aware localization/SNF gate and exact checker for complete proposed boundaries |
| [INTRINSIC_A2_BOUNDARY_GATE.md](INTRINSIC_A2_BOUNDARY_GATE.md) | Adjunction/Noether reconstruction, pole-vector ramification gate, and intrinsic dicritical depth obstruction |
| [PLANE_BOUNDARY_EXCLUSION.md](PLANE_BOUNDARY_EXCLUSION.md) | Conditional smooth-target residue calculation and finite-flat conductor-packet inequality; the former residue-immersion claim for singular target curves is corrected by the cubic cusp audit |
| [LOG_BOUNDARY_COMPILER.md](LOG_BOUNDARY_COMPILER.md) | Certified branch scales to regular toroidal proximity graphs, complete boundary matrices, differents, and conductors |
| [FRONTIER_LOG_SCALE_AUDIT.md](FRONTIER_LOG_SCALE_AUDIT.md) | Fixed-completion replay of the `(72,108)` rays, the Wronskian-forced `E3∩E4` cluster, all five plane-return root-partition fans, the two 23-component terminal packages, their corrected `X^2` ramification, and the finite residue-cover split |
| [FRONTIER_CLOSING_ATTACKS.md](FRONTIER_CLOSING_ATTACKS.md) | Map-decorated boundary package, completed Case-2 composition sieve, and finite harmonic-cover, Pluecker, and log-Chern closure attacks |
| [POISSON_SQUARE_RIGIDITY.md](POISSON_SQUARE_RIGIDITY.md) | Reduced classification and exact eight-prime embedded filtration of the three-layer `[P,Q]=X^2` box |
| [NEXT_DEGREE_FRONTIER.md](NEXT_DEGREE_FRONTIER.md) | Deterministic 125--150 candidate-table regression and ranked worklist |
| [SEARCH_POLICY.md](SEARCH_POLICY.md) | Consequences for future JC(2) computation |
| [CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md](CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md) | Exact arbitrary-degree classification of every normalized support with at most six nonlinear monomial occurrences and the affine-normalized support lower bound seven |
| [AFFINE_SUPPORT_NEWTON_BRIDGE.md](AFFINE_SUPPORT_NEWTON_BRIDGE.md) | Exact obstruction to a coarse Newton/support bridge and the Kummer-character gate for the live `(75,125)` terminal block |
| [F2_75_125_DERIVATION.md](F2_75_125_DERIVATION.md) | Corrected \([t,z]=-z\) recurrence, complete `r=3` B0 tail, exact fixed-endpoint substitution, and tridiagonal unit elimination through layer `29`; the coupled block starts exactly at layer `28` |
<!-- status-consumer: PF2ER1 64378dad616fc3f2 -->
| [F2_MODIFIED_LAURENT_FAMILY.md](F2_MODIFIED_LAURENT_FAMILY.md) | Conditional all-`r` modified-series compiler, exact 14/22-function windows and Fitting residues, projected top-band unit obstruction, and the genuine terminal `A_(2r)` theorem |
| [F2_MODIFIED_CHART_BRIDGE.md](F2_MODIFIED_CHART_BRIDGE.md) | Corner-derived `gamma=2` chart, exact binomial-jet source ranks, polynomial lift of the negative tail, corrected tangent resonances, and length-27/resultant proof that every `r=3` literal polynomial-projection branch is empty |
<!-- status-consumer: PF2MCB1 6ff13314e0090f52 -->
<!-- status-consumer: PF2BH1 dcd3e54be59f32de -->
| [F2_BOUNDARY_HANDOFF.md](F2_BOUNDARY_HANDOFF.md) | Retained four-stratum contact census, failed contact-to-ramification surrogate, and updated handoff to the exact Kummer/target rows |
<!-- status-consumer: PF2KO1 c3a129906d2f75d2 -->
| [F2_KUMMER_ORBIT_TRANSFER.md](F2_KUMMER_ORBIT_TRANSFER.md) | Exact transfer around nonzero fifth-root orbits, exclusion of zero-root strata, and reduction to one principal chain or two copies |
<!-- status-consumer: PF2TR1 bb41ccb3d135dbf2 -->
<!-- status-consumer: PF2GC1 33dbc5ff48b5d064 -->
| [F2_TERMINAL_RESIDUE_COVER.md](F2_TERMINAL_RESIDUE_COVER.md) | Target ray `(5,2)`, degree floors `6`/same-target `12`, three forced node attachments, genus-25 `A_6` Galois closure, purity/different ledger, geometric `A_6` versus arithmetic `S_6`, and trivial target-fixed deck group |
| [F2_A6_SIMPLE_SPECTATOR_GLUING.md](F2_A6_SIMPLE_SPECTATOR_GLUING.md) | Conditional spectator gluing: six genus-zero `S_7` one-cycle classes, exact Kummer/terminal order-five comparison, and one inertia-supported degree-11 `S_11` class in the rational-source fivefold model |
| [UNIBRANCH_SPECTATOR_COUNTERMODELS.md](UNIBRANCH_SPECTATOR_COUNTERMODELS.md) | Universal finite-free unibranch packets with an étale spectator; exact refutation of a purely local exclusion and isolation of the global \(\mathbb A^2\)-open obstruction |
| [cas/README.md](cas/README.md) | Replay commands, hashes, and independent checker |
| [AUDIT_COMPLETION.md](AUDIT_COMPLETION.md) | Objective-by-objective completion and residual limitations |
| [Weighted tangent suspension](../extended-geometry/WEIGHTED_TANGENT_SUSPENSION.md) | Exact Poisson-square and weighted-Wronskian bridge from the weighted JC(3) model |

Current scoped conclusion:

> Externally reduced and locally reproduced: subject to the exact published
> minimal/standard normal-form reduction, a hypothetical plane Keller
> counterexample has larger coordinate degree at least 125.
>
> Independently and unconditionally: the canonical finite normalization of
> every plane Keller map is a finite free cover of \(\mathbb A^2\).  Its
> missing-boundary primes freely generate the normalization's class group.
> Orevkov's three-sheeted theorem excludes geometric degrees two and three,
> so the global finite-normalization classification starts at degree four.
> In degree four, Orevkov's budget leaves exactly a one-boundary \(3+1\)
> jump packet or a ramified-plus-unramified two-boundary packet; the lone
> cusp without a \(2+2\) self-collision is excluded by monodromy.
>
> In fixed tangent-to-identity coordinates, every exact support with at most
> six nonlinear monomial occurrences is certified invertible in arbitrary
> degree.  Consequently any noninvertible plane Keller map has minimum
> support at least seven across all tangent-to-identity affine normalizations.
>
> For the first conditional coordinate-degree row `(75,125)`, the selected F2
> terminal valuation has `(e,f)=(1,6)`.  Any global realization therefore has
> geometric degree at least six; two distinct double-root packets over the
> same target divisor force at least twelve.  The valuation is centered at
> target infinity, so no affine-sheet increment applies, and purity requires
> a separate affine ramification row.  Three interior preimages of the target
> toric nodes already force three boundary-attachment points.  Its geometric
> residue monodromy is `A_6`, so global monodromy has `A_6` as a nonabelian
> simple composition factor.

This does not prove JC(2), require both degrees to be at least 125, or show
that any candidate at or beyond 125 exists.  The finite-flatness theorem
removes the module-theoretic obstruction but leaves the surface boundary
classification open.
