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
| [PAIR_72_108_REPRODUCTION.md](PAIR_72_108_REPRODUCTION.md) | Newton-to-coefficient reconstruction, exact ideals, certificates, and division audit |
| [WEIGHTED_WRONSKIAN_FIRST_BLOCK.md](WEIGHTED_WRONSKIAN_FIRST_BLOCK.md) | Hyperelliptic/de Rham interpretation and residual-scaling quotient of the audited first block |
| [SUPERELLIPTIC_DERHAM_ENGINE.md](SUPERELLIPTIC_DERHAM_ENGINE.md) | Reusable character-wise Hermite reduction, exact implementation, and frontier experiment design |
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

This does not prove JC(2), require both degrees to be at least 125, or show
that any candidate at or beyond 125 exists.  The finite-flatness theorem
removes the module-theoretic obstruction but leaves the surface boundary
classification open.
