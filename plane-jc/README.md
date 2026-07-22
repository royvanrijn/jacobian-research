# Plane Jacobian constraint program

This directory contains the audit and local reproduction of the current plane
degree frontier.  It is separate from the repository's three-dimensional
counterexample construction.

| Document | Purpose |
| --- | --- |
| [PROVENANCE.md](PROVENANCE.md) | Exact Zenodo/arXiv versions, files, licenses, and hashes |
| [DEGREE_FRONTIER_125.md](DEGREE_FRONTIER_125.md) | Theorem scope, reduction chain, and historical frontier |
| [PAIR_72_108_REPRODUCTION.md](PAIR_72_108_REPRODUCTION.md) | Newton-to-coefficient reconstruction, exact ideals, certificates, and division audit |
| [WEIGHTED_WRONSKIAN_FIRST_BLOCK.md](WEIGHTED_WRONSKIAN_FIRST_BLOCK.md) | Hyperelliptic/de Rham interpretation and residual-scaling quotient of the audited first block |
| [SUPERELLIPTIC_DERHAM_ENGINE.md](SUPERELLIPTIC_DERHAM_ENGINE.md) | Reusable character-wise Hermite reduction, exact implementation, and frontier experiment design |
| [NEWTON_BOUNDARY_DICTIONARY.md](NEWTON_BOUNDARY_DICTIONARY.md) | Qualified comparison with boundary/valuation language |
| [BOUNDARY_LATTICE_PREFILTER.md](BOUNDARY_LATTICE_PREFILTER.md) | Chart-aware localization/SNF gate and exact checker for complete proposed boundaries |
| [NEXT_DEGREE_FRONTIER.md](NEXT_DEGREE_FRONTIER.md) | Deterministic 125--150 candidate-table regression and ranked worklist |
| [SEARCH_POLICY.md](SEARCH_POLICY.md) | Consequences for future JC(2) computation |
| [cas/README.md](cas/README.md) | Replay commands, hashes, and independent checker |
| [AUDIT_COMPLETION.md](AUDIT_COMPLETION.md) | Objective-by-objective completion and residual limitations |
| [Weighted tangent suspension](../extended-geometry/WEIGHTED_TANGENT_SUSPENSION.md) | Exact Poisson-square and weighted-Wronskian bridge from the weighted JC(3) model |

Current scoped conclusion:

> Externally reduced and locally reproduced: subject to the exact published
> minimal/standard normal-form reduction, a hypothetical plane Keller
> counterexample has larger coordinate degree at least 125.

This does not prove JC(2), require both degrees to be at least 125, or show
that any candidate at or beyond 125 exists.
