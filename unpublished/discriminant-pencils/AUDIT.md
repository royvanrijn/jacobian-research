# Proof audit for the polynomial tangent-pencil note

Edition: revised source and PDF dated 17 November 2026. This record describes the present standalone note; it does not claim that material removed from earlier versions remains in the paper.

## Current proof and contribution

The note gives a self-contained contact-stratum proof of a classical dual-curve specialization, with explicit thin-family genericity and an infinity calculation. Priority for the cusp and node counts is not established. Yuran, arXiv:2204.04626v1, Theorems 1.3 and 1.16, assumes a Newton polygon containing 5 times the standard triangle. The graph Y-H(W) has vertical width one, so that sufficient hypothesis does not cover this family. Yuran's conditional counts assume contact exclusions; the note proves those exclusions directly. This primary text was inspected during revision.

| Risk | Argument in the current paper | Scope |
|---|---|---|
| Infinity | Basepoint-free homogenization has one preimage over Z=0; S/T has nonzero linear term there. | Proves smoothness at the sole point at infinity for each exact-degree H. No relative no-escape theorem is asserted. |
| Higher one-branch singularity | The derivative is H second derivative times (1,r); contact three gives nonzero jet determinant. | Contact at least four lies in the forbidden (4) stratum. |
| Several branches | Distinct unramified parameters r,u have tangent determinant u-r. | Ramified-plus-other and triple fibers lie in (3,2) and (2,2,2). |
| Genericity | Factor H-ell=M_mu Q and take the closure of its image in coefficient space. | Each bad-image closure has dimension at most n-2, strictly below n-1. |
| Counting | Exclude all other singularities first; use squarefree H second derivative for n-2 cusps, then genus-delta. | Genus counts nodes; it does not exclude additional singularities. |
| Small example | H=W^3 gives (s,t)=(3r^2,2r^3), 27t^2=4s^3 and a smooth point at infinity. | An exact illustration of n=3, not proof for arbitrary n. |

The current paper contains no Jacobian-map motivation, admissible endpoint-slice argument, or relative equal-image escape argument. Earlier receipts describing those are superseded for this edition.

## Partial formalization and review limits

The pinned project and build command are linked in README.md. The written dimension-of-image argument, complete normalization/degree theorem, scheme-theoretic cusp/node criteria, and genus step remain outside the partial Lean development. No fresh Lean build or axiom audit was performed during this revision. A prior statement that the project uses no additional mathematical axioms is not a new transitive axiom check.

The revised PDF was compiled and visually inspected; details and source/PDF hashes are in verification.json. A PDF-only isolated internal AI reading of this revised edition is recorded with its exact hash and five reader answers in verification.json. The source-aware AI revision is internal work and is not external peer review. Novelty comparison beyond the explicitly inspected Yuran statements remains open, so this note stays parked as a possible companion or appendix.
