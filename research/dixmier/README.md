# Dixmier and Weyl algebras

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Formal Weyl-algebra lifts, finite scalar mates and centralizers.

**Authored assessment:** Formal lifting and restricted scalar classifications are established internally. Regular closed Newton support and algebraic order-(3,4) scalar pairs are excluded. Actual scalar centralizer fields are finite, and exact positional Long inputs have a global Ore presentation and conditional finite-completion criteria. Weighted Euler rank reduction is obstructed at every level by a normalized trace. Actual bounded rank-two Long realization and DC1/DC2 remain open; external review is separate.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| Formal quantum lifting and the finite-degree realization criterion | Internal result | proved | [Source](dixmier-formal-lifting/PROOF.md) |
| Scalar mate classifications and finite-order obstructions | Internal result | proved | [Source](dixmier-formal-lifting/ALGEBRAIC_INFINITY_LIMITS.md) |
| Spectral centralizer fields and finite quantum completion | Internal result | proved | [Source](dixmier-formal-lifting/FINITE_CENTRALIZERS_AND_COMPLETION.md) |
| Weighted Euler reduction and the rank-two descent obstruction | Internal result | proved | [Source](dixmier-formal-lifting/EULER_REDUCTION_TRACE_OBSTRUCTION.md) |

**Computed obligations:** 14 unresolved; 1 resolved within their exact scope.


Scope limits:

- All-order formal lifting has polynomial coefficients, but supplies no common finite spatial-degree bound. Evaluation of an unbounded h-adic series at one is unavailable.
- A proved bounded-realization criterion would produce a finite nonsurjective Weyl endomorphism from the specified noninvertible classical frame and common degree bound; that bound is not supplied.
- The pure-form scalar finite-mate classification and regular closed-support obstruction concern actual finite operators. Coefficients above the closed cone, negative-v terms, nonpositional first images and nonscalar carriers remain outside these results; DC1 and DC2 remain unresolved.
- The F1 cotangent transfer supplies finite nonsurjective Weyl endomorphisms in ranks at least three, with generator total degree at most twelve at rank three. This route does not provide a rank-two endomorphism or settle DC1/DC2.
- Actual scalar joint-centralizer fields are automatically finite, and their degrees divide every certified positive spectral rank; finite field degree still does not imply centralizer module finiteness.
- For exact positional Long R, one actual second-image mate forces polynomial joint centralizer; the two actual individual mates suffice for finite full completion. Their construction and a bounded normalization fixing R remain open.
- For regular exact-R Long deformations the closed Newton support cone is excluded, with no generic momentum-order cutoff. If generic S has a-order three, its prescribed T must have order at least five.
- Weighted Euler reduction of the rank-three cotangent example has a normalized trace at every scalar level, so neither the whole reduction, its finite matrices nor its full corners give a Weyl algebra.

Tracked obligations (computed states):

- **mathematical / open:** [A bounded full lift of the rank-two Long map, and DC1/DC2, remain open.](README.md) — GAP-BOUNDED-FULL-LIFT
- **mathematical / open:** [Negative-v coefficients, support above QPN2 closed cone, other infinity faces, deformed first images and nonscalar carriers remain open; order(3,4) is excluded but DC2 is not settled.](README.md) — GAP-DIXMIER-2
- **mathematical / open:** [No actual Long finite carrier or new Long radical certificate is supplied; bounded rank-two realization remains open.](README.md) — GAP-ACTUAL-LONG-CENTRALIZER
- **mathematical / open:** [Construct one actual finite commuting deformation of the Long positions together with actual generic individual mates and finitely many certified positive spectral ranks in exact finite Weyl polarizations whose gcd is one or two. All certificates must concern this same pair.](dixmier-formal-lifting/SPECTRAL_RANK_AND_FINITE_CENTRALIZER_FIELDS.md) — DIX-GAP-LONG-SPECTRAL-CERTIFICATE
- **mathematical / open:** [Construct actual finite S and actual generic finite D,T for the same exact R=X(2-3XQ), with Long classical S, [R,S]=0,[S,T]=h,[R,T]=0,[D,R]=h,[D,S]=0. The mate curvature may be nonzero; QEP1 removes it finitely.](dixmier-formal-lifting/EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md) — DIX-GAP-EXACT-R-INDIVIDUAL-MATES

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
