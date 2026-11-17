# Dixmier and Weyl algebras

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Formal Weyl lifts, scalar mate obstructions, finite centralizers and exact quantum division spectra.

**Authored assessment:** Formal lifting and scoped scalar results are established internally. Every regular polynomial scalar Long deformation is excluded; remaining exact-R candidates need a negative Laurent tail escaping the quadratic critical-face test. Algebraic orders (3,4) and (3,5) are excluded, and a generic cubic Long repair needs partner order at least seven. The fixed homogeneous Euler/cotangent cubic family has exact scalar spectrum 5chZ, excluding its proposed compatibility throughout the actual division ring. Finite bounded Long lifts and DC1/DC2 remain open; external review is separate.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| Formal quantum lifting and the finite-degree realization criterion | Internal result | proved | [Source](dixmier-formal-lifting/PROOF.md) |
| Scalar mate classifications and finite-order obstructions | Internal result | proved | [Source](dixmier-formal-lifting/ALGEBRAIC_INFINITY_LIMITS.md) |
| Spectral centralizer fields and finite quantum completion | Internal result | proved | [Source](dixmier-formal-lifting/FINITE_CENTRALIZERS_AND_COMPLETION.md) |
| Weighted Euler reduction and the rank-two descent obstruction | Internal result | proved | [Source](dixmier-formal-lifting/EULER_REDUCTION_TRACE_OBSTRUCTION.md) |
| The homogeneous Euler spectrum and cubic cover obstruction | Internal result | proved | [Source](dixmier-formal-lifting/QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md) |

**Computed obligations:** 16 unresolved; 1 resolved within their exact scope.


Scope limits:

- All-order formal lifting has polynomial coefficients, but supplies no common finite spatial-degree bound. Evaluation of an unbounded h-adic series at one is unavailable.
- A proved bounded-realization criterion would produce a finite nonsurjective Weyl endomorphism from the specified noninvertible classical frame and common degree bound; that bound is not supplied.
- The F1 cotangent transfer supplies finite nonsurjective Weyl endomorphisms in ranks at least three, with generator total degree at most twelve at rank three. This route does not provide a rank-two endomorphism or settle DC1/DC2.
- Actual scalar joint-centralizer fields are automatically finite, and their degrees divide every certified positive spectral rank; finite field degree still does not imply centralizer module finiteness.
- For exact positional Long R, one actual second-image mate forces polynomial joint centralizer; the two actual individual mates suffice for finite full completion. Their construction and a bounded normalization fixing R remain open.
- Weighted Euler reduction of the rank-three cotangent example has a normalized trace at every scalar level, so neither the whole reduction, its finite matrices nor its full corners give a Weyl algebra.
- QCL1 excludes every regular finite polynomial v,a deformation of Long scalar S, including arbitrary mixed positive support and generic momentum order. An exact-R finite scalar candidate needs an actual negative-v tail and, after critical normalization, a negative weighted coefficient or a critical face of order>=3.
- Algebraic scalar exact order pairs(3,4) and(3,5) are excluded. A regular generic cubic Long first image requires partner order>=7; no all-cubic-mate theorem is proved.
- Actual internal radical-cover towers and their quantum division automorphism orbit cannot be bounded full Long lifts, even with parameter poles. This uses only commutative-residue Hensel lifting and Long geometric cubic.
- The actual cubic Weyl model, binary-dihedral realization, classical Long quintic square and individual quantum covers are valid. Their displayed compatible fifth root is excluded, so the conditional fourth-image identities do not construct an actual Long tuple.
- For every rational homogeneous source Euler function in the fixed scalar cotangent cubic transport, the actual scalar eigenvalue set is exactly 5chZ. This excludes any Euler pair with canonical step ch and every rational target fifth-root equation in the entire Weyl division ring, without regularity or residue assumptions.
- The next bounded construction must change this quantum source or cubic transport, or use a different route. Formal completions may have additional eigenoperators; no theorem places every bounded Long lift in the excluded family.

Tracked obligations (computed states):

- **mathematical / open:** [A bounded full lift of the rank-two Long map, and DC1/DC2, remain open.](README.md) — GAP-BOUNDED-FULL-LIFT
- **mathematical / open:** [Negative-weight or higher-order Laurent critical faces, deformed first images and nonscalar carriers remain open. All regular polynomial scalar Long deformations and algebraic orders(3,4),(3,5) are excluded; DC1/DC2 are unresolved.](README.md) — GAP-DIXMIER-2
- **mathematical / open:** [No actual Long finite carrier or new Long radical certificate is supplied; bounded rank-two realization remains open.](README.md) — GAP-ACTUAL-LONG-CENTRALIZER
- **mathematical / open:** [Construct one actual finite commuting deformation of the Long positions together with actual generic individual mates and finitely many certified positive spectral ranks in exact finite Weyl polarizations whose gcd is one or two. All certificates must concern this same pair.](dixmier-formal-lifting/SPECTRAL_RANK_AND_FINITE_CENTRALIZER_FIELDS.md) — DIX-GAP-LONG-SPECTRAL-CERTIFICATE
- **mathematical / open:** [Construct actual finite S and actual generic finite D,T for the same exact R=X(2-3XQ), with Long classical S, [R,S]=0,[S,T]=h,[R,T]=0,[D,R]=h,[D,S]=0. The mate curvature may be nonzero; QEP1 removes it finitely.](dixmier-formal-lifting/EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md) — DIX-GAP-EXACT-R-INDIVIDUAL-MATES
- **mathematical / open:** [For every partner order coprime to three and all lower flows, justify admission of the normalized rank-three spectral connection and single-valued meromorphic scalar U,V reconstruction over every finite time, including resonant apparent poles, degenerate monodromy and gauge normalization. Miwa nonresonance and Malgrange meromorphic connection coefficients alone do not discharge this application. No claim that a modern suitable theorem is absent is made.](dixmier-formal-lifting/CUBIC_QUINTIC_STRING_OBSTRUCTION.md) — DIX-GAP-CUBIC-ISOMONODROMY-RECONSTRUCTION
- **mathematical / open:** [Construct actual quantum Long images outside the excluded homogeneous Euler source and scalar cotangent cubic transport, and prove finite global admission. Rational source Euler functions, target reorderings, mixed division denominators and parameter poles cannot repair that entire fixed family. No theorem normalizes every bounded Long lift into it.](dixmier-formal-lifting/QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md) — DIX-GAP-NONHOMOGENEOUS-COVER-ADMISSION

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
