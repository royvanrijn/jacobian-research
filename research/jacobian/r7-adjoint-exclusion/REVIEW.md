# Proof-obligation audit

14 November 2026. This is an author-side audit, not an independent external
review. The finite checker supports the named algebraic and numerical
steps only. Its positive two-contact-family control and tampered-certificate
controls guard against two specific false exclusions.

| Obligation | Evidence and boundary |
|---|---|
| Exact original scope | PROOF Section1 retains degree6, sole ramified affine prime e2,f1 and one of the two nonimmersed coefficient classes. It never replaces this by a complete one-prime omitted ledger. |
| Existing positivity/genus input | JC2POS1--3 are explicit dependencies. No positivity on an arbitrary further SNC blowup is asserted. |
| Adjoint pencil | PROOF Section2 uses the complete h0=2 system, integer Hodge bound M^2<1, connected Stein fibres and rational base to obtain M=F. H.F1 is excluded by degree-one image and Q-nefness, not by an unsupported rationality assumption. |
| Fixed components | Only after H.F2 is established is H.Z0 used. Minimality gives K.Z>=0, forcing rho8 and Z0. |
| Infinity residue degree | Borisov Theorem3.5(1) is newly imported at the target-P2 scope. One f2 prime from the older numerical list is removed; it is not miscounted as two primes. |
| Two-sheet cusp | The standard normal cusp double cover is A2, hence requires two crepant exceptional curves. W must also be singular over infinity. The complete rho8 boundary permits only two contracted curves in this alternative, a contradiction. |
| Three-sheet cusp | The standard S3 cubic cusp cover is smooth. Separated node/A5 local orbits and etale remaining factors make the affine normalization smooth. All three contracted curves are above the common infinity point. |
| Source finality | SOURCE_AUDIT identifies the exact assertion inside the primary proof of Borisov Theorem3.8. It is an external input, not proved by our checker or inferred from relative minimality. Separate source choices for individual curves suffice. |
| Model transport | Q.C starts at1 on a final type1 curve and drops by2v at each smooth vertical blowdown. Final nefness forces every v0, so the two type1 curves stay smooth -1 curves on Y. No SNC preservation for arbitrary curves is assumed. |
| Horizontal poles | H.F2 and verticality of both type1 curves give the complete positive-integer budget sum m_j*d_j2. Conic injectivity leaves one m2,d1 section, rather than a finite support search. |
| Double-line case | Field-degree conservation gives a cubic base map. Its poles are omitted by the whole A2 source, making the source pencil parameter polynomial after one Mobius change. Surjectivity and a critical value contradict the nonsingular target-linear component. |
| Generic parabola chart | The unique infinity section and degree-one x coordinate identify the whole generic finite-x fibre with A1. After shrinking the base, the chart is regular there, including unramified omitted points. No Jacobian from a birational parameterization is treated as a Keller determinant. |
| One/two contact coverage | The generic affine determinant is a quadratic in x. The sole e2 ramified prime has coefficient one and is horizontal (image degree10 exceeds a fibre degree2). Therefore its projection degree is exactly one or two. |
| One-contact family | Reduced slope denominator degree4 has a Bezout-coprime numerator. The generic-x resultant prevents cancellation for every constant A, including exceptional A. Function degrees14 and16 both exceed6. |
| All Mobius involutions | A degree-two rational map on P1 has a Mobius deck involution. The simple-pole bound rules out those moving infinity; all others are 2h-t. No arbitrary rational involution is silently replaced by an affine one. |
| Characteristic-zero exclusion | Exact coefficient identities and independent Fraction Bezout multiplication show no common zero of the two selected coefficient polynomials over C. No modular or floating-point promotion is used. |
| Both classes | Certificate inputs match both retained classification rows exactly, including their rational coefficients. |
| Remaining programme | Immersed r7, possible additional ramified rows, r5, nonnormal slices and JC2 are explicitly left unresolved. No new coordinate-degree bound follows. |

The written deduction supplies a contradiction for both coefficient classes
relative to the cited published source-finality input. Independent
end-to-end review should begin with that import and its model scope, then
the generic parabola chart. If either fails, the full exclusion must be
reopened; the rank-eight pencil result and algebraic architecture
exclusions are separately retained and do not depend on source finality.
