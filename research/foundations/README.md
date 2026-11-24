# Formal verification

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Formal developments and characteristic-zero Hopf/Image classifications: complete even-start purity via marked Davenport-Zannier pairs, Fuss-Catalan geometric counts and exact field descent; all-height safety and sharp degree/prime thresholds.

**Authored assessment:** Formal scope varies by project. GMC2 is conditional on two named axiom interfaces; binary GVC and full generic discriminant geometry remain partial.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| GVC counterexamples and binary supporting algebra | Partial result | proved | [Source](formal/gvc/README.md) |
| Binary Gaussian moments in Lean | Internal result | conditional | [Source](formal/gmc2/README.md) |
| Characteristic-zero finite etale Keller construction | Internal result | proved | [Source](formal/finite-etale-keller/README.md) |
| Tangent-pencil algebra in Lean | Partial result | proved | [Source](formal/discriminant-pencils/README.md) |
| Primitive quadratic norms and binary correlations in Lean | Partial result | proved | [Source](formal/quadratic-norms/README.md) |
| Complete complex quadratic-height classification and exact finite moment ideals | Internal result | proved | [Source](extended-geometry/HOPF_LIFT_CLASSIFICATION.md) |
| Complete quadratic-height Hopf purity and safety classification | Internal result | proved | [Source](extended-geometry/QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md) |
| Odd cubic-height rational reconstruction and sharp low-degree safety | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md) |
| Complete safety classification in the three-height quintic Image box | Internal result | proved | [Source](extended-geometry/THREE_HEIGHT_QUINTIC_IMAGE_SAFETY.md) |
| Complete one-pole odd-cubic purity and safety classification | Internal result | proved | [Source](extended-geometry/ONE_POLE_ODD_CUBIC_CLASSIFICATION.md) |
| Exact/logarithmic purity reduction and one-double-root prime obstruction | Internal result | proved | [Source](extended-geometry/MIXED_CENTER_LOGARITHMIC_SEPARATION.md) |
| Odd-start rationality, full zero-center safety and a quintic exclusion | Internal result | proved | [Source](extended-geometry/ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md) |
| Odd-cubic septic safety and the sharp odd unsafe threshold nine | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_DEGREE_SEVEN_SAFETY.md) |
| Odd-cubic degree11 safety and the sharp prime threshold thirteen | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_DEGREE_ELEVEN_SAFETY.md) |
| Canonical polynomial descent and reversible purity transport for odd-cubic phases | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_CANONICAL_POLYNOMIAL_DESCENT.md) |
| All-height safety classification by the two fixed original markers Z and Z squared | Internal result | proved | [Source](extended-geometry/ALL_START_HEIGHT_MARKER_SAFETY.md) |
| Full odd-cubic degree17 safety by canonical quotient and literal pole obstruction | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_DEGREE_SEVENTEEN_SAFETY.md) |
| Sharp unsafe degree 2a plus 4 and unique equality family at every starting height | Internal result | proved | [Source](extended-geometry/ALL_START_HEIGHT_SHARP_UNSAFE_DEGREE.md) |
| Full odd-cubic degree19 safety with exact root-orbit and integer certificates | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_DEGREE_NINETEEN_SAFETY.md) |
| Sharp prime threshold thirteen across every two-height starting profile | Internal result | proved | [Source](extended-geometry/TWO_HEIGHT_MINIMUM_UNSAFE_PRIME_THIRTEEN.md) |
| Admissible all-height pullbacks and an explicit obstruction to cubic polynomial reduction | Internal result | proved | [Source](extended-geometry/HIGHER_START_CUBIC_REDUCTION_COUNTEREXAMPLES.md) |
| Complete even-start purity and Fuss-Catalan enumeration | Internal result | proved | [Source](extended-geometry/CONSTANT_ENDPOINT_DAVENPORT_ZANNIER_CLASSIFICATION.md) |
| Canonical rational composition normal form for every pure odd-start unit phase | Internal result | proved | [Source](extended-geometry/ODD_START_RATIONAL_COMPOSITION_NORMAL_FORM.md) |
| Marked gcd recognition and finite purity at every odd starting height | Internal result | proved | [Source](extended-geometry/ODD_START_MARKED_GCD_AND_FINITE_PURITY.md) |
| Uniform finite purity cutoff and reduced moment locus for every starting height | Internal result | proved | [Source](extended-geometry/ALL_START_TWO_HEIGHT_EFFECTIVE_PURITY_PREFIX.md) |
| Lean verification of Mathieu kernels on quadratic complements | Partial result | stale | [Source](formal/mathieu-quadratic/README.md) |
| All-field fixed-unit-marker counterexamples on Cayley-plane projector products | Internal result | proved | [Source](extended-geometry/CAYLEY_PLANE_PROJECTOR_PRODUCT_COUNTEREXAMPLE.md) |
| Compact unitary EIV products over all fields | Internal result | proved | [Source](extended-geometry/ALBERT_EIV_AND_REDUCED_NORM_ONE_PRODUCT_COUNTEREXAMPLES.md) |
| Reduced Albert norm-one product counterexamples over all fields | Internal result | proved | [Source](extended-geometry/ALBERT_EIV_AND_REDUCED_NORM_ONE_PRODUCT_COUNTEREXAMPLES.md) |

**Computed obligations:** 11 unresolved; 1 resolved within their exact scope.


Scope limits:

- Formal verification is scoped by exact declarations, transitive assumptions and the paper bridge. A paper theorem does not imply that its full formalization is complete.
- GMC2 retains two explicit mathematical axiom interfaces, shown as open formalization obligations. Its propagated formal claim remains conditional until discharged; its written mathematical theorem has separate proof support.
- Binary GVC and generic tangent-pencil geometry are partial formalizations. The focused finite-etale certificate excludes the broader supplied-presentation and analytic counting scopes. No Lean build was rerun for this registry migration.
- The retained two-pair Image cubic proof covers exact bihomogeneous bidegree (3,3), with sharp mixed cutoff and minimum failing balanced degree four. Independent review and formal coverage are separate obligations.
- The entire complex polynomial winding-one quadratic-height numerator class is classified. On the nonzero constant-coefficient branch, its ordinary full moment ideal is exactly the finite coefficient ideal of (C-xCprime)^2D+C^3. The higher Hopf residual gap is discharged at this scope; general profiles and nilpotence indices are separate.
- The full polynomial winding-one quadratic-height Hopf family, including arbitrary linear height, is now classified over every characteristic-zero field. Its only unsafe branch has even least balanced degree and a fixed marker surviving at every power. This strengthens the earlier B=0 field classification while preserving its separate localized ordinary-ideal statement; arbitrary quintic and septic phases remain open.
- In the entire winding-one odd cubic-height class t*x^-1(C+t^2D), pure phases have an automatic rational center on the nonzero-constant branch. Every pure phase of least degree<=5 has an explicit one-sided safe normal form; the retained degree6 example makes the class threshold sharp. Finite Bezout identities follow the all-order rationality reduction and do not provide a universal moment-prefix cutoff.
- The full unit-center one-pole odd-cubic class is now classified, including constant and linear boundaries, two explicit nonconstant unsafe branches and fixed original markers. It is exactly the rational-primitive dx/C pure class. All nonconstant squarefree centers are excluded in every degree; the remaining multi-root repeated-center class has an explicit open obligation. Restricted odd prime degrees are exactly3mod10 greater than3, without a global prime-degree conclusion.
- All unit-center odd-cubic purity now separates into rational and logarithmic inverse-branch conditions; the logarithmic component cannot be silently discarded. A center with exactly one double root and otherwise simple roots has no odd prime unsafe degree. General triple/multiple repeated-center classification remains open.
- Rational inverse coordinate, height scale and amplitude are forced for every pure odd-start t^a(C+t^2D)/x. The full zero-center branch is safe and determined by just its second and fourth moments via a symbolic Z[a] determinant. Starting height3 is safe through degree5, and any pure unit-center phase with finite nonzero u at infinity has composite degree. These restrictions preserve the broader quintic/septic gaps.
- The entire odd-cubic winding-one class is safe at every odd least degree<=7; the retained F9 makes9 the exact least odd unsafe degree in this class. Every two-root unit center is excluded in all degrees. Prime29,31 and higher-degree multi-root centers remain research targets; this does not settle arbitrary quintics or septics.
- Every pure odd-cubic phase of exact least degree11 is safe; there is no unit-center phase at that degree. The retained F13 and all lower prime exclusions give exact least PRIME unsafe degree13 in this full coordinate class, while arbitrary Image phases remain outside scope.
- Every nonconstant pure unit-center odd-cubic phase has a unique normalized polynomial descent over its original coefficient field to a genuine rational-return quotient. Admissible polynomial covers preserve purity both ways and every composed x-only marker moment. Polynomiality of BOTH coefficients remains essential; the remaining quotient classification is open.
- Across every integer starting height a>=0 in the exact two-term profile, ALL-ORDER PURITY plus eventual vanishing for the two fixed original markers Z,Z^2 is equivalent to safety for every fixed original marker and to the explicit normal forms. Even a needs onlyZ. This is a safety theorem under purity, not a finite-prefix test or a full pure-unit classification.
- There is no unit-center pure odd-cubic phase of exact least degree17; all zero-center pure phases at that degree are safe. The claim uses the full canonical polynomial quotient and an exact rational-identity obstruction, not a moment prefix or a global degree17 Image classification.
- The minimum unsafe degree2a+4 is now exact at EVERY integer starting height a>=0, with unique normalized equality pair C=(1+lambda x)^(a+1),D=-(1+lambda x)^(a+3). A local common-root multiplicity forces the lower bound; the all-k witness was already retained. This settles the former a3 degree6--9 safety gap.
- Pure odd-cubic phases of exact least degrees17 and19 are safe, with no unit-center phases at either degree. Degree19 uses complete root-orbit rational normal forms and exact integer ideal certificates. Prime13 remains the sharp minimum in this full odd-cubic class; the individual next uncovered primes include29,31, while prime23 is retained.
- The least unsafe PRIME sphere/radial-primitive balanced degree is13 across EVERY integer two-height starting profile a>=0. Every nontrivial even-start pure phase has EVEN degree; height3degree11 reduces to three excluded a1 profiles. This strengthens the earlier a1 prime threshold without extending the separate degree17/19 a1 exclusions to all heights.
- Higher odd starting heights do not universally reduce polynomially to a1: explicit pure primitive nonpower degree15 and nonconstant-Phi prime59 examples violate the necessary all-root divisibility. The generic all-height admissible pullback instead preserves all polynomial-composed marker moments and proves purity equivalence in both directions.
- All Phi1 profiles are in marked DZ bijection; for EVEN starts this completes purity classification. At even a>=2 the unsafe degree spectrum is(a+2){2,3,...}; at a0 only4. Geometric unit classes at degP=r have Fuss-Catalan count1/r*binom(ar,r-1), with exact canonical marked Galois descent. Explicit geometric classification of nonconstant-Phi odd-start pure families remains open; exact individual recognition and degree-only finite prefixes are now supplied separately.
- Every nontrivial Phi1 geometric family now has a radial-primitive NONPOWER original lift over every field extension, by exact T/L factor orders and gcd(a,ar+1)=1. Thus the Fuss-Catalan count is also a count of primitive nonpower Image families under the specified source-scaling equivalence.
- All pure ODD-start unit phases now have a unique primitive rational-composition normal form, with exact finite C/D polynomiality and generic-fiber primitivity tests. This does not impose a polynomial source or polynomial outer coefficients, or enumerate all admissible triples.
- Every odd-start unit-center two-height input now has an exact MARKED primitive-gcd purity criterion, valid before rationality of the selected square root or inverse coordinate. The full common gcd is removed from BOTH polynomials before the nonzero residual-resultant cutoff; an off-base common factor alone is explicitly insufficient.
- For every integer starting heighta>=0 and least sphere degree<=d, the first B_d=3 floor((d+5)^2/4)(floor((d+3)^2/4)-1)-4 moments now determine all-order purity in the exact two-height class. This gives equality of RADICAL moment ideals for each fixed coefficient box, not ordinary ideal equality or a nilpotence exponent. Delayed examples prove that degree dependence is necessary.

Tracked obligations (computed states):

- **formalization / open:** [Formalize the complete binary GVC proof, including the moving-envelope argument.](formal/gvc/README.md) — FORMAL-BINARY-GVC-COMPLETE
- **formalization / open:** [Formalize the full generic tangent-pencil geometric theorem beyond the checked tangent, contact and collision algebra.](formal/discriminant-pencils/README.md) — FORMAL-DISCRIMINANT-COMPLETE
- **verification / open:** [Respect the named imported theorem inputs and the analytic statements outside the focused finite-etale certificate.](README.md) — GAP-FOUNDATIONS-2
- **mathematical / open:** [Classify the remaining pure rational-return odd-cubic quotients and admissible polynomial covers. In the regular winding-one class, pure phases of exact least degree 17 or 19 are safe and have no unit-center cases. General higher quotients and covers remain open, including the individual odd-cubic questions at 29 and 31. The unrestricted degree-17 and degree-19 Image questions remain open. Exact membership recognition and degree-bounded purity prefixes are now settled by the separate gcd/resultant theorem; this question asks for geometric classification.](extended-geometry/ODD_CUBIC_DEGREE_NINETEEN_SAFETY.md) — GAP-HOPF-ODD-CUBIC-MULTIROOT-CENTERS
- **mathematical / open:** [Classify which varying normalized rational endpoint coordinates carry nonconstant-Phi pure unit centers at odd starting heights, and give a geometric classification of the resulting families. For a supplied endpoint the all-height atlas and exact reconstruction are already proved.](extended-geometry/ODD_START_RATIONAL_COMPOSITION_NORMAL_FORM.md) — GAP-HOPF-ODD-START-UNIT-CENTER-CLASSIFICATION

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
