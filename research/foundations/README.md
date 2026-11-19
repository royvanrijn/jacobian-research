# Formal verification

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Precisely scoped Lean developments supporting published papers and current manuscripts.

**Authored assessment:** Formal scope varies by project. GMC2 is conditional on two named axiom interfaces; binary GVC and full generic discriminant geometry remain partial.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| GVC counterexamples and binary supporting algebra | Partial result | proved | [Source](formal/gvc/README.md) |
| Binary Gaussian moments in Lean | Internal result | conditional | [Source](formal/gmc2/README.md) |
| Characteristic-zero finite etale Keller construction | Internal result | proved | [Source](formal/finite-etale-keller/README.md) |
| Tangent-pencil algebra in Lean | Partial result | proved | [Source](formal/discriminant-pencils/README.md) |
| Complete complex quadratic-height classification and exact finite moment ideals | Internal result | proved | [Source](extended-geometry/HOPF_LIFT_CLASSIFICATION.md) |
| Complete quadratic-height Hopf purity and safety classification | Internal result | proved | [Source](extended-geometry/QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md) |
| Odd cubic-height rational reconstruction and sharp low-degree safety | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md) |
| Complete one-pole odd-cubic purity and safety classification | Internal result | proved | [Source](extended-geometry/ONE_POLE_ODD_CUBIC_CLASSIFICATION.md) |
| Exact/logarithmic purity reduction and one-double-root prime obstruction | Internal result | proved | [Source](extended-geometry/MIXED_CENTER_LOGARITHMIC_SEPARATION.md) |
| Odd-start rationality, full zero-center safety and a quintic exclusion | Internal result | proved | [Source](extended-geometry/ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md) |
| Odd-cubic septic safety and the sharp odd unsafe threshold nine | Internal result | proved | [Source](extended-geometry/ODD_CUBIC_DEGREE_SEVEN_SAFETY.md) |

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
- The entire odd-cubic winding-one class is safe at every odd least degree<=7; the retained F9 makes9 the exact least odd unsafe degree in this class. Every two-root unit center is excluded in all degrees. Prime11 and higher-degree multi-root centers remain research targets; this does not settle arbitrary quintics or septics.

Tracked obligations (computed states):

- **formalization / open:** [Formalize the complete binary GVC proof, including the moving-envelope argument.](formal/gvc/README.md) — FORMAL-BINARY-GVC-COMPLETE
- **formalization / open:** [Formalize the full generic tangent-pencil geometric theorem beyond the checked tangent, contact and collision algebra.](formal/discriminant-pencils/README.md) — FORMAL-DISCRIMINANT-COMPLETE
- **verification / open:** [Respect the named imported theorem inputs and the analytic statements outside the focused finite-etale certificate.](README.md) — GAP-FOUNDATIONS-2
- **mathematical / open:** [Classify the remaining odd-cubic centers with at least three geometric roots and repeated factors; determine the residual prime11 case beyond the new sharp odd threshold9.](extended-geometry/ODD_CUBIC_DEGREE_SEVEN_SAFETY.md) — GAP-HOPF-ODD-CUBIC-MULTIROOT-CENTERS
- **mathematical / open:** [Classify pure unit-center odd-start phases beyond the exact known inverse-coordinate fibers, especially zero/pole behavior at infinity and the remaining low-prime cases.](extended-geometry/ODD_START_HEIGHT_RATIONALITY_AND_QUINTIC_SAFETY.md) — GAP-HOPF-ODD-START-UNIT-CENTER-CLASSIFICATION

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
