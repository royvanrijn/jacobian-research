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

**Computed obligations:** 9 unresolved; 1 resolved within their exact scope.


Scope limits:

- Formal verification is scoped by exact declarations, transitive assumptions and the paper bridge. A paper theorem does not imply that its full formalization is complete.
- GMC2 retains two explicit mathematical axiom interfaces, shown as open formalization obligations. Its propagated formal claim remains conditional until discharged; its written mathematical theorem has separate proof support.
- Binary GVC and generic tangent-pencil geometry are partial formalizations. The focused finite-etale certificate excludes the broader supplied-presentation and analytic counting scopes. No Lean build was rerun for this registry migration.
- The retained two-pair Image cubic proof covers exact bihomogeneous bidegree (3,3), with sharp mixed cutoff and minimum failing balanced degree four. Independent review and formal coverage are separate obligations.
- The entire complex polynomial winding-one quadratic-height numerator class is classified. On the nonzero constant-coefficient branch, its ordinary full moment ideal is exactly the finite coefficient ideal of (C-xCprime)^2D+C^3. The higher Hopf residual gap is discharged at this scope; general profiles and nilpotence indices are separate.
- The full polynomial winding-one quadratic-height Hopf family, including arbitrary linear height, is now classified over every characteristic-zero field. Its only unsafe branch has even least balanced degree and a fixed marker surviving at every power. This strengthens the earlier B=0 field classification while preserving its separate localized ordinary-ideal statement; arbitrary quintic and septic phases remain open.
- In the entire winding-one odd cubic-height class t*x^-1(C+t^2D), pure phases have an automatic rational center on the nonzero-constant branch. Every pure phase of least degree<=5 has an explicit one-sided safe normal form; the retained degree6 example makes the class threshold sharp. Finite Bezout identities follow the all-order rationality reduction and do not provide a universal moment-prefix cutoff.

Tracked obligations (computed states):

- **formalization / open:** [Formalize the complete binary GVC proof, including the moving-envelope argument.](formal/gvc/README.md) — FORMAL-BINARY-GVC-COMPLETE
- **formalization / open:** [Formalize the full generic tangent-pencil geometric theorem beyond the checked tangent, contact and collision algebra.](formal/discriminant-pencils/README.md) — FORMAL-DISCRIMINANT-COMPLETE
- **verification / open:** [Respect the named imported theorem inputs and the analytic statements outside the focused finite-etale certificate.](README.md) — GAP-FOUNDATIONS-2

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
