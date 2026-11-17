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

**Computed obligations:** 10 unresolved; 0 resolved within their exact scope.


Scope limits:

- Formal verification is scoped by exact declarations, transitive assumptions and the paper bridge. A paper theorem does not imply that its full formalization is complete.
- GMC2 retains two explicit mathematical axiom interfaces, shown as open formalization obligations. Its propagated formal claim remains conditional until discharged; its written mathematical theorem has separate proof support.
- Binary GVC and generic tangent-pencil geometry are partial formalizations. The focused finite-etale certificate excludes the broader supplied-presentation and analytic counting scopes. No Lean build was rerun for this registry migration.
- The retained two-pair Image cubic proof covers exact bihomogeneous bidegree (3,3), with sharp mixed cutoff and minimum failing balanced degree four. Independent review and formal coverage are separate obligations.

Tracked obligations (computed states):

- **formalization / open:** [Complete the binary GVC formal theorem and the full generic tangent-pencil geometric theorem.](README.md) — FORMAL-BINARY-GVC-COMPLETE
- **verification / open:** [Respect the named imported theorem inputs and the analytic statements outside the focused finite-etale certificate.](README.md) — GAP-FOUNDATIONS-2

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
