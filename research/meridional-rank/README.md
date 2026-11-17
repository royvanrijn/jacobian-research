# Meridional rank

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Knot-group representations, branched-cover bounds and explicit sharp families.

**Status:** Explicit sharp examples, infinite families and lower bounds are established. The full conjecture and global novelty claims remain open.

| Current result | Status | Main source |
|---|---|---|
| Twenty named equalities and established-family scope | Internal result | [Source](representations/ORIGINAL_SIXTEEN_KNOWN_GRAPH_COVER_FAMILY.tex) |
| Infinite periodic families and connected sums | Internal result | [Source](representations/PERIODIC_TWIST_FAMILIES.tex) |
| Representation and branched-cover lower bounds | Internal result | [Source](representations/PROOF.tex) |

Open questions:

- Full Meridional Rank Conjecture.
- Sharper lower bounds beyond known deficient representation/cover methods.
- Precise novelty and existing-family comparison for the four other sharp seeds and infinite periodic families.

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses.

The [current verification record](representations/graph_cover_family_scope_verification.json)
records an internal replay of eight current finite-certificate entry points and
[the compact sharp-five named-diagram audit](representations/verify_diagram_audit.py).
Run a listed checker from the repository root with
`python3 research/meridional-rank/representations/CHECKER.py`, replacing
`CHECKER.py` with its recorded filename. These checks cover their exact retained
inputs and controls; the uniform theorems also require the written arguments
and explicit classical imports.

The restored sharp-five, branched-cover, parabolic and first-witness subgroup
certificates preserve their earlier source/input identities. Their original
receipts report historical executions; only the nine entry points in the current
verification record were replayed during this cleanup. The original workbook,
complete 512-matrix run, broader 243/293-knot campaigns and removed companion
notes remain in Git history. In particular, `verify_integer_lifts.py` retains
helper routines but its original full-campaign entry point requires those
undistributed inputs. Historical receipt hashes were not replaced to imply a
fresh run. Formal verification, independent replay and external review remain
separate from these internal checks.
