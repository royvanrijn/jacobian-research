# Elliptic-curve research instructions

This directory inherits the repository rules. [MATH_STATUS.json](../MATH_STATUS.json)
is the sole mathematical-status authority. The programme remains open for
theorem-directed work.

## Required preflight

Read the [programme map](README.md), [method memory](../KNOWLEDGE_BASE.md),
[shared runtime](notes/SHARED_RESEARCH_RUNTIME.md) and the relevant canonical
claim. For rank-jump work also read the
[structural reassessment](notes/RANK_JUMP_REASSESSMENT_2026-09-05.md) and
[implemented performance improvements](notes/V3_FUTURE_SEARCH_PERFORMANCE_2026-09-08.md).

Curve302's [alternative MW17 parent](notes/CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md)
is complete. Use its equation and basis loader before reconstructing them.
The current rank targets and mechanism boundary are `OP-EC-NEXT` and
`OP-EC-RANK-JUMP-MECHANISM-20260910`; do not duplicate their chronology here.

## Search and proof discipline

- Large rank, parameter, point, class-group and descent campaigns require an
  explicit mathematical gate, declared limits, checkpoints and a certificate plan.
- Exactly verified independent points prove an unconditional lower bound.
  Exact rank requires matching unconditional lower and upper bounds.
- Only an unconditional certified upper bound below a target mathematically
  excludes a fibre. Incomplete or conditional descent, scores and timeouts
  can schedule bounded work but cannot supply an exclusion.
- Keep candidate incidence, conditional point visibility and global cover
  solubility separate. Prove the generic geometry of a deformation first.
  Preserving the cubic field does not preserve the Mordell–Weil group.
- State the certified subgroup and base change in every jump count. Promoting
  anchor directions to the generic subgroup consumes those quotient directions.
- Deduplicate PGL2/Weierstrass-equivalent fibrations before subgroup
  transversality tests; Curve398's duplicate MW16 pair is a retained regression.
- Keep known-record equations, parameters, points, ranks, target j-invariants
  and jump labels out of prospective A1/MW16 selection and execution. Masked
  or retrospective controls remain separate and are not new rank discoveries.
- Normalize only when the arithmetic endpoint needs it. Raw pointed search
  needs neither global minimality nor factorization. Preserve both coordinate
  maps and exact transports where their bounded boxes differ.
- Reuse exact arithmetic contexts, finite-reduction signatures, sealed
  landscapes and compatible chart receipts. New bases require new landscape
  verification; missing cached data remain missing.
- Before calculating, search the active claim and lesson index for the object
  and method. A matching proof, certificate, or sealed landscape is the
  baseline; record the different input or theorem gate before repeating work.
- Conductor is certified through minimal-model/local reduction data, not a
  discriminant radical. Numerical heights and scores retain their stated scope.
- Preserve raw checkpoints under ignored local-artifact paths and compact
  certificates under `../artifacts/generated-results/elliptic-curves/`.
  Keep versions, inputs, limits, timeout semantics and failed runs.
- A certificate input hash identifies the bytes used at its production time.
  Never retag a frozen certificate to match changed source code; produce a
  successor certificate or retain an explicit historical boundary instead.
- Add reusable lessons to `../knowledge/lessons.json`, then regenerate navigation.
  An old runbook is not authorization to restart its completed campaign.

## Hygiene passes

- Keep hygiene changes small: repair a stale pointer, duplicate status narrative,
  or missing provenance link at a time.
- For an already reviewed source, stop when its status entry, source hash,
  theorem boundary, and artifact/review record agree. No identified defect means
  no source rewrite.
- Do not condense a canonical proof or certificate note unless a specific
  navigation or status defect requires it; preserve its exact scope and replay
  boundary.
- Use static validation for organisation-only changes.  A replay is separate
  mathematical work and must be explicitly in scope.

For K3 construction work follow the [K3 instructions](../elkies-k3/AGENTS.md).
