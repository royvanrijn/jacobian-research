# AGENTS.md

Research repository. Keep changes reproducible, claims fail-closed, and navigation concise.

## Research posture

All programmes are open for theorem-directed breakthrough work. Large searches,
descents, neighbour enumerations, specialization sweeps, and comparable
compute-heavy campaigns must still be explicitly in scope, mathematically
motivated, checkpointed, and reproducible.

## Authority order

1. `MATH_STATUS.json` — mathematical status.
2. Canonical proof/source notes referenced by that file.
3. Generated certificates under `artifacts/generated-results/`.
4. Exploratory notes.
5. `archive/` — historical context only.

`STATUS.md` is generated. Regenerate it through the repository status renderer; never edit it by hand.

## Editing discipline

- Preserve proof artifacts and replay inputs.
- Do not infer a theorem from a bounded search or heuristic score.
- Keep `UNKNOWN` as `UNKNOWN` until an exact certificate closes it.
- Archive or replace superseded handoffs with short navigation/tombstone records instead of maintaining parallel status narratives.
- Keep README pages short and link to canonical proofs.
- Prefer narrow, cheap checks. Run expensive whole-suite or research calculations only when mathematically necessary and explicitly in scope.
- Do not rewrite unrelated active proof documents merely for stylistic consistency.

## Repository hygiene

When a result supersedes an old route or handoff:

1. update the active README/navigation surface;
2. keep the current theorem in `MATH_STATUS.json` and its canonical proof note;
3. move historical narrative to `archive/`, or retain it at the old path only when current certificates/ledgers depend on that path and index it as historical;
4. retain scripts/certificates when they are useful regressions;
5. avoid deleting evidence just because it is no longer operational.

The pre-streamlining repository state is pinned in [`archive/STREAMLINING_2026-09-03.md`](archive/STREAMLINING_2026-09-03.md).

For K3-specific work, follow [`elkies-k3/AGENTS.md`](elkies-k3/AGENTS.md). For elliptic-curve work, follow [`elliptic-curves/AGENTS.md`](elliptic-curves/AGENTS.md).
