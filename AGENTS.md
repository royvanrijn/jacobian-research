# AGENTS.md

Research repository. Keep changes reproducible, claims fail-closed, and navigation concise.

## Current compute priority

**Prime gaps.** The elliptic-K3 / high-rank elliptic-curve programme is paused as of 2026-09-03. Do not start expensive K3, elliptic-neighbour, descent, specialization, or broad lattice searches unless that programme is explicitly resumed.

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
- Archive or replace superseded handoffs with short tombstones instead of maintaining parallel status narratives.
- Keep README pages short and link to canonical proofs.
- Prefer narrow, cheap checks. Run expensive whole-suite or research calculations only when mathematically necessary and explicitly in scope.
- Do not rewrite unrelated active branches merely for stylistic consistency.

## Repository hygiene

When a result supersedes an old route or handoff:

1. update the active README/navigation surface;
2. keep the current theorem in `MATH_STATUS.json` and its canonical proof note;
3. move historical narrative to `archive/` or leave a concise tombstone pointing at immutable Git history;
4. retain scripts/certificates when they are useful regressions;
5. avoid deleting evidence just because it is no longer operational.

For K3-specific work, follow [`elkies-k3/AGENTS.md`](elkies-k3/AGENTS.md). While the programme is paused, maintenance and small exact verification are fine; new compute campaigns are not.
