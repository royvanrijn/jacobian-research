# AGENTS.md

Research repository. Keep changes reproducible, claims fail-closed, and navigation concise.

## Research posture

The active programme is elliptic curves and their supporting K3 constructions.
Other projects are preserved under `research/archive/non-elliptic/`; their old
runbooks and TODOs are historical, not an active work queue. Large searches,
descents, neighbour enumerations, specialization sweeps, and comparable
compute-heavy campaigns must still be explicitly in scope, mathematically
motivated, checkpointed, and reproducible.

## Authority order

1. `research/MATH_STATUS.json` — mathematical status.
2. Canonical proof/source notes referenced by that file.
3. Generated certificates under `research/artifacts/generated-results/`.
4. Exploratory notes.
5. `research/archive/` — historical context. Archived programmes retain canonical
   evidence through the live ledger, with `programme_status: archived`; this flag
   is separate from mathematical state.

`research/STATUS.md` is generated. Regenerate it through the repository status
renderer; never edit it by hand.

## Find prior work before computing

Start with `research/README.md`, `research/KNOWLEDGE_BASE.md` and the generated
catalogue in `research/index/`. Search the current tree before designing a
calculation, including historical failure reasons:

```sh
python3 research/scripts/research.py search "your topic"
python3 research/scripts/research.py show CLAIM-ID
python3 research/scripts/research.py routes --area elliptic-curves
python3 research/scripts/research.py work --area elkies-k3
```

Read full scopes, replacements and retained checkpoints. Dated runbooks,
historical inventory totals and `LIVE_STATUS.json` snapshots are not current
instructions or evidence that a process is running. Missing local artifacts
must not trigger automatic reconstruction.

When work changes an algorithm or invalidates an approach, update the sourced
record in `research/knowledge/lessons.json`: applicability, implementation,
failed approach, exact boundary and the condition for revisiting it. Preserve
the underlying experiment. `make render-navigation` regenerates status and
indexes; `make check-navigation` checks metadata and links without research
calculations. Generated navigation is not a second mathematical authority.

Use `research/knowledge/WORK_LEDGER.md` for unknowns and suggested next gates.
Maintain `research/knowledge/work_items.json` with the canonical claim, next
step, completion evidence and prerequisites; review affected scope fingerprints
after a result changes. Work items are unscheduled proposals. Every inherited
retrospective checkbox is preserved in `research/knowledge/legacy_work_review.json`;
do not recreate its completed subsets as new tasks. `research.py resources`
also exposes K3 process mechanisms and the curve inventory. Use `--history`
to include other archived programmes. `show ID` can read any retained claim.

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
2. keep the current theorem in `research/MATH_STATUS.json` and its canonical proof note;
3. move historical narrative to `research/archive/`, or retain it at the old path only when current certificates/ledgers depend on that path and index it as historical;
4. retain scripts/certificates when they are useful regressions;
5. avoid deleting evidence just because it is no longer operational.

The pre-streamlining repository state is pinned in [`research/archive/STREAMLINING_2026-09-03.md`](research/archive/STREAMLINING_2026-09-03.md).

For K3-specific work, follow [`research/elkies-k3/AGENTS.md`](research/elkies-k3/AGENTS.md). For elliptic-curve work, follow [`research/elliptic-curves/AGENTS.md`](research/elliptic-curves/AGENTS.md).
