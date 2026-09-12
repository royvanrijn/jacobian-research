# Elliptic-curve cleanup — 12 September 2026

The other projects are archived. Active navigation, search, method memory and
work lists now focus on elliptic curves and the K3 constructions supporting
them. The archival move and the active partial-result review are complete.
**The repository has not received a complete independent proof and replay audit.**

This is a dated maintenance report. Current mathematics belongs to
[MATH_STATUS.json](MATH_STATUS.json) and its canonical sources.

| Completed work | Evidence |
|---|---|
| Moved 2,319 files without changing their bytes; preserved 901 archived claims and their mathematical states | [Archive manifest and navigation](archive/non-elliptic/README.md) |
| Reviewed all 77 active EC/K3 partial results: exact surviving gate, later results, checker purpose and input boundary | [Source review](knowledge/PARTIAL_REVIEW.md) |
| Made 55 active algorithmic lessons searchable, with sources, failed approaches and reopening conditions | [Method memory](knowledge/ALGORITHMS.md) |
| Condensed the 4,394-line K3 script map, 592-line CAS map and 1,227-line obsolete foundry shortlist; preserved the originals | [K3 scripts](elkies-k3/scripts/README.md) · [CAS map](elliptic-curves/cas/README.md) · [Foundry history](elkies-k3/LATTICE_FOUNDRY_EQUATION_FIRST_SHORTLIST_2026-09-02.md) |
| Recounted and compared every active claim against the September4 audit baseline and cleanup-start revision | [Snapshot and per-entry differences](archive/repository-cleanup-2026-09-12/EC_SCOPE_RECONCILIATION.json) |

The snapshot contains 616 EC/K3 claims: 515 proved, 77 partial, 3 open, 3 parked,
16 mathematically archived and 2 falsified. These are recorded states, not new
assurance granted by cleanup. Concurrent discoveries remain in the comparison.

## Results that agents should reuse

Curve302's alternative MW17 parent and saturated generic basis are complete.
Curve398's equivalent MW16 presentations, the q8/q12 endpoint and the direct
alternate11952 equation already have canonical sources. The maps now lead to
those results instead of sending agents back into construction recovery.

The legacy Mestre/Fermigier rank13 entry point now rejects execution. Use the
coherent-label rank11 generic lower-bound result and its explicit historical
provenance check. Historical source hashes were preserved.

The q323 candidate5887 lift already exists and gives the wrong 6A1 child.
NS0031's formal branch remains a local theorem, while its exact rational
marking is excluded. Determinant1236 already has its candidate cover and fibre
evaluations; only the particular CM branch-orbit identification remains open.
These boundaries and reusable methods are in the source review and lessons.

The CRT experiment's completed 2,560-row bounded miss is distinct from the
unrun 1,536-row family commitment. An old “running” handoff was corrected, and
the replacement worker/analyzer input mismatch is explicit. Missing successor
data do not justify repeating the completed experiment.

The T-first planner now represents unproved rational-point existence as
`null`, retaining `false` for an exact exclusion. Its classifications, order
and empty equation handoff are unchanged.

## What is still unfinished

Of the 50 inherited checklist items in the active scope, 20 have explicit
completion records and **30 remain unfinished**. The other 27 belong to the
archived projects. Archiving did not complete them. See the
[full reconciliation](knowledge/LEGACY_WORK_REVIEW.md) and current
[work ledger](knowledge/WORK_LEDGER.md).

The remaining items include broader theorem/dependency/checker-purpose reviews,
sufficient portable inputs for every external or historical replay, independent
witnesses for high-risk claims, and explicitly unscheduled long replays.
Historical H3 intermediate inputs and amended-runtime experiment chunks are
still missing at some required paths. Each inspected absence is named in its
source review; optional local evidence receipts are distinguished from maintained
proof sources. Existing theorems were not downgraded merely because a local
cache is absent, and absent evidence was not replaced by an invented result.

## Validation and compute

`make check`, strict partial-review coverage, Markdown links, archive preservation,
process-ledger validation, inventory-summary validation and all 46 maintenance
tests passed. The default check uses system Python and launches no CAS research.
The historical H3 prefix passed its 14 script and 14 artifact hash checks.

The only mathematical recalculation in this pass was the bounded coherent-label
control: 63 finite groups and 3,324 elements in 1.95 seconds. The T-first
metadata rerender changed only the point-existence typing described above.
No large search, descent, neighbour enumeration or census was rerun.

Future navigation checks require a fresh source review for every active partial
result and reject stale claim/source fingerprints. A passing navigation check
still does not certify every theorem or every historical replay mode.
