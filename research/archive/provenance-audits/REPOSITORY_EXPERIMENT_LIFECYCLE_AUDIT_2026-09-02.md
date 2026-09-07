# Repository experiment-lifecycle audit — 2026-09-02

This is a maintenance and provenance audit, not a mathematical status source.
`MATH_STATUS.json` remains the sole authority.

## Scope and method

The pass covered all 6,824 tracked files, including 933 Markdown files and
1,064 typed status entries.  It classified 503 active programs whose names
identify them as searches, scans, probes, experiments, diagnostics, batches,
sweeps, or enumerations.  Classification used:

- repository and nested `AGENTS.md` instructions;
- status roles (`canonical_source`, `checker`, `software_lock`, and
  `consumers`);
- Makefile and reproduction-catalogue entry points;
- active-document references and the Elkies--K3 process/success-path ledgers;
- explicit bounded/partial/proved language in experiment notes; and
- tracked backup names and archive placement.

This was a lifecycle audit, not a claim that every experimental calculation
was independently recomputed.  Negative experiments that remain cited,
status-locked, or useful as route diagnostics were retained.  A bounded miss
was never promoted to a proof.

## Corrections

1. `HC4_FINITE_FIELD_SEARCH.md` still named simultaneous cubic/sextic
   interaction as the next experiment.  The later support-free chain
   `HC4T31`, `HC4T21`, `HC4T11`, `HC4TC1` closes the coordinate
   cubic--quartic--sextic chart.  The note now identifies quintic/higher layers
   and non-coordinate coisotropic embeddings as the live boundary.
2. Several 2026-08-20 through 2026-08-22 Elkies--K3 documents still used
   “current frontier”, “active lifting target”, or model-recovery language.
   They now identify themselves as historical snapshots and point to the
   completed q12/orbit5867 rootless endpoint and the residual 2-Selmer gate.
3. The active K3 README's late “Next strategic gate” heading contradicted its
   own current-priority section.  It is now a completed route handoff.

## Archived material

Eighteen unreferenced early K3 experiment notes or superseded branch audits were moved to
[`archive/elkies-k3/early-experiments/`](../elkies-k3/early-experiments/).
At the time of the move they had no mathematical-status role, software lock,
status consumer, or inbound active-document reference.  Their original paths
and archived-byte hashes are pinned in that directory's `MANIFEST.tsv`.

Three unreferenced Fermigier search/diagnostic notes were moved to
[`archive/elliptic-curves/notes/`](../elliptic-curves/notes/).  Two record
bounded misses and one records a superseded pre-descent priority list.  They
remain computational provenance, not active search instructions; their paths
and hashes were appended to the elliptic-curve archive manifest.

Eleven timestamped q24 `*.bak-*` source snapshots were moved out of the active
script directory into
[`elkies-k3/scripts/archive/q24-bak-20260823/`](../../elkies-k3/scripts/archive/q24-bak-20260823/).
A second manifest records those moves.  The default repository-hygiene audit
now rejects future tracked backup snapshots outside an explicit archive.

## Literature and public-source check

- Mondello's characteristic-two plane result remains arXiv v1, submitted
  2026-07-29; the repository's explicit `v1` attribution is current:
  [arXiv:2608.02634](https://arxiv.org/abs/2608.02634).
- Padurariu--Saia is currently arXiv v2 (2025-10-28).  The repository uses the
  version-independent URL and separately pins the level-474 source data:
  [arXiv:2509.25368](https://arxiv.org/abs/2509.25368).
- The public rank leaderboard lists curve 302 as rank at least 31, consistent
  with `ECR31`; no rank-at-least-32 result was found in the checked source:
  [ICARM curve 302](https://elliptic-rank.icarm.cloud/curve/302).
- The Pasten--Salgado theorem used for non-thin rank jumps remains the
  published 2024 theorem cited by the canonical note:
  [DOI 10.1007/s00229-024-01554-2](https://doi.org/10.1007/s00229-024-01554-2).

No mathematical-status entry was changed solely from this literature pass.

## Verification

The maintenance pass used:

```bash
make check
python3 elkies-k3/scripts/analyze_process_ledger.py --check-document
python3 elkies-k3/scripts/success-path/verify_ledger.py
awk -F '\t' 'NR>1 {print $3 "  " $2}' \
  archive/elkies-k3/early-experiments/MANIFEST.tsv | sha256sum -c -
awk -F '\t' 'NR>1 {print $3 "  " $2}' \
  archive/elliptic-curves/MANIFEST.tsv | sha256sum -c -
awk -F '\t' 'NR>1 {print $3 "  " $2}' \
  elkies-k3/scripts/archive/q24-bak-20260823/MANIFEST.tsv | sha256sum -c -
```

Long Sage, Singular, Macaulay2, Magma, Julia, Lean, LaTeX, and expensive
symbolic replays were not rerun because no proof, certificate, or experiment
output changed.
