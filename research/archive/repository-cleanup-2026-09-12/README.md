# Repository cleanup — 12 September 2026

This archive preserves the navigation state before the current cleanup.
[MANIFEST.json](MANIFEST.json) records the base revision, original paths,
preserved paths and SHA-256 hashes. The snapshots were copied from the working
tree and include contemporaneous principal28 research updates. Their priorities,
inventory counts and claims of pending work are historical.

The cleanup replaces long READMEs with short programme maps, generates a
complete claim/source catalogue, and extracts sourced methods into
[the algorithmic memory](../../knowledge/ALGORITHMS.md). The authority's
failed-route records are exposed in [a generated index](../../knowledge/FAILED_ROUTES.md).
No mathematical state was promoted by this cleanup.

The manifest also binds the [LND status baseline](LND_STATUS_BEFORE_SOURCE_REVIEW.json)
and [GVC/Gaussian/Danielewski status baseline](GVC_GAUSSIAN_DANIELEWSKI_STATUS_BEFORE_SOURCE_REVIEW.json),
plus the [K3 bisection baseline](K3_BISECTION_STATUS_BEFORE_SOURCE_REVIEW.json),
from before the source-level corrections. The [partial review](../../knowledge/PARTIAL_REVIEW.md)
states its actual coverage; a passing navigation check does not complete it.

The follow-up [work ledger](../../knowledge/WORK_LEDGER.md) recovers all 77
unchecked retrospective items with explicit destinations and preserved original
wording. The September 10 high-rank synthesis is also preserved here; its old
path now points to completed campaigns and current canonical evidence. Its
inverse-parent construction recipe is retained in the algorithmic memory.

The old [OP-EC-NEXT chronology](OP-EC-NEXT.before.json) is retained in full.
Its current scope keeps the open target and narrowing references, corrects the
obsolete “unknown302 parent” description, and points inventory totals to the
generated curve database. All prior narrowing dependencies remain in the ledger.

The long replay guide moved to [replay/CATALOGUE.md](../../replay/CATALOGUE.md).
Only its heading, historical notice and relative Markdown links changed; commands
still run from `research/`. The two canonical H3 q24 sources and the document's
consumer bindings moved with it. The old [REPRODUCE.md](../../REPRODUCE.md) path
now gives the short guide.

Snapshots use `.txt` to preserve original Markdown bytes and relative-link text
without presenting those old links as active navigation. Search them with
`python3 research/scripts/research.py search "topic" --history` from the root.

Proof notes, mathematical checkers, certificates, point packets, generated curve
data and checkpoint inputs retain their paths. Historical artifacts are not
deleted to make an active programme look simpler. Concurrent research work is
kept separate from this navigation migration.

Verify the preserved bytes and replay relocation without computations:

```sh
python3 research/scripts/audit_cleanup_snapshot.py
make check-navigation
```
