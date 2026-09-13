# Verified foundry additions to the curve ledger

> **Historical fixed-ledger snapshot.** The high-rank foundry and its
> downstream conductor interface are stopped. Their queue and controller
> details below preserve the provenance of this snapshot; they do not authorize
> factorization, a launch, a resume, or an automatic successor. Use the
> [current inventory](../INVENTORY.md) and [elliptic-curve programme](../README.md)
> for current work.

The first fixed foundry publication snapshot adds **30 curves with certified
lower bounds at least22** and strengthens the existing `08f72` fibre at
`1245/2519` from23 to25. These additions are distinct overQ from each other
and from the prior291-row inventory. Equal-j twists are compared by exact
rational isomorphism, not merged by j alone. No worldwide novelty or exact
rank is asserted.

The new curves comprise nine lower bounds22, fourteen23, two24, four25 and
one26. Four come from commissioning and26 from the production runs at this
cutoff. The strongest is family `07ca9` at `-1508/909`, with26 independent
rational points. The earlier commissioning fibre `07ca9` at `-943/412` is
also included with lower bound25. The existing rank27 fibre at `877/781`
is not counted as a new discovery.

The full ledger grows from291 to**321 curves**; the main README highlights
**274**, including the four retained below22 structural examples. New curves
retain their exact native equations, point lists and portable proofs. Their
conductors and minimal-model height metrics remain `UNKNOWN`/uncomputed;
no heuristic conductor screen is promoted to an exact value. The publication
therefore contains195 exact conductors and126 unresolved, using the existing
third conductor snapshot. The then-running factorization pass kept its original
roster; these30 additions had a separate frozen successor queue.

## Evidence and selection

The immutable [selection snapshot](../../artifacts/generated-results/elliptic-curves/foundry_curve_ledger_snapshot_v1.json)
embeds each selected job export, its source path/hash, its previous bound when
applicable, and the pre-update curve models used for deduplication. Its SHA-256
is `52be94faff04dd66cf19c6a2696414b0cb96cb0659e80520ccdbba5b85827628`.

Selection reads sealed `PASS_CERTIFIED_SEARCH` exports from the commissioning
and production directories, requires matching packet hashes and two-replay
receipts, then retains the strongest bound per rational isomorphism class.
Only bounds at least22 are admitted. Lower-rank ordinary discoveries, failed
jobs and unresolved generic seeds stay in their original search archives.
Equal-bound continuations do not replace existing ledger entries. The list of
source paths is fixed before replay; later discoveries need another snapshot.

The [publisher/checker](../cas/refresh_foundry_curve_ledger.py) freshly replays
all31 selected packets through both existing finite-group replay paths,
verifies the native MW17 prefix and rational2-torsion
exclusion, and checks pairwise/baseline rational nonisomorphism. No point
search, prospective selection change or conductor factorization occurs.
All rank values are certified subgroup lower bounds. Existing source packets
and their complete adaptive histories remain untouched.

Publication replay uses a fresh RAM-backed arithmetic cache to avoid thousands
of synchronous temporary fact-file writes. This changes storage only: complete
finite quotients and all exact point/reduction checks are still recomputed.

For the stopped experiment and its retained discoveries, see the
[foundry note](HIGH_RANK_SEARCH_FOUNDRY_2026-09-09.md). Automatic job exports
never rewrite this fixed publication snapshot or the mathematical-status register.

## Retained conductor interface

The [successor queue](../cas/run_foundry_conductor_queue_v1.py) freezes these30
new curves only, ordered by decreasing certified lower bound and then ID. Its
foreground Sage model-transport preflight passed all30. Its former detached
waiting-controller condition is retained as historical protocol, not an
automatic action.

The successor uses the same four workers,1,800-second build cap,120-second
independent replay cap and2GiB per process. Original equations and all older
proofs remain unchanged. Partial results stay `UNKNOWN`; engineering/replay
failures require review. If the predecessor stops early, this queue waits
rather than silently bypassing it. This completed interface has no live
controller or relaunch path; publication remains a separate replayed
snapshot/index/render step.
