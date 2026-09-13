# Broad-rank curve ledger and completed conductor pass

The [completed broad campaign](BROAD_RANK_SEARCH_2026-09-10.md) contributes
**30 distinct curves with certified lower bound above 22** to the prior
415-curve ledger. Fresh replay verifies 25 lower bounds of 23, four of 24,
and one of 25. All 30 are pairwise nonisomorphic over Q and absent from the
pinned 415-curve baseline.

At this ledger's publication snapshot, the [main README](../../../README.md#elliptic-curve-inventory)
and [complete inventory](../INVENTORY.md) contained **445 curves**, with
**282 exact conductors and 163 unresolved**; the README displayed 398 under
its then-current rank and structural-example selection rule. Current totals
are generated in the linked inventory.

The portable [curve snapshot](../../artifacts/generated-results/elliptic-curves/broad_rank_ledger_snapshot_v1.json)
retains all selected point packets, normalized parents, final states, source
bindings and baseline equations. The [fresh replay](../../artifacts/generated-results/elliptic-curves/broad_rank_ledger_replay_v1.json)
checks both existing finite-group certificate paths, all point memberships,
rational 2-torsion exclusion, specialized parent equations and generic prefixes.
The original frozen campaign, its runtime and all replay inputs remain intact.
These are subgroup lower bounds; exact ranks and worldwide novelty are unknown.

All 30 curves completed the separate
[bounded conductor queue](../cas/run_broad_rank_conductors.py), ordered by
certified lower bound. Four detached workers each had a 1,800-second build
limit, a 120-second independent replay limit and a 2 GiB memory limit. Exact
denominator-clearing transports passed the zero-search Sage preflight for every
input. The existing conductor regression suite passed all five tests, including
rational transport and partial-result handling. No point search was repeated.
The pass finished on September 12 at 01:46 CEST after 62.36 minutes; all workers
have exited.

The [completed conductor snapshot](../../artifacts/generated-results/elliptic-curves/broad_rank_conductor_snapshot_v2.json)
contains **26 exact conductors and four partial `UNKNOWN` results**. All 30
pass [fresh publication replay](../../artifacts/generated-results/elliptic-curves/broad_rank_conductor_replay_v2.json):
exact model transport, prime certificates and agreement of Sage's generic Tate
algorithm with PARI local reduction, without factor discovery. The four partial
results reached the 30-minute build limit and retain certified conductor divisor
and upper bounds. No additional conductor pass has been launched.

| Unresolved parent | Parameter | Certified rank lower bound |
| --- | --- | ---: |
| 074d9 | 361/942 | 23 |
| 07ca9 | 345/286 | 23 |
| 08f72 | -514/143 | 24 |
| 11952 | 277/226 | 23 |

The preserved [first snapshot](../../artifacts/generated-results/elliptic-curves/broad_rank_conductor_snapshot_v1.json)
and its [replay](../../artifacts/generated-results/elliptic-curves/broad_rank_conductor_replay_v1.json)
contain the first 11 exact certificates, including the rank-at-least-25 curve
`broad-11952-7bb187bc9254e81c6283` at parameter `921/653`. That historical cutoff
had 267 exact conductors and 178 unresolved across the inventory; the completed
snapshot adds 15 exact conductors and preserves those earlier certificates.
Minimal-model height and discriminant columns remain uncomputed for all 30 new
rows. No conductor-record claim is made.

The immutable snapshots are already frozen. The
[rank publisher](../cas/publish_broad_rank_results.py),
[conductor publisher](../cas/publish_broad_rank_conductors_v2.py),
[inventory renderer](../cas/render_main_readme_curves.py), and
[status renderer](../../scripts/render_status.py) are retained validation
interfaces; they are not instructions to refreeze, rerun, or extend the
completed queue.
