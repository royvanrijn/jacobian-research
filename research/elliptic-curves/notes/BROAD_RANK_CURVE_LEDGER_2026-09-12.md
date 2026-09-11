# Broad-rank curve ledger and first conductor cutoff

The [completed broad campaign](BROAD_RANK_SEARCH_2026-09-10.md) contributes
**30 distinct curves with certified lower bound above 22** to the prior
415-curve ledger. Fresh replay verifies 25 lower bounds of 23, four of 24,
and one of 25. All 30 are pairwise nonisomorphic over Q and absent from the
pinned 415-curve baseline.

The [main README](../../../README.md#elliptic-curve-inventory) and
[complete inventory](../INVENTORY.md) now contain **445 curves**, with
**267 exact conductors and 178 unresolved**. The README displays 398 curves
under its existing rank and structural-example selection rule.

The portable [curve snapshot](../../artifacts/generated-results/elliptic-curves/broad_rank_ledger_snapshot_v1.json)
retains all selected point packets, normalized parents, final states, source
bindings and baseline equations. The [fresh replay](../../artifacts/generated-results/elliptic-curves/broad_rank_ledger_replay_v1.json)
checks both existing finite-group certificate paths, all point memberships,
rational 2-torsion exclusion, specialized parent equations and generic prefixes.
The original frozen campaign, its runtime and all replay inputs remain intact.
These are subgroup lower bounds; exact ranks and worldwide novelty are unknown.

All 30 curves have entered the separate
[bounded conductor queue](../cas/run_broad_rank_conductors.py), ordered by
certified lower bound. Four detached workers each have a 1,800-second build
limit, a 120-second independent replay limit and a 2 GiB memory limit. Exact
denominator-clearing transports passed the zero-search Sage preflight for every
input. The existing conductor regression suite passed all five tests, including
rational transport and partial-result handling. No point search is repeated.

The [first conductor snapshot](../../artifacts/generated-results/elliptic-curves/broad_rank_conductor_snapshot_v1.json)
freezes the first 11 completed certificates, all exact, including the
rank-at-least-25 curve `broad-11952-7bb187bc9254e81c6283` at parameter `921/653`.
Its [fresh publication replay](../../artifacts/generated-results/elliptic-curves/broad_rank_conductor_replay_v1.json)
verifies exact model transport, prime certificates and agreement of Sage's
generic Tate algorithm with PARI local reduction, without factor discovery.
Only these 11 conductor results enter this ledger cutoff. The other 19 remain
unknown in the ledger while the live queue continues; later completions require
another fixed, replayed publication snapshot. Minimal-model height and
discriminant columns remain uncomputed for all 30 new rows. No conductor-record
claim is made.

From the repository root, replay the publications and check the generated views:

```sh
python3 research/elliptic-curves/cas/publish_broad_rank_results.py check
~/.local/bin/sage -python research/elliptic-curves/cas/publish_broad_rank_conductors.py check
python3 research/elliptic-curves/cas/render_main_readme_curves.py --check
python3 research/scripts/render_status.py --check
```

The immutable snapshots are already frozen; do not rerun `freeze` over them.
The two publishers' `index` commands bind the replayed snapshots to the proved
entries in `research/MATH_STATUS.json`; the existing inventory renderer then
regenerates the JSON, CSV, curve pages and README table.

Inspect the continuing queue with:

```sh
python3 research/elliptic-curves/cas/run_broad_rank_conductors.py status
```
