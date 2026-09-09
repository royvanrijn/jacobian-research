# Curve inventory refresh — September 9

The [main README](../../../README.md#elliptic-curve-inventory) and
[expanded inventory](../INVENTORY.md) now index **321 distinct curves**,
with **195 exact conductors** and **126 unresolved conductors**, using the third
replayed conductor snapshot plus the [foundry additions](FOUNDRY_CURVE_LEDGER_2026-09-09.md). Rank entries
are certified lower bounds, not exact ranks.

The selected updates are:

| Source | Inventory change |
| --- | --- |
| [Curve52](CURVE52_SHORT_V3_SEARCH_2026-09-09.md), [curve113](CURVE113_PRODUCTIVE_V3_SEARCH_2026-09-09.md), [curve116](CURVE116_PRODUCTIVE_V3_SEARCH_2026-09-09.md), [curve200](CURVE200_PRODUCTIVE_V3_SEARCH_2026-09-09.md) | Four existing lower bounds rise from 25 to 26; the new points are transported to the existing certified minimal models. |
| [First fresh R17 cohort](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md) | Six new curves with lower bounds 24, 23, 17, 24, 25, 17 in family order 074d9, 07ca9, 08234, 08f72, 103b2, 11952. The 08f72 conductor is exact; three others have certified partial bounds. |
| [Second fresh R17 cohort](SECOND_FRESH6_SEED_COHORT_2026-09-09.md) | Six curves with updated lower bounds 20, 20, 22, 17, 26, 22 in the same family order. |
| [Lower-height cohort](LOWHEIGHT_FRESH6_SEED_COHORT_2026-09-09.md) | Six previously omitted certified curves, using their strongest retained subgroup packets. |
| [Frozen sixty-fibre panel](R17_SIXTY_SEED_COMPLEMENT_PANEL_2026-09-09.md) | All60 independently replayed curves retained, including one lower bound27, three26 and four25. |
| [First foundry ledger snapshot](FOUNDRY_CURVE_LEDGER_2026-09-09.md) | Thirty further distinct curves at lower bounds22–26; the existing08f72 fibre1245/2519 strengthens23→25. All31 selected packets replay. |
| [Bifibration conic](DET1092_BIFIBRATION_CONIC_SEEDS_2026-09-09.md) | One exported determinant1092 specialization of orbit47755 at conic parameter u=0, with lower bound 18. |

The [earlier supplement](INVENTORY_SEED_SUPPLEMENT_2026-09-08.md) remains
included. All 120 additions to the original 201-row inventory retain source
models and uncomputed minimal-model height/discriminant fields. An exact
conductor does not supply those missing metrics. Bounded misses remain
unresolved, and the inventory makes no current public novelty or record claim.

## Main README curation

The main README highlights **274 of321 curves**: all270 with certified lower
bound at least22, plus four structural examples below22. These are the
orbit8044 seed factory's exported `000000` example, the uniform progression's
`n0` example, the orbit47755 alternate-fibration seed, and the small
equation-derived conic control. Their individual pages link to canonical proofs.
Other exported members of a represented family need not all occupy the front page.

The deterministic [display rule](../cas/render_main_readme_curves.py) also retains
any below22 curve with an **exact** conductor at or below the smallest conductor
reported at that rank or higher in the [pinned ICARM snapshot](../data/icarm_current.json).
This is a conservative editorial benchmark, not a claim about current world records.
No such additional exception is certified in this snapshot. Of47 hidden curves,
sixteen now have exact conductors above that benchmark and31 remain unresolved;
hiding an unresolved curve is **not** a high-conductor classification.
They remain in the complete inventory, JSON/CSV and individual
pages, with points and certificates untouched. Future verified conductor updates
can automatically bring qualifying curves back into the main table.

Both inventory tables bold the smallest displayed `log N`, naive height,
Faltings height and `log abs(Δ)` separately within each certified-rank group.
Ties at the displayed two-decimal precision are all bold; missing values do
not compete. The main README recomputes minima after its row filter, so its
below22 highlights can differ from the full inventory. These are within-table
comparisons, not exact-rank or public-record assertions.

## Bounded missing-conductor pass

The frozen [runner](../cas/run_inventory_conductor_pass.py) covers all161 conductors
missing at launch, rank-descending, with four detached workers, a90-second build
cap and120-second independent replay cap per curve, and2GiB per process. It reuses
retained factors and checkpoints exact local arithmetic before factoring the
remaining discriminant. Factorization timeouts preserve only proved bounds;
engineering or replay failures stop dispatch. No point search is performed.

The immutable [first snapshot](../../artifacts/generated-results/elliptic-curves/inventory291_conductor_snapshot_v1.json)
contains12 independently replayed certificates, seven exact and five partial.
Its SHA-256 is `7eddf86858dc59a8f676879e9863e1a556ebaac576e1838be7a812e92b2accac`.
The new rank-at-least27 panel fibre at877/781 has exact conductor

```text
277381620390410218245942534536512005669171519187199792818529590968773776639717845818054857352552455956407733220880597259130200
```

Thus its natural log conductor is288.84336. This is not a conductor-record claim.
Sage's generic Tate algorithm (using a degree-one number-field presentation ofQ)
and PARI `elllocalred` independently agree on each retained conductor exponent.
Exact PARI prime certificates, discriminant products and conductor products are
replayed without factor discovery. Minimal and nonminimal presentations of the
conductor11 regression both pass. The live pass may advance beyond this immutable
snapshot; later results require a new replayed snapshot before publication.

The [second snapshot](../../artifacts/generated-results/elliptic-curves/inventory291_conductor_snapshot_v2.json)
publishes all128 sealed V1 results:54 exact conductors and74 certified partial
bounds. Its SHA-256 is
`9621c25652553978de0332426d8ada196388d85a45b429a99afbc063987e4003`.
All128 certificates pass fresh portable Sage replay without factor discovery.
This adds47 exact conductors beyond the seven in the first snapshot, bringing
the full ledger to184 exact and107 unresolved. Snapshot version2 refers to
this publication revision, not to results of the V2 continuation worker.

The [third snapshot](../../artifacts/generated-results/elliptic-curves/inventory291_conductor_snapshot_v3.json)
merges the latest sealed result for each of all161 original cases, including
the completed33-case continuation and the first three exact results of the
longer pass (`new-20260906-123`, `new-20260906-142`, `new-20260906-151`).
It contains65 exact conductors and96 certified partial bounds, all passing
fresh portable replay with exact model-transport checks where applicable.
Its SHA-256 is
`2348fdb738c0a17dce0c76616df16e7ed789fcc45d3e03d4089c209b186510d9`.
This adds eleven exact conductors to the second snapshot's publication,
bringing the ledger to195 exact and96 unresolved. It is a fixed cutoff;
live checkpoints and later completions are not silently included.

### Preserved V1 stop and V2 continuation

The post-reboot audit found **128 sealed results:54 exact and74 partial**.
All frozen inputs and all retained result/binding hashes matched, including
the duplicate generated copies. There were no unsealed conductor checkpoints.
The125-case controller summary was stale. The actual stop was an engineering
guard, `integral model required`, on `det1092-bifibration-47755-u0`, not evidence
of lost arithmetic data. The other33 cases had no sealed result.

The [versioned V2 continuation](../cas/run_inventory_conductor_resume_v2.py)
preserves V1's sources, protocol, failed-case logs and128 sealed certificates.
It freezes only those33 unfinished inputs, in their original order. For each
rational equation, let `d` be the least common multiple of its coefficient
denominators and replace `a_i` by `d^i*a_i`. The exact map
`(x,y) -> (d^2*x,d^3*y)` gives an integral equation overQ with the same conductor;
this is not a global minimal-model computation. Original equations and rank
certificates in the ledger are not replaced.

The [V2 wrapper](../cas/inventory_conductor_worker_v2.py) verifies coefficient
identities, discriminant scaling by `d^12`, and rational isomorphism before
both build and replay. It then delegates local arithmetic to the unchanged V1
worker. Two rational presentations of the conductor11 regression build and
independently replay exactly; corrupted transport metadata is rejected.
The zero-search Sage preflight passed all33 actual continuation models before
launch. Four detached workers retain the same90-second build,120-second replay
and2GiB per-process caps. Completed V1 cases are not repeated. The published
ledger includes all128 retained V1 results and the completed33-case
continuation through the third immutable snapshot.

```sh
python3 research/elliptic-curves/cas/run_inventory_conductor_resume_v2.py status
```

## Reproduction

### Longer bounded factorization pass

V1 plus its rational-model V2 continuation completed all161 cases, with62
exact conductors and99 partial results, each independently replayed. The full
available inventory at that point therefore had192 exact conductors. The
third publication snapshot adds three more from the longer pass; live results
are not silently indexed.

At the user's request, the separate
[long V3 runner](../cas/run_inventory_conductor_long_v3.py) freezes all99
unresolved cases with **1,800 seconds per build**, up from90, and120 seconds
per independent replay. It retains four workers and2GiB per process. The
maximum allotted worker time corresponds to roughly13.3 hours at four-way
concurrency, plus launch/replay overhead; early completions can shorten it.
This is one bounded pass, not a promise to factor every remaining integer.

All saved certified primes and prior factor hints are reused. Exact pairwise
GCDs of the99 residual cofactors provide proper factor splits for two cases;
these splits are hints, not primality certificates. Scheduling is frozen by
increasing cofactor digit count, then decreasing certified rank and curve ID.
The local Tate/PARI algorithms, primality proof requirements and exact replay
remain unchanged. Ninety of the99 residual cofactors have78–149 digits; the
largest two have13,569 and22,693 digits. Extra runtime alone need not settle
those large cases, and nonminimal-model factors can contribute to their size.

Shared-split unit tests and the rational-model arithmetic regressions pass.
The foreground Sage preflight verifies all99 transported models and checks
that every saved prime and residual cofactor is retained. The pass then runs
detached, preserving all V1/V2 sources, inputs, certificates and logs.
Timeouts keep certified partial bounds and `UNKNOWN`; engineering or replay
failures stop dispatch for review. No point search or rank change is involved.

```sh
python3 research/elliptic-curves/cas/run_inventory_conductor_long_v3.py status
```

### Published snapshot replay

The [refresh manifest](../data/research_curve_refresh.json) selects proved
`MATH_STATUS.json` entries and pins saved result packets by SHA-256. The
[loader](../cas/research_curve_refresh.py) checks source bindings, curve
identities, rational nonisomorphism for additions, and exact point transports
for existing rows. It replays the saved finite reduction certificates and
2-torsion exclusions. Conductor values come from the pinned arithmetic audit,
including matching independent replay evidence for exact results.

```sh
python3 research/elliptic-curves/cas/render_main_readme_curves.py
python3 research/elliptic-curves/cas/render_main_readme_curves.py --check
~/.local/bin/sage -python research/elliptic-curves/cas/index_inventory_conductor_merged.py check
python3 research/elliptic-curves/cas/run_inventory_conductor_long_v3.py status
```

These commands synchronize the README, expanded inventory, JSON/CSV downloads,
and individual curve pages. They perform no point searches, specialization
sweeps, or conductor factorizations. Historical certificates remain preserved.
Computation checkpoints are automatic; README publication is not. A new
immutable snapshot must be replayed, registered in `MATH_STATUS.json`, indexed
in the refresh manifest, and rendered before newly computed values appear.
