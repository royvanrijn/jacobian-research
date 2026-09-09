# Curve inventory refresh — September 9

The [main README](../../../README.md#elliptic-curve-inventory) and
[expanded inventory](../INVENTORY.md) now include **225 distinct curves**,
with **130 exact conductors** and **95 unresolved conductors**. Rank entries
are certified lower bounds, not exact ranks.

The selected updates are:

| Source | Inventory change |
| --- | --- |
| [Curve52](CURVE52_SHORT_V3_SEARCH_2026-09-09.md), [curve113](CURVE113_PRODUCTIVE_V3_SEARCH_2026-09-09.md), [curve116](CURVE116_PRODUCTIVE_V3_SEARCH_2026-09-09.md), [curve200](CURVE200_PRODUCTIVE_V3_SEARCH_2026-09-09.md) | Four existing lower bounds rise from 25 to 26; the new points are transported to the existing certified minimal models. |
| [First fresh R17 cohort](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md) | Six new curves with lower bounds 24, 23, 17, 24, 25, 17 in family order 074d9, 07ca9, 08234, 08f72, 103b2, 11952. The 08f72 conductor is exact; three others have certified partial bounds. |
| [Second fresh R17 cohort](SECOND_FRESH6_SEED_COHORT_2026-09-09.md) | Six new curves with lower bounds 18, 18, 18, 17, 26, 22 in the same family order. |
| [Bifibration conic](DET1092_BIFIBRATION_CONIC_SEEDS_2026-09-09.md) | One exported determinant1092 specialization of orbit47755 at conic parameter u=0, with lower bound 18. |

The [earlier supplement](INVENTORY_SEED_SUPPLEMENT_2026-09-08.md) remains
included. All 24 additions to the original 201-row inventory retain source
models and uncomputed minimal-model height/discriminant fields. An exact
conductor does not supply those missing metrics. Bounded misses remain
unresolved, and the inventory makes no current public novelty or record claim.

## Reproduction

The [refresh manifest](../data/research_curve_refresh.json) selects proved
`MATH_STATUS.json` entries and pins saved result packets by SHA-256. The
[loader](../cas/research_curve_refresh.py) checks source bindings, curve
identities, rational nonisomorphism for additions, and exact point transports
for existing rows. It replays the saved finite reduction certificates and
2-torsion exclusions. Conductor values come from the pinned arithmetic audit,
including its matching independent replay digest for the exact result.

```sh
python3 research/elliptic-curves/cas/render_main_readme_curves.py
python3 research/elliptic-curves/cas/render_main_readme_curves.py --check
```

These commands synchronize the README, expanded inventory, JSON/CSV downloads,
and individual curve pages. They perform no point searches, specialization
sweeps, or conductor factorizations. Historical certificates remain preserved.
