# Explicit seed supplement to the curve inventory

The [inventory](../INVENTORY.md) adds **11 distinct curves** to the original
201-row discovery inventory. The supplement lists the exported determinant1092
seed examples certified in the current mathematical-status ledger:

| Source | Additional distinct curves | Certified rank lower bounds |
|---|---:|---|
| [Small conic seed](DET1092_SMALL_CONIC_SEED_2026-09-08.md) | 1 | 18 |
| [Funnel first seeds](DET1092_FUNNEL_FIRST_SEEDS_2026-09-08.md) | 3 | 18, 21, 18 |
| [Orbit8044 factory](ORBIT8044_SEED_FACTORY_2026-09-08.md) | 6 | 18 each |
| [Progression, exported n=0 example](DET1092_CONIC_SEED_PROGRESSION_2026-09-08.md) | 1 | 18 |

The factory's `u=0` packet is the small conic seed and appears once, with its
alias retained. Exact rational j-invariants select comparison buckets; exact
rational isomorphism decides duplicates, so twists are not merged merely for
having the same j. The progression uses the original conic slope coordinate,
which differs from the factory coordinate.

These entries retain their certified source equations and rational points.
Global minimality, minimal discriminants and ICARM-style heights have not been
certified for this supplement, so those fields are null and display as dashes.
Their conductors remain UNKNOWN. The original 201 rows and their 129 exact
conductors retain their existing arithmetic evidence. Rank21 is the certified
completed epoch at reduced parameter `s=1926/2699`; later search progress does
not change the table without a selected certificate.

This is a finite discovery inventory, not a census of every public calibration
curve, experimental control, or member of an infinite family. It asserts no
literature-wide novelty or exact rank. The canonical proof notes above and
`MATH_STATUS.json` remain the mathematical authorities.

## Regeneration

The [supplement manifest](../data/research_curve_supplement.json) pins each
equation packet and independent proof by SHA-256 and names its proved ledger
entry. The renderer checks those bindings, point identities, proof ranks and
deduplication. It does not rerun point searches or conductor factorization.
The original [201-row metric and conductor replays](INVENTORY201_TABLE_AND_CONDUCTORS_2026-09-07.md)
remain available separately.

```sh
python3 research/elliptic-curves/cas/render_main_readme_curves.py
python3 research/elliptic-curves/cas/render_main_readme_curves.py --check
```

Both commands cover the main README, expanded inventory, JSON/CSV downloads
and individual curve pages. To add a later certified packet, extend the
manifest with its proof binding and retain the previous certificates.
