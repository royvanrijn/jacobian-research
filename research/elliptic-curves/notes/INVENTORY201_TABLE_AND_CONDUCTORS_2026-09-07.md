# Main README inventory and continuing conductor calculations

The [elliptic-curve inventory](../INVENTORY.md) contains all201
research curves in an expandable ICARM-style table. The [full inventory page](../INVENTORY.md) shows the same generated table expanded. Each row links to a complete minimal
equation, exact discriminant, transported rational points, rank provenance and
the available conductor proof. [JSON](../data/research_curves/database.json)
and [CSV](../data/research_curves/database.csv) exports retain full numbers.

## Models, points and numerical columns

The [metric certificate](../../artifacts/generated-results/elliptic-curves/inventory201_table_metrics_v1.json)
checks all201 global minimal models and4,796 exact point transports. Local
minimality is independently checked at2,3 and all prime divisors of the invariant
GCD, which is at most10^7. Outside that set, an invariant is a unit and a smaller
integral model is impossible. Original rank certificates remain attached to the
original models; exact isomorphisms transfer independence to the displayed models.
The rank28 entry is a public-point reproduction, and seven rows match ICARM.

Columns follow the [ICARM table](https://elliptic-rank.icarm.cloud/curves):
a-invariants, certified rank lower bound, log conductor, naive height, Faltings
height and log absolute minimal discriminant. All logarithms are natural.
Naive height is `log max(abs(c4)^3,c6^2)` on the minimal model. The Faltings column
uses ICARM's actual period-area formula, `-1/2 log(area)`, on that model. This
matches Sage's `stable=False` convention; it must not silently be replaced by
Sage's different stable-height option for additive curves. The numerical columns
are approximations, checked at96 and160 bits and displayed to two decimal places.
All seven existing public matches reproduce ICARM's minimal discriminant and
height values within10^-9.

Coefficients are clipped after14 characters in the table, as on ICARM; linked
curve pages and JSON retain every digit. Unknown exact conductors appear as a
dash. Certified divisors, upper bounds and partial prime lists are kept on the
curve pages and in separately named JSON fields. Partial lists are explicitly
unsuitable as complete submission lists.

## Continued factorization — complete bounded pass

The [201-row conductor certificates](../../artifacts/generated-results/elliptic-curves/inventory201_conductors_v2/summary.json)
now prove **129 exact conductors**, up from29: **100 additional completions**.
The other72 remain UNKNOWN. Both previously public rows600 and619 now have
local exact conductor proofs. All129 exact values exceed their recorded
rank-threshold benchmark; there is no new conductor record from this pass.
The current database, main README and expanded inventory select this supplement;
previous snapshots and the original four submission certificates are preserved.

The continuation pass is frozen in
`artifacts/local/elliptic-curves/conductor-inventory-continuation-v2/factor_protocol.json`.
It covers all172 rows lacking a local exact-conductor proof after the earlier
29 completions, including the two already-public rows whose conductor had not
been proved locally. Each receives at most15 seconds of PARI partial factoring
and30 seconds of deterministic ECM, with four workers and per-attempt checkpoints.
Previous factors are imported as candidates and checked by product; exact prime
certificates and two local conductor implementations are required before a value
appears in the README's log N column. No point search or rank change is involved.

The earlier [24-curve screen](NEW_CURVE_CONDUCTOR_RECORD_SCREEN_2026-09-07.md)
and all its certificates remain reproducible as a historical snapshot.

## Regeneration

```sh
cd research
sage -python elliptic-curves/cas/build_inventory201_table_metrics.sage --check
sage -python elliptic-curves/cas/certify_inventory201_conductors_v2.sage --check
python3 elliptic-curves/cas/local_conductor_database.py
python3 elliptic-curves/cas/render_main_readme_curves.py
python3 elliptic-curves/cas/render_main_readme_curves.py --check
```

The main README wrapper copies the canonical inventory into an expandable marked section. The inventory renderer produces the full table and curve data.
It preserves the research introduction, navigation and status-consumer markers.
The downloadable data and201 individual pages are generated from the same
certified sources as the table.

During the pass, the repository was relocated under `research/`. The interrupted
worker log and all completed checkpoints were retained, and the pass resumed
without restarting completed curves. Certificate-relative paths still use the
research root; the renderer locates the Git root and keeps the main table there.
The duplicate generated table in `research/README.md` is replaced by a link.

To select a later certified batch and regenerate both tables:

```sh
cd research
python3 elliptic-curves/cas/local_conductor_database.py --index --summary artifacts/generated-results/elliptic-curves/inventory201_conductors_v2/summary.json
python3 elliptic-curves/cas/render_main_readme_curves.py
```

Further factorization needs a new bounded protocol and versioned certificates;
the completed pass preserves all partial factors and exhausted attempts. A
compatibility symlink to the existing root `.python-version` preserves the
research status ledger’s runtime lock after relocation.
