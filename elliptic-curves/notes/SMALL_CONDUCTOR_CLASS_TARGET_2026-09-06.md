# MW16 at 3/17: reducing the class-group quotient toward 16

Mathematical status: `MATH_STATUS.json`, entry `EC-SMALL-CONDUCTOR-CLASS-TARGET-20260906`.

The current certified **GRH-conditional class-2-rank upper bound is 1142**.
The unconditional curve-rank lower bound remains **22**. Exact rank remains unknown.

The [preceding descent study](SMALL_CONDUCTOR_DESCENT_SHORTCUT_2026-09-06.md) establishes
the cubic field, exact local correction `rank <= dim Sel_2 <= g+7`, proved even
Selmer parity, and the interval-certified generating cutoff 37,638 under GRH.
Here `g=dim Cl(K)[2]`. A certified `g<=16` suffices for rank 22 under the same
assumption. This criterion need not be attainable if the actual class-2-rank exceeds 16.

## Authorized continuation

The user explicitly set the goal of reducing this quotient to 16. Each wave fixes
a finite target list, candidate region, smoothness bound and single-worker resource
limits before execution. Prime ideals are selected from free columns after exact
supported row reduction, excluding unsuitable ramified or index-dividing primes.
The Hessian-reduced index-prime lattices and previously searched regions are replayed.
Target checkpoints permit bounded resumes. A wave stops early only at dimension 16.

Each saved witness is reconstructed and its principal ideal is factored exactly.
The checker verifies the nonmonic norm identity including the fixed square factor.
Relations may use primes up to the declared smoothness cutoff, never beyond the
inherited complete 400,000 factor base. All coordinates above 37,638 are retained
until exact elimination cancels them. The supported rank is independently checked
as `rank(all relation rows) - rank(their outside projection)`.

Only successful witnesses are required for the rank upper bound. Candidate counts
and population digests describe the worker run; rejected values are not replayed
and no exhaustive-search or smoothness-completeness assertion is made. Earlier
full scalar sieve replays remain preserved as implementation calibration.

Near-root strips use the three real roots of the norm form, transformed into each
target lattice. Exact algebraic arithmetic floors each slope times `2^96`.
At each permitted positive denominator v, the search tries the integer center
`floor(v*slope_scaled/2^96)` and its two neighbors. All accepted integer
coordinates, norm values and principal ideals are checked exactly; no general
smoothness or completeness claim follows from this choice of region.

## Audited waves

| Wave | Targets completed | Region | Smooth bound | Candidate occurrences | Relation occurrences | Independent supported gain | Remaining dimension | Worker seconds |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| [box 1](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_001_v1.json) | 512 | box 128 | 37638 | 9891662 | 824 | 588 | 2291 | 130.5 |
| [box 2](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_002_v1.json) | 512 | box 128 | 37638 | 9901917 | 807 | 434 | 1857 | 128.32 |
| [box 3](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_003_v1.json) | 512 | box 128 | 50000 | 9918007 | 1248 | 346 | 1511 | 143.13 |
| [box 4](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_004_v1.json) | 512 | box 128 | 50000 | 9847135 | 1181 | 345 | 1166 | 140.91 |
| [strip 1](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_strip_wave_001_v1.json) | 64 | v <= 1024 | 50000 | 358109 | 55 | 24 | 1142 | 4.09 |

Starting dimension was **2,879**; the audited supported gain is **1737**.
**1126 further independent supported rows** suffice for the target.
The current coarse curve-rank upper bound is **1148 under GRH**.
A remaining dimension is an upper bound, not a proven number of actual class-group
directions. It gives neither a new rational point nor an algebraic-rank parity claim.

## Replay

The checker first replays the inherited proof chain, then every new principal
relation, target selection, matrix transition and independent rank identity.

```bash
sage -python elliptic-curves/cas/pursue_small_conductor_class_target_strips.sage check --wave 1
```

Protocols, successful witnesses, logs and supervisor outcomes are retained under
`artifacts/local/elliptic-curves/small-conductor-class-target-v1/` and
`artifacts/local/elliptic-curves/small-conductor-class-target-strips-v1/`.
Every linked generated certificate pins its protocol, source and witness chunks.
No new portable-archive replay is claimed until its separate report is available.
