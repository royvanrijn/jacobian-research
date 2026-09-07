# Current conductor comparison for the201-curve inventory

The630-entry ICARM snapshot contains seven exact Q-isomorphism matches to the
201 local research curves. The194 unmatched curves were screened; **24 now have
exact conductors**, including all13 selected priorities. **None of these24 is
a conductor record at its certified rank threshold.** The170 other unsubmitted
curves remain UNKNOWN. A large upper bound alone does not exclude a record.

All four submitted curves627–630 are now exact; see the
[complete prime lists](../../artifacts/generated-results/elliptic-curves/submitted628_conductor_v2/PRIMES_TO_PASTE.md).
Including626 and the24 new completions, the combined201-row local view contains
29 exact conductors,171 unknown conductors and one publicly reported value.
The extra unknown is an already-public curve without a recorded conductor.

## Closest solved unsubmitted candidates

Rank means the certified lower bound, and the benchmark cohort is public rank
at least that threshold. These are the best solved candidates, also the smallest
upper bounds in the initial unsubmitted roster. Unknown curves with larger bounds
could still have smaller actual conductors; no full-inventory exclusion is claimed.

| Rank at least | Local curve | Family / parameter | Exact N digits | N / recorded minimum | Benchmark |
|---|---|---|---:|---:|---|
| 22 | `new-20260906-181` | `11952`, `-8/27` | 78 | 3.73206e+06 | [ICARM#625](https://elliptic-rank.icarm.cloud/curve/625) |
| 23 | `new-20260906-162` | `103b2`, `48/7` | 87 | 9.76705e+08 | [ICARM#624](https://elliptic-rank.icarm.cloud/curve/624) |
| 24 | `new-20260906-129` | `07ca9`, `561/26` | 95 | 6.26381e+08 | [ICARM#623](https://elliptic-rank.icarm.cloud/curve/623) |
| 25 | `new-20260906-54` | `103b2`, `-538/249` | 100 | 1.96157e+09 | [ICARM#622](https://elliptic-rank.icarm.cloud/curve/622) |
| 26 | `new-20260906-92` | `a1-fibration-02`, `3161/432` | 105 | 1.28677e+10 | [ICARM#542](https://elliptic-rank.icarm.cloud/curve/542) |
| 27 | `new-20260906-41` | `11952`, `-2448/11` | 123 | 1.37802e+17 | [ICARM#614](https://elliptic-rank.icarm.cloud/curve/614) |

The public benchmark conductors are hash-pinned reported values, not independently
proved universal records. Missing public conductor values remain missing.

## Results and local database

- [24 exact conductor rows, CSV](../../artifacts/generated-results/elliptic-curves/new_curve_conductor_screen_v1/exact_conductors.csv)
- [Three remaining rank27 prime lists and submission links](../../artifacts/generated-results/elliptic-curves/new_curve_conductor_screen_v1/NEW_RANK27_PRIMES.md)
- [194 individual proof references and comparison summary](../../artifacts/generated-results/elliptic-curves/new_curve_conductor_screen_v1/summary.json)
- [Full201 local-bound audit](../../artifacts/generated-results/elliptic-curves/inventory201_conductor_bounds_v1.json)
- [Queryable conductor manifest](../data/conductor_screen_current.json)

```python
from local_conductor_database import load_conductor_inventory
rows = load_conductor_inventory()
```

Each row retains its original equation, points, rank certificate and publication
provenance. The added `conductor_information` exposes EXACT, UNKNOWN or REPORTED,
certificate paths, exact bad-prime lists where available, and proved conductor
divisors/upper bounds for unresolved new curves. These data support retrospective
search design; existing prospective protocols retain their frozen inputs.

## Computation and proof boundary

The finite local audit covers all201 curves, including the14 absent from the old
187-curve audit. It checks every discriminant prime through10000 plus2,3 using
Sage generic Tate and independent PARI conductor exponents. The unfactored residue
bounds the remaining conductor without assuming squarefreeness. No initial bound
beats the fresh recorded minimum.

The discovery protocol screens residual primality for all194 unpublished curves
and selects the two smallest bounds at each rank22–27 plus the third rank27 curve:
13 distinct priorities. Each gets a15-second PARI pass and at most275 deterministic
ECM attempts within120 seconds, with two workers and per-attempt checkpoints.
Eleven complete within this schedule. The two remaining rank25 candidates43 and54
complete with CADO-NFS in104.63 and444.90 seconds respectively (see raw measured
results for precise times), each within1800 seconds and confined to eight CPUs.
Raw inputs, factor products, attempts, protocols, NFS relations and checkpoints
are retained in `artifacts/local/elliptic-curves/conductor-record-screen-v1/`.

Every asserted prime is certified exactly, with its certificate root checked.
Factored and unresolved pieces reconstruct the displayed integral discriminant.
Two local implementations agree on all processed conductor exponents. Their product
divides the true conductor; multiplying it by the unprocessed discriminant cofactor
bounds the true conductor above. Cofactor1 proves an exact conductor and complete
bad-prime list. Local minimization handles nonminimal input equations. No new rank
claim or point search is involved.

## Replay

```sh
sage -python elliptic-curves/cas/audit_inventory201_conductor_bounds_v1.sage check
sage -python elliptic-curves/cas/certify_new_curve_conductor_screen.sage --check
python3 elliptic-curves/cas/local_conductor_database.py
python3 elliptic-curves/cas/index_submitted_conductors_v2.py --check
python3 elliptic-curves/cas/refresh_icarm_local_database.py --check
python3 elliptic-curves/cas/replay_icarm626_controls.py
```

The arithmetic replays, current database loader, mathematical-status audit and161
local links in the edited documents pass. The repository-wide Markdown check also
reports an existing unrelated broken reference from
`archive/elliptic-curves/EXTENDED_CACHE_RETENTION_GATE_2026-09-06_PRE_COMPLETION.md`
to `OUTER131072_TRIAL_2026-09-06.md`; that historical document was not changed.
