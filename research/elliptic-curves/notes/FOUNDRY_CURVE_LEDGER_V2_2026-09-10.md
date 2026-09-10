# Second foundry curve-ledger and conductor cutoff

The second fixed foundry snapshot adds **94 distinct curves** of certified
lower bound at least22 to the prior321-curve ledger. The complete inventory is
now **415 curves:256 exact conductors and159 unresolved**. The main README
retains its concise rank-at-least22 presentation plus the existing structural
exceptions; the full inventory retains every row.

| Certified lower bound | New curves |
| ---: | ---: |
| 22 | 39 |
| 23 | 25 |
| 24 | 23 |
| 25 | 6 |
| 26 | 1 |

The [curve snapshot](../../artifacts/generated-results/elliptic-curves/foundry_curve_ledger_snapshot_v2.json)
pins180 sealed eligible endpoints from `high-rank-foundry-v3` and selects the
strongest packet in each exact rational-isomorphism class relative to the
321-curve baseline. Fresh replay verifies every selected point subgroup twice,
including rational2-torsion exclusion and pairwise/baseline nonisomorphism.

The simultaneous [conductor snapshot](../../artifacts/generated-results/elliptic-curves/foundry_v3_conductor_snapshot_v1.json)
contains the eight selected curves whose autonomous conductor jobs had already
finished and independently replayed at the cutoff. Six conductors are exact;
two remain `UNKNOWN` with certified divisor and upper bounds. Later source and
conductor completions remain in the live queues for a future snapshot.

Ranks are lower bounds, distinctness is relative to the pinned ledger, and no
exact-rank, worldwide-novelty or conductor-record claim is made.

```sh
~/.local/bin/sage -python research/elliptic-curves/cas/refresh_foundry_curve_ledger_v2.py check
python3 research/elliptic-curves/cas/render_main_readme_curves.py --check
```
