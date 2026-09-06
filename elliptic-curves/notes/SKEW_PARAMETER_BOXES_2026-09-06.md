# Equal-area skew boxes improve a coefficient bound but add no high-rank curve

A new rectangular parameter population completes 40,733,526 primitive
addresses, 2,048 retained prime-score extensions, and 368 point-search boxes
on eight fixed candidates. The exact lower bounds are six at 17 and two at
18, unchanged modulo 3 and 5. No curve enters the 101-curve high-rank inventory.
The [aggregate certificate](../../artifacts/generated-results/elliptic-curves/skew_r17_experiment_v1.json)
binds the geometry, selection, trace caches and point results.

## The omitted box-shape choice

For homogeneous coefficients
`A_h(n,d) = sum A_i n^i d^(8-i)` and
`B_h(n,d) = sum B_i n^i d^(12-i)`, the rectangle
`|n| <= H 2^k`, `1 <= d <= H 2^-k` gives exact triangle bounds

```
|A_h| <= H^8 U_A(k),    U_A(k) = sum |A_i| 2^(k(2i-8)),
|B_h| <= H^12 U_B(k),   U_B(k) = sum |B_i| 2^(k(2i-12)).
```

Consequently `max(4|A_h|^3,27|B_h|^2)` is at most
`H^24 max(4 U_A(k)^3,27 U_B(k)^2)`. All terms and the minimum over
`k = -4,...,4` are exact rationals. The [bound certificate](../../artifacts/generated-results/elliptic-curves/r17_parameter_box_skew_v1.json)
checks all six families. Four choose the original square. Only `08234` and
`08f72` improve, by factors greater than 2^54 and 2^9 respectively in this
weighted bound. This is a worst-case coefficient bound, not minimal
arithmetic height or a rank predictor.

At `H = 32768`, the new rectangles are:

| Family | Numerator bound | Denominator bound | Chosen k |
| --- | ---: | ---: | ---: |
| `08234` | 262144 | 4096 | 3 |
| `08f72` | 65536 | 16384 | 1 |

Both preserve `N*D = 32768^2`. Primitive counts in denominator slices need
separate arithmetic and are not inferred from equal area. The known
`08234`/published-R17 base equivalence remains in force; this is a changed
parameter region, not a new generic family.

## Frozen execution

The two signs and the same previously fixed denominator residues modulo 64
are retained. The four unchanged-family square boxes are not rerun. Before
any large call, all 2,634 primitive scores and complete orders in four small
rectangles with the same aspect ratios are checked against canonical traces.
The unchanged periodic scanner then retains 512 rows per signed slice.
All 2,048 returned 562-prime S1 scores replay.

Of the retained addresses, 1,381 have `|n| > 32768`; 654 overlap the original
S1 pool. Trace caching against both earlier pools reuses 656 rows and computes
1,392 fresh rows. Sixteen independent character sums and the existing cost
gate pass before bulk extension. There are 8,321,376 fresh and 3,921,568 reused
extension traces. Selection uses unchanged S1 through 32,749 and excludes
validation primes through 65,521 even from ties.

Exactly four candidates per family are selected only where `|n| > 32768`.
Their generic-17-only policy uses all 43/49 generic maximum parity classes,
height 125,000 and ten seconds per chart. Every one of the 368 boxes completes.
All histories, rational maps and full-cloud point proofs replay. The eight
curves are distinct and unmatched in the pinned 593-equation catalogue and
520 prior measured address-equations, but their lower bounds do not meet the
inventory threshold of 22.

The scan/replay takes 14.303/2.010 seconds; trace extension/replay takes
134.443/10.482 seconds. Point attempts and history/cloud verification take
193.539/19.776 seconds.

## The improved bound did not substantially shrink selected equations

The preceding eight square-box candidates in these same two families have
normalized coefficient sizes 268–338 bits, median 294.5. The new eight have
272–308 bits, median 294. The aggregate records which normalizations have a
proved global-minimality witness and leaves other minimality claims unresolved.

Thus the large worst-case bound improvement does not establish a comparable
reduction for the score-selected candidates. It neither explains their low
point bounds nor demonstrates that skew selection is superior. No larger
rectangle or automatic refill follows this bounded result.

The [recent generic twist-contraction criterion](../rank-jump/GENERIC_SUBGROUP_FORCES_TWIST_CONTRACTION.md)
also does not supply a replacement high-rank selector: its relative
contractions are already forced by generic sections on all three tested
controls, including the observed-zero control. Its stated missing implication
is rational solubility inside the common Selmer quotient. Those retrospective
results were not inserted into the new score policy.

The separate [zero/infinity audit](COMPACT_ENDPOINT_AUDIT_2026-09-06.md)
checks another explicit omission from the compact parameter boxes.
