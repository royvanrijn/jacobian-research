# Extreme-anchored MW>=18 covers (2026-09-04)

<!-- status-consumer: EC-K3-R17-EXTREME-ANCHORED-MW18-COVERS 4763babcbf5d923c -->

## Status

The MW18 search lane is restarted around exact high-jump anchors.  The old
first-cover Nagao run remains a regression, not the default search.

The complete smooth rootless rational-bisection frames on the native
`norm12-orbit-07ca9` and `norm12-orbit-08234` charts were tested at eight
refreshed `+10/+11/+12` controls: 312,960 bisection--fibre pairs in total.
Every modular nonresidue is retained as an exact finite-field obstruction and
every survivor is reconstructed over `QQ`.

The exact result is eight refreshed extreme-anchored covers:

| curve | displayed jump over MW17 | split covers | nonzero anchored covers | anchored quotient span |
|---:|---:|---:|---:|---:|
| 543 | 12 | 0 | 0 | 0 |
| 544 | 11 | 0 | 0 | 0 |
| 545 | 11 | 2 | 2 | 2 |
| 531 | 11 | 3 | 3 | 3 |
| 534 | 11 | 2 | 2 | 2 |
| 535 | 11 | 0 | 0 | 0 |
| 536 | 11 | 1 | 1 | 1 |
| 537 | 10 | 0 | 0 | 0 |

Thus the attractive `+12` premise does **not** hold in this complete
bisection layer: curve 543 has no split.  The positive refreshed anchors are
all `+11`.  The historical published rank-at-least-28 control at
`t=-9529/5471` supplies one further exact `+11` anchor, `orbit-15a68`.

## Exact cover roster

| anchor | cover | published priority | maximum branch-quadratic coefficient bits |
|---|---|---:|---:|
| curve 545 | `07ca9-orbit-08c1e` | 27,252 | 1,078 |
| curve 545 | `07ca9-orbit-1d516` | 8,947 | 1,120 |
| curve 531 | `08234-orbit-0a9bf` | 29,425 | 1,087 |
| curve 531 | `08234-orbit-12f61` | 24,492 | 1,095 |
| curve 531 | `08234-orbit-1293d` | 35,877 | 1,116 |
| curve 534 | `08234-orbit-13d7a` | 24,637 | 1,084 |
| curve 534 | `08234-orbit-1a371` | 6,662 | 1,108 |
| curve 536 | `08234-orbit-19188` | 38,128 | 1,051 |
| historical rank-28 | `orbit-15a68` | 35,421 | small published chart |

For each refreshed cover the certificate verifies:

- the transported trace word is integral, distinct, and has norm ten;
- `u^2=q(t)` is a smooth rational conic and the lifted section identities
  hold over `QQ[t,u]/(u^2-q(t))`;
- the two branches at the declared extreme fibre lie on the fibre and add to
  the trace;
- the positive branch has a nonzero class modulo the specialized generic
  MW17 subgroup;
- the conic parameter is normalized by `t(0)=t0` and `dt/dr(0)=1`.

Numerical canonical heights only propose public-group relations.  Every
retained relation is checked by exact elliptic-curve addition.  The exact
certificate is
[`../artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json`](../artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json),
SHA-256
`7b94bf9b3dec377ce118c64e4455c6e2f089e88108cdabc9b3c5ca16c3a94b27`.

## Rank accounting

Each listed cover has generic Mordell--Weil rank **at least** 18.  Exact
generic rank 18 is not proved.  On a `+11` anchor, one exceptional direction
becomes generic and ten displayed exceptional directions remain at the
anchored specialization.  Reaching rank 32 from a generic-rank-at-least-18
cover therefore needs at most fourteen further directions, four beyond the
ten already displayed at the anchor.

The proposed `+12` accounting—eleven remaining displayed directions and only
three further ones—remains correct conditionally, but no such split exists in
the two complete frames tested here.

## Uniform Nagao triage

All nine exact covers were searched with one protocol:

- primitive rational `r` with numerator and denominator at most 1,000;
- the same three prime blocks through 197;
- per-height-bucket retention `32,16,8` with bucket width 100;
- at most 100 final records per cover.

This scores 1,216,768 parameters per cover and 10,950,912 in total.  It leaves
178 heuristic survivors.  Ordering covers by their best final Nagao score
gives:

| order | anchor / cover | survivors | best `r` | best score |
|---:|---|---:|---:|---:|
| 1 | curve 536 / `08234-orbit-19188` | 25 | `405/124` | 16.688242839829 |
| 2 | curve 531 / `08234-orbit-12f61` | 20 | `-901/737` | 16.512548206412 |
| 3 | curve 531 / `08234-orbit-0a9bf` | 26 | `808/403` | 16.351297229128 |
| 4 | curve 534 / `08234-orbit-1a371` | 21 | `24/901` | 16.055093926971 |
| 5 | curve 534 / `08234-orbit-13d7a` | 19 | `-701/399` | 15.694291160633 |
| 6 | historical rank-28 / `orbit-15a68` | 20 | `-401/294` | 14.823366947598 |
| 7 | curve 545 / `07ca9-orbit-08c1e` | 19 | `-82/101` | 14.342330260311 |
| 8 | curve 531 / `08234-orbit-1293d` | 18 | `704/71` | 14.291996919480 |
| 9 | curve 545 / `07ca9-orbit-1d516` | 10 | `-501/89` | 11.168867285069 |

The compiled campaign is
[`../artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-nagao-h1000-summary-v1.json`](../artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-nagao-h1000-summary-v1.json),
SHA-256
`314fa5f91fd3ee86b5cbc0a25d60e81d3e1ddd3bf0eb3be0f29622b217bb9d45`.
The individual ledgers preserve every finalist and exact base parameter.

This table is a search ordering only and proves no superiority of one cover
beyond this bounded protocol.  The subsequent exact audit specialized all
178 finalists and certified eighteen independent points on every fibre.  The
residual-Selmer and point-search continuation is paused; see
[`R17_EXTREME_ANCHORED_MW18_CONTINUATION_HANDOFF_2026-09-04.md`](R17_EXTREME_ANCHORED_MW18_CONTINUATION_HANDOFF_2026-09-04.md).

## Replay

```bash
sage -python elkies-k3/scripts/certify_r17_extreme_anchored_mw18_covers.sage

python3 elkies-k3/scripts/search_r17_extreme_anchored_mw18_nagao.py \
  --curve-id 536 --cover-label 08234-orbit-19188 \
  --numerator-bound 1000 --denominator-bound 1000 \
  --height-bucket-width 100 --finalists 100 \
  --output artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-nagao-curve536-orbit19188-h1000-v1.json

python3 elkies-k3/scripts/summarize_r17_extreme_anchored_mw18_nagao.py --check
sage -python elkies-k3/scripts/specialize_r17_extreme_anchored_mw18_finalists.sage --check --no-resume
python3 -m unittest elliptic-curves/tests/test_r17_extreme_anchored_mw18.py
```

The full exact replay is checkpointed under `artifacts/local/` for speed, but
the local cache is not a proof artifact.  `--check` on the exact compiler
rebuilds without trusting that cache.
