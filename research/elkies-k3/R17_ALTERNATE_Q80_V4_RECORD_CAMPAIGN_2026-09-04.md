# Native alternate-Q80 record campaign from rank-one V4 bases

<!-- status-consumer: EC-K3-R17-NORM12-11952-V4-RANK-ONE-RECORD-CAMPAIGN db2be7ab9f65ce57 -->

## Exact bounded result

The seventeen exact rank-one bases in the alternate-Q80 V4 shortlist now
have primitive free generators and explicit maps to the alternate parameter
`u`.  The preparation distinguishes the product quotient

```text
w^2 = q_i(u) q_j(u)
```

from the actual V4 base.  It parametrizes the first conic through the stored
rational intersection point, constructs the pointed quartic for the V4 base,
identifies its Jacobian by the exact degree-two isogeny from the product
quotient, and saturates the image at 2.  On all seventeen bases:

```text
product-quotient visible-point saturation index     1
isogeny-image saturation index on the V4 base       2
V4-base rational torsion order                       2
primitive free-generator canonical height      90--154 approximately
```

The first native record pass ranks the multiples with

```text
2 <= |n| <= 8
```

using complete periods of `n -> u(nP) mod p` for the three disjoint blocks

```text
131,137,151 ; 157,167,173 ; 181,191,193.
```

At every good base period, the alternate-surface trace contribution and the
number of native bisection squareclasses that split modulo `p` are centered
and population-standardized separately.  Integers are ordered by the weakest
combined block, followed by the weakest splitting and trace blocks.  This is
a heuristic ordering only.  Exactly the best integer on each of the seventeen
bases is then materialized over `QQ`.

For every selected rational `u`, the complete set of 39,147 native alternate
bisection values is passed through all nine modular square filters before an
exact rational-square test.  The pass reduces

```text
17 * 39,147 = 665,499 possible exact tests
```

to `1,772` exact tests, or about `0.266%`.  Every exact square is one of the
two defining characters.  There are no branch collisions and no additional
native split points in this bounded pass.

All seventeen generic alternate-Q80 sections and the two defining character
sections are nevertheless specialized at every selected fibre.  Exact
finite-quotient certificates at relation prime 2 use between 13 and 17 good
reduction primes and prove all nineteen displayed points independent.  Thus
the result is:

```text
selected exact fibres                         17
certified rank lower bound 19                 17
additional native bisection splits             0
rank at least 20 promotions                     0
rank at least 24 promotions                     0
rank at least 28 promotions                     0
rank at least 32 promotions                     0
```

The generator/map artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-rank-one-bases-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-rank-one-bases-v1.json),
SHA-256
`ebb52552ecd76cc0d2334dbfd36364095a23816b4014beec53aa72aa44977733`.
The campaign artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-record-campaign-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-record-campaign-v1.json),
SHA-256
`06a5569ee51fa27ce1818779ec9554457e0200ff1ca7b13ec416af06a2c103d0`.

## Promotion and descent gates

The campaign records the required thresholds without weakening them:

| certified rank lower bound | action |
|---:|---|
| 20 | native mechanism control |
| 24 | useful score calibration |
| 28 | serious alternate-Q80 positive control |
| 32 | record candidate |

For the current generic rank-19 subgroup, a completed residual 2-Selmer
dimension below 13 rejects rank 32.  If a separately certified generic
twentieth section is added, the rejection threshold becomes 12.  No fibre in
this pass reached rank 32, so residual descent was not triggered.  In
particular, no incomplete descent or heuristic score authorized an
unrestricted point search.

## Replay

```bash
sage -python \
  elkies-k3/scripts/prepare_r17_norm12_11952_alternate_v4_rank_one_bases.sage \
  --jobs 2 --timeout 180

sage -python \
  elkies-k3/scripts/prepare_r17_norm12_11952_alternate_v4_rank_one_bases.sage \
  --check

sage -python \
  elkies-k3/scripts/run_r17_norm12_11952_alternate_v4_record_campaign.sage \
  --base-limit 17 --n-bound 8 --minimum-abs-n 2 --exact-per-base 1

sage -python \
  elkies-k3/scripts/run_r17_norm12_11952_alternate_v4_record_campaign.sage \
  --check
```

## Claim boundary

The base generators, maps, selected rational points, 39,147-value modular
and exact split tests, specialized sections, and displayed rank lower bounds
are exact.  The negative split result is exhaustive only for the seventeen
selected multiples.  It is not a rank upper bound for any fibre or a negative
statement about other multiples.  The periodic score is not a rank or Selmer
calculation, and the missing native rank-20 calibration remains open.
