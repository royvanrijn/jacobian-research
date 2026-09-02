# Fermigier mixed-small-prime CRT--Gauss search

## Closed search

The second Fermigier CRT--Gauss search uses the exact local targets

```text
5^3, 7^5, 11^4, 13^3, 17^4, 19^4, 23^2, 29^3, 31^2
```

in the canonical adapter coordinate `u=s/2`.  It exhausts both projective
charts for every four-prime subset, compresses complete singular Hensel trees
to maximal p-adic balls, and enumerates a radius-4 box in every exact
two-dimensional Gauss-reduced kernel.  The ICARM-282 effective profile
`(11,13,23,31)` additionally receives radius 64; the known parameter
`u=11671/42` occurs at reduced coordinates `(-52,3)`.

At projective height at most `100000`, the closed population contains:

- 126 mixed four-prime profiles;
- 36,352 projective CRT classes;
- 10,062,080 bounded lattice vectors;
- 963,684 unique positive parameters;
- SHA-256 `b96a68752dfbd7cdaf4030fc4c4656adedb9300a3a6772a0baa18eb514da1e75`
  for the canonical parameter stream.

The known curve 282 is recovered and excluded from novel selection.  The
previous frozen parameter census removes a further 48 parameters.

## Exact gates and result

Three exact lanes select 381 fibers: smallest known radical upper proxy,
lowest height, and largest known powerful part.  Before any accidental-point
score, an exact mod-5 finite-reduction gate measures the twelve specialized
generic section differences:

| certified subset rank | fibers |
| ---: | ---: |
| 9 | 1 |
| 10 | 6 |
| 11 | 37 |
| 12 | 337 |

Only the 337 rank-12 survivors continue.  Trial division through 2000 selects
32 radical-best fibers.  Each receives the bounded exact quartic slice

```text
x = +/-2u + n,  |n| <= 5000.
```

Thirty of the 32 fibers have at least one additional quartic point, but exact
finite-reduction certificates give only rank at least 13 at best: 12 fibers
reach 13 and 20 remain at 12.  The search therefore does not meet either the
rank-21/conductor target or the rank-at-least-30 target.

This negative result exposes an important selection effect.  The strongest
small-prime radical profiles are rich in accidental quartic points, but those
points mostly span one extra Mordell--Weil direction; without the retention
gate, the smoothest fibers also frequently sacrifice generic directions.

No exact conductor was computed because PARI/GP was unavailable in this
workspace.  That does not affect the population or rank certificates, but no
candidate can pass the conductor promotion gate.  The integer-offset point
slice is bounded and is not a complete quartic point search.

## Denominator-aware follow-up

All twelve rank-13 fibers were subsequently searched on the reduced rational
offset slice

```text
x = +/-2u + n/d,  2 <= d <= 64,  |n| <= 5000.
```

The implementation clears the normalized quartic by the exact square
`B^6*d^4`, applies only lossless modular square filters, and finishes with
integer `isqrt`.  Every hit is then replayed on the raw quartic, mapped to the
canonical curve, and included in a new finite-reduction certificate.

The scan found 135 genuinely new quartic abscissas.  Ten fibers remained at
rank at least 13, while two gained one certified direction:

| adapter `u` | new denominator points | certified rank lower bound |
| ---: | ---: | ---: |
| `1785/4` | 14 | 14 |
| `539` | 6 | 14 |

A deeper pass on exactly those two fibers used `65 <= d <= 256` and
`|n| <= 20000`.  It found four further points at `u=1785/4` and none at
`u=539`; neither increased the certified rank beyond 14.  This supports moving
to a general rational-point engine or a different slice geometry rather than
merely increasing the same denominator box.

## Reproduction

```bash
PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/search_fermigier_mixed_small_prime_crt_gauss.py

PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/search_fermigier_mixed_small_prime_crt_gauss.py --check

PYTHONPATH=elliptic-curves python3 \
  elliptic-curves/scripts/search_fermigier_denominator_offsets.py --check
```

Pinned artifact:

`artifacts/generated-results/elliptic-curves/fermigier_mixed_small_prime_crt_gauss_v1.json`

`artifacts/generated-results/elliptic-curves/fermigier_denominator_offsets_v1.json`
