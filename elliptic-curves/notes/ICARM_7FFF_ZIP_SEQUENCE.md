# ICARM `7fff-zip` rank-19/20/21 sequence

## Result

ICARM curves 281, 282, 285 and 286 were submitted by `7fff-zip` on
2026-08-22, with reported rank lower bounds 19, 20, 21 and 21.  The
dependency-free checker
[`analyze_icarm_7fff_zip_sequence.py`](../cas/analyze_icarm_7fff_zip_sequence.py)
independently verifies all 81 displayed rational points and proves each
displayed set independent using exact products of good-reduction quotients
`E(F_p)/2E(F_p)`.  It also proves trivial rational torsion for each curve.

The most important consequence is curve 285:

```text
y^2 = x^3 - x^2
      - 108391804584990603814796450120*x
      + 13755098120758219348428060610562384253731136

rank >= 21
N = 1746512979501225597635228204059275382729246622846831247192424343764785450672
log(N) = 173.25150319151186...
```

Thus its 21-point lower bound is independently certified and its
ICARM-reported conductor lies well below the strict `182.72` threshold.  The
current lightweight replay does not independently run Tate's algorithm or
global minimalization, so the conductor component remains attributed to the
public database pending a PARI/Sage replay.

## Exact comparison

| curve | submitted (UTC) | rank certified here | `log(N)` | bad primes | repeated valuations | integral `x` |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 281 | 05:30:32 | 19 | 159.072016 | 10 | 6 | 7 |
| 282 | 16:53:59 | 20 | 185.165304 | 12 | 8 | 1 |
| 285 | 18:12:14 | 21 | 173.251503 | 6 | 5 | 2 |
| 286 | 18:12:50 | 21 | 177.880553 | 10 | 5 | 5 |

The discriminant factorizations, reconstructed exactly from the reported bad
primes, are:

```text
281: 2^14 3^11 11^2 23^2 79^2 463^2 61651
     * 6857983864531 * 743660058281821
     * 69519711615691712465655885193

282: 2^10 3^13 5^2 7^4 11^4 13^3 23^2 31^2 599
     * 544533839 * 11178415412615466013
     * 417544577316007236807560709864792021583751

285: 2^10 3^2 7^5 17^4 59^2
     * 5182408071918843462574265904842838694420448840522572897950269262939709

286: 2^13 3^8 5^8 19^4 193^2 313 * 30133 * 145063
     * 2052244018292853705583
     * 193026453955572059145432582176976218143
```

## What the pattern supports

The strongest common fingerprint is conductor engineering: all four
discriminants concentrate large valuations at several small primes and leave
only a few exponent-one factors.  Curve 285 is the extreme outcome, with five
of its six bad primes repeated.  This is consistent with a search that rewards
high rank while forcing or selecting repeated discriminant roots at small
primes.

The timing is also informative.  The sequence rose from 19 to 20 to 21 in
under thirteen hours, and the two rank-21 submissions arrived 36 seconds
apart.  That strongly suggests a common automated search/certification
pipeline rather than four unrelated hand constructions.  Timing alone does
not identify its family.

## What the data rule out or do not yet show

- Every pair has a different exact `j`-invariant.  No two are quadratic
  twists or alternate models of one geometric elliptic curve.
- Coefficient height and conductor do not vary monotonically with the reported
  rank.  There is no credible one-parameter interpolation visible from the
  four coefficients alone.
- Apart from the inevitable denominator root `1`, the rational-point
  denominator fingerprints have only sporadic pairwise overlap and no common
  four-curve core.  The submitted bases therefore do not expose a transported
  section basis in these coordinates.
- The public records contain no construction commentary.  Membership in a
  Mestre, Nagao, or K3 family is not established by this audit.

The next structural test should be family recognition, not coefficient
regression: compare the exact `j`-values against the repository's corrected
K3/Mestre maps, and if a specialization is found, transport the generic
sections and measure the exceptional quotient directions.

## Reproduction

```bash
python3 elliptic-curves/cas/analyze_icarm_7fff_zip_sequence.py
```

Pinned input and output:

- `artifacts/generated-results/elliptic-curves/icarm_7fff_zip_281_282_285_286.json`
- `artifacts/generated-results/elliptic-curves/icarm_7fff_zip_sequence_analysis.json`

This proves the rank lower bounds from the displayed points.  It does not
prove exact ranks, reveal the submitter's construction, or independently
recompute the four conductors.
