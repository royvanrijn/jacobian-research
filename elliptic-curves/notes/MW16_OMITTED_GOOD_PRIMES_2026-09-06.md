# Good-prime terms omitted by the displayed MW16 models

**Exact family-wide scaling classification and finite score diagnostic pass.
The frozen narrow and broad higher-parameter searches retain their original
scores, finalists and point budgets.**

The projective caches deliberately omit primes where the displayed short
equation is singular. This agrees with their declared score and all saved
cache replays. It can nevertheless omit a good-reduction trace of the
underlying curve when the displayed equation admits a removable scale.
Unlike the earlier R17 result, some MW16 cells do restore good reduction.

## Complete support in the active prime range

Let A_h(n,d), B_h(n,d) be the five atlas binary forms of degrees 8 and 12,
with gcd(n,d)=1. For a prime p>=5, a removable short-model scale requires
p^4 dividing A_h and p^6 dividing B_h. The homogeneous resultant R has
integer Bezout identities expressing R times powers of n and d in the ideal
(A_h,B_h). At least one of n,d is a p-adic unit, so the necessary scaling
condition implies v_p(R)>=4.

All five exact 20-by-20 Sylvester determinants agree with independent Sage
polynomial resultants. Trial division through 131071 leaves precisely 33
family/prime pairs with p>=5 and resultant valuation at least four. All are
at most 101. The large residual cofactors remain unfactored; bounded trial
division still excludes their containing a prime in the active range.

Complete coefficient residue trees classify all 33 pairs, using the affine
coordinate t=n/d when d is a unit and the infinity coordinate s=d/n otherwise.
Each excluded branch fails a necessary congruence. Each admitted ball has
exact polynomial coefficients divisible by p^4 and p^6. Independent branch
partitions and Horner polynomial composition replay all witnesses. There are
27 exclusions and the following six eligible pairs:

| Family | Prime | All first-scale balls |
|---|---:|---|
| a1-fibration-01 | 5 | t=12 modulo 25 |
| a1-fibration-01 | 13 | t=6,10,11 modulo 13 |
| a1-fibration-02 | 13 | t=1 modulo 13 |
| a1-fibration-03 | 13 | t=0 modulo 13 |
| a1-fibration-04 | 13 | t=5,12 modulo 13, or s=0 modulo 13 |
| a1-fibration-05 | 13 | t=2 modulo 13 |

On each of these ten balls, division by p^4,p^6 gives integral coefficient
polynomials. A separate bounded tree excludes a second scale everywhere.
Independent enumeration of 5,023 whole-ring residues replays all ten
exclusions by depth three. At p>=5, a nonsingular short equation with no
further p^4,p^6 division is minimal at p. The discriminant of the divided
model therefore decides intrinsic good reduction on each next-digit cell.

Good reduction is restored exactly in these cells:

| Family | Prime | Good cells after scaling |
|---|---:|---|
| a1-fibration-01 | 5 | t=12,37,62,112 modulo 125 |
| a1-fibration-01 | 13 | t=6 modulo 13; t=10 modulo 13 except 127 modulo 169; t=11 modulo 13 except 167 modulo 169 |
| a1-fibration-04 | 13 | t=5 modulo 13 except 57 modulo 169; t=12 modulo 13 except 155 modulo 169; s=0 modulo 13 |

The other first-scale balls remain intrinsically bad. These statements cover
every primitive parameter in the five displayed families for primes
5 through 131071. They do not assert global minimality, conductors, or a
classification at 2,3 or beyond that range.

## Effect on the frozen narrow scalar roster

The separate diagnostic reads all 10,240 frozen higher-band scalar models,
without point outcomes or catalogue labels. It divides removable scales
exactly and checks every restored trace both by a character sum and by
enumerating all finite-field point pairs. Fifty-one reduced models are
checked. There are **904 candidates with restored good terms: 874 terms at
13 and 39 at 5**, with nine candidates contributing both.

Adding only those terms, using the existing short-band quantization, leaves
55 of the original sixty finalists in the diagnostic top-six lists. The
five changes are distributed as follows:

| Band | Family | Candidates with restored terms | Shared finalists out of six |
|---|---|---:|---:|
| 16384<H<=65536 | a1-fibration-01 | 253 | 5 |
| 16384<H<=65536 | a1-fibration-04 | 205 | 5 |
| 65536<H<=262144 | a1-fibration-01 | 236 | 4 |
| 65536<H<=262144 | a1-fibration-04 | 210 | 5 |

The other six band/family groups have no restored terms and retain all six
finalists. This diagnostic does not recover candidates lost before scalar
selection, show improved rank prediction, or authorize changing the running
experiments. A future corrected score can use local data modulo 125 and 169;
the score cannot be repaired by a table indexed only by t modulo p.

The [local correction table builder](../cas/build_mw16_local_score_corrections.py)
now covers all 1,660 cells across the five families on the two projective
rings: 150 cells modulo 125 and 182 modulo 169 per family. Affine cells use
n/d when d is a unit; infinity cells use d/n in pZ otherwise. The exact
polynomial classification proves each correction constant on its whole
cell. All 10,240 higher scalar models reproduce the table lookup, including
all 913 restored terms. Build and read-only replay pass in
[`mw16_local_score_corrections_v1.json`](../../artifacts/generated-results/elliptic-curves/mw16_local_score_corrections_v1.json).
The five binary correction files total 15,160 bytes and replay byte for byte.
A separate [corrected scanner](../cas/newfamily/scan_corrected_mw16_annulus.cpp)
passes three regression groups covering all five families, signed parameters,
both local charts, repeated periods, annulus boundaries, exact ties and
malformed inputs. It adds the restored terms before heap retention.

Its [finite benchmark](../cas/benchmark_corrected_mw16_annulus.py) passes:
twenty complete real-cache signed frames precede one untouched positive
family01 slice with 6,379,630 primitive addresses in 16384<H<=65536. The
denominator residue is 25 modulo 256 in the scanner's zero-based convention;
the original outer residue, narrow residue and sixteen broad residues are
excluded. Exactly 4,096 corrected-score candidates are retained, with
independent cache-component and full-model local-scaling checks. One worker,
120 seconds per call and a 45-second full-slice cost gate apply. The
[controller](../cas/finish_corrected_mw16_benchmark.py) also requires all
fifteen base/local binary hashes to match their upstream exact encoding
proofs. No retry, automatic wider scan or point search is included.

The complete slice takes **30.2566 seconds**, below the 45-second gate.
All twenty complete real-cache frames and top-seven orderings pass, and all
4,096 retained full-slice scores agree with independent short/extended cache
components plus exact local scaling and finite-field point counts. The
[benchmark certificate](../../artifacts/generated-results/elliptic-curves/corrected_mw16_annulus_benchmark_v1.json)
and full read-only replay pass. This clears the implementation and cost gate
for a separately frozen corrected-score campaign in new territory; it does
not itself produce a certified curve.

A separately frozen [corrected higher-population campaign](CORRECTED_MW16_HIGHER_POPULATION_2026-09-06.md)
now uses this gate for 320 slices, with one exact benchmark reuse and 319
untouched slices. Its later 10,240-scalar and sixty-point-fibre budgets stay
fixed, and both earlier experiments retain their original evidence.

The [broader early-stage experiment](BROAD_MW16_HIGHER_POPULATION_2026-09-06.md)
continues with its frozen score. Its comparison with the preserved narrow
trial therefore still compares initial coverage under the same score.

## Reproducible evidence

- [Resultant audit](../cas/audit_mw16_scaling_prime_support.py) and
  [independent polynomial check](../cas/verify_mw16_scaling_prime_support.sage):
  [`mw16_scaling_prime_support_v1.json`](../../artifacts/generated-results/elliptic-curves/mw16_scaling_prime_support_v1.json).
- [First-scale classification](../cas/classify_mw16_score_prime_scalings.py)
  and [independent tree check](../cas/verify_mw16_score_prime_scalings.py):
  [`mw16_score_prime_scalings_v1.json`](../../artifacts/generated-results/elliptic-curves/mw16_score_prime_scalings_v1.json).
- [Second-scale/reduction audit](../cas/audit_mw16_postscale_reduction.py)
  and [whole-ring check](../cas/verify_mw16_postscale_reduction.py):
  [`mw16_postscale_reduction_v1.json`](../../artifacts/generated-results/elliptic-curves/mw16_postscale_reduction_v1.json),
  [`mw16_postscale_reduction_replay_v1.json`](../../artifacts/generated-results/elliptic-curves/mw16_postscale_reduction_replay_v1.json).
- [Frozen-roster diagnostic](../cas/audit_higher_mw16_omitted_good_primes.py):
  [`higher_mw16_omitted_good_primes_v1.json`](../../artifacts/generated-results/elliptic-curves/higher_mw16_omitted_good_primes_v1.json).

Build and read-only replay pass for all stages. Support and residue-tree
checks are local independent replays; no standalone isolated archive is
claimed for this supplement.
