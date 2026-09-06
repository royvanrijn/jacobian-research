# Full score range on the remaining compact R17 fibrations

**All six R17 models are now proved minimal at every prime from 5 through
131071 except the separately classified prime 13. The five fibrations other
than `11952` pass all fifteen extended-cache cost cases and their replay.
The fixed 14,740-table build is active. No new parameter or point campaign
has started from these tables.**

The user requests broader initial populations while retaining the working
selection stages. The [completed broader MW16 trial](BROAD_MW16_HIGHER_POPULATION_2026-09-06.md)
supplies a finite comparison in that direction. A separate machinery gap
remains on the [six compact R17 fibrations](COMPACT_SIX_R17_ATLAS_2026-09-05.md):
`11952` has a complete extended projective trace cache, while `103b2`,
`074d9`, `07ca9`, `08234` and `08f72` have short-prime projective caches
followed by scalar extension only after short-score retention. Their existing
high-parameter trials therefore do not implement the same longer scoring
before retention used in the recent `11952` and MW16 trials.

These are existing fibrations with proved generic sections, not new generic
groups. The intended use of new tables is to score previously unscanned
parameter territory on them. No target parameter, catalogue equation, public
point, jump label or rank-conditioned local residue enters this preparation.

## Prime range theorem

The [earlier theorem](R17_SMALL_PRIME_MINIMALITY_2026-09-06.md) excludes a
removable scale for all 990 family/prime pairs with 5 <= p <= 997, p != 13.
For a primitive integer parameter pair (n,d), a removable short-model scale
at p >= 5 requires p^4 dividing A_h(n,d) and p^6 dividing B_h(n,d).
The integral homogeneous resultant Bezout identity then forces
v_p(Res(A_h,B_h)) >= 4: one of n,d is a p-adic unit, and the two projective
Bezout identities cover the respective charts.

The new audit recomputes all six 20-by-20 Sylvester determinants and checks
every one of the 12,083 primes in 1009 through 131071. Only one additional
factor occurs: prime 1213 in family `074d9`, with valuation one. Thus every
one of the additional 72,498 family/prime pairs has resultant valuation
below four and admits no removable scale. Combined with the earlier theorem,
this proves minimality at **73,488 family/prime pairs** through 131071,
excluding 13. This statement concerns nonsingular specializations; the
divisibility exclusion itself holds for every primitive parameter pair.

Independent Sage polynomial resultants and a separately enumerated complete
prime interval agree with all determinants, valuations and residual cofactors.
Both generation and read-only replay pass:

- [Exact range certificate](../../artifacts/generated-results/elliptic-curves/r17_extended_score_prime_minimality_v1.json)
- [Independent Sage replay](../../artifacts/generated-results/elliptic-curves/r17_extended_score_prime_minimality_sage_replay_v1.json)
- [Generator and protocol](../cas/extend_r17_score_prime_minimality.py)
- [Independent verifier](../cas/verify_r17_extended_score_prime_minimality.sage)

The separately proved [13-scaling classification](R17_INTEGRAL_13_PARAMETER_CHARTS_2026-09-06.md)
is unchanged: scaled models remain bad at 13. Hence the MW16 restored-good-prime
defect does not recur in these R17 score intervals. Primes 2 and 3 and primes
above 131071 remain outside this theorem. Large residual cofactors remain
unfactored and complete global prime support remains UNKNOWN. This audit
establishes neither a conductor nor a new rank bound.

## Bounded table gate

The [remaining-five benchmark](../cas/benchmark_r17_remaining_extended_projective_tables.py)
uses the first, middle and last selection primes in 4099 through 32749,
equally on each of the five families lacking the full extended cache. Each
case computes all p+1 projective residues, verifies every discriminant/Hasse
frame and checks five independent full character sums. Existing `11952`
tables are not recomputed. This is fifteen GP calls, one worker, twenty
seconds per case and 300 seconds total, with immutable raw transcripts.
Replay makes no new GP calls. Each family must project below 1800 serial
seconds by the same worst measured cost per projective residue. Projection
is a scheduling estimate, not a mathematical runtime bound.

The primes are 4099, 17749 and 32749. All fifteen tables and 75 independent
character sums pass, followed by read-only replay. Per-family projected
serial times range from 1641.91 to 1686.38 seconds, all below the fixed
1800-second threshold. The complete
[benchmark certificate](../../artifacts/generated-results/elliptic-curves/r17_remaining_extended_projective_benchmark_v1.json)
retains individual timings and exact table hashes.

The separately frozen [full-cache builder](../cas/build_r17_remaining_extended_projective_caches.py)
now covers all 2,948 primes from 4099 through 32749 on each of the five
families: 14,740 complete projective tables. It reuses exactly the fifteen
benchmark transcripts and makes at most 14,725 new GP calls, with five
workers, twenty seconds per call, eighty-table checkpoints and a 7,200-second
build cap. Every table gets the same full frame and five independent character
sum checks. A separate 3,600-second read-only replay follows. The
[controller](../cas/finish_r17_remaining_projective_caches.py) stops on failure
or censoring without a retry. These tables contain 264,948,100 projective
residues and require 73,700 independent character sums per full pass.

The user explicitly requests filling this five-fibration gap. After the cache
and binary-score checks pass, the next finite protocol will select fresh
high-height parameter slices, apply all selection primes before retention,
and retain bounded scalar and generic17-only point stages. It must exclude
the earlier square and skew-rectangle scan slices explicitly; parameter
freshness and rational-isomorphism novelty remain separate checks. No
compact-height revisit, adaptive point exposure or record-conditioned
candidate selection is part of this preparation.
