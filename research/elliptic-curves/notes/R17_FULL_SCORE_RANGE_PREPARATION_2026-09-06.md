# Full score range on the remaining compact R17 fibrations

**The five queued R17 search/downstream controllers were cancelled before
any parameter scan or point exposure by the latest user instruction. Their
frozen protocols and waiting ledgers are preserved. The already-started cache
construction, read-only replay and byte verification have finished successfully. The next priority is the
[retained MW16 score-stratum comparison](RETAINED_MW16_SCORE_STRATA_2026-09-06.md),
after the corrected MW16 experiment finishes unchanged.**

The user requests broader initial populations while retaining the working
selection stages. The [completed broader MW16 trial](BROAD_MW16_HIGHER_POPULATION_2026-09-06.md)
supplies a finite comparison in that direction. A separate machinery gap
has now been filled on the [six compact R17 fibrations](COMPACT_SIX_R17_ATLAS_2026-09-05.md):
`103b2`, `074d9`, `07ca9`, `08234` and `08f72` now join `11952` with complete
extended projective trace caches. Earlier searches on these five used
short-prime retention followed by scalar extension. Their existing
high-parameter trials therefore do not implement the same longer scoring
before retention used in the recent `11952` and MW16 trials.

These are existing fibrations with proved generic sections. The new tables
are preserved for possible future work; their proposed new-parameter
sweep is cancelled. No target parameter, catalogue equation, public
point, jump label or rank-conditioned local residue enters this preparation.

The standalone [Sage export](../../artifacts/generated-results/elliptic-curves/remaining_five_r17_fibrations.sage)
provides all five exact Q(t) equations, their 85 generic sections and their
transported generic height forms. Running it checks every rational-function
point identity and positive definiteness of all five supplied forms, each
with determinant 948. Direct Sage execution and exact extraction replay pass.
The helper `fibre(family, parameter)` returns an integral rational fibre and
the transported seventeen points; specialized independence must still be
certified. The [export manifest](../../artifacts/generated-results/elliptic-curves/remaining_five_r17_fibrations_export_v1.json)
pins the original atlas and [generator](../cas/export_remaining_five_r17_fibrations.py).

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

The full cache build and its read-only replay now pass for all 14,740 tables.
The ten binary encodings replay exactly, and the independent compiled reader
agrees on both score components for all 4,831 saved fixtures: 9,662 checks.
See the [cache certificate](../../artifacts/generated-results/elliptic-curves/r17_remaining_extended_projective_caches_v1.json)
and [binary/replay certificate](../../artifacts/generated-results/elliptic-curves/r17_remaining_joint_binary_caches_v1.json).
These are arithmetic and encoding results; they launch no new fibre search.

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

The [binary encoder](../cas/encode_r17_remaining_joint_caches.py) and its
[waiting controller](../cas/finish_r17_remaining_binary_caches.py) follow
successful table replay. All ten short/extended cache files are checked byte
for byte. The independently compiled retained-list reader must reproduce
both score components and good-prime counts for every one of the 4,831
previously scalar-scored outer candidates on these five families: 9,662
component comparisons. These fixtures are existing arithmetic checks, not
new candidate selection. The explicit reader height limit is 524288.

## Cancelled-before-execution fresh-fibre campaign

The earlier request to fill this five-fibration gap produced the following
frozen protocol. Its scan and downstream controllers have now been stopped
before execution; the following allocation is preserved historical scope. The separate
[scanner](../cas/scan_r17_remaining_higher_annuli.py) and
[controller](../cas/finish_r17_remaining_higher_scan.py) freeze the following
before new scores or point outcomes exist:

| Stage | Fixed allocation |
|---|---:|
| Families | Five, excluding `11952` |
| Height bands | 32768 < H <= 131072; 131072 < H <= 524288 |
| New signed denominator slices | 320 |
| Primitive parameter addresses | 3,059,808,912 |
| Retention per slice | 4,096 |
| Complete-score survivors | 1,310,720 |
| Fresh scalar candidates | 10,240 |
| Distinct prospective point finalists | 60 |
| Maximum generic17 point boxes | 2,652 |

For each family, sign and band, a fixed SHA256 seed chooses sixteen residues
of the prescribed parity, modulo 1024 and 16384 respectively. Every slice
excludes the earlier outer131072 residue modulo 1024, and the applicable
skew-rectangle residue modulo 64. Earlier square scans through 32768 fall
below the inner cut. These exclusions are stored explicitly with every
slice; the complete new slices are mutually disjoint within each family.
Parameter freshness does not establish equation novelty.

Every address receives all 3,510 selection primes through 32749 before
retention. All 320 complete signed actual-modulus frames and top-seven
orderings must agree with the independent retained-list reader. The first
full slice in each band has a 45-second gate and is reused once within this
campaign. Four workers, 120 seconds per call, a 7,200-second main cap and
1,800-second replay cap apply. All retained component scores and exact
primitive counts must replay; immutable files checkpoint individual slices.

The [scalar stage](../cas/score_r17_remaining_higher.py) chooses 1,024
within-roster Q-distinct equations per band/family, applies fresh traces
through 65521, and selects six per group. It retains the existing four-worker,
20-second per-case, 2,400-second total and first-twenty cost gates, with
20,480 independent character sums. Disjoint primes 65537 through 131071
validate the frozen finalists without changing them.

The [point controller](../cas/finish_remaining60_r17_points.py) starts only
after those proofs pass. Every fibre starts from its seventeen generic
sections. Four families use their 43 recorded generic parity labels; `08f72`
uses its 49 labels. Twelve finalists per family therefore give exactly
2,652 possible boxes at height 125000 and ten seconds per chart. All sixty
map files precede every point attempt. Two workers, 600 seconds per fibre
and a provisional 28-direction stop apply. There is no adaptive wave,
refill or retry.

Exact histories, generic transports, rational maps and full-cloud independence
modulo 2, 3 and 5 precede any catalogue access. The
[finalizer](../cas/finalize_remaining60_r17_results.py) additionally waits for
the corrected MW16 cohort before comparing against all 1,345 prior equations
and the 593 pinned catalogue equations. The
[portable controller](../cas/finish_remaining60_r17_portable.py) then packages
and runs 182 isolated point-only proof stages. No known-record parameter,
equation, point, rank, j-invariant or jump label enters prospective selection
or execution. No rank, record or universal-novelty result is asserted by
these queued protocols.
