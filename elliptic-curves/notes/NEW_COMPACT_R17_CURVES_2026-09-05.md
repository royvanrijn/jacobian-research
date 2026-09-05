# Newly certified compact R17 specializations

The subsequent [six-family compact-atlas experiment](NEW_COMPACT_ATLAS_CURVES_2026-09-05.md)
adds six distinct curves, including rank at least 25. This note retains the
certificates and scope of the earlier fifteen-curve experiment.

Fifteen pairwise nonisomorphic prospective curves have unconditional rank lower
bounds 22–24 and no rational-isomorphism match in the pinned 584-curve ICARM
snapshot. These are actual point certificates, not Nagao estimates. They do
not yet meet the rank-at-least-28 near-record or rank-at-least-32 record target.
There is no exact-rank, conductor-record, or universal novelty claim.

| Published R17 parameter | Certified rank at least | Certified gain beyond the specialized seventeen |
|---|---:|---:|
| `33/119` | 24 | 7 |
| `-695/97` | 24 | 7 |
| `44/35` | 23 | 6 |
| `-3/148` | 22 | 5 |
| `-70/61` | 22 | 5 |
| `8/39` | 22 | 5 |
| `-2300/843` | 24 | 7 |
| `-1264/1047` | 23 | 6 |
| `2291/2392` | 23 | 6 |
| `-3895/3749` | 23 | 6 |
| `-1193/1560` | 23 | 6 |
| `1229/894` | 23 | 6 |
| `129/70` | 23 | 6 |
| `1348/1431` | 22 | 5 |
| `-7540/2317` | 22 | 5 |

The exact models, every independent point, family transports and finite
quotient matrices are in the [height-256 certificates](../../artifacts/generated-results/elliptic-curves/compact_r17_new_curves_v1.json)
and [height-4096 certificate](../../artifacts/generated-results/elliptic-curves/compact_r17_wide_new_curves_v1.json),
with eight further curves in the [top-64 interim export](../../artifacts/generated-results/elliptic-curves/compact_r17_top64_interim_curves_v1.json).
The [largest-initial-gain follow-up certificate](../../artifacts/generated-results/elliptic-curves/compact_r17_largest_gain_curve_v1.json)
adds `-7540/2317` at rank at least 22.
Every curve in that interim export already has all 43 initial chart records;
the export label refers to the unfinished cohort at export time.
The parameter refers to the literal published R17 model, not a native alternate
Q80 coordinate. The first seventeen points are exactly the transported generic
sections, checked with their common rational Weierstrass scale and sign.

For example, a convenient integral presentation of the `33/119` curve is

```text
y^2 + x*y = x^3
 - 9941757705488943928466323874475426901458080*x
 + 12011979816231939408913238860529441116586016252018764063900132352.
```

From the short model in its certificate, set `x=X-1/12`, `y=Y-x/2`.
This transports all 24 certified points. This displayed integral equation
does not require a conductor computation or a minimality assertion.

## Exact rank argument

The discovery worker admits points using the shared incremental finite
quotients. A separate Sage-free checker recomputes point membership and the
older standalone finite-group quotient construction. For each listed curve,
the columns of the retained product of `E(F_p)/2E(F_p)` have full rank equal
to the displayed lower bound. A separately checked modular irreducibility
witness for the short cubic proves `E(Q)[2]=0`.

An integral relation must therefore have every coefficient even. Dividing
by two still gives a relation because there is no rational 2-torsion.
Infinite descent excludes a nonzero relation. Neither numerical heights,
analytic rank, a Selmer dimension, nor completeness of point search is used.
The two implementations share finite group-law primitives; this is not an
external or formal verification claim.

All fifteen lower bounds can be replayed without Sage or the discovery caches:

```sh
python3 elliptic-curves/cas/certify_compact_r17_candidates.py --check \
  artifacts/generated-results/elliptic-curves/compact_r17_new_curves_v1.json
python3 elliptic-curves/cas/certify_compact_r17_candidates.py --check \
  artifacts/generated-results/elliptic-curves/compact_r17_wide_new_curves_v1.json
python3 elliptic-curves/cas/certify_compact_r17_candidates.py --check \
  artifacts/generated-results/elliptic-curves/compact_r17_top64_interim_curves_v1.json
python3 elliptic-curves/cas/certify_compact_r17_candidates.py --check \
  artifacts/generated-results/elliptic-curves/compact_r17_largest_gain_curve_v1.json
```

The checker also reconstructs each generic specialization and its point
transport, and repeats the exact snapshot comparison. Cross-certificate
`j`-invariants are distinct. A partial search can still yield a valid rank
certificate: the certificate proves independence of its explicit points,
not completeness of the worker or an upper bound on the curve.

## Discovery and novelty boundary

The height-256 experiment evaluated 79,791 signed reduced rational parameters
in the compact published family, used prime scores through 997 and then
4093, and froze sixteen finalists. Searches use the retained 43 generic deep
parities, specialized numerical height geometry only for representative
selection, and the common exact pointed-quartic sieve at height 100,000.
Worker and chart time limits, source hashes, protocol amendments, failed
attempts and exact witnesses are retained locally.

The separate height-4096 experiment evaluated 20,400,078 signed primitive
parameters. Its short-prime cutoff was subsequently found to discard the
three published rank-25/26/27 controls inside that box. The original run is
preserved. A separate full-prime population evaluates every prime through
4093 before retaining any parameter. Its initial recovery of known curves
414 and 417 to ranks 25 and 26 is calibration, not discovery of those curves.
The top-64 continuation has its own fixed roster, worker limits and checkpoints.

A separately frozen height-16384 expansion scores 326,397,350 parameters by
the same full-prime rule. The known rank-28 curve at `-9529/5471` ranks third;
known curves are removed from prospective claims by the exact post-freeze
comparison. Its bounded point measurements are a separate continuation.

## Adaptive and ambiguity checks

The two first rank-24 curves, `33/119` and `-695/97`, received a separately
frozen sparse adaptive pilot. Each enlarged rank-24 lattice supplied 1,232
parity proposals: the 43 generic deep masks plus zero, crossed with every
weight-one/two word in the seven discovered quotient coordinates. Numerical
specialized height/CVP norms scheduled 128 charts per curve, at height
100,000 and three seconds per chart. Neither first round increased the
certified rank, so the declared second-round condition did not trigger.

All 256 chart records and rank admissions replay without another sieve.
The [replay manifest](../../artifacts/generated-results/elliptic-curves/compact_r17_adaptive_replay_v1.json)
and [portable transcripts](../../artifacts/generated-results/elliptic-curves/compact_r17_adaptive_witnesses_v1.zip)
retain the complete experiment. The preceding failed startup used a mistaken
string/rational type comparison; it stopped before searching and is retained
locally as adaptive-v1, separately from the executed adaptive-v2 cohort.

The 1,307 and 758 distinct ambiguous points, respectively, were then audited
modulo 3 and 5 at all usable odd primes through 997. Both complete finite
column ranks are 24 for both curves. The [first](../../artifacts/generated-results/elliptic-curves/compact_r17_adaptive_ambiguities_0_v1.json)
and [second](../../artifacts/generated-results/elliptic-curves/compact_r17_adaptive_ambiguities_1_v1.json)
certificates replay with `audit_compact_r17_ambiguous.py --check FILE`.
No extra lower bound emerged. None of these finite failures is a dependence
proof or an upper bound; a quotient escape from a nonsaturated subgroup would
not by itself add to rational rank either.

The broader 301-class follow-up recovers the known rank-28 control from its
initial rank-26 search subgroup after eleven charts. The
[independent recovery certificate](../../artifacts/generated-results/elliptic-curves/compact_r17_blind_rank28_recovery_v1.json)
explicitly matches ICARM curve 11; it is excluded from the fifteen new
prospective curves. All three new rank-24 fibres and four score-selected fibres
with initial ranks 18/19 completed their broader 301-chart follow-ups without
an additional certified direction. Two rank-24 workers needed separately
declared 300-second continuations of their remaining frozen charts. Every
chart and admission replayed. Their retained point clouds also gave no higher
finite column rank modulo 3 or 5 through prime 997. These remain bounded results.

The height-16384 cohort finished with 58 fresh measurements, three reused
measurements and three known curves. Its unique strongest fresh initial
measurement, `-7540/2317` at rank at least 21, then completed a separate
301-chart follow-up and gained one certified direction. Its rank-at-least-22
certificate is included above; all 301 chart and admission records replayed.

A coverage audit found that the original four-second chart runs often stopped
before completing their height-100,000 boxes. A separate continuation on the
four highest-scoring fresh zero-gain fibres completed all 172 original boxes.
None added a certified direction beyond the seventeen generic sections. The
[tail replay manifest](../../artifacts/generated-results/elliptic-curves/compact_r17_tail_replay_v1.json)
and [portable witnesses](../../artifacts/generated-results/elliptic-curves/compact_r17_tail_witnesses_v1.zip)
retain the prefix-to-tail identities and completed coverage. This does not
prove any of those curves has exact rank seventeen.

## Pinned novelty comparison

Novelty means **no rational-isomorphism match in this pinned snapshot**:

```text
source: https://elliptic-rank.icarm.cloud/database.json
count: 584
SHA-256: 7e80549befa11a07422a3960967f4cd80264d8675cb3e0a99f0c9c5afb340f72
```

The certificates retain all 584 comparison equations. Equality of `j` alone
is not the test: the rational fourth/sixth-power invariant scaling is checked
as well. The public equation comparison occurs after selection is frozen;
public exceptional points are not search inputs. Absence from this database
does not prove absence from unpublished work or every historical search.
This research used the ICARM leaderboard, supported by NSF Grant DMS 2425401.

See the [machinery audit](ELLIPTIC_BREAKTHROUGH_AUDIT_2026-09-05.md) for the
selection and cache failures, broader programme gaps and continuation limits.
