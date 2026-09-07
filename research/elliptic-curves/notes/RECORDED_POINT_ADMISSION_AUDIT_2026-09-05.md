# Recorded points can outrun the independence admission budget

The blind compact-MW16 search had already found enough points to prove rank
at least **26** on the subsequently identified ICARM 542 control, although its
worker reported only 25. The finite admission primes stopped at 251. Adding
prime **257** certifies a 26th independent direction from the retained cloud,
without using any public point. This is a certification loss, not a missing
point or a new curve.

The subsequent exact audit covers **202 retained search transcripts** and
**57,599 point occurrences up to sign within each transcript**. These are not
202 distinct curves or independent experiments. All 202 rank certificates
passed separate replay. Only the known 542 transcript gains a direction;
the 32 distinct new curves retain their earlier certified lower bounds, with
three at least 25. Neither a new rank-28 nor a new rank-32 curve is certified.

The [portable audit manifest](../../artifacts/generated-results/elliptic-curves/recorded_mod2_admission_audit_v1.json)
indexes every input hash, point-cloud proof, original lower bound and audited
lower bound. Its accompanying ZIP contains the standalone point proofs,
checker sources, protocols and successful replay logs. The original search
snapshots remain preserved separately; their completion claims are not
required for these rational-point lower bounds.

## The exact control and failure

The control is `a1-fibration-04`, compact parameter `-1905/52`, in the fixed
height-4096 MW16 roster. Its 43 initial charts completed their declared plans
before the catalogue comparison. The stored cloud has 865 distinct points
up to sign, including the original 25 independent points.

Exhaustive finite quotient computations give the following ranks of this
cloud's mod-2 image:

| Prime prefix | Available quotient dimensions | Cloud image rank |
| --- | ---: | ---: |
| through 211 | 26 | 25 |
| through 251 | 32 | 25 |
| through 257 | 33 | 26 |
| through 997 | 122 | 26 |

Having more available quotient rows than the desired rank does not ensure
those rows separate the rational points. The worker's `AMBIGUOUS_FINITE_REDUCTIONS`
label correctly left the lower bound at 25, but it was incorrect to interpret
that measurement as evidence that a 26th direction had not been found.

The replayable [rank-26 proof](../../artifacts/generated-results/elliptic-curves/recorded_rank25_mod2_recertification_0_v1.json)
selects columns `0..24,27` from the retained cloud. Its last point, on the
short model in that artifact, is

```
x = 4515295075368258616073647635034/1494413283
y = 1204771734849394056900512987780621008244331750/11117936687759
```

All coordinates satisfy the curve equation over Q. Full-rank finite mod-2
images, together with the separately checked absence of rational 2-torsion,
prove independence: a primitive integral relation would have some odd
coefficient and contradict those images. No numerical height, analytic rank,
search score, or completeness assumption enters this lower bound.

## Why raw public-point visibility gave an incomplete diagnosis

An explicitly retrospective audit transported all 26 public control points
and checked both signs in all 43 frozen charts. All 2,236 displayed
representatives were outside the selected boxes. The escaping public point
at zero-based index 17 had minimum displayed height 187,924,541,015.

That does not describe all translates by the discovered subgroup. A separate
floating-CVP proposal step chose one translate for each sign and chart.
Independent rational group arithmetic replayed all 86 resulting words and
their exact chart maps: **38 were visible and recorded**, and 48 outside the
boxes. None was visible in completed coverage but missing from the output.
No CVP optimality is claimed.

One recorded translate has chart coordinate `[-257,253]`, height 257, and
point

```
x = 10265254988536785272073124/363
y = -6133007986093269833545773910345650750/1331
```

Appending it to the original 25-point subgroup also certifies rank 26. The
[initial visibility artifact](../../artifacts/generated-results/elliptic-curves/curve542_initial_visibility_v1.json),
[86 translations](../../artifacts/generated-results/elliptic-curves/curve542_translated_visibility_v1.json),
and [independent replay](../../artifacts/generated-results/elliptic-curves/curve542_translated_visibility_replay_v1.json)
retain the exact maps and group words. These oracle calculations never fed
the prospective selector, chart centres, worker, or bulk rank audit.

## Scope of the wider retained-cloud check

Two frozen rosters used all good odd primes at most 1,000, with early stopping
only at rank 32. Each job had a 300-second build and 180-second replay limit,
1.5 GiB memory cap, and at most two concurrent jobs.

- 50 retained compact/prospective transcripts with original lower bounds at
  least 22: 37,900 points counted within their transcripts. Only the known
  control improves, from 25 to 26.
- 152 retained transcripts with original lower bounds 16–21 and at least
  one additional stored point up to sign: 19,699 point occurrences. None
  improves its lower bound.

Previous checkpoint copies and the then-running newest rank-25 follow-up
were excluded. Stopped source searches were permitted as explicitly
hash-pinned snapshots; the new certificates do not turn them into complete
searches. Exact point membership and quotient tables were recomputed, the
independent subset extracted, and its standalone proof replayed. A null
change bounds only what these finite reductions certify, not the rank of
the curve or even the rational span of its whole point cloud.

That follow-up subsequently completed all 301 planned charts and passed
exact chart/admission replay. A separate certificate checks its 3,739-point
cloud through prime 997, still at least rank 25. Its
[diagnostic bundle](../../artifacts/generated-results/elliptic-curves/rank25_followup_diagnostics_v2.json)
also retains a fixed sample of 6,144 new specialized parity classes across
the three new rank-25 curves. All exact norm/parity checks passed, but all
three samples fell below the declared 5% median-norm improvement gate;
they triggered no additional point searches.

## A separate cache cost and a compatible option for future workers

The newest rank-25 curve's 217-chart checkpoint contains 178,112 immutable
`finite-field/point-mod2` facts, occupying 93,546,992 bytes even as compact
JSON. It contains only 44 distinct finite quotient facts. Persisting every
point mask is unnecessary when exact point coordinates and the complete
verified quotient table are retained.

The optional `QuotientOnlyReductionCache` retains those quotient witnesses
and uses a bounded in-memory point-mask cache. A fixed control validation
compared **19,895 point-prime masks** with the independently replayed 865-point
proof: all agreed, with only 23 stored quotient facts, zero stored point
facts and 512 in-memory point entries. Tests additionally cover eviction,
nonintegral points reducing to infinity, a general Weierstrass model,
off-curve rejection and state replay using the original cache.

This is an optional implementation for newly frozen workers. Existing cache
sources, worker protocols and old certificates remain unchanged. No speedup
factor is asserted from these storage measurements.

## Reproduction

From the repository root, or after extracting the proof bundle:

```sh
python3 elliptic-curves/cas/audit_recorded_point_mod2_rank_v2.py --check \
  artifacts/generated-results/elliptic-curves/recorded_high_rank_mod2_recertification_47_v2.json
python3 elliptic-curves/cas/audit_curve542_initial_visibility.py --check \
  artifacts/generated-results/elliptic-curves/curve542_initial_visibility_v1.json
python3 -m unittest discover -s elliptic-curves/tests -p test_quotient_only_reduction.py
```

The first command also passed from an isolated extraction of the bundle.
To audit another retained snapshot, supply its exact SHA-256 with
`--input-sha256`, a new output path and an explicit prime bound. Keep the
original transcript and its supervision outcome. A larger prime allowance
is a finite certification policy, not a theorem that all independent points
will be recognized.
