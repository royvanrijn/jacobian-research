# Inventory188 / ICARM619: M28 continuation complete, no certified gain

The native R17 fibre `11952` at `110314/102227` remains at certified lower
bound28 after900 completed invocations. Independent replay passes. The
retained cloud has473 distinct recorded points and finite ranks28 modulo2,3,5.
There are no timeouts. This is not an exact-rank result or absence of a29th direction.

## Seed and bounded experiment

Preparation reconstructs the native seventeen by the compact atlas and exact
scale156. It rechecks the retained17→27 history and the later fixed49-chart
27→28 recovery. The known28th direction is seed input, not a discovery in this
experiment. The curve is the already-published ICARM619 rank-at-least28 example.
The two histories provide nine productive generic masks:
51446,63347,69511,76282,83323,83747,86925,91511,96095.
Their generic minima are checked with both integer and original rational CVP.

The M28 landscape scores18432 extension classes, selects450 centres and takes
59.414 seconds. Both quartic-minimized and factor-free coordinate maps run on
each selected centre. The whole search takes637.509 seconds with peak RSS
363581440 bytes. Independent replay takes187.306 seconds, peak319451136 bytes.
It regenerates the full landscape using the original rational CVP solver,
both coordinate maps and the invocation schedule, exact point witnesses,
and complete-cloud certificates.

Frozen limits: height125000,10 seconds per invocation,4096 actual invocations,
target32, at most five adaptive epochs, one worker, and separate
7200-second/3-GiB search and replay supervisors. No gain triggers a second epoch.
The tested parent bank is the explicit historical productive subset, not all
generic parity classes. No conductor record or general success rate is inferred.

## Evidence

- [Sealed result, full seed and completion records](../../artifacts/generated-results/elliptic-curves/curve188_productive_v3_v1/result.json).
- [Native M28 preparation](../cas/prepare_curve188_productive.py).
- [Bounded runner and independent replay](../cas/run_productive_seed_v3_v2.py).
- `artifacts/local/elliptic-curves/curve188-v3-preparation-v1/`.
- `artifacts/local/elliptic-curves/curve188-productive-v3-discovery-v1/`,
  including the frozen source archive, plus sibling supervision directories.

## Separate admission diagnostic

The exact classification is complete. All445 sign-normalized nonbasis points
are integer combinations of the28 seed points: every final relation has
multiplier1. Classification finishes in526.371 seconds, peak RSS299515904
bytes, under the600-second/1-GiB limit. All decisions finish in2–4 steps,
below the eight-step cap. No new point-box search is performed.

A separate Fraction-arithmetic verifier checks every relation, reconstructs
coverage directly from the retained cloud, and freshly checks the28-point
independence certificate. Thus the retained point subgroup has exact rank28;
indeed these retained points add no generator to the seed subgroup. This
closes the finite-fingerprint admission gap for this cloud only. It gives no
upper bound for the elliptic curve and does not classify other searches.

- [Full span certificate and445 explicit relations](../../artifacts/generated-results/elliptic-curves/curve188_retained_cloud_span_v1.json).
- [Classifier driver](../cas/classify_productive_cloud.py).
- [Independent full-cloud verifier](../cas/verify_productive_full_cloud_span.py).
- Raw frame and per-point receipts:
  `artifacts/local/elliptic-curves/curve188-cloud-classification-v1/`.

## New parent-class follow-up

A separate follow-up uses all36 distinct pairwise XORs of the nine productive
masks, excluding zero and the original masks. Their exact generic minima are
independently replayed in1.637 seconds. The new classes have generic norm4
(one),6 (seven),8 (fifteen), and10 (thirteen). Every shell fits within the existing
V3 quota of16 anchors, so all36 are retained without increasing the search
height or changing those quotas. They are untested classes derived from
productive ones; their productivity is not assumed.

[Preparation](../cas/prepare_pairwise_parent_bank.py) binds the class derivation,
M28 seed and generic CVPs. [The parent-subset runner](../cas/run_parent_seed_v3.py)
retains the same search, point admission and replay mechanisms while accepting
the explicitly labelled new parent bank. The run is under
`artifacts/local/elliptic-curves/curve188-pairwise-v3-discovery-v1/`, with separate
7200-second/3-GiB search/replay bounds and4096 actual invocations.

The follow-up completes3454 invocations on1727 distinct centres, with no
timeouts and no certified gain: M28→M28. The initial centre set is disjoint
up to sign from the previous450-centre policy on the same ordered basis;
the exact comparison is retained in `initial-centre-disjointness.json`.
The new landscape takes239.795 seconds. Search completes in2387.337 seconds,
peak RSS590807040 bytes. Independent replay passes in751.599 seconds,
peak RSS524881920 bytes, and the
[sealed result](../../artifacts/generated-results/elliptic-curves/curve188_pairwise_v3_v1/result.json)
is exported. Complete-cloud finite ranks modulo2,3,5 are all28.

The earlier445 exact dependence relations apply only to the original900-chart
cloud; they do not classify this new cloud. This second bounded miss supplies
neither an exact rank for the curve nor absence of a29th direction. Its point
searches account for84.79% of search wall time, motivating better scheduling
of future parent/seed candidates before another chart-budget expansion.
