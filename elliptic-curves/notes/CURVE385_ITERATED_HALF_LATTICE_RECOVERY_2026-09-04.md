# Curve 385 iterated half-lattice recovery (2026-09-04)

## Result

The first quotient-aware iteration succeeds decisively. Starting with only the
specialized generic `MW17` and the points found by the previously frozen
generic-deepest 43 charts, exact finite-reduction certificates select three
new independent directions. All ten nonbasis points from those charts are
integral combinations of the resulting basis, so the discovered subgroup is
exactly `M20`; its saturation index *inside the discovered group* is one.

The next round searches the `43(2^3-1)=301` lifted parity classes involving a
nonzero new quotient word. All 301 pointed quartics are freshly minimized and
reduced and all bounded searches complete. The exact blind classifier finds
nine further independent directions:

```text
M17  -- old deep 43 -->  M20  -- 301 quotient-bit lifts -->  M29.
```

Only after the blind ledger was frozen did the verifier load the 29 displayed
public points. The blind and public bases have mutually inverse integral
coordinate matrices, both of determinant `1`. Hence the blindly recovered
`M29` is exactly the displayed public rank-29 subgroup, and all twelve public
directions modulo generic `MW17` were recovered. Equivalently, if the fourth
direction from the earlier generic/specialized union is inserted first, the
iteration recovers the eight remaining displayed directions.

This is an exact subgroup-recovery statement, not an exact-rank theorem. The
displayed public group has rank at least 29; it is not proved to be all of
`E(Q)`.

## Quotient-weight structure

The first enlargement is primitive inside the discovered group: the retained
`M20` basis is the generic `M17` basis followed by three new generators, and
the exact discovered-group classifier performs no finite-index enlargement.
Thus the runner's chosen coordinates give the literal splitting

```text
M20/2M20 = M17/2M17 + F2^3.
```

The old generic-deepest set `D` has 43 classes.  Including its already-searched
zero-word slice, the structured search space is `D x F2^3`; the fresh round is
`D x (F2^3 - {0})`.  Exact integral coordinates in the final blind `M29` basis
give the following quotient ranks:

| quotient words | charts | distinct iteration discoveries | rank over `M20` |
|---|---:|---:|---:|
| weight `1` | 129 | 125 | 7 |
| weight at most `2` | 258 | 261 | 9 |
| all nonzero words | 301 | 287 | 9 |

All nine basis-extension events have quotient weight one or two.  Consequently
the weight-at-most-two subcylinder already generates the full discovered
`M29/M20`; the 43 weight-three charts have zero marginal rank gain.

Hamming weight is basis-dependent.  The complete census of the 28 unordered
bases of `F2^3` finds weight-at-most-two quotient rank nine for 20 bases and
rank eight for eight bases.  Weight one alone has rank nine for one basis and
rank seven for the natural runner basis.  This is exact posthoc structure in a
completed blind ledger, not a prospective success probability.

The compact replay is
`curve385_quotient_weight_profile_v1.json`, generated without loading the
public rank-29 fixture.

## Frozen sparse rank-32 protocol

The next search is no longer the monolithic `43(2^12-1)=176085` round.  The
separately frozen protocol starts each discovered lattice state with all 516
natural weight-one charts and then the remaining 2,838 natural weight-two
charts.  Any exact rank or finite-index enlargement causes the height lattice
and quotient complement to be recomputed and the search to restart at weight
one.  A certified rank of at least 32 stops the campaign successfully.

If both natural stages miss, two SHA-256-derived, precommitted quotient bases
are available before the natural weight-three stage.  Their physical quotient
words are deduplicated exactly against earlier stages.  Timeouts, PARI
failures, unclassified points, stage limits, and completed sparse misses remain
fail-closed and provide no rank upper bound.

The primary natural-basis campaign has now completed: all 3,354 planned charts
are accounted for by 3,116 fresh completed searches and 238 exact prior-chart
skips, with zero timeouts or PARI failures.  Both exact classifiers leave the
rank-29 basis unchanged.  This bounded no-growth result and its complete
compressed ledger are recorded in
[`CURVE385_SPARSE_QUOTIENT_RANK32_PRIMARY_2026-09-04.md`](CURVE385_SPARSE_QUOTIENT_RANK32_PRIMARY_2026-09-04.md).
Rank at least 32 remains open.

## Frozen search

The blind runner reads only
`half_lattice_search_ablation_rank29_holdout_blind_v1.json`, whose prior search
did not contain public exceptional-point coordinates. It does not import or
read the public 29-point fixture.

For `M20`, the canonical-height Gram matrix is recomputed at 110 decimal
digits. Every old deep class is lifted through all eight new-bit words, but the
already searched zero word is omitted. The 301 operative representatives are
the shortest vectors for the height form rounded at scale `10^6`, ordered by
their actual decimal canonical depths. The complete priority order agrees at
the audit scale `10^5`; 145 representatives differ, so representative identity
remains numerical rather than exact lattice evidence.

Each quartic uses one PARI `hyperellminimalmodel`, `hyperellred`, and
`hyperellratpoints` call with:

- reduced-coordinate height bound `100000`;
- 15-second wall limit per quartic;
- 1 GB GP stack;
- no retry.

The round takes 187.4 wall seconds on the recorded host. It has 301 completed
searches, zero timeouts, and zero PARI failures. There are 310 finite point
occurrences on 118 charts and 296 distinct blind discoveries after adjoining
the ten initial points. Every discovery has an exact integral relation in the
final blind basis. The nine independent additions first occur on charts with
priority numbers `43,111,184,189,213,218,220,221,235` (listed in search order,
not in the classifier's point-height order).

## Stop boundary

The run does **not** claim bounded-search stability at `M29`. Once rank 29 is
reached, the same rule would request

```text
43(2^12-1) = 176085
```

new-bit lifts. This exceeds the predeclared four-bit / 688-total-lift first
campaign limit, so the runner stops with
`STOPPED_AT_DECLARED_LIFT_LIMIT`. A bounded miss would not prove exact rank in
any event.

## Replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage \
  --height-bound 100000 \
  --timeout-seconds 15 \
  --stack-bytes 1000000000 \
  --max-quotient-bits 4 \
  --max-planned-lifts 688

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/verify_curve385_iterated_half_lattice_search.sage

python3 elliptic-curves/cas/analyze_curve385_quotient_weight_profile.py --check
python3 elliptic-curves/cas/build_curve385_sparse_quotient_rank32_protocol.py --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_curve385_sparse_quotient_rank32_search.sage \
  --plan-only --max-stage 6

python3 -m unittest -v \
  elliptic-curves/tests/test_curve385_iterated_half_lattice_search.py \
  elliptic-curves/tests/test_curve385_sparse_quotient_rank32.py
```

The blind artifact SHA-256 is
`356001898f738f607d984e081663a015825e11de0c606d35055af156eb2d7502`.
The post-freeze verification artifact SHA-256 is
`b281556f5d08250f67b69b2c62a640ac17ba4d03325e4402e85c7d60882c3ae5`.
The exact quotient-weight profile SHA-256 is
`c321d1b40d9e5fc77ebff64e5d6584feeab5f503b13eadda4f6d524d0e38162a`.
The prospective sparse protocol definition hash is
`5723679da2907e036095f90376cdabde457a4f7ba5bc284ad4a4ca3edea1aa37`;
its whole-file SHA-256 is
`2c9150f50f305b8aa3763590cd5e81c4d7e121f9373177827780789ce472834f`.
