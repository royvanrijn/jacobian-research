# Fermigier rank-20 relative 2-descent experiment

Status: **experiment / external CAS task**, not a rank claim.

## Question

For the certified rank-at-least-20 Fermigier specialization at adapter parameter
`u = 28917/20`, does a 2-descent relative to the known 20-dimensional subgroup
leave any Selmer classes that can be converted into a new Mordell--Weil
direction?

The pinned near-miss is the curve represented in the repository by
`artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json`.
Its exact global minimal model is

```text
[1, 1, 1,
 -4437412060110743641525245114305,
 3586842216822165612930264910099076801587288127]
```

The existing certificate proves only `rank >= 20`; the selected subgroup is not
claimed saturated.

## Reproduction

First generate the Magma job from repository data:

```bash
python3 elliptic-curves/cas/build_fermigier_rank20_relative_descent.py \
  --manifest artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json \
  --output artifacts/local/elliptic-curves/fermigier_rank20_relative_descent.m
```

Then run it with Magma 2.29 or later and retain the complete stdout/stderr log:

```bash
timeout 48h magma -b \
  artifacts/local/elliptic-curves/fermigier_rank20_relative_descent.m \
  > artifacts/local/elliptic-curves/fermigier_rank20_relative_descent.log 2>&1
```

The generator reconstructs the same searched quartic point cloud from the
pinned abscissas, selects the exact 20 indices stored in the manifest, verifies
all selected points on the canonical model, and emits those points as
`RemoveGens` for Magma's `TwoDescent`.

## Required output

The generated job prints machine-readable records beginning with
`R20DESCENT|`. Preserve at least:

- Magma version and exact elliptic-curve coefficients;
- the 20 known points and the fact that they lie on the curve;
- number of residual 2-covers returned by `TwoDescent`;
- the full Cassels--Tate pairing matrix on those residual covers, when the
  computation reaches that stage;
- number of 4-covers returned for every residual 2-cover in the pairing radical;
- bounded `PointsQI` results for each returned 4-cover;
- wall-clock outcome (`complete`, `timeout`, or `error`).

Raw logs stay under `artifacts/local/elliptic-curves/`. A compact deterministic
summary may be promoted to `artifacts/generated-results/elliptic-curves/` only
after manual inspection.

## Classification

Classify a completed run into exactly one of the following buckets.

### R0 — no residual 2-cover

`TwoDescent` returns no cover after quotienting by the 20 supplied generators
(and torsion). This is strong evidence that this 2-descent route is exhausted.
It is **not** by itself an exact-rank claim unless Magma's returned data are
converted into a rigorous upper bound with all hypotheses recorded.

### R1 — residual classes, nonzero Cassels--Tate obstruction

Residual 2-covers exist and at least one pairing entry is nonzero. Those
classes can be explained by nontrivial 2-primary Tate--Shafarevich information;
do not count them as rational points.

### R2 — pairing-radical classes, no rational 4-cover point found

A residual class lies in the Cassels--Tate radical and `FourDescent` succeeds,
but the bounded `PointsQI` search finds no rational point. This is the main
**continue descent / increase search** bucket, not a negative result.

### R3 — candidate new rational point

A 4-cover produces a rational point that maps to the original elliptic curve.
This is only a candidate until the mapped point is checked exactly, shown to
escape the certified rank-20 subgroup, and the enlarged set is independently
certified. If successful, build a separate rank-at-least-21 candidate artifact
and re-check the exact conductor of its global minimal model.

### R4 — incomplete

Timeout, memory exhaustion, Magma exception, unsupported parameter/API, or an
unresolved local-solubility computation. Record the last completed stage. This
bucket has no mathematical interpretation.

## Promotion gate for an R3 result

A newly mapped point must pass all of these before any status change:

1. exact on-curve verification;
2. exact finite-quotient escape from the existing subgroup;
3. a new independence certificate for 21 points;
4. saturation attempt / index information where feasible;
5. replay on the global minimal model;
6. exact conductor replay (the current curve already satisfies `log(N)<182.72`);
7. independent implementation before a public record-style claim.

Until that gate is passed, `MATH_STATUS.json` must not be changed.

## Why this experiment is prioritized

The current fixed-fiber point search is already broad enough that simply
raising another height box has diminishing value. Relative descent asks a
different arithmetic question: whether the 2-Selmer information left after the
known rank-20 subgroup contains any unresolved direction. The Cassels--Tate and
4-descent stages distinguish a plausible Mordell--Weil direction from a
Tate--Shafarevich obstruction instead of treating residual Selmer dimension as
rank.
