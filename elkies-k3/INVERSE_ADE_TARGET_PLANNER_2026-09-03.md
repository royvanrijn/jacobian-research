# Inverse-ADE target-planner benchmark

Date: 2026-09-03.

## Outcome

Inverse ADE is now implemented as a deterministic, target-conditioned
candidate generator, but the blind benchmark **does not pass** the requested
constructive-algorithm gate.

The planner withholds each historical isotropic line and retains only:

- the parent rank-15 core and fixed rank-two bridge;
- the good prime;
- the exact parent completion-root lines prescribed to survive;
- the desired child ADE metric.

It first intersects the finite-field quadric with the survivor linear space,
then enumerates a fixed sparse-plus-dense projective parametrization.  Exact
death incidences are checked before constructing a candidate core.  Only
incidence survivors are completed with the bridge and classified against the
desired child root metric.  There is no call to Sage's random Kneser-line
generator in this planner.

The complete fixed-rule results are:

| corridor | withheld transitions recovered | materialized completions | historical raw-neighbour baseline | gate |
|---|---:|---:|---:|---|
| H3 | 7/8 | 3,182 | at least 42,300 | final rootless edge missed |
| NS0024 | 3/3 | 2,609 | 7,477 | only `7477/2609 = 2.87...x` |
| Q80 | 11/12 | 3,798 | 72,528 | final rootless edge missed |

The H3 and Q80 totals include the declared cap of 3,000 materialized
completions on the missed terminal edge.  The Q80 baseline combines 42,300
raw prefix neighbours with 30,228 defect-directed suffix neighbours.  The H3
baseline is only a lower bound because the 42,300-candidate fixed-rule run
missed before a separate search supplied the stored path.

Thus the intermediate prescribed-root transitions are often extremely
selective: before their terminal edges, H3 uses 182 materialized completions
and Q80 uses 798.  But the actual target is rootless.  On a terminal
`3A1 -> rootless` or analogous transition there is no surviving or born root
to impose a positive equality.  The condition becomes only

```text
isotropic line
+ nonorthogonality to every old root
+ emptiness of every affine birth shell.
```

The current deterministic parametrization has no constructive rule for that
last simultaneous emptiness condition.  This is exactly where it falls back
to classifying many completions.  The benchmark therefore confirms H0l as an
exact *predicate* and shows that it is not yet an effective *solver*.

## What is exact and what is experimental

The following are exact for every tested candidate:

1. the survivor linear subspace over `GF(p)`;
2. the isotropy and complete parent-root survival/death incidence;
3. the completed child root enumeration and ADE metric;
4. the absence of random Kneser-line sampling;
5. the candidate accounting and the 3,000-completion stopping boundary.

The sparse-plus-dense projective sequence is bounded and is not a complete
enumeration of the constrained quadric.  A bounded miss is not a
nonexistence theorem.  The benchmark also accepts an abstract child ADE
metric; it does not require the withheld intermediate core isometry class by
default.  Since even this weaker gate misses the two terminal transitions,
no endpoint-isometry speedup is claimed.

The exact affine-CVP pre-materialization oracle remains
[`certify_ns0024_inverse_ade_mutation.sage`](scripts/certify_ns0024_inverse_ade_mutation.sage).
Applying that oracle independently to hundreds of cells and affine layers was
slower than direct completion classification in this experiment.  A genuine
next algorithm must batch those affine queries or invert them into a finite
set of forbidden line strata; moving the same work into repeated CVP calls is
not a speedup.

## Foundry-wide decision

The entire current foundry was audited before launching any route search:

| item | count |
|---|---:|
| rootless target frames | 138 |
| target NS classes | 33 |
| source-target route pairs | 936 |
| pairs with both positive frame Gram matrices | 936 |
| pairs with complete inverse-ADE planner inputs | 0 |
| curated source-equation attempts | 6 |
| certified characteristic-zero source equations among them | 0 |

Every one of the 936 route rows lacks a common marked NS basis/source `U`, a
rank-15-core/rank-two-bridge decomposition, graph glue, a prime plan,
prescribed survival/birth templates, and an elliptic-neighbour or relative-`U`
transport.  Positive source and target frame Grams do not supply those data.
All 936 routes remain `NOT_YET_ENUMERATED` in the foundry authority.

The foundry-wide action is therefore deliberately fail-closed: do not launch
936 unconstrained searches after the blind controls failed their recovery and
`10x` gates.  More importantly, none of the six curated source searches has a
characteristic-zero equation yet.  Even a successful core-Kneser plan would
not produce the requested explicit family: core neighbours are not
elliptic-neighbour pencils, rational maps, or arithmetic descent data.

This does not refute the target-planner programme.  It identifies two concrete
missing constructions:

1. a batched inverse for the terminal affine birth-emptiness conditions;
2. a compiler from a planned frame/core path to marked primitive `U`
   embeddings and equation-level elliptic neighbours.

The standard Kneser line parametrization and asymptotic neighbour statistics
are from Gaetan Chenevier,
[*Statistics for Kneser p-neighbors*](https://arxiv.org/abs/2104.06846).
Visible-root filtering in the unimodular setting is developed further in
[*Unimodular Hunting*](https://arxiv.org/abs/2410.18788).  The graph-glue
completion layer uses Nikulin's
[*Integral symmetric bilinear forms and some of their applications*](https://www.mathnet.ru/eng/im1677).
None of those sources supplies the missing affine-shell target solver.

## Replay

```bash
sage -python elkies-k3/scripts/plan_inverse_ade_targets.sage \
  --maximum-parameter-support 3 \
  --dense-probes 150000 \
  --max-materialized 3000

python3 elkies-k3/scripts/audit_inverse_ade_foundry_readiness.py

sage -python elkies-k3/scripts/plan_inverse_ade_targets.sage --check
python3 elkies-k3/scripts/audit_inverse_ade_foundry_readiness.py --check
```

The generated records are
[`elkies-k3-inverse-ade-target-planner-benchmark-v1.json`](../artifacts/generated-results/elkies-k3-inverse-ade-target-planner-benchmark-v1.json)
and
[`elkies-k3-inverse-ade-foundry-readiness-v1.json`](../artifacts/generated-results/elkies-k3-inverse-ade-foundry-readiness-v1.json).
No `MATH_STATUS.json` entry is changed: this is a bounded algorithmic
experiment and readiness audit, not a new existence or equation theorem.
