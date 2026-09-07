# Inverse-ADE target-planner benchmark

Date: 2026-09-03.

<!-- status-consumer: EC-K3-INVERSE-ADE-PROJECTIVE-BIRTH-STRATA b4a7edb452e6dcc7 -->

## Outcome

The recovery-first v2 benchmark passes the requested gate

```text
recover every withheld historical edge.
```

It does so by strengthening the two terminal rootless fixtures with a marked
target core in the parent rational quadratic space. The historical isotropic
line is still absent as an input field, but it is mathematically determined by
the marked parent--target intersection. This is intentionally stronger than
the ADE-only v1 fixture and must not be advertised as an ADE-only speedup.

| corridor | withheld transitions recovered | materialized completions | historical raw-neighbour baseline | speed comparison |
|---|---:|---:|---:|---|
| H3 | 8/8 | 183 | at least 42,300 | deferred: terminal input strengthened |
| NS0024 | 3/3 | 2,609 | 7,477 | `7477/2609 = 2.87...x` |
| Q80 | 12/12 | 799 | 72,528 | deferred: terminal input strengthened |

Thus the correctness gate is closed, while the `10x` gate remains open. H3
and Q80 materialization counts are not comparable with their historical
ADE-only baselines because their terminal target marking was added. NS0024
retains the old inputs and remains far below an order-of-magnitude reduction.

## Why the terminal H3 and Q80 lines were difficult

The v1 terminal fixture prescribed a rootless child. It therefore had no
surviving root to impose as a positive modular equality. On both missed edges
the equality matrix has rank zero and its kernel is the full 15-dimensional
core space. This is the predicate with no pruning power:

```text
prescribed surviving-parent-root equalities.
```

The nonzero death inequalities still reject some lines, but only after the
zero-rank equality stage has left the full quadric. The exact replay of the
first 3,000 v1 materialized rejects gives:

| terminal edge | isotropic proposals needed | rejected by death incidence | materialized rejects | child core has roots | core rootless, completion regrows graph-glue roots |
|---|---:|---:|---:|---:|---:|
| H3 `3A1 -> rootless`, `p=11` | 3,997 | 997 | 3,000 | 2,068 | 932 |
| Q80 `2A1 -> rootless`, `p=29` | 3,218 | 218 | 3,000 | 2,521 | 479 |

Every one of the 3,000 incidence survivors fails: most already have a root in
the child core, and every core-rootless remainder acquires a norm-two witness
in a nonzero graph-glue coset. The hidden successful line has exactly the same
root survival/death signature as all 3,000 rejects.

The marked target supplies an independent low-norm intersection fingerprint.
Counting unoriented parent shell lines that survive in the parent--target
intersection gives:

| terminal edge | norm 4 | norm 6 | norm 8 | v1 rejects sharing all three counts |
|---|---:|---:|---:|---:|
| H3 | 9 | 133 | 921 | 2/3,000 |
| Q80 | 3 | 74 | 490 | 0/3,000 |

These batched modular counts have real pruning power, unlike the empty
survivor-equality predicate. They are necessary target-isometry data, not a
sufficient abstract-ADE classifier: the two H3 count matches are still
rootful after completion.

## Stronger pre-materialization constraint

Let `K` be the parent core, `K'` the marked target core inside `K tensor QQ`,
and suppose they are good `p`-neighbours. Then

```text
I = K intersect K',
K'/I = Z/p,
p*K' subset K.
```

For any marked target basis row `v` outside `K`, the nonzero residue

```text
p*v mod p*K
```

spans the unique neighbour line. All nonintegral target basis rows yield the
same projective residue. The planner now derives that line before constructing
a child, checks its old-root incidences and norm-4/6/8 parent--target overlap,
then materializes one child and independently verifies the rootless graph
completion. The terminal H3 and Q80 edges each use one proposal and one
materialization.

This closes historical recovery without pretending that the marked target
solves the harder problem posed by v1: a marked target core in the parent
rational space is nearly equivalent to the neighbour line itself.  Theorem
H0l.2 now supplies the previously missing target-free logical predicate.  It
eliminates the affine variable by projecting scaled dual-lattice shells to
explicit birth points on the isotropic quadric; rootlessness is the complement
of those points and the old-root hyperplane sections.  The exact construction
and exhaustive ternary controls are recorded in
[`INVERSE_ADE_PROJECTIVE_BIRTH_STRATA_2026-09-03.md`](INVERSE_ADE_PROJECTIVE_BIRTH_STRATA_2026-09-03.md).

This is not yet a rank-15 speed result.  A full scaled shell can have the same
asymptotic size as the quadric, so a practical planner must choose among shell
expansion, orbit compression, and lazy affine CVP.  The v2 materialization
counts above therefore remain historical benchmark data, not measurements of
the new compiler.

## What is exact and what remains experimental

The following are exact in v2:

1. the v1 replay and the split of all 3,000 terminal rejects;
2. the zero rank of the terminal survivor-equality systems;
3. recovery of the projective neighbour line from the marked target quotient;
4. the norm-4/6/8 parent--target overlap counts;
5. the independently materialized child core and completed root metric;
6. recovery of every stored H3, NS0024, and Q80 historical transition.

The bounded sparse-plus-dense generator used away from marked terminal edges
is not a complete quadric enumeration. The v2 result proves neither a
complexity improvement nor a general solver from abstract ADE data. These are
rank-15 core Kneser moves, not elliptic-neighbour pencils, rational maps, or
arithmetic descent data.

## Foundry-wide decision

Foundry-wide deployment remains fail-closed. All 936 bulk source--target route
rows lack the original bridge/glue and survival/birth data, and now also lack
the marked target core in the source rational space that closes the terminal
benchmark. A separately curated H3 `A1/MW16 -> NS0001-F001` positive-control
route is planner-ready and has an exact characteristic-zero source equation;
it is not one of the 936 bulk rows.

| item | count |
|---|---:|
| rootless target frames | 138 |
| target NS classes | 33 |
| source-target route pairs | 936 |
| pairs with both positive frame Gram matrices | 936 |
| bulk pairs with complete v2 planner inputs | 0 |
| separately curated planner-ready controls | 1 |
| curated source-equation attempts | 6 |
| certified characteristic-zero source equations after curated addition | 1 |

The recovery gate is therefore necessary but not permission to launch 936
searches. Run the single certified positive control, but keep the bulk foundry
closed. The no-birth logic no longer needs marked parent--target cores; the
remaining data task is to reconstruct compatible source markings and the
core/bridge/graph inputs that define the scaled shells.  The remaining
algorithmic task is to compress those projected shells at rank 15.  The `10x`
benchmark comes only after both boundaries are fixed.

The standard Kneser line parametrization and asymptotic neighbour statistics
are from Gaetan Chenevier,
[*Statistics for Kneser p-neighbors*](https://arxiv.org/abs/2104.06846).
Visible-root filtering in the unimodular setting is developed further in
[*Unimodular Hunting*](https://arxiv.org/abs/2410.18788). The graph-glue
completion layer uses Nikulin's
[*Integral symmetric bilinear forms and some of their applications*](https://www.mathnet.ru/eng/im1677).
Those sources provide the standard neighbour and visible-root framework; the
target-free scaled-shell elimination used here is Theorem H0l.2.

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

The v1 miss record is preserved at
[`elkies-k3-inverse-ade-target-planner-benchmark-v1.json`](../artifacts/generated-results/elkies-k3-inverse-ade-target-planner-benchmark-v1.json).
The v2 recovery record and readiness audit are
[`elkies-k3-inverse-ade-target-planner-benchmark-v2.json`](../artifacts/generated-results/elkies-k3-inverse-ade-target-planner-benchmark-v2.json)
and
[`elkies-k3-inverse-ade-foundry-readiness-v2.json`](../artifacts/generated-results/elkies-k3-inverse-ade-foundry-readiness-v2.json).
No `MATH_STATUS.json` entry changes: this is an exact retrospective algorithm
benchmark and data-readiness audit, not a new existence or equation theorem.
