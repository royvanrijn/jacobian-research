# Marked-`U` realization planner (2026-09-03)

<!-- status-consumer: EC-K3-R17-NONCYCLIC-4A1-DIRECT-EQUATION f657620e07f8f3f0 -->

## Outcome

The elliptic-fibration realization search is now a separate program from the
core inverse-ADE planner:

- `plan_inverse_ade_targets.sage` proposes positive-definite core mutations;
- `plan_marked_u_realizations.sage` accepts only an explicit marked
  `(NS,U,W)` and searches for literal primitive `U'` embeddings in that same
  Neron--Severi lattice.

This separation is operational, not just terminology. A target Gram matrix or
a core-planner hit cannot enter the realization graph unless the input also
contains a unimodular source splitting. This is the fail-closed rule exposed
by NS0024.

## Planner contract

The source basis convention is

```text
U = <F,F+O>,                 Gram(U)=J=[[0,1],[1,0]],
NS = U orthogonal_sum W(-1).
```

An input provides the ambient `NS` Gram matrix, two rows for the ordered source
`U`, and an integral frame basis. The planner verifies the Gram identities,
orthogonality, positivity, and determinant-one full transport before searching.

The finite search box names the four physical intersections

```text
(d,s,t,z) = (F.F', F.O', O.F', O.O').
```

They are enumerated lexicographically. The corresponding splitting-basis cross
matrix and positive projection Gram are

```text
A   = [[d,   d+s],
       [d+t, d+s+t+z]],

G_A = A^t J A-J.
```

Representation-independent determinant, saturation-index, Smith-primary, and
`det(A)` parity constraints run before a source-frame shell is enumerated.
Every surviving ordered representation of `G_A` constructs a literal `U'`,
its ambient fibre and pseudo-zero, a unimodular orthogonal-frame transport, and
an exact ADE signature. A desired target may be either an integral frame Gram
or only a root rank/ADE type.

There are two representation modes:

- `exact` uses complete PARI norm shells and gives completeness in the declared
  four-dimensional box unless an explicit representation cap is set;
- `catalog` replays supplied literal representations and always records that it
  makes no global completeness claim.

Nefness and zero effectivity are not inferred from root rank. They pass only
through hash-pinned exact evidence whose JSON fields or literal classes bind to
the constructed candidate. A lattice hit without both gates is retained as
pending, never as a realized fibration.

Each retained result keeps the physical incidence coordinates, projection
norm/coordinate size, declared Riemann--Roch ambient size, and any supplied
equation coefficient/bit metrics as an equation-facing cost record.

A minimal input has this shape (matrix specifications may instead use a
whitespace-delimited `path` plus `sha256`):

```json
{
  "source": {
    "label": "rank-four syntax control",
    "ns_gram": [
      [0, 1, 0, 0], [1, 0, 0, 0],
      [0, 0, -12, -16], [0, 0, -16, -24]
    ],
    "u_basis_in_ns": [[1, 0, 0, 0], [0, 1, 0, 0]],
    "frame_basis_in_ns": [[0, 0, 1, 0], [0, 0, 0, 1]]
  },
  "target": {},
  "intersection_box": {
    "F_dot_F_prime": [2],
    "F_dot_O_prime": [1],
    "O_dot_F_prime": [1],
    "O_dot_O_prime": [0]
  },
  "prime_local_bridge_constraints": {
    "relative_rank": 2,
    "saturation_index": 1,
    "saturated_bridge_determinant": 32
  },
  "representations": {"mode": "exact"},
  "physical_gates": {"required": false}
}
```

The displayed rank-four source is only a syntax/control example; K3 use
supplies the full rank-19 `NS` and rank-17 frame. Evidence artifacts are repository-relative
and must name their exact schema, status, and SHA-256.

## Controls

The exact control certificate is
[`elkies-k3-marked-u-realization-planner-controls-v1.json`](../artifacts/generated-results/elkies-k3-marked-u-realization-planner-controls-v1.json).

| order | source | request | result |
|---:|---|---|---|
| 1 | published R17 | alternate-Q80 frame | ten degree-two copies; `norm12-orbit-11952` is cheapest by the retained equation metrics |
| 2 | published R17 | `4A1` plus the noncyclic bridge profile | the new primitive nef `4A1/MW13` marking with saturated bridge `Z/4+Z/8` |
| 3 | that new `4A1/MW13` | root rank zero | the reverse degree-two hop to the published rootless frame, with the published equation zero restored literally |
| 4 | first nonhistorical foundry rootless target (`NS0002-F007`) | attempt only after the marked controls | blocked before enumeration because its bulk route has no common marked `NS`, source `U`, or relative-`U` lift |

The order-three reverse request contains only `root_rank=0` and
`ADE=rootless`; no target frame Gram is supplied.  The resulting selection is
therefore a target-free rootless hop rather than a replay forced to the known
R17 Gram.  Its equation-level inverse is now completed in
[`R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md`](R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md).

The certificate also contains a small complete exact-shell control. Its
negative bridge determinant profile is rejected by the square-index constraint
with zero representations enumerated.

## Existing R17 end-to-end milestone replay

A second invocation supplies no target Gram. It asks for a rootless rank-17
frame distinct from the source inside

```text
B = {
  old fibre degree <= 2,
  old-zero degree <= 1,
  new-zero/old-fibre degree <= 1,
  equation coefficient L1 <= 13,
  coordinate input bits <= 3482
}.
```

Within the exact 43-member published-R17 genus-one catalog, the planner selects
`norm12-orbit-11952`. The already independent downstream artifacts then check
the characteristic-zero `24 I1` equation, a saturated rational section basis
of rank 17, and arithmetic generic rank 17 over `QQ`.

This closes the requested end-to-end *control* using an existing certified
route. It does not claim that arbitrary new planner hits can yet be compiled
automatically, and it does not make the 936 bulk foundry rows planner-ready.

## Replay

```bash
sage -python \
  elkies-k3/scripts/certify_marked_u_realization_planner_controls.sage

sage -python \
  elkies-k3/scripts/certify_marked_u_realization_planner_controls.sage --check
```

For a standalone configuration:

```bash
sage -python elkies-k3/scripts/plan_marked_u_realizations.sage \
  path/to/marked-u-input.json \
  --output artifacts/generated-results/my-marked-u-plan.json
```

## Boundary

The relative-`U` theorem and the older exact source certificates remain the
mathematical authority. This planner is the reproducible realization and
routing layer around them. No theorem status is upgraded from a catalog hit,
and no core graph edge is silently treated as an elliptic-fibration edge.
