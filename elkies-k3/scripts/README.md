# Elkies K3 script map

This directory is the executable history of the elliptic-neighbour reconstruction.
It contains current proof replays, reusable compiler code, active searches, regression
fixtures, and some historical attacks that remain in the root because other notes or
launchers still refer to them.

This audit was made against repository commit
`4eac04a442a132b696a85384f08b81569870a940`. Every file in `scripts/` and
`scripts/archive/` was enumerated. Current and ambiguous entry points were inspected
against their code headers, outputs, certificates, and the notes that consume them;
the archive is classified by historical programme rather than promoted file-by-file.

## Status vocabulary

| status | meaning |
|---|---|
| **ACTIVE_PROOF** | Replays a claim used by a current note or pinned certificate. |
| **ACTIVE_COMPILER** | Reusable exact equation/Riemann--Roch/Jacobian infrastructure. |
| **ACTIVE_SEARCH** | Current bounded exploration. Its output is not a theorem until separately certified. |
| **REGRESSION** | Reproduces a specialization or previously solved compiler case. |
| **HISTORICAL_DIAGNOSTIC** | Preserves a failed or superseded route and the evidence that rejected it. |
| **ARCHIVED_SNAPSHOT** | Historical source copy. Never use it as the current proof entry point. |

Location alone is not authority. A root-level script is authoritative only when a current
note or pinned certificate names it. Conversely, a failed script is worth retaining when
it records a bounded negative result, a normalization bug, or a useful local model.

## Current proof and compiler entry points

### Generic neighbour and compiler infrastructure

- `exact_neighbor_examples.sage`
- `run_exact_neighbor_engine.sage`
- `verify_exact_neighbor_engine.sage`
- `elliptic_neighbor_compiler_field_generic.sage`
- `verify_elliptic_neighbor_compiler_field_generic.sage`

These are reusable infrastructure rather than a claim that one particular route is best.

### H3 source-family recovery

The source is the H3 `E7+E8/MW2` family on the level-474 `H21 cap H92` component.
The following scripts are the current source-side proof chain:

- `classify_kumar_e7e8_anchors.sage`
- `deconstruct_x0679_quotients.sage`
- `verify_h21_h92_level474_branch.sage`
- `factor_h21_h92_level474_modp.sage`
- `reconstruct_h21_h92_level474_qq.sage`
- `export_h3_level474_source_family.sage`
- `prove_h3_level474_rational_points.m`
- `verify_h3_noncm_q6_source_anchor.sage`
- `verify_h21_q6_section_descent.sage`
- `lift_h21_p1_modular.sage`
- `verify_h92_section_descent.sage`

### Exact H3 q6 and q8 equation route

The characteristic-zero equation route is exact through `D13/MW4`:

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4.
```

Authoritative replays:

- `derive_h92_q6_child_q8_marking.sage`
- `derive_h92_q6_child_q8_physical_root_target.sage`
- `derive_h92_q6_child_q8_corrected2cover_qq.sage`
- `certify_h92_q6_child_jacobian.sage`

The corrected q8 route deliberately supersedes the degree-46, `true1600`, and
`corrected1278` experiments in `archive/`. See
[`../H3_Q8_REAUDIT_2026-08-22.md`](../H3_Q8_REAUDIT_2026-08-22.md).

### H3 selected degree-two corridor

The following scripts certify one complete lattice/chamber corridor from H3 to the
rootless rank-17 frame:

- `analyze_h3_first_q6_chamber.sage`
- `classify_h3_q6_child_q8_orbits.sage`
- `analyze_h3_d13_q4_chamber.sage`
- `analyze_h3_rank_growing_degree2_chain.sage`
- `analyze_h3_mw10_to_rootless_chambers.sage`
- `verify_h3_d13_to_mw17_path.sage`
- `verify_rank17_to_h3_reverse_transport.sage`

The corridor is

```text
H3 E7+E8/MW2
 --q6 --> E8+E6/MW3
 --q8 --> D13/MW4
 --q24--> D12/MW5
 --q6 --> A11/MW6
 --q8 --> 2A5/MW7
 --q4 --> 3A3/MW8
 --q4 --> A3+2A2/MW10
 --q4 --> 5A1/MW12
 --q4 --> 4A1/MW13
 --q4 --> 3A1/MW14
 --q4 --> 2A1/MW15
 --q4 --> A1/MW16
 --q6 --> rootless/MW17.
```

This is a **selected certified corridor**, not a shortest-path or equation-cost
optimality theorem. At D13, all stated proper presentations through `q=23` were checked
and the first rank growth occurs at `q=24`; three q24 orbits lead to `D12/MW5`, and
orbit 85 was selected. The later continuation follows deterministic first hits from
that one root-adapted frame. Lateral moves, larger-q exits, the other q24 children, and
alternative multi-step corridors can still be easier at equation level.

The last checker supplies the previously missing positive-frame isometry from
the corridor's rootless endpoint to pinned `rank17_gram.txt`, inverts the full
H3-to-R17 transport, and exports every stage basis in both H3 and pinned-R17
coordinates. It also retains the exact bridge between the dominant D13
lattice marking and the distinct component-nef D13 equation marking; that
bridge changes the embedded `U`.

### Current H3 q24 equation frontier

- `close_h92_q8_q24_by_q6_translation.sage` — exact marking/NS bridge.
- `recover_h92_q24_exact_by_qq_trace_interpolation.sage` — exact interpolation
  machinery for the transported q24 horizontal.

These are **ACTIVE_SEARCH** until an exact q24 Riemann--Roch pencil, quartic/Jacobian,
and `D12/MW5` fibre classification are pinned. They should not yet be cited as a
completed characteristic-zero D12 child.

### Q80 compiler and regression route

The generic Q80 lattice corridor and the exact CM24 characteristic-zero shadow are
separate claims. Current entry points include:

- `verify_q80_to_rootless_path.sage`
- `trace_q80_candidate1_marked_transport.sage`
- `recover_q80_final_q6_via_basis_sections.sage`
- `certify_q80_final_q6_char0_rr_from_basis.sage`
- `compile_q80_final_q6_char0_child.sage`
- `verify_q80_lowq_cm24_equations.sage`

The CM24 child `4A2+A3+A5/MW2` is a regression/specialization shadow; it is not the
generic rootless `MW17` equation.

### Backtracks and comparison fibrations

These remain useful and reproducible, but they are not the H3 source construction:

- `verify_e6_neighbor_chain.sage`
- `verify_rank17_h8_split.sage`
- `analyze_rank17_h8_q9_fibers.sage`
- `verify_humbert8_d9e7_two_neighbor.sage`
- `verify_picard20_mw1_path.sage`
- `recover_mw1_*_glue.sage`

Use the names **Low-q MW2 Backtrack**, **E6 Backtrack**, **H2 Symmetry Comparison**,
and **H2 Minimal-MW Comparison** from
[`../CONSTRUCTION_ROUTES.md`](../CONSTRUCTION_ROUTES.md).

## Historical root-level launchers retained intentionally

The old E6/MW3 attack was run through scripts such as

- `start_e6_attack.sh`
- `run_e6_mw3_probe.py`
- `run_mw3_local_probe.py`
- `run_iv_mod_p2.py`

They are retained because they preserve exact bounded negative work and are referenced by
historical notes. They are **HISTORICAL_DIAGNOSTIC**, not source-construction entry
points. The guessed split E6 chart was never obtained by transporting the certified
neighbour chain, and its old survivors fail the missing off-diagonal height gate.

## Archive policy

[`archive/README.md`](archive/README.md) classifies the imported snapshots by failed or
superseded programme. Rules:

1. An archived script is never the authoritative proof source.
2. If the same filename exists in the root, use the root copy; the archived copy is an
   older snapshot unless a note explicitly says otherwise.
3. Keep unique failed scripts when they explain why an approach was abandoned.
4. Delete only demonstrated byte-identical copies or accidental empty files.
5. New failed experiments should be archived together with a short outcome note, not
   dropped into an unlabelled filename pile.

The 2026-08-23 audit removed five byte-identical ` (1)` copies and made no broad moves.
Moving historical root scripts without repairing every old command would make the record
less reproducible, not cleaner.

## Where to read the history

- [`../SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](../SCRIPT_ROUTE_AND_FAILURE_LEDGER.md) —
  route chronology, selection reasons, rejected assumptions, and remaining gaps.
- [`../CONSTRUCTION_ROUTES.md`](../CONSTRUCTION_ROUTES.md) — named geometric routes.
- [`../H3_Q8_REAUDIT_2026-08-22.md`](../H3_Q8_REAUDIT_2026-08-22.md) — exact q8 bug
  diagnosis and repair.
- [`../Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](../Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md)
  — Q80 terminal marking/RR closeout.

## Rule for future additions

Every new script should state near its header:

```text
status: ACTIVE_PROOF | ACTIVE_COMPILER | ACTIVE_SEARCH | REGRESSION | HISTORICAL_DIAGNOSTIC
claim: the exact claim or bounded search it supports
inputs: pinned files/certificates
outputs: generated artifact path
supersedes/superseded-by: optional script or note
```

A successful search result becomes a proof entry point only after an independent replay
checks its exact divisor, chamber/nefness, equation identity, and claimed fibre or
Mordell--Weil data.
