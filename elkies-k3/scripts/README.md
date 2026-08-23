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

The selected route is fixed. Do not reopen q32, native q24 suffixes, generic q24
section Hensel lifting, or easy orbit42 zero-pole searches as the active path.

Active q24 corridor scripts:

- `certify_h92_q24_equation_d13_to_pinned_r17.sage` — **ACTIVE_PROOF** for the
  current-equation D13 marking, q24/orbit85 D12 marking, pinned-R17 suffix, and
  exact backward-manifest transition/profile lock. This is still a lattice/NS
  boundary, not an equation-level child.
- `replay_h92_q24_d12_discovery_and_theta.sage` — **ACTIVE_PROOF** for the exact
  D12 fibre structure `D24 = Theta + sum m_i C_i`; use this Theta/component
  construction to build the q24 RR pencil.
- `probe_h92_q24_d12_component_valuation_rr_modp.sage` — **REGRESSION** for the
  modular q24/orbit85 component-valuation RR replay. It records the 56-dimensional
  geometric ambient, smooth-collision rank 48, component-cover resolved rank 6,
  `h0=2`, quartic degree 4, and D12 fibre data.
- `extract_h92_q24_d12_modp_signature.sage` — **REGRESSION** for the modular
  q24/orbit85 resolved-RR signature: ambient 56, collision rank 48, resolved rank
  6, kernel 2, quartic degree 4, D12 fibre data. New outputs retain the embedded
  component compiler status and remain `CANDIDATE_*` unless an actual
  equation-level terminal PASS exists; do not treat the modular signature alone
  as a D12 equation certificate.
- `run_h92_q24_orbit85_modp_stack.py` — **ACTIVE_PROOF_AID** to replay the
  fixed q24/orbit85 modular construction stack for canonical split-centre
  `C10` primes. It runs the preflight, affine component graph, effective D13
  transport, component-valuation RR compiler, and compact D12 signature
  extractor in order; it is for CRT/rational-reconstruction data collection,
  not Hensel lifting or neighbor search.
- `reconstruct_h92_q24_orbit85_rr_kernel_crt.sage` — **ACTIVE_PROOF_AID** to
  CRT/rational-reconstruct the canonical q24/orbit85 `2 x 56` resolved-RR
  kernel, binary quartic, and Jacobian data from modular signatures with a
  held-out prime check. A passing output is a QQ kernel candidate, not yet the
  full exact D12 equation certificate.
- `reconstruct_h92_q24_orbit85_compact_rows_lll.sage` — **ACTIVE_PROOF_AID** to
  recover the same q24/orbit85 kernel projectively in compact intrinsic
  `(A, BZ, C)` coordinates using CRT plus LLL and a held-out prime. This is the
  fallback when scalar RREF rational reconstruction is denominator-heavy.
- `extend_h92_q24_orbit85_crt_precision.py` — **ACTIVE_PROOF_AID** to add fresh
  q24/orbit85 modular signatures beyond the archived small-prime pool. It
  builds the q24 direct bridge, I9* resolution, canonical `C10` stack, and then
  retries scalar CRT plus compact-row LLL after each accepted prime.
- `build_h92_q24_orbit85_exact_construction_manifest.sage` — **ACTIVE_PROOF** for
  the q24/orbit85 exact-construction contract. It binds the exact Theta divisor,
  component thresholds, modular 2D RR kernel, quartic/Jacobian regression, and
  D12 gates into
  `artifacts/local/elkies-k3/q24-orbit85-exact-construction-manifest.json`.

Orbit42/D12 profile rule:

- Do not use parity/minimum-`P.O` shortcuts for D12 correction classes. The
  R17-directed orbit42 class is the spinor class with correction `3`,
  `P.O=3`, denominator degree `3`, and no extra fibre twist. Live q24/D12
  profile helpers should use the exact D12 discriminant-class lookup table.
- `extract_h92_q24_orbit42_current_equation_bridge.sage` — **ACTIVE_PROOF** for
  the corrected orbit42 divisor in the current-equation D12 frame:
  `mw=(-1,0,-1,-1,0)`, height 7, correction 3, `P.O=3`, fibre twist 0.
- `recover_h92_q24_orbit42_current_equation_section_modp.sage`,
  `recover_h92_q24_d12_orbit42_section_modp.sage`,
  `archive/recover_h92_q24_pointed_zero_pole_sections.sage`, and
  `archive/recover_h92_q24_r17_a11_zero_pole_sections.sage` —
  **HISTORICAL_DIAGNOSTIC**.
  They prove the old easy zero-pole route does not generate the selected
  orbit42 target.
- `audit_h92_q24_orbit42_explicit_multisections.sage` — **REGRESSION** for the
  NS/Picard-level multisection span audit. The target is in the span, but the
  artifact does not compute the fibrewise Abel-Jacobi functions or the A11
  equation.
- `compile_h92_q24_orbit42_a11_chord_modp.sage` — **REGRESSION** only when fed an
  actual P42 section artifact. The preferred exact route is now the direct
  resolved-RR compiler for the degree-two divisor `D42`, not a zero-pole P42
  search.
- `recover_h92_q24_exact_by_qq_trace_interpolation.sage` — **HISTORICAL_DIAGNOSTIC**
  unless it is refactored to serve the Theta/component construction directly.

Archived q24 dead ends:

- `archive/lift_h92_q24_direct_hensel.sage`
- `archive/exactify_h92_q24_from_padic_srr.sage`

These record why p-adic direct lifting/SRR was abandoned. Do not increase Hensel
precision or reconstruct the generic q24 section as the next step. The remaining
frontier is the exact characteristic-zero resolved-RR compiler from the explicit
Theta/component construction to the q24/orbit85 `D12/MW5` child, followed by the
direct exact resolved-RR compiler for the current-equation orbit42 divisor
`D42 = O + P + V`.

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
