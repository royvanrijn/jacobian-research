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

<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-MOD131-HORIZONTAL 2249c509c1217d7c -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-HORIZONTAL cd7f46981d77f99b -->

## Local Sage 10.9 installation

On this machine the pinned SageMath 10.9 build is installed from conda-forge at

```text
/home/royvanrijn/.local/share/jacobian-sage-10.9
```

It is isolated from the repository and system packages.  The launcher in this
conda build does not accept the historical `sage -python` spelling; invoke Sage
Python scripts with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python script.sage
```

or prepend its `bin` directory to `PATH` and use that environment's `python`.

## Status vocabulary

| status | meaning |
|---|---|
| **ACTIVE_PROOF** | Replays a claim used by a current note or pinned certificate. |
| **ACTIVE_COMPILER** | Reusable exact equation/Riemann--Roch/Jacobian infrastructure. |
| **ACTIVE_SEARCH** | Current bounded exploration. Its output is not a theorem until separately certified. |
| **UNPROMOTED_RESULT** | An exact-looking output exists, but the canonical note, locks, or `MATH_STATUS.json` have not admitted it yet. |
| **REGRESSION** | Reproduces a specialization or previously solved compiler case. |
| **HISTORICAL_DIAGNOSTIC** | Preserves a failed or superseded route and the evidence that rejected it. |
| **ARCHIVED_SNAPSHOT** | Historical source copy. Never use it as the current proof entry point. |

Location alone is not authority. A root-level script is authoritative only when a current
note or pinned certificate names it. Conversely, a failed script is worth retaining when
it records a bounded negative result, a normalization bug, or a useful local model.

## Quick decision guide

- For the story and lessons, read `../ELKIES_K3_PROCESS_ATLAS.md`.
- For pinned commands and hashes, use `success-path/ledger.json`.
- Equation proofs currently stop at the C8-pointed `2A3+2A1/MW9` child; the
  marked lattice route reaches R17.
- `lift_h92_q24_orbit42_resolved_rr_qq.sage` and
  `certify_h92_q24_orbit42_a11_equation_marking.sage` are the exact equation
  and marking proofs for that child.
- Work next on q8/orbit376 from the exact C8-pointed `2A3+2A1/MW9` child,
  followed preferably by q12/orbit5867, with q12/orbit4484 retained as the
  certified fallback. The A11 q8/orbit12, physical q4/orbit208,
  q4/orbit1584, and q4/orbit164 equation edges are closed. Do not restart the
  withdrawn q6/orbit1307, q323, changed-zero, zero-pole, point-transport, or
  halving searches unless an exact gate on the selected suffix requires it.

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

The former q4/orbit230--q6/orbit1315 and q6/orbit1307 promotions are withdrawn
as equation-cost targets.  Exact physical-chamber audits show respectively a
pseudo-zero after the q4 return and a non-section component-10 zero after the
q6 Weyl repair.  The relevant route and correction scripts are:

- `search_h92_a5a5_zero_changing_loops.sage` — **ACTIVE_SEARCH** for the
  exhaustive 283-first-edge/558-zero q4/q6 loop ranking;
- `certify_h92_a5a5_q6o1307_promoted_route.sage` — **ACTIVE_PROOF** for the
  earlier abstract q6/orbit1307 splice, retained as an exact historical
  comparator only;
- `audit_h92_a5a5_q6o1307_physical_nef.sage` — exact physical-I6 Weyl repair;
  it preserves P1307, rejects component 10 as a zero, and exposes physical
  degree-one components 3, 5, and 9 with estimated RR profile `9 -> 3 -> 2`;
- `export_h92_a5a5_physical_source_marking.sage` — exports the exact
  component-9-zero `2A5` equation frame with all suffix/reverse target fibres;
- `export_h92_a5a5_zero_loop_returned_marking.sage` — exports any exact
  changed-zero `2A5/MW7` marking with its equation-A11 transport;
- `search_h92_d13_zero_changing_d12_presentations.sage --mode a5` — evaluates
  second zero changes using inherited-explicit curve degrees and per-edge
  horizontal-cost floors;
- `certify_h92_a5a5_q6o3372_q6o2052_promoted_route.sage --variant q230` —
  **ACTIVE_PROOF** for q4/orbit230, q4 return, q6/orbit1315, q4 return, q4
  exit, the exact current-`3A3` landing, and pinned R17 endpoint;
- `build_h92_a11_route_optimization_handoff.py` — the compact machine handoff
  generator consumed by
  [`../A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md`](../A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md).

The historical 4,199 and 10,334 totals are not equation-realizable as stated,
and the old q104/13,518 comparator is also non-nef in the physical I6 chamber.
`certify_h92_a5a5_direct_physical_q10.sage` reduces it exactly to q10 and
certifies the canonical current-3A3 landing.
`certify_h92_a5a5_direct_physical_q10_promoted_route.sage` composes that landing
to pinned R17 and promotes the 4,471 operational target (5,071 with the older
affine-omitting curve convention).  This is the current lifting target while
lower-score lateral searches continue.

The promoted replay target now begins one stage earlier:

- `search_h92_d13_zero_changing_d12_presentations.sage --mode d13` — exact
  search over the compact q4/q6/q8 equation-D13 frontier, all degree-one old
  D13 component zeros, exact return fibres, and the fixed current-D12 exit;
- `certify_h92_d13_q4o11_promoted_route.sage` — **ACTIVE_PROOF** for the
  q4/orbit11, q4-return, q24-exit splice, exact current-D12 landing, and full
  pinned-R17 endpoint;
- `search_h92_d13_zero_changing_d12_presentations.sage --mode a11` — the
  companion first-q8 audit.  It scans all 1,119 declared-nef candidates and
  ranks inherited-explicit costs; no presentation beats orbit12. Twenty-four
  nonprimitive-root re-zeroings are handled by saturated unimodular frames
  retaining the embedded simple-root lattice.
- `export_h92_first_q8_source_marking.sage` and
  `export_h92_first_q8_zero_loop_returned_marking.sage` — export the exact
  equation-explicit E8+E6 source and its q4/orbit11 changed-zero return state;
- `certify_h92_first_q8_q4o11_promoted_route.sage` — **ACTIVE_PROOF** for the
  q4/orbit11, q4-return, q4-exit replacement of the first q8, its exact
  equation-D13 landing, and the full pinned-R17 endpoint.

The D13 splice scores 25,323 against the measured-RR-calibrated direct q24
score 27,885. It rejoins the exact current D12 basis, keeps A11 q8/orbit12,
then uses the promoted q4/orbit230 and q6/orbit1315 double-zero suffix.

The earlier first q8 now has a cheaper exact replay: q4/orbit11, q4 return,
q4 exit scores 3,961 against 5,802 for direct q8. A complete second q4/q6/q8
zero-loop layer from the returned E8+E6 marking has 38 exact presentations and
no winner. Widened q10 degree-two, q6/q9/q12 degree-three, and
q8/q12/q16 degree-four searches likewise have no winner. The promoted full H3 route is recorded in
`elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json`.
The exact changed-zero D13 landing is exported separately and its 73-presentation
D12 crossover scan rejects carrying that zero into the canonical D13 splice.

### H3 q24 closeout and orbit42 frontier

The direct q24 equation edge is closed.  For a new replay, prefer the promoted
D13 q4/orbit11 zero-changing presentation above. Do not reopen
q32, native q24 suffixes, generic q24 section Hensel lifting, or easy orbit42
zero-pole searches as the active path.

Pinned q24 closeout scripts:

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

- `success-path/` — **ACTIVE_ROUTE_LEDGER** containing the version-locked
  launcher/index for every successful equation-lift stage, the exact artifact
  statuses and hashes, the complete pending route to pinned R17, and the
  separately labelled shortcut audits.  Canonical implementations stay in
  this directory and are not copied into the ledger directory.
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
- `preflight_h92_q24_orbit42_component_valuation_qq.sage` and
  `map_h92_q24_orbit42_i8star_physical_components_qq.sage` —
  **ACTIVE_EXACT_PREREQUISITES** for the resolved-RR compiler.  They pin the
  corrected divisor profile and the two exact physical spinor orientations.
- `recover_h92_q24_orbit42_zero_pole_smallprime.sage`,
  `scan_h92_q24_orbit42_zero_pole_model_modp.sage`, and
  `lift_h92_q24_orbit42_zero_pole_sections_qq.sage` — **EXACT_BOUNDARY_AUDIT**.
  Together they reconstruct eighteen exact rational identity-class zero-pole
  sections.  The modular scan is only a seed-selection aid; the final QQ
  identities are exact.
- `lift_h92_q24_orbit42_spinor_zero_pole_sections_qq.sage` —
  **EXACT_CONSTRUCTION_AID** for the remaining opposite spinor-class pair.
  It solves the cubic-`x`, constant-`y` cancellation branch directly over QQ
  and pins its mod-100003 residue.  Together with the preceding audit this
  represents the full twenty-point zero-pole shell exactly, but does not
  construct the resolved RR kernel or the A11 child.
- `construct_h92_q24_orbit42_exact_section_candidates_qq.sage` —
  **EXACT_POINT_AID_WITH_MODULAR_MARKING_BOUNDARY**.  Exact QQ(u) group law
  combines the 18+2 shell into four projective `(9,9,3)` points with exact
  Weierstrass identities.  Their four-way orbit42 identification uses the
  pinned mod-100003 shell isometry; resolved RR must still select the physical
  orientation.
- `lift_h92_q24_orbit42_resolved_rr_qq.sage` — **ACTIVE_PROOF**. It certifies
  the exact weighted resolved-RR two-plane, quartic, minimized Jacobian and
  A11 fibre classification.
- `certify_h92_q24_orbit42_a11_equation_marking.sage` —
  **ACTIVE_PROOF_WITH_GOOD_REDUCTION_MARKING_BOUNDARY**. The exact
  identity-shell degree fingerprint at p=100003 selects orbit64/mapping7 in
  the C10 orientation, with orbit65/mapping6 as spinor conjugate.
- `lift_h92_q24_a11_q8_residual_resolved_hensel.sage` and
  `derive_h92_q24_a11_q8_difference_qq.sage` — **EXACT_CONSTRUCTION**.  The
  first lifts the regular six-variable component-3 residual chart; the second
  uses fraction-free group law to derive the selected `(16,24,6)` q8
  horizontal without the former 36-variable lift.
- `lift_h92_q24_a11_q8_resolved_rr_qq.sage` — **ACTIVE_PROOF**.  It certifies
  the exact 14-to-2 recurrence RR plane, quartic, globally minimal
  `2I6+12I1` Jacobian, Euler number 24, and `2A5/MW7` classification without
  a Groebner basis.
- `certify_h92_q24_a11_q8_equation_marking_qq.sage` — **ACTIVE_PROOF**.  The
  nodal-cubic sign identity selects `old_A11_component_9` as the exact child
  zero, transports the affine component as a `(4,6,0)` section, attaches the
  two physical A5 chains, and verifies determinant-`-1` NS transports in both
  directions.
- `certify_h92_q4o208_physical_q4o1584_rr_qq.sage` and
  `certify_h92_q4o1584_equation_marking_qq.sage` — **ACTIVE_PROOF**. They
  certify the exact q4/orbit1584 RR plane and `D4+A3+3A1/MW7` Jacobian, then
  point the second-I6-affine zero and the old-`C0` branch required next.
- `certify_h92_q4o1584_physical_q4o164_rr_qq.sage` and
  `certify_h92_q4o164_c8_equation_marking_qq.sage` — **ACTIVE_PROOF**. They
  recover the q4/orbit164 pencil by rational branch-value interpolation,
  certify the `2A3+2A1/MW9` Jacobian, and identify old `C8` as its exact zero.
- `recover_h92_q4o164_odd_i4_sections_modp.sage` and
  `lift_h92_q4o164_odd_i4_sections_qq.sage` — **EXACT/MODULAR FRONTIER
  DIAGNOSTIC**. They enumerate odd finite-`I4` polynomial sections without a
  Groebner basis, Newton-lift the regular branches, and verify every retained
  `QQ(t)` section by literal substitution. The optional bounded signed-sum
  audit replays with
  `--precision 180 --rational-sum-scan --output artifacts/generated-results/elkies-k3-h3-q4o164-odd-i4-rational-sum-scan-p41.json`;
  it filters 364 words at three fibres, interpolates 22 survivors in the
  marked q8 degree window, and finds no exact rational sum. This is a bounded
  construction audit, not a non-existence theorem for the missing ninth
  direction.
- `recover_h92_q4o164_zero_node_sections_modp.sage` and
  `lift_h92_q4o164_zero_node_sections_qq.sage` — **EXACT/MODULAR FRONTIER
  DIAGNOSTIC**. The exhaustive degree-four zero-node scan at `p=23` has one
  inverse pair; both lift over `QQ`, but specialization quotient tests keep
  them in the known rank-eight shell.
- `probe_h92_q4o164_inherited_p1_abel_trace_modp.sage` — **MODULAR q8
  CONSTRUCTION FRONTIER**. It transports the exact inherited `P1` curve
  through the three certified q4 models using their degree-one pointed
  quartic maps, obtaining degrees `3 -> 6 -> 7`. Fibrewise Abel reduction uses
  only the unique `7 x 8` kernel in `L(8O)`. At `p=131`, `--interpolate`
  reconstructs the full trace section with degrees `x=(32,28)` and
  `y=(48,42)` from 122 good fibres, with 61/31 holdout fibres and an exact
  modular Weierstrass identity. Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/probe_h92_q4o164_inherited_p1_abel_trace_modp.sage \
    --interpolate \
    --output artifacts/local/elkies-k3/q4o164-inherited-p1-abel-trace-section-mod131.json
  ```

  The trace is not the q8/orbit376 horizontal modulo the already named
  degree-one sections: their marked tails leave the exact correction
  `(0,0,0,0,-1,-1,2,0,0)`. This remains a modular construction, not an
  equation certificate.
- `audit_h92_q4o164_integral_basis_height_gram.sage` — **EXACT HEIGHT/MARKING
  AUDIT**. Raw node incidence had mislabeled the infinity-`I4` component of
  `B7`. Multiplying every basis section and pair sum by four clears all
  component groups; compact pole growth then gives the corrected height Gram
  with determinant `459/8`, rather than relying on unresolved node labels.
  A finite positive-definite enumeration finds 16 integral embeddings in the
  C8-pointed marked MW9 lattice, eight compatible with the valid
  `B0,...,B6` profiles. All eight contain the q8 residual but with different
  basis words. Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/audit_h92_q4o164_integral_basis_height_gram.sage
  ```

  The old coarse profile `[1,0,0,2]` for `B7` and the resulting height `3`
  for `N=2*B0+B5+B7` are withdrawn; exact pole growth gives `h(N)=13/4` and
  an odd infinity-`I4` label.
- `identify_h92_q4o164_q8_horizontal_mod131.sage` — **EXACT MODULAR q8
  HORIZONTAL**. It combines each of those eight embeddings with the complete
  inherited-`P1` Abel trace and the exact saturated relation
  `3*C8opp=-2*B0-3*B1-4*B2+3*B3-2*B4+2*B5-B6-2*B7`. Exactly one candidate
  has the certified q8 pole profile `x=(12,8)`, `y=(18,12)`, `P.O=4`, giving

  ```text
  H=T-C8opp-B0+2*B1+B2-3*B3-B4-2*B5+B7.
  ```

  Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/identify_h92_q4o164_q8_horizontal_mod131.sage
  ```

  Pass `--prime p` to replay any exact good-prime Abel-trace artifact.  For
  larger CRT primes, the trace probe accepts `--good-fibre-limit 100`; this
  retains 91 interpolation fibres and nine independent holdouts without
  scanning all of `GF(p)`.
- `reconstruct_h92_q4o164_q8_horizontal_crt_qq.sage` — **EXACT QQ q8
  HORIZONTAL**. It combines the monic-normalized compact coordinates of `H`
  over 22 good primes. Independent coefficient reconstruction exposes the
  large common scales; simultaneous projective LLL on the 22- and 32-entry
  coordinate vectors recovers primitive vectors of at most 363 and 526 bits
  from the 566-bit CRT modulus. It accepts the result only after exact
  substitution in the compact Weierstrass equation and reduction back to
  every input prime. This targets the much smaller `(12,8)/(18,12)` q8
  section directly rather than first reconstructing the `(32,28)/(47,42)`
  Abel trace. Replay with:

  ```bash
  ~/.local/share/jacobian-sage-10.9/bin/python \
    elkies-k3/scripts/reconstruct_h92_q4o164_q8_horizontal_crt_qq.sage
  ```

  Terminal status is `PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT`. The exact
  fourfold pole degree is 172, giving canonical height 11.  The horizontal is
  not yet the q8 child equation; the remaining gate is the resolved q8 RR
  calculation.
- `build_h92_q4o164_all_zero_q8_known_horizontal_audit.py` — **BOUNDED
  NEGATIVE ROUTING AUDIT**. For each of the eight certified effective q4/o164
  origins, a q8/degree-two shell was generated with
  `search_root_adapted_weyl_neighbors.sage --q 8 --degree 2
  --mw-vector-cap 10000` and filtered with
  `score_h92_marked_frontier_equation_cost.sage
  --require-known-horizontal`. The aggregate covers 119,220 primitive
  candidates and finds zero horizontals in the corresponding explicit-section
  subgroup. Replay the compact aggregate after the individual shell/score
  files are present with:

  ```bash
  python3 elkies-k3/scripts/build_h92_q4o164_all_zero_q8_known_horizontal_audit.py
  ```

  This is exact for the stored capped shells, not a global non-existence
  theorem. Together with the 364-word odd-I4 rational-sum scan, it supports
  the inherited-`P1` Abel-reduction construction target above.
- `certify_h92_marked_route_to_pinned_r17.sage` — **ACTIVE_ROUTE_PROOF**. It
  certifies the q323-free marked q4/q4/q8/q12 suffix and terminal integral
  isometry to pinned R17. It is not an equation certificate for q8/orbit376 or
  q12/orbit5867 (or fallback q12/orbit4484) and does not supply the endpoint
  section/rank package.
- `transport_h92_q24_a11_degree_one_shell_qq.sage` —
  **EXACT_CONSTRUCTION_AID**.  It transports the two degree-one identity
  curves and the degree-one spinor curve by Möbius inversion and exact binary
  quartic covariants, and verifies the pointed opposite relation over QQ.
- `audit_h92_q24_a11_missing_direction_alternatives.sage` —
  **EXACT_NEGATIVE_CONSTRUCTION_AUDIT**.  It proves that all three exact
  degree-one transports remain in the known fifth-coordinate-zero parent
  hyperplane, that q8 orbit2162 only reverses the required child coordinate,
  and that the smallest parent carrier is the D12 `P.O=4` vector
  `(0,0,0,0,1)`.
- `audit_h92_a11_quintic_bridge_zero_mismatch.sage` —
  **EXACT_REJECTION_CERTIFICATE**.  It replays the selected `R3`-zero frame
  and rejects the former quintic shortcut, which had mixed it with `A0`-zero
  coordinates.  The compatible degrees are 46 and 4, and the claimed word
  does not equal `M`.
- `audit_h92_a11_explicit_aj_carriers.sage` —
  **EXACT_ALTERNATIVE-CONSTRUCTION_AUDIT**.  It scans all stored explicit
  (-2)-classes in the selected marking.  No positive single carrier works;
  the first positive subset uses degree-40 and degree-44 traces, so it is a
  fallback rather than the active route.
- `score_h92_a11_equation_cost_neighbors.sage` and
  `certify_h92_a11_equation_cost_orbit849.sage` —
  **EXACT_ROUTE-PLANNING/AUDIT**.  The first scores all declared neighbours;
  the second applies the stronger nefness and marked-U gates.  Low-score
  orbits 849 and 591 fail nefness; passing lateral candidates do not yet have
  the required 2A5-to-pinned-R17 continuation.
- `probe_h92_q24_a11_close_p24_quintic_modp.sage` —
  **MODULAR_RESOLUTION_GATE**.  It composes both pointed quartics using only
  univariate arithmetic and rejects the naive unresolved restriction:
  q24 degree 14 and A11 degrees 39/41 replace the compatible
  strict-transform degree 46.  Both tangent signs vanish in the raw chord
  discriminant but are cancelled from the normalized quartic, locating the
  issue at the resolved base locus.
- `preflight_h92_q24_o12_p42_exact_q6_points.sage` and
  `run_h92_q24_orbit42_fast_parallel.py` — **EXACT_NEGATIVE_AUDIT** for the
  rejected q6-point transport.  The named equation-D13 coordinate conversion
  is exact, but `O12` and `P42` have q6 degrees `435` and `703`; the runner
  stops cleanly before the archived transport/orientation experiments.
- `analyze_h92_q24_orbit42_identity_halving.sage` and
  `recover_h92_q24_orbit42_by_identity_halving_qq.sage` —
  **SHORTCUT_BOUNDARY_AUDIT**.  The first proves the exact identity-shell
  doubling relation.  The second uses the exact QQ points but only a modular
  shell marking/halving census; its four rational degree-three candidates
  have degree-18 squarefree chord branches and no A11 fibre.  This is not a
  characteristic-zero non-existence theorem.
- `recover_h92_q24_exact_by_qq_trace_interpolation.sage` — **HISTORICAL_DIAGNOSTIC**
  unless it is refactored to serve the Theta/component construction directly.

Archived q24 dead ends:

- `archive/lift_h92_q24_direct_hensel.sage`
- `archive/exactify_h92_q24_from_padic_srr.sage`

These record why p-adic direct lifting/SRR was abandoned. Do not increase Hensel
precision or reconstruct the generic q24 section as the next step. The q24
`D12/MW5`, orbit42 `A11/MW6`, q8 orbit12 `2A5/MW7`, physical q4/orbit208
`3A3/MW8`, q4/orbit1584 `D4+A3+3A1/MW7`, and q4/orbit164
`2A3+2A1/MW9` children are closed. Continue from the C8-pointed child to
q8/orbit376.

### Q80 compiler and regression route

The generic Q80 lattice corridor and the exact CM24 characteristic-zero shadow are
separate claims. Current entry points include:

- `verify_q80_to_rootless_path.sage`
- `trace_q80_candidate1_marked_transport.sage`
- `score_h3_d12_q80_crossovers.sage` — exact transport of all eleven retained
  Q80 fibre classes into the current equation-side H3 D12 frame, with the
  negative crossover-cost audit documented in
  [`../H3_D12_Q80_CROSSOVER_AUDIT_2026-08-24.md`](../H3_D12_Q80_CROSSOVER_AUDIT_2026-08-24.md)
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

### q12/o4484 low-pole section compiler

The direct final q12 edge remains the certified lattice route, but its horizontal
section has an exact cheaper group-law compilation plan:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_h92_q12o4484_p0_section_word.sage \
  --marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --q12-cost artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap10000-rootless-equation-cost.json \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o4484-four-p0-section-word.json
```

This exhausts the marked norm-at-most-four shell, imposes the physical simple
and affine component gate, rules out words of length at most three modulo the
explicit subgroup, and selects four genuine `P.O=0` chamber sections with
q4/o164 parent degrees `3,2,3,2`. It is an exact MW/group-law compiler target;
the four characteristic-zero section equations are still to be constructed.

### Rootless-exit low-pole compiler frontier

The same construction is compared across both retained q12 rootless exits and
all four retained q20 exits by:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_h92_rootless_p0_section_word_frontier.sage \
  --marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --cost artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q12d2-cap10000-rootless-equation-cost.json \
  --cost artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-4a1-p1229-q20d2-cap10000-rootless-equation-cost.json \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-rootless-p0-section-word-frontier.json
```

The exact marked quotient search promotes q12/orbit5867 as an optional final
compiler target: its four `P.O=0` branches have q4/o164 parent degrees
`3,2,1,2`, versus total parent degree 10 for orbit4484. Optimizing the same
word by parent `a-b` gives values `2,2,1,1` (sum 6, maximum 2), versus minimum
sum 8 for orbit4484. A separate exact marked-route certificate pins the
orbit5867 rootless child to R17. Exact meet-in-the-middle searches also rule
out five- and six-section physical polynomial words with parent `a-b` total
at most 5, so the four-section total 6 is stable against those short detours.
The same artifact records a three-section
`P.O=(0,0,1)` fallback and the negative length-at-most-four results for all four
q20 exits. These remain compiler plans until the new sections are constructed
over characteristic zero.

The optional orbit5867 endpoint and its pinned-route composition replay with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_degree2_candidate.sage \
  --source-marking artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-marking.json \
  --source-frame artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-p1229-frame.txt \
  --fibre 6,2,-16,18,2,5,-11,2,13,11,-6,-41,6,5,-14,0,-3,-1,-2 \
  --candidate-label q12o5867_after_q8o376 --target pinned_R17 \
  --frame-output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-frame.txt \
  --output artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-certificate.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_route_to_pinned_r17.sage \
  --edge artifacts/generated-results/elkies-k3-h3-q4o208-physical-q4o1584-lateral-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-certificate.json \
  --edge artifacts/generated-results/elkies-k3-h3-q4o164-q8o376-q12o5867-rootless-certificate.json \
  --source q4o208_equation_physical_3A3_C5_zero \
  --status PASS_EXACT_Q323_FREE_Q4O1584_Q4O164_Q8O376_Q12O5867_PINNED_R17_ROUTE \
  --output artifacts/generated-results/elkies-k3-h3-q4o208-q323-free-q4o1584-q4o164-q8o376-q12o5867-pinned-r17-route-certificate.json
```

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
