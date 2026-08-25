# Elkies K3 / neighbor programme

This directory contains the K3-surface and elliptic-neighbor work used to search for high-rank elliptic fibrations and explicit equation routes.

For the chronological account of rank exchange, cancellations, corrections,
equation lifts, and reusable methods, start with
[`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md). For the
general rank, determinant, fibration, specialization, and lift statements,
see [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
For the
human-readable map of how the reverse searches, H3 source construction, H2
comparison, Q80 route, and CM24 scaffold fit together, see
[`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md). The exact lossless reverse
transport from the pinned recovered 17-by-17 lattice to H3, including the q8
marking distinction, is documented in
[`RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md`](RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md).
For executable entry points and the research history behind
failed/superseded scripts, see [`scripts/README.md`](scripts/README.md) and
[`SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](SCRIPT_ROUTE_AND_FAILURE_LEDGER.md).

## One-minute status

- Think of one K3 surface with several coordinate systems. A neighbour changes
  the elliptic fibration, moving rank between fibre roots and sections.
- The generic budget is always `root rank + MW rank = 17`.
- Equations are certified through the component-9-zero `2A5/MW7` child; the marked
  lattice route continues exactly to the pinned rootless `MW17` endpoint.
- The q4/orbit230 equation and q4 return are exact, but the formerly promoted
  q6/orbit1315 continuation is withdrawn as an equation-cost target: its cheap
  score used a chamber pseudo-zero rather than the effective P230 section.
- The older q6/orbit1307 horizontal survives physical Weyl reduction, but its
  advertised component-10 zero does not: after correction component 10 has
  degree zero, while components 3, 5, and 9 have degree one.  Its 10,334-point
  continuation is therefore also withdrawn pending a physical-zero rerank;
  see
  [`../artifacts/local/elkies-k3/q24-2a5-q6o1307-physical-nef-audit.json`](../artifacts/local/elkies-k3/q24-2a5-q6o1307-physical-nef-audit.json).
- The active replacement is physical q4/orbit208: `P.O=0`, exact RR dimensions
  `4 -> 2 -> 2`, quartic degree 4, and a minimal `3I4+12I1` Jacobian.  Its full
  lattice route lands on canonical `3A3/MW8` and then pinned R17; only effective
  C5 equation pointing and the full old-curve marking remain.
- The separately certified first-q8 and D13 lattice detours remain available,
  but combined totals using the withdrawn 4,199 suffix are not valid compiler
  comparisons.  See
  [`A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md`](A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md).

Use `R17` for the recovered rootless endpoint and `H3 source` for the level-474
Kumar polarization. Low-q, E6, H2 and Q80/CM24 are comparison or regression
routes, not substitutes for the selected H3 marking. See
[`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md).

## Current priority

The **H3 source family** is the primary characteristic-zero starting point. The currently selected **H3 degree-two corridor** is one certified route from that source to the pinned recovered rootless/MW17 frame; the endpoint identification and complete inverse NS transport are now exact. The route is not proved shortest, globally optimal, or cheapest to compile. Its first four neighbours are exact at equation level:

```text
H3 E7+E8/MW2
 --q6--> E8+E6/MW3
 --q8--> D13/MW4
 --q24 orbit85--> D12/MW5
 --q6 orbit42--> A11/MW6
 --q8 orbit12--> 2A5/MW7.
```

The q8 repair found two independent bugs in the previous child-side compiler: a binary-quartic 2-cover multiplier was applied twice, and the q-frame CRT normalizer omitted the `Dx` factor of `x(S)=Nx/Dx`. After both corrections the q8 Riemann--Roch problem collapses to an exact `13 -> 2` characteristic-zero intersection and its Jacobian has one `I9*` fibre plus nine `I1`, hence `D13/MW4` exactly.

The selected `D13/MW4 --q24 orbit85--> D12/MW5` equation arrow is now closed. Its exact q24 horizontal section is exported as `artifacts/local/elkies-k3/q8-q24-horizontal-section-qq.json` with status `PASS_EXACT_Q24_HORIZONTAL_SECTION`, and the resolved component-valuation RR lift is exported as `artifacts/local/elkies-k3/q24-d13-to-d12-component-valuation-qq.json` with status `PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR`. The RR profile is `ambient=56`, `collision=48`, `post=8`, `resolved=6`, `kernel=2`; the compiled quartic has degree four and its Jacobian has root rank `12`, root determinant `4`, Euler number `24`, and MW rank `5`.

At lattice/planning level, the candidate D13 detour is

```text
D13 --q4 orbit11--> D5+D9 --q4--> D13(new zero) --q24--> current D12.
```

The landing is the exact stored D12 full basis and the composite endpoint is
the pinned R17 full basis. Its historical planning score is 25,323 against
the calibrated direct score 27,885. Because the q4/orbit230 return exposed a
physical-chamber pseudo-zero, this changed-zero score is not treated as an
equation target until its own return equation identifies an effective zero.
A complete scan of all 1,119 declared-nef A11 q8
candidates, including all 24 torsion-root presentations via saturated adapted
frames, found no cheaper inherited-explicit operational-cost zero loop.

For a replay beginning one stage earlier, the corresponding exact lattice
detour is

```text
E8+E6 --q4 orbit11, zero=old_E8E6_component_1--> A2+D5+E7
      --q4--> E8+E6(new zero) --q4--> equation D13.
```

Its historical planning score is 3,961 versus 5,802 for the direct q8. The
D13 landing and pinned R17 endpoint
are identified by full integral unimodular bases, not ADE/MW labels. An exact
second-loop q4/q6/q8 scan finds no further improvement from the returned
E8+E6 marking. Widened q10 degree-two, q6/q9/q12 degree-three, and
q8/q12/q16 degree-four first-edge scans also have no winner. The machine certificate is
[`../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json).
The changed-zero D13 landing was also searched directly toward current D12;
its degree-two q60 target and all 73 exact q4/q6/q8 crossover presentations
are operationally more expensive than the canonical D13 continuation. As
with the D13 detour, its changed-zero compiler score awaits a physical return
chamber audit.

The following formerly promoted lattice splice is no longer an equation-cost
target. It was

```text
2A5 --q4 orbit230--> A1+A4+A5 --q4--> 2A5(new zero)
     --q6 orbit1315--> 3A2+A3 --q4--> 2A5(second new zero) --q4--> 3A3,
```

using `old_A11_component_10` and then `old_A5A5_component_1` as explicit
zeros, then resume the existing suffix. Its marked fibres and determinant-one
transports to stored `3A3` and pinned R17 remain exact. However, the 4,199
score is withdrawn: the stored return zero has `P.O=26` and differs from the
effective P230 section by vertical roots; q6/orbit1315 meets effective P230 in
degree 54. Thus this is an exact lattice path, not a compiler-cost promotion.

For orbit42, use `mw=(-1,0,-1,-1,0)`, height `7`, correction `3`,
`P.O=3`, fibre twist `0`, and `D42=O+P+V`. The physical I8* marking and
complete twenty-point zero-pole shell are exact. Four exact section candidates
are also known. The resolved-RR artifact certifies the A11 child, and the
shell-degree marking binds it to equation-side orbit64/mapping7 in the chosen
C10 orientation.  The next q8 lift is also exact: its resolved RR profile is
`ambient=14`, `collision rank=12`, `h0=2`; its quartic Jacobian has
`2I6+12I1`, and old A11 component 9 is the exact child zero. See
[`ORBIT42_EQUATION_LIFT.md`](ORBIT42_EQUATION_LIFT.md)
and the version-locked [`scripts/success-path/`](scripts/success-path/) ledger.

Useful negative lessons: the fast q6 transport has degrees `435/703`; the
rational-halving candidates have no A11 chord; and the old correction-one,
`P.O=2`, easy zero-pole and archived Hensel paths are not the selected lift.
These close shortcuts, not every possible construction.

Marking warning: the pinned dominant D13 representative used by the lattice corridor and the component-nef D13 representative used by the q8 equation compiler are not the same stored 17-by-17 frame. Their exact full-NS bridge changes the embedded `U`; use the reverse-transport ledger rather than substituting one frame for the other by ADE/MW label.

The **Q80 Low-q Compiler Route** is now a completed secondary compiler route through its terminal characteristic-zero specialization shadow. Its generic lattice corridor reaches rootless/MW17, while the exact CM24 characteristic-zero terminal equation is `4A2+A3+A5/MW2`. The former final-q6 MW-marking problem is closed; see [`Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md).

## H3 q8 exact certificate

Reproduce the corrected marking with

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_marking.sage
```

Expected geometry:

```text
MW coordinate = (-2,-2,0)
height = 24
O-intersection = 10
smooth collision degree = 10
II*, IV* = identity component
```

Then reproduce the complete second neighbour with

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --representative component-nef \
  --output artifacts/local/elkies-k3/q8-target-component-nef-audit.json
cmp artifacts/local/elkies-k3/q8-target-component-nef-audit.json \
  elkies-k3/data/fibrations/h3_q8_component_nef_physical_root_target.json
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_corrected2cover_qq.sage
```

Expected endpoint:

```text
Q8QQRR|ambient=13|rows=11|rank=11|kernel=2|...
Q8QQQUARTIC|degree=4
Q8QQCHILD|finite=[(1,[2,3,15],'I9*'),(9,[0,0,1],'I1')]|infinity=((0,0,0),'smooth')|root_rank=13|root_euler=24|root_det=4|MW_rank=4|status=PASS_EXACT_CORRECTED_Q8_D13_CHILD
```

The detailed diagnosis and supersession boundary is in [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md).

## Q80 status: generic route and exact CM24 characteristic-zero shadow

The generic Q80 alternate route is certified to a new rootless `MW17` frame:

```text
E6+D5+A3/MW3
 --q4--> D9+A4/MW4
 --q4--> D7+D5/MW5
 --q6--> D7+D4/MW6
 --q4--> A6+A4/MW7
 --q4--> A6+A3/MW8
 --q6--> A4+A2+A1/MW10
 --q4--> A3+A2/MW12
 --q4--> 4A1/MW13
 --q4--> A1/MW16
 --q6--> rootless/MW17.
```

Every retained new divisor from `D7+D5/MW5` onward has chamber-reduced old-fibre degree `2`.

The old suffix labels `7774`, `1938`, `6855`, and `candidate1` are search-result identifiers rather than intrinsic names. New notes should use the canonical stage names in [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md), while [`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md) records the exact historical search lineage and the successful terminal-marking reconstruction.

The former characteristic-zero obstruction at the `A1/MW16` parent is closed. The selected generic terminal divisor is replayed with exact unimodular neighbour transports, the exact final horizontal is reconstructed as a difference of easier high-incidence MW sections, and the historical CM24 `P2-P3` point is used only to identify the correct difference modulo `73`.

The exact final resolved Riemann--Roch pencil over `QQ(sqrt(-3))` has

```text
ambient = 4
A4 quotient rank = 1
A5 quotient rank = 1
condition rank = 2
kernel = 2
h0(D) = 2.
```

Its exact degree-four binary quartic compiles to a characteristic-zero CM24 child with

```text
finite fibres = 4 I3 + I4 + I6 + 2 I1
infinity = smooth
root lattice = 4A2+A3+A5
root rank = 16
root determinant = 1944
MW rank = 2
Euler number = 24.
```

The late specialized sequence is therefore exact in characteristic zero:

```text
2A6+3A1/MW3
 --q6--> A5+2A4+2A1/MW3
 --q4--> 2A4+2A3+A1/MW3
 --q4--> A1+2A3+2D4/MW3
 --q4--> A1+A2+A3+A4+A5/MW3
 --q6--> 4A2+A3+A5/MW2.
```

This does **not** identify the specialized child with the generic endpoint. The generic terminal child remains the independently certified rootless/MW17 frame; CM24 has higher Picard rank and acquires extra algebraic classes.

Reproduce the terminal reconstruction with

```bash
sage elkies-k3/scripts/trace_q80_candidate1_marked_transport.sage
sage elkies-k3/scripts/recover_q80_final_q6_via_basis_sections.sage
sage elkies-k3/scripts/certify_q80_final_q6_char0_rr_from_basis.sage
sage elkies-k3/scripts/compile_q80_final_q6_char0_child.sage
```

The final pinned certificate/model are

```text
data/fibrations/q80-final-q6-char0/Q80_CHAR0_FINAL_Q6_CERTIFICATE.md
data/fibrations/q80-final-q6-char0/q80_char0_final_q6_child.sage
```

The historical GF73 final certificate remains useful as a regression:

```text
data/fibrations/kumar_q80_final_q6_cm24_equation_gf73.txt
```

and the machine-readable stage ledger is

```text
data/fibrations/kumar_q80_cm24_equation_progress.tsv.
```

## Reusable equation/compiler results

The combined H3/Q80 work now gives these durable rules:

1. **Binary-quartic covariant maps are 2-covering maps.** Differences of transported quartic points must not be interpreted as primitive MW differences without the factor-of-two check.
2. **Clear the full rational expression before CRT normalization.** For the corrected H3 q8 frame, the `Nx` residue is `Ny*Dx/(h*Dy)`; omitting `Dx` leaves a hidden vertical pole.
3. **Specialize the actual divisor before searching equations.** Generic `P.O`, MW height, twist, and vertical support can change dramatically at CM24.
4. **Connected ADE corrections compile as resolved quotient-line/module conditions**, not one independent row per listed exceptional component.
5. `root_component_data()` may return arbitrary integral root-lattice bases; use discriminant groups rather than assuming Cartan coordinates.
6. For an A3 correction `(-2,-1,-1)`, the exact local module is the middle-component double-vanishing condition.
7. For a D4 correction `(-1,0,-1,-1)`, the ramified-chart outer-complement condition is the deterministic quotient residue `c=0`.
8. In the final A5 correction `(-1,0,-1,-1,0)`, the exact quotient line is the `+/-4` residue pair for the two horizontal signs.
9. **A singular modular marking need not be lifted directly.** Reconstruct easier exact sections first and identify the required section by exact group-law difference plus modular regression.
10. **Keep generic and specialized endpoints separate.** A characteristic-zero equation on CM24 certifies the specialization shadow; it does not replace the generic rootless lattice certificate.
11. **Optimize equation paths separately from lattice paths.** The first rank-growing or smallest-q neighbour need not minimize pole order, coefficient growth, RR dimension, or marking complexity.

The field-generic exact module compatibility layer remains in

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage.
```

## Start here

- [`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md) — short chronology, rank flow, cancellations, and reusable lessons.
- [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md) — proved deductions, conditional lift theorem, and open navigation conjectures.
- [`AGENTS.md`](AGENTS.md) — compact handoff rules, hints, and known dead ends.
- [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md) — human-readable geometric route map and naming scheme.
- [`scripts/README.md`](scripts/README.md) — current proof/compiler/search entry points and archive policy.
- [`SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](SCRIPT_ROUTE_AND_FAILURE_LEDGER.md) — exact path-selection reasons, failed lanes, supersessions, and remaining gaps.
- [`Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md) — completed Q80 terminal characteristic-zero reconstruction and supersession boundary.
- [`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md) — exact Q80 suffix search lineage and successful marking recovery.
- [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md) — Q80 CM24 equation ledger, updated with the characteristic-zero closeout.
- [`H3_Q8_CURRENT_FRONTIER.md`](H3_Q8_CURRENT_FRONTIER.md) — concise H3 q8 exact frontier.
- [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md) — two-bug diagnosis and exact H3 q8 repair.
- [`RESEARCH_UPDATE_2026-08-22.md`](RESEARCH_UPDATE_2026-08-22.md) — historical 2026-08-22 repository-wide checkpoint; later Q80 status is superseded by the closeout note above.
- [`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md) — source recovery, H2/H3 distinction, and exact H3 lattice corridor.
- [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv) — machine-readable Q80 CM24 stage summary.
- [`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md) — Q80 rootless lattice certificate and historical detailed ledger.

## Next strategic gate

For **H3**, the active lifting target after the exact A11 q8/orbit12 equation is
the physical degree-two q4/orbit208 pencil from component-9-zero `2A5/MW7` to
canonical current `3A3/MW8`.  Its special member is the already-explicit I4
cycle `old_zero + P1229 + C10 + C8`; hence `P.O=0` and the RR ambient
dimension is exactly 4.  The characteristic-zero RR and Jacobian compilation now
passes exactly: dimensions `4 -> 2 -> 2`, quartic degree 4, minimal degrees
`(8,12,24)`, fibres `3I4 + 12I1`, smooth infinity and Euler number 24.  Complete
physical-component, all-section, and finite-horizontal-wall gates pass.  The
remaining equation gate is to attach the effective C5 point/sign and full
old-curve marking.  Every NS transport
is unimodular in both directions, and the endpoint is exactly pinned R17.
Its operational score is -1,412 (gross positive burden 1,388), down from the
superseded q10 score 4,471.  See
[`../artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-promoted-route-certificate.json).
The exact compiled equation artifact is
[`../artifacts/local/elkies-k3/q24-2a5-physical-q4o208-rr-qq.json`](../artifacts/local/elkies-k3/q24-2a5-physical-q4o208-rr-qq.json).
The former 4,199, 10,334, and 13,518 targets are withdrawn: respectively they
used a pseudo-zero, a non-section component-10 zero, and a non-nef q104
representative.  The q10 certificates remain exact provenance, but q10 is no
longer the lifting target.

For **Q80**, the terminal marking/RR/equation gate is closed. Treat the direct two-parameter final-section resultants, digit-by-digit `73`-adic lift, and local-73 singularity probes as diagnostics only. Any new Q80 work should start from the exact certificates above rather than reopening the final q6 reconstruction.
