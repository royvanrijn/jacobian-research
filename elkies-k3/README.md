# Elkies K3 / neighbor programme

This directory contains the K3-surface and elliptic-neighbor work used to search for high-rank elliptic fibrations and explicit equation routes.

For the human-readable map of how the reverse searches, H3 source construction, H2 comparison, Q80 route, and CM24 scaffold fit together, start with [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md). The exact lossless reverse transport from the pinned recovered 17-by-17 lattice to H3, including the q8 marking distinction, is documented in [`RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md`](RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md). For executable entry points and the research history behind failed/superseded scripts, see [`scripts/README.md`](scripts/README.md) and [`SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](SCRIPT_ROUTE_AND_FAILURE_LEDGER.md).

Route discipline: use `R17` for the recovered rootless lattice endpoint,
`H3 source` for the level-474 Kumar polarization, and `next equation child` for the
model produced by the active neighbor compiler.  The current compiler input is
the exact `D13/MW4` parent, its operation is the selected q24 pencil, and its
immediate output is the `D12/MW5` child.  The reverse Low-q MW2/E6 backtracks
and the Q80 CM24 regression endpoint are different routes and must not be used
as substitutes for that parent or child.  See the target vocabulary and
decision ledger in [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md).

## Current priority

The **H3 source family** is the primary characteristic-zero starting point. The currently selected **H3 degree-two corridor** is one certified route from that source to the pinned recovered rootless/MW17 frame; the endpoint identification and complete inverse NS transport are now exact. The route is not proved shortest, globally optimal, or cheapest to compile. Its first three neighbours are exact at equation level:

```text
H3 E7+E8/MW2
 --q6--> E8+E6/MW3
 --q8--> D13/MW4
 --q24 orbit85--> D12/MW5.
```

The q8 repair found two independent bugs in the previous child-side compiler: a binary-quartic 2-cover multiplier was applied twice, and the q-frame CRT normalizer omitted the `Dx` factor of `x(S)=Nx/Dx`. After both corrections the q8 Riemann--Roch problem collapses to an exact `13 -> 2` characteristic-zero intersection and its Jacobian has one `I9*` fibre plus nine `I1`, hence `D13/MW4` exactly.

The selected `D13/MW4 --q24 orbit85--> D12/MW5` equation arrow is now closed. Its exact q24 horizontal section is exported as `artifacts/local/elkies-k3/q8-q24-horizontal-section-qq.json` with status `PASS_EXACT_Q24_HORIZONTAL_SECTION`, and the resolved component-valuation RR lift is exported as `artifacts/local/elkies-k3/q24-d13-to-d12-component-valuation-qq.json` with status `PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR`. The RR profile is `ambient=56`, `collision=48`, `post=8`, `resolved=6`, `kernel=2`; the compiled quartic has degree four and its Jacobian has root rank `12`, root determinant `4`, Euler number `24`, and MW rank `5`.

The active continuation is now the selected `D12/MW5 --q6 orbit42--> A11/MW6` equation arrow. Use the exact orbit42 profile from the backward manifest: `mw=(-1,0,-1,-1,0)`, height `7`, correction `3`, `P.O=3`, fibre twist `0`, and divisor `D42=O+P+V`. Do not reuse the old correction-one, `P.O=2`, or zero-pole search assumptions.

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

For **H3**, continue the selected `D12/MW5 --q6 orbit42--> A11/MW6` equation reconstruction toward R17 as the nearest certified continuation. The completed q24 D12 equation lift should be treated as the parent for this gate; do not reopen q24 search, the historical degree-46 trace interpolation path, `true1600`, or hand-built `corrected1278` q8 compilers except as diagnostics. Do not call the full suffix optimal.

For **Q80**, the terminal marking/RR/equation gate is closed. Treat the direct two-parameter final-section resultants, digit-by-digit `73`-adic lift, and local-73 singularity probes as diagnostics only. Any new Q80 work should start from the exact certificates above rather than reopening the final q6 reconstruction.
