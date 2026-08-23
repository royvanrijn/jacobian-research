# Elkies K3 / neighbor programme

This directory contains the K3-surface and elliptic-neighbor work used to search for high-rank elliptic fibrations and explicit equation routes.

For the human-readable map of how the reverse searches, H3 source construction, H2 comparison, Q80 route, and CM24 scaffold fit together, start with [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md).

## Current priority

The **H3 Source Route** is the primary characteristic-zero equation route. Both first neighbours are now exact:

```text
H3 E7+E8/MW2
 --q6--> E8+E6/MW3
 --q8--> D13/MW4.
```

The q8 repair found two independent bugs in the previous child-side compiler: a binary-quartic 2-cover multiplier was applied twice, and the q-frame CRT normalizer omitted the `Dx` factor of `x(S)=Nx/Dx`. After both corrections the q8 Riemann--Roch problem collapses to an exact `13 -> 2` characteristic-zero intersection and its Jacobian has one `I9*` fibre plus nine `I1`, hence `D13/MW4` exactly.

The **Q80 Low-q Compiler Route** remains a strong secondary route and an independent compiler scaffold. Its generic lattice corridor reaches rootless/MW17. The characteristic-zero equation pipeline has now progressed through the penultimate `A1/MW16` child; the final q6 problem is preserving/recovering the intended MW marking after the preceding quartic-to-Jacobian conversions. The complete CM24 equation corridor remains the regression scaffold for that final lift.

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

## Q80 status: generic route and CM24 scaffold

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

The old suffix labels `7774`, `1938`, `6855`, and `candidate1` are search-result identifiers rather than intrinsic names. New notes should use the canonical stage names in [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md), while [`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md) records the exact historical search lineage.

The characteristic-zero pipeline now reaches the `A1/MW16` parent of the final q6. The current obstruction is that the preceding equation compiler did not persist enough of the child MW marking through the binary-quartic/Jacobian conversions. The preferred repair is therefore to backtrack the exact selected neighbor transports and preserve the required rational points, not to rediscover the final section from its singular CM24 specialization.

The entire CM24 equation-development corridor is algebraized through the final q6. The late specialized sequence is

```text
2A6+3A1/MW3
 --q6--> A5+2A4+2A1/MW3
 --q4--> 2A4+2A3+A1/MW3
 --q4--> A1+2A3+2D4/MW3
 --q4--> A1+A2+A3+A4+A5/MW3
 --q6--> 4A2+A3+A5/MW2.
```

The final q6 simultaneously passes an independent regression that its **generic** child is rootless/MW17. Thus the CM24 corridor is a complete equation/compiler scaffold for the generic rootless route, not the generic rootless equation itself.

The pinned final CM24 equation certificate is

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

The field-generic exact module compatibility layer remains in

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage.
```

## Start here

- [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md) — canonical human-readable route map and naming scheme.
- [`RESEARCH_UPDATE_2026-08-22.md`](RESEARCH_UPDATE_2026-08-22.md) — current repository-wide K3 status checkpoint.
- [`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md) — source recovery, H2/H3 distinction, and exact H3 lattice route.
- [`H3_Q8_CURRENT_FRONTIER.md`](H3_Q8_CURRENT_FRONTIER.md) — concise H3 q8 exact frontier.
- [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md) — two-bug diagnosis and exact repair.
- [`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md) — exact Q80 suffix search lineage and final-marking backtrack plan.
- [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md) — complete Q80 CM24 equation ledger.
- [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv) — machine-readable complete CM24 stage summary.
- [`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md) — Q80 rootless lattice certificate and historical detailed ledger.

## Next strategic gate

For **H3**, continue from the exact `D13/MW4` child along the certified H3 Source Route toward the rootless/high-rank target. Do not return to the historical degree-46, `true1600`, or hand-built `corrected1278` q8 compilers except as diagnostics.

For **Q80**, preserve/recover the exact MW marking on the characteristic-zero `A1/MW16` parent, then execute the selected terminal `q6` divisor. Use the CM24 `P2-P3` section and quotient residues as regression checks only; they are specialization data, not the generic definition of the final horizontal.
