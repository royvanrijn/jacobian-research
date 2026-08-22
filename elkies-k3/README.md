# Elkies K3 / neighbor programme

This directory contains the K3-surface and elliptic-neighbor work used to search for high-rank elliptic fibrations and explicit equation routes.

## Current priority

The first H3 q6 neighbour remains exact, but the **second H3 q8 equation-level hop is paused** after a point-to-Mordell--Weil marking audit found a hard height contradiction. The H3 lattice route remains valid; what is withdrawn is the bridge from one exact rational child section to its claimed q8 MW coordinate.

The Q80 programme is therefore the **live equation-construction route**. This is an execution-priority change, not a claim that Q80 is intrinsically superior to H3.

## Q80 status: CM24 corridor complete

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

As of 2026-08-22 evening, the **entire CM24 equation-development corridor is also algebraized**, through the final q6. The specialized sequence is not rootless: specialization repeatedly changes the horizontal section and ADE type. The late CM24 stages are

```text
2A6+3A1/MW3
 --q6_7774--> A5+2A4+2A1/MW3
 --q4_1938--> 2A4+2A3+A1/MW3
 --q4_6855--> A1+2A3+2D4/MW3
 --q4 candidate1--> A1+A2+A3+A4+A5/MW3
 --final q6--> 4A2+A3+A5/MW2.
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

The completed Q80 corridor established several reusable rules:

1. **Specialize the actual divisor before searching equations.** Generic `P.O`, MW height, twist, and vertical support can change dramatically at CM24.
2. **Connected ADE corrections compile as resolved quotient-line/module conditions**, not one independent row per listed exceptional component.
3. `root_component_data()` may return arbitrary integral root-lattice bases; use discriminant groups rather than assuming Cartan coordinates.
4. For an A3 correction `(-2,-1,-1)`, the exact local module is the middle-component double-vanishing condition.
5. For a D4 correction `(-1,0,-1,-1)`, the ramified-chart outer-complement condition is the deterministic quotient residue `c=0`.
6. In the final A5 correction `(-1,0,-1,-1,0)`, the exact quotient line is the `+/-4` residue pair for the two horizontal signs.

The field-generic exact module compatibility layer remains in

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage.
```

## Start here

- [`RESEARCH_UPDATE_2026-08-22.md`](RESEARCH_UPDATE_2026-08-22.md) — current repository-wide K3 status.
- [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md) — complete Q80 CM24 equation ledger.
- [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv) — machine-readable complete CM24 stage summary.
- [`Q80_LOW_Q_ALTERNATE_2026-08-22.md`](Q80_LOW_Q_ALTERNATE_2026-08-22.md) — generic low-q route plus specialization/equation summary.
- [`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md) — canonical Q80 rootless lattice certificate.
- [`data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv`](data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv) — exact generic geometry of the alternate route.
- [`scripts/verify_q80_new_lowq_rootless_geometry.py`](scripts/verify_q80_new_lowq_rootless_geometry.py) — exact replay of the generic alternate route.
- [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md) — current H3 q8 trust boundary.

## Next strategic gate

Q80 no longer needs additional CM24 neighbour discovery or local-module search. The next Q80 problem is the **generic characteristic-zero lift** from orbit 1222 onward: recover the generic horizontal sections and resolved quotient data, lift the modular identities, control fields of definition, and verify that the final seventeen sections live over the intended characteristic-zero specialization.

Do not infer generic rootlessness from the CM24 endpoint: the special final child has root rank 16. Conversely, do not interpret that special root rank as a failure of the generic route: the generic final q6 is independently certified rootless/MW17.

For H3, do not resume q8 local rank searches until the exact rational child points have been independently matched back to the pinned MW lattice by canonical heights, pairings, zero intersections, and reducible-fibre corrections.
