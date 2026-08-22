# Elkies K3 / neighbor programme

This directory contains the K3-surface and elliptic-neighbor work used to search for high-rank elliptic fibrations and explicit equation routes.

## Current priority

The first H3 q6 neighbour remains exact, but the **second H3 q8 equation-level hop is paused** after a new point-to-Mordell--Weil marking audit found a hard contradiction. The H3 lattice route itself remains valid; what is withdrawn is the bridge from one exact rational child section to the claimed q8 MW coordinate.

The Q80 programme is therefore the **live equation-construction route** for now. This is an execution-priority change, not a claim that Q80 is mathematically superior to H3.

### H3 q8 audit status

The exact q6 child is still certified as `E8+E6/MW3`, with MW height Gram

```text
[[8/3,1/3,-1],
 [1/3,8/3,1],
 [-1,1,46]].
```

The q8 marking script assigns its constructed rational section the coordinate `(-2,-2,0)`, hence height `24`. But that same rational section meets the standard zero section transversely at 46 smooth fibres. Shioda's height formula, even with a deliberately loose `E8+E6` correction bound, forces height at least `60`. The point-to-MW bridge is therefore invalid.

See:

- [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md) — authoritative re-audit, retractions, and next gate.
- [`H3_Q8_CURRENT_FRONTIER.md`](H3_Q8_CURRENT_FRONTIER.md) — compact active H3 q8 state.
- [`scripts/audit_h92_q6_child_q8_marking_height.sage`](scripts/audit_h92_q6_child_q8_marking_height.sage) — exact regression reproducing the height contradiction.
- [`scripts/audit_h92_q8_representative_selection.sage`](scripts/audit_h92_q8_representative_selection.sage) — reconciles the degree-18 dominant and degree-16 classifier-nef source representatives.

The source-side degree-16 lattice class remains re-authorized: it is exactly the old-fibre root reduction of the independently certified classifier-nef q8 representative. What is **not** authorized is the later hand-converted `q6^8` local RR compiler built around that class.

## Q80 live route

The Q80 generic route is certified to a new rootless `MW17` frame. Equation work has advanced through the CM24 low-q corridor, including q6_7774 and q4_1938, and continues along the retained degree-two neighbours toward the rootless endpoint.

The important reusable equation lesson from q6_7774 and q4_1938 is that a connected vertical ADE divisor must be compiled as a **single resolved quotient-line condition**. Treating each listed exceptional component as an independent evaluation row can overconstrain the Riemann--Roch pencil.

## Start here

- [`RESEARCH_UPDATE_2026-08-22.md`](RESEARCH_UPDATE_2026-08-22.md) — concise repository-wide K3 status.
- [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md) — current H3 q8 trust boundary.
- [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md) — detailed Q80 equation-level ledger.
- [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv) — machine-readable CM24 stage summary.
- [`Q80_LOW_Q_ALTERNATE_2026-08-22.md`](Q80_LOW_Q_ALTERNATE_2026-08-22.md) — current low-q route summary and discovery history.
- [`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md) — canonical Q80 rootless lattice certificate.
- [`data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv`](data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv) — exact generic geometry of the completed alternate corridor.
- [`scripts/verify_q80_new_lowq_rootless_geometry.py`](scripts/verify_q80_new_lowq_rootless_geometry.py) — exact replay of the alternate generic route and chamber geometry.

## Generic alternate Q80 route

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
 --q6--> rootless/MW17
```

Do not resume H3 q8 local rank searches until the exact rational child points have been independently matched back to the pinned MW lattice by canonical heights, pairings, zero intersections, and reducible-fibre corrections.
