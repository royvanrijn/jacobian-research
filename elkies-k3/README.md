# Elkies K3 / neighbor programme

This directory contains the K3-surface and elliptic-neighbor work used to search for high-rank elliptic fibrations and explicit equation routes.

## Current priority

The corrected **H3 source-polarization route remains primary**. The Q80 programme is a secondary/fallback route: its generic lattice search is closed, and current Q80 work is equation-level algebraization of the already-certified low-q corridor.

As of 2026-08-22 evening, the Q80 fallback has:

- a complete exact generic route from `E6+D5+A3/MW3` to a new rootless `MW17` frame;
- old-fibre degree `2` at every retained new hop from `D7+D5/MW5` to rootless;
- exact CM24 characteristic-zero equations through the q6 escape and orbit 424;
- an exact modular CM24 equation for orbit 1222, together with a proof that this `2A6+3A1/MW3` fibration is the same GF(73) fibration as the old pinned third-q12 CM24 child;
- an exact modular CM24 equation for q6_7774, with special child `A5+2A4+2A1/MW3` and root data `(15,74,600)`;
- an exact modular CM24 equation for q4_1938, with special child `2A4+2A3+A1/MW3` and root data `(15,66,800)`;
- the q4_1938 horizontal reconstructed as `-P1+P2+2P3`, height `12/5`, profile `(0,1,4,4,3)`, `P.O=1`;
- the live CM24 equation frontier advanced to q4_6855.

The important reusable equation lesson from q6_7774 and q4_1938 is that a connected vertical ADE divisor must be compiled as a **single resolved quotient-line condition**. Treating each listed exceptional component as an independent evaluation row can overconstrain the Riemann--Roch pencil.

## Start here

- [`RESEARCH_UPDATE_2026-08-22.md`](RESEARCH_UPDATE_2026-08-22.md) — concise current status and priorities.
- [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md) — detailed equation-level ledger, including passed vs pending gates.
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
 --q6--> D7+D4/MW6                  escape
 --q4--> A6+A4/MW7                  orbit 424
 --q4--> A6+A3/MW8                  orbit 1222
 --q6--> A4+A2+A1/MW10              7774
 --q4--> A3+A2/MW12                 1938
 --q4--> 4A1/MW13                   6855
 --q4--> A1/MW16                    candidate 1
 --q6--> rootless/MW17
```

The alternate route does **not** replace the corrected H3 route merely because its lattice and CM24 modular equation corridor is progressing well. The remaining strategic gate is a credible generic characteristic-zero equation path through the later Q80 neighbors.
