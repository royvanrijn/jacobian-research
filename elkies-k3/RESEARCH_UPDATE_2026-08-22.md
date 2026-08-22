# Research update — 2026-08-22

## Q80 secondary route: complete all-degree-two rootless corridor

The Q80 secondary/fallback route is now certified all the way to a **new rootless MW17 frame**. It remains secondary to the corrected H3 source-polarization route, but it is no longer merely a partial lattice corridor.

The exact retained continuation from the common `D7+D5/MW5` source is

```text
D7+D5/MW5
 --q6 (2,3)--> D7+D4/MW6
 --q4 (2,2)--> A6+A4/MW7
 --q4 (2,2)--> A6+A3/MW8
 --q6 (2,3)--> A4+A2+A1/MW10
 --q4 (2,2)--> A3+A2/MW12
 --q4 (2,2)--> 4A1/MW13
 --q4 (2,2)--> A1/MW16
 --q6 (2,3)--> rootless/MW17.
```

The main structural result is that **all eight retained new divisors chamber-reduce to old-fibre degree two**:

```text
degree_distribution = ((2, 8),)
P.O distribution:
  1 : 3 moves
  2 : 3 moves
  3 : 1 move
  4 : 1 move
```

The corrected effective-section normalization also improves the earlier generic vertical-complexity records: orbit 424 has one-fibre/three-component vertical correction with L1 `3` rather than `12`, and orbit 1222 has one-fibre/two-component correction with L1 `2` rather than `9`. The old values came from arbitrary shortest root-coset representatives rather than the unique shortest section effective in the old fibration chamber.

The later retained branch is

```text
A6+A3/MW8
 --q6--> A4+A2+A1/MW10       [7774]
 --q4--> A3+A2/MW12           [1938]
 --q4--> 4A1/MW13             [6855]
 --q4--> A1/MW16              [candidate 1]
 --q6--> rootless/MW17.
```

The new `4A1/MW13` and `A1/MW16` frames are not isometric to their canonical-route counterparts. The final rootless q6 was found directly in the new A1 class by decomposing the rank-16 MW search by A1 root pairing and using an exact rational LDL / Fincke--Pohst shell streamer. The successful final shell has A1 pairing `p=1`.

Machine-readable results and replay/search scripts are:

```text
data/fibrations/kumar_q80_new_lowq_rootless_path.tsv
data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv
data/fibrations/kumar_q80_new_lowq_rootless_final_q6.txt
data/fibrations/kumar_q80_a6a3_q6_chamber_scores.tsv
data/fibrations/kumar_q80_7774_q4_rank5_scores.tsv
data/fibrations/kumar_q80_1938_q4_4a1_scores.tsv
scripts/verify_q80_new_lowq_rootless_geometry.py
scripts/search_q80_new_lowq_final_q6_rootless.py
```

The first q6 and orbit-424 moves still have explicit CM24 characteristic-zero binary-quartic equations as recorded in [`Q80_LOW_Q_ALTERNATE_2026-08-22.md`](Q80_LOW_Q_ALTERNATE_2026-08-22.md). The remaining degree-two hops are now exact lattice/equation-geometry targets for the generic neighbour compiler.

## Reusable compiler update

The Q80 third-q12 module work also identified that the exact neighbour
compiler's finite-quotient stack was unnecessarily restricted to `QQ`.
A compatibility layer now makes the local quotient/module intersection
machinery field-generic while preserving the historical `QQ` default:

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage
```

This supports exact number-field module intersections such as the Q80 CM24
compositum `QQ(sqrt(-6),sqrt(-3))`.

## Current priority

The Q80 lattice search is closed for this fallback: it now has a complete exact all-degree-two route to rootless/MW17. Do not spend more time extending Q80 neighbour shells unless needed for a specific equation obstruction.

If Q80 is pursued further, the priority is to algebraize the retained degree-two hops after `A6+A3/MW8` with the field-generic neighbour/module machinery and determine whether they give a materially faster equation-level rootless model than the corrected H3 route.

The H3 route remains the primary source-polarization route.

