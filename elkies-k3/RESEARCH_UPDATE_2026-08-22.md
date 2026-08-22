# Research update — 2026-08-22

## H3 q8 re-audit changes the active equation priority

The first H3 neighbour remains exact:

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3.
```

Its complete resolved RR cover, two-dimensional q6 pencil, minimal child equation, fibre classification, MW height lattice, and Neron--Severi discriminant check remain re-authorized.

The second q8 equation-level hop is **paused**. A regression audit proves that the exact rational child section used by `derive_h92_q6_child_q8_marking.sage` cannot have the MW coordinate assigned to it. The claimed coordinate `(-2,-2,0)` has height `24`, while the same rational section has `S.O=46`; even a deliberately loose reducible-fibre correction bound forces height at least `60`.

The authoritative H3 q8 checkpoint is [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md). The H3 lattice q8 representative remains valid; the point-to-MW bridge and dependent local compiler do not.

## Q80 CM24 equation corridor is complete

The Q80 generic alternate route remains certified all the way to a new rootless `MW17` frame:

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

All eight retained divisors chamber-reduce to old-fibre degree `2`.

The **full CM24 equation corridor is now exact through the final q6**. The specialization sequence after orbit 1222 is

| stage | CM24 child | special horizontal | equation status |
|---|---|---|---|
| orbit 1222 | `2A6+3A1/MW3` | `P.O=1`, height `25/8` | exact GF(73), spliced to old q12 child |
| q6_7774 | `A5+2A4+2A1/MW3` | `P3`, `P.O=0`, height `8/7` | exact GF(73) |
| q4_1938 | `2A4+2A3+A1/MW3` | `-P1+P2+2P3`, height `12/5` | exact GF(73) |
| q4_6855 | `A1+2A3+2D4/MW3` | `2P1`, height `3/5` | exact GF(73) |
| q4 candidate 1 | `A1+A2+A3+A4+A5/MW3` | `-P3`, height `3/4` | exact GF(73) |
| final q6 | `4A2+A3+A5/MW2` | `P2-P3`, height `1` | exact GF(73) |

The final special q6 has fibre signature

```text
4 I3 + I4 + I6 + 2 I1
= 4A2 + A3 + A5,
root data = (16,66,1944),
MW = 2.
```

Its exact modular construction has six node-compatible horizontal sections (three `+/-` pairs). A multi-specialization filter reduces the A5 quotient search to eight symbolic cases. Exactly two symmetry-related cases hit the lattice target; the pinned one has

```text
H = P2-P3
H_X = 15 + 29W + 9W^2 + 59W^3 + 4W^4
H_Y = 7 + 61W + 7W^2 + 3W^3 + 14W^4 + 57W^5 + 8W^6
whole A4 at I5 W=25, chord value 59
A5 at I6 W=47, quotient residue 69 = -4
kernel = ((1,0,50,37),(0,1,16,33)).
```

The opposite horizontal sign uses residue `+4` and gives the same quartic/Jacobian.

The final equation certificate is

```text
data/fibrations/kumar_q80_final_q6_cm24_equation_gf73.txt
```

and the complete stage ledger is

```text
data/fibrations/kumar_q80_cm24_equation_progress.tsv.
```

## Important interpretation

The final q6 passes two different statements simultaneously:

```text
generic child: rootless/MW17
CM24 special child: 4A2+A3+A5/MW2.
```

There is no contradiction. CM24 is a rank-18 specialization and repeatedly changes the horizontal section, fibre correction, and root system. The completed CM24 corridor is an exact **equation/compiler scaffold** for the generic route; it is not itself the desired characteristic-zero rootless equation.

## Reusable compiler/math results

The completed Q80 corridor gives a useful library of local rules:

1. specialize the actual divisor before equation search;
2. connected ADE corrections are quotient-line/module conditions, not independent component rows;
3. component labels must be discriminant-group based when root bases are non-Cartan;
4. A3 `(-2,-1,-1)` is the middle-component double-vanishing module;
5. D4 `(-1,0,-1,-1)` is the ramified-chart outer-complement row with residue `0`;
6. the final A5 `(-1,0,-1,-1,0)` resolves to the `+/-4` quotient pair for the two horizontal signs.

The field-generic compatibility layer remains:

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage.
```

## Current execution order

1. **Q80:** stop CM24 neighbour/module discovery; begin the generic characteristic-zero lift from orbit 1222 onward, using the now-complete modular corridor as the exact target/regression scaffold.
2. Control the generic fields of definition, horizontal sections, and resolved quotient lines at each later hop; verify the final seventeen sections on the intended characteristic-zero specialization before any rank claim.
3. **H3:** independently repair the q8 rational-point/MW bridge by canonical heights, pairings, zero intersections, and fibre corrections.
4. Do not resume `true1600`, `corrected1278`, or further H3 q8 local jet screens until that bridge is certified.
