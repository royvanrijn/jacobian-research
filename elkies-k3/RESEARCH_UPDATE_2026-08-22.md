# Research update — 2026-08-22

## Q80 secondary route: lattice corridor complete, CM24 equations now through q4_1938

The Q80 secondary/fallback route is certified all the way to a **new rootless MW17 frame**. It remains secondary to the corrected H3 source-polarization route, but the equation-level picture has advanced substantially during 2026-08-22.

The exact retained generic continuation from the common `D7+D5/MW5` source is

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

All eight retained new divisors chamber-reduce to old-fibre degree `2`. The corrected effective-section normalization remains the generic source of truth for `P.O`, fibre twist, and vertical coefficients.

Machine-readable generic results remain in

```text
data/fibrations/kumar_q80_new_lowq_rootless_path.tsv
data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv
data/fibrations/kumar_q80_new_lowq_rootless_final_q6.txt
data/fibrations/kumar_q80_a6a3_q6_chamber_scores.tsv
data/fibrations/kumar_q80_7774_q4_rank5_scores.tsv
data/fibrations/kumar_q80_1938_q4_4a1_scores.tsv
```

with replay/search scripts

```text
scripts/verify_q80_new_lowq_rootless_geometry.py
scripts/search_q80_new_lowq_final_q6_rootless.py
```

## CM24 equation corridor

The equation-level CM24 corridor is now certified through **q4_1938**:

| stage | CM24 child | status |
|---|---|---|
| q6 escape | `D8+D6+2A1/MW2` | exact over `QQ(sqrt(-6))` |
| orbit 424 | `2A7/MW4` | exact over `QQ(sqrt(-6))` |
| orbit 1222 | `2A6+3A1/MW3` | exact over GF(73); same fibration as old q12 CM24 child |
| q6_7774 | `A5+2A4+2A1/MW3` | exact over GF(73) |
| q4_1938 | `2A4+2A3+A1/MW3` | exact over GF(73) |
| q4_6855 | TBD | next live gate |

The detailed status is in [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md), with machine-readable stage data in [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv).

## Orbit 1222: the CM24 splice is proved

The orbit-1222 child is `2A6+3A1/MW3` at CM24. The successful Riemann--Roch construction uses the canonical seven-dimensional raw marked-chord ambient, not the earlier guessed `qsat` subspaces.

The resulting model is exactly the same GF(73) elliptic fibration as the old pinned third-q12 CM24 child. Two PGL2/Weierstrass identifications exist, differing by the old model's base involution

```text
V -> (V+28)/(23V+72).
```

This closes the previously open “same root signature vs same fibration” question.

## q6_7774: specialization-first correction

The generic q6_7774 horizontal has `P.O=3` and height `219/28`, but its actual CM24 specialization is completely different:

```text
P.O = 0
height = 8/7
horizontal = old pinned P3
fiber twist = 2
vertical support = one A1 plus one connected A6.
```

The rank-18 CM24 lattice forces the special child

```text
A5 + 2A4 + 2A1
root data = (15,74,600)
MW = 3.
```

A key reusable lesson came from the failed first module: imposing one independent row per exceptional component of the connected A6 overconstrains the pencil. The correct object is a **single resolved quotient-line condition for the connected ADE divisor**.

Scanning the 73 quotient lines in the resolved A6 row span gives exactly four symmetry-related target hits. A pinned representative has

```text
sign = +1
I7 = 6
I2 = 5
c7 = 58
c2 = 62
m = T(V-6)(V-5) - (4V+64).
```

The pinned q6_7774 CM24 Weierstrass model is

```text
A(T) =
46*T^8 + 5*T^7 + 16*T^6 + 44*T^5 + 6*T^4
+ 13*T^3 + T^2 + T,

B(T) =
54*T^12 + 58*T^11 + 48*T^10 + 16*T^9 + 42*T^8
+ 67*T^7 + 25*T^6 + 19*T^5 + 27*T^4 + 45*T^3
+ 61*T^2 + 44*T + 49.
```

Its reducible fibres are `2I2+2I5+I6`, giving `2A1+2A4+A5/MW3` exactly.

## q4_1938: equation now pinned

The actual q4_1938 specialization is comparatively stable:

```text
D.F = 2
D.O = 0
P.O = 1
MW height = 12/5
fiber twist = 1
vertical support = one A4
vertical coefficients = (-1,-1,-1,0).
```

The CM24 child is predicted by the lattice as

```text
2A4 + 2A3 + A1
root data = (15,66,800)
MW = 3.
```

The parent MW basis is pinned with profiles

```text
P1 = (0,0,3,2,4)
P2 = (0,1,1,1,3)
P3 = (1,1,3,0,2)
```

and the q4 horizontal is

```text
H = -P1 + P2 + 2P3
coordinates = (-1,1,2)
profile = (0,1,4,4,3)
height = 12/5
P.O = 1.
```

Equation-side node-constrained searches recover a single `+/-` pair for `H`. A pinned sign has

```text
X_H =
(60 + 40T + 34T^2 + 58T^3 + 48T^4 + 30T^5 + 3T^6)
/ (19 + 52T + T^2),

Y_H =
(24 + 15T + 57T^3 + 28T^4 + 65T^5 + 15T^6)
/ (56 + 57T + 5T^2 + T^3).
```

The connected-A4 scan checked

```text
2 signs x 2 I5 fibres x 73 residues = 292 cases.
```

Exactly four cases reduce to quartics, and all four have the required `(15,66,800)/MW3` child. The pinned representative is

```text
sign = +1
I5 = 14
residue = 34
evaluation row = (1,14,50,43,46)
kernel = ((1,0,51,40,23),(0,1,36,27,48)).
```

Thus q4_1938 is now **closed at the CM24 equation level**.

## Reusable compiler/math updates

Two durable lessons should be carried into later neighbors:

1. **Specialize the actual divisor before searching sections.** Generic `P.O`, MW height, and vertical support can change drastically at CM24, as q6_7774 demonstrated.
2. **Compile connected ADE corrections as quotient-line/module conditions.** Independent exceptional-component evaluations can produce spurious extra rank.

A second implementation lesson from q4_1938 is that component-group labeling must be basis-independent. `root_component_data()` can return an arbitrary integral basis of an `A_n` root sublattice; use the Smith discriminant group and a minimal endpoint-weight class rather than assuming a Cartan/simple-root basis.

The earlier field-generic compatibility layer remains relevant:

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage
```

## Current priority

The Q80 lattice search is closed. The **live Q80 equation gate is q4_6855**.

The next sequence is:

1. transport the actual q4_6855 divisor through the full rank-18 CM24 lattice;
2. chamber-reduce it and identify its true specialized horizontal/vertical decomposition;
3. predict the CM24 child root data before equation search;
4. reconstruct the horizontal on the pinned q4_1938 equation;
5. compile any connected ADE vertical correction using the quotient-line rule;
6. continue to q4 candidate 1 and the final rootless q6.

Do not spend more time on new Q80 lattice shell searches unless an equation obstruction demands it. The H3 route remains the primary source-polarization route, and the strategic Q80 question is still whether this modular/equation corridor can be lifted to a materially simpler generic characteristic-zero rootless construction.
