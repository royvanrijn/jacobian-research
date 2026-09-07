# Q80 low-q alternate corridor — 2026-08-22

## Status

The Q80 low-q alternate route is complete at the generic lattice level and reaches a new rootless `MW17` frame. Every retained new neighbour from `D7+D5/MW5` through rootless has chamber-reduced old-fibre degree two.

The **CM24 equation corridor is also now complete through the final q6**. There is no remaining CM24 neighbour/module search. The remaining Q80 problem is the generic characteristic-zero lift.

The generic and CM24 final endpoints differ:

```text
generic final q6: rootless/MW17
CM24 final q6:    4A2+A3+A5/MW2.
```

This is expected specialization behaviour in the rank-18 CM24 Neron--Severi lattice and does not weaken the generic rootless certificate.

For the detailed equation ledger use [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md). Machine-readable status is in [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv).

## 1. Complete generic low-q route

```text
E6+D5+A3/MW3
 --q4 (2,2)--> D9+A4/MW4
 --q4 (2,2)--> D7+D5/MW5
 --q6 (2,3)--> D7+D4/MW6                 [escape]
 --q4 (2,2)--> A6+A4/MW7                 [orbit 424]
 --q4 (2,2)--> A6+A3/MW8                 [orbit 1222]
 --q6 (2,3)--> A4+A2+A1/MW10             [7774]
 --q4 (2,2)--> A3+A2/MW12                [1938]
 --q4 (2,2)--> 4A1/MW13                  [6855]
 --q4 (2,2)--> A1/MW16                   [candidate 1]
 --q6 (2,3)--> rootless/MW17.
```

Exact vectors and generic chamber geometry are pinned in

```text
data/fibrations/kumar_q80_lowq_alternate_prefix.tsv
data/fibrations/kumar_q80_new_lowq_rootless_path.tsv
data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv
data/fibrations/kumar_q80_new_lowq_rootless_final_q6.txt
```

with replay in

```text
scripts/verify_q80_new_lowq_rootless_geometry.py.
```

The alternate `4A1/MW13` and `A1/MW16` frames are not integrally isometric to the canonical-route frames.

## 2. Generic equation geometry

| step | move | child | D.F | P.O | MW height | twist | generic vertical |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | escape q6 | `D7+D4/MW6` | 2 | 2 | `8` | 1 | 1 fibre / 5 components |
| 2 | orbit424 q4 | `A6+A4/MW7` | 2 | 1 | `13/4` | 1 | 1 fibre / 3 components |
| 3 | orbit1222 q4 | `A6+A3/MW8` | 2 | 1 | `122/35` | 1 | 1 fibre / 2 components |
| 4 | q6_7774 | `A4+A2+A1/MW10` | 2 | 3 | `219/28` | 0 | 1 fibre / 2 components |
| 5 | q4_1938 | `A3+A2/MW12` | 2 | 1 | `47/10` | 1 | 1 fibre / 2 components |
| 6 | q4_6855 | `4A1/MW13` | 2 | 2 | `19/3` | 0 | none |
| 7 | q4 candidate 1 | `A1/MW16` | 2 | 2 | `6` | 0 | none |
| 8 | final q6 | rootless/MW17 | 2 | 4 | `23/2` | -1 | none |

The effective-section rule is: chamber-reduce first, enumerate exact shortest lifts in the root coset, choose the unique lift nonnegative on the old chamber, then compute the vertical decomposition.

## 3. Complete CM24 equation corridor

The late specialized chain is

```text
2A6+3A1/MW3
 --q6_7774--> A5+2A4+2A1/MW3
 --q4_1938--> 2A4+2A3+A1/MW3
 --q4_6855--> A1+2A3+2D4/MW3
 --q4 candidate1--> A1+A2+A3+A4+A5/MW3
 --final q6--> 4A2+A3+A5/MW2.
```

The corresponding specialized horizontals are

```text
q6_7774:       P3,              P.O=0, height=8/7
q4_1938:      -P1+P2+2P3,      P.O=1, height=12/5
q4_6855:       2P1,             P.O=0, height=3/5
q4 candidate1: -P3,             P.O=0, height=3/4
final q6:      P2-P3,           P.O=0, height=1.
```

This repeated collapse is why generic `P.O` and vertical data must never be assumed to survive specialization.

## 4. Equation/compiler lessons

### Connected ADE quotient rule

A connected vertical ADE divisor is one resolved quotient/module condition, not one independent condition per exceptional component. q6_7774 was the decisive counterexample: componentwise A6 evaluation overconstrained the pencil, while a single connected quotient-line scan recovered the exact child.

### q4_6855 A3 middle-double rule

At CM24 q4_6855 has

```text
H=2P1
vertical A3=(-2,-1,-1).
```

The A3 multiplicity pattern compiles deterministically as a node-value plus pure-base first-order condition. The resulting child is

```text
A1+2A3+2D4/MW3
root_data=(15,74,512).
```

Certificate:

```text
data/fibrations/kumar_q80_q4_6855_cm24_equation_gf73.txt
```

### q4 candidate-1 D4 outer-complement rule

Candidate 1 specializes to

```text
H=-P3
vertical A3=(-1,-1,-1)
vertical D4=(-1,0,-1,-1).
```

The ramified I0* chart predicts D4 quotient residue `c=0`. All eight target-compatible horizontal cases use `c=0`, independently validating the rule. The child is

```text
A1+A2+A3+A4+A5/MW3
root_data=(15,70,720).
```

Certificate:

```text
data/fibrations/kumar_q80_q4_a1_candidate1_cm24_equation_gf73.txt
```

### final q6 A5 quotient

The final specialized divisor is

```text
H=P2-P3
profile=(0,2,2,0,4)
vertical A4=(-1,-1,-1,-1)
vertical A5=(-1,0,-1,-1,0).
```

The minimal parent has fibres

```text
I2@60, I3@23, I4@24, I5@25, I6@47.
```

Six polynomial horizontals occur. The A5 quotient search reduces to eight symbolic survivors and exactly two symmetry-related target hits. The pinned one has residue `-4`; the opposite horizontal sign has `+4` and gives the same quartic/Jacobian.

The child fibres are

```text
4 I3 + I4 + I6 + 2 I1,
```

hence

```text
4A2+A3+A5/MW2
root_data=(16,66,1944).
```

Certificate:

```text
data/fibrations/kumar_q80_final_q6_cm24_equation_gf73.txt
```

## 5. Earlier CM24 stages

The q6 escape and orbit 424 have exact characteristic-zero equations over `QQ(sqrt(-6))`. Orbit 1222 has an exact GF(73) equation and is proved to be the same fibration as the old pinned third-q12 CM24 child. q6_7774 and q4_1938 have their own pinned GF(73) certificates:

```text
data/fibrations/kumar_q80_q6_7774_cm24_weierstrass_gf73.txt
data/fibrations/kumar_q80_q4_1938_cm24_equation_gf73.txt
```

## 6. Reusable implementation results

The corridor also established:

- field-generic exact quotient/module intersections via `elliptic_neighbor_compiler_field_generic.sage`;
- discriminant-group component labeling when root-component bases are arbitrary;
- specialization-first section recovery;
- exact node-constrained polynomial-section reconstruction;
- cheap finite-specialization prefilters before expensive symbolic quartic classification.

## 7. Next strategic problem

The Q80 CM24 corridor is closed. Do not resume broad Q80 shell or modular local-module search without a specific obstruction.

The next Q80 task is the **generic characteristic-zero lift from orbit 1222 onward**:

1. recover the generic horizontal sections and fields of definition;
2. lift the CM24 resolved quotient conditions to characteristic zero;
3. construct and minimize each generic neighbour equation;
4. track fibre components, sections, and Galois fields through the chain;
5. verify the final seventeen independent sections on the intended characteristic-zero specialization before any rank claim.

The completed CM24 corridor is now the regression scaffold for that lift.
