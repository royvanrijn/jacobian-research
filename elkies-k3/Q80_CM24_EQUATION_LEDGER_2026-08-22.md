# Q80 CM24 equation ledger — 2026-08-22 evening

## Status

The generic Q80 low-q route is certified to a new rootless `MW17` frame, and the **entire CM24 equation corridor is now certified through the final q6**. No CM24 neighbour or local-module gate remains open.

The final generic and specialized endpoints are different and both are exact:

```text
generic final q6 child: rootless/MW17
CM24 final q6 child:    4A2+A3+A5/MW2
```

CM24 is a rank-18 specialization and repeatedly changes horizontal sections, fibre corrections, and ADE type. It is the equation/compiler scaffold for the generic route, not the generic rootless equation itself.

## Complete stage table

| stage | generic child | CM24 child | special horizontal | special module | status |
|---|---|---|---|---|---|
| q6 escape | `D7+D4/MW6` | `D8+D6+2A1/MW2` | old q12 section | E6 jet | exact `QQ(sqrt(-6))` |
| orbit 424 | `A6+A4/MW7` | `2A7/MW4` | 2-torsion | degree-two chord | exact `QQ(sqrt(-6))` |
| orbit 1222 | `A6+A3/MW8` | `2A6+3A1/MW3` | `P.O=1`, `25/8` | raw saturated A7 chord | exact GF(73), old-q12 splice proved |
| q6_7774 | `A4+A2+A1/MW10` | `A5+2A4+2A1/MW3` | `P3`, `P.O=0`, `8/7` | A1 + connected A6 quotient | exact GF(73) |
| q4_1938 | `A3+A2/MW12` | `2A4+2A3+A1/MW3` | `-P1+P2+2P3`, `12/5` | smooth saturation + A4 quotient | exact GF(73) |
| q4_6855 | `4A1/MW13` | `A1+2A3+2D4/MW3` | `2P1`, `3/5` | A3 middle-double | exact GF(73) |
| q4 candidate 1 | `A1/MW16` | `A1+A2+A3+A4+A5/MW3` | `-P3`, `3/4` | whole A3 + D4 `c=0` | exact GF(73) |
| final q6 | rootless/MW17 | `4A2+A3+A5/MW2` | `P2-P3`, `1` | whole A4 + A5 `c=+/-4` | exact GF(73) |

Machine-readable status:

```text
data/fibrations/kumar_q80_cm24_equation_progress.tsv
```

## Pinned equation certificates

```text
data/fibrations/kumar_q80_q6_7774_cm24_weierstrass_gf73.txt
data/fibrations/kumar_q80_q4_1938_cm24_equation_gf73.txt
data/fibrations/kumar_q80_q4_6855_cm24_equation_gf73.txt
data/fibrations/kumar_q80_q4_a1_candidate1_cm24_equation_gf73.txt
data/fibrations/kumar_q80_final_q6_cm24_equation_gf73.txt
```

The orbit-1222 child is additionally proved to be the same GF(73) fibration as the old pinned third-q12 CM24 child.

## Late-stage specialization data

### q4_6855

```text
H=2P1
P.O=0
height=3/5
twist=2
vertical A3=(-2,-1,-1)
child=A1+2A3+2D4/MW3
root_data=(15,74,512)
```

The A3 pattern is the exact middle-component double-vanishing module; no residue scan is needed.

### q4 candidate 1

```text
H=-P3
P.O=0
height=3/4
twist=2
vertical A3=(-1,-1,-1)
vertical D4=(-1,0,-1,-1)
child=A1+A2+A3+A4+A5/MW3
root_data=(15,70,720)
```

The D4 ramified-chart outer-complement rule gives quotient residue `c=0`; all eight admissible horizontal candidates that hit the target use `c=0`.

### final q6

```text
H=P2-P3
profile=(0,2,2,0,4)
P.O=0
height=1
twist=2
vertical A4=(-1,-1,-1,-1)
vertical A5=(-1,0,-1,-1,0)
```

The minimal candidate-1 parent has reducible fibres

```text
I2@60, I3@23, I4@24, I5@25, I6@47.
```

The horizontal hits exactly I3, I4, and I6. Six polynomial candidates occur. The A5 quotient prefilter reduces `6x73` cases to eight symbolic survivors. Exactly two symmetry-related cases hit the target. The pinned case has

```text
H_X=15+29W+9W^2+59W^3+4W^4
H_Y=7+61W+7W^2+3W^3+14W^4+57W^5+8W^6
A4 root=25, chord value=59
A5 root=47, residue=69=-4
kernel=((1,0,50,37),(0,1,16,33)).
```

Its fibres are

```text
4 I3 + I4 + I6 + 2 I1
```

so the child is exactly

```text
4A2+A3+A5/MW2
root_data=(16,66,1944).
```

The opposite horizontal sign uses residue `+4` and gives the same quartic/Jacobian.

## Reusable rules established

1. Specialize the actual divisor before equation search.
2. Connected ADE support is one resolved quotient/module condition, not independent component rows.
3. Use discriminant groups when root-component bases are non-Cartan.
4. A3 `(-2,-1,-1)` is a middle-double condition.
5. D4 `(-1,0,-1,-1)` is the outer-complement `c=0` condition.
6. Final A5 `(-1,0,-1,-1,0)` gives the `+/-4` quotient pair.

## Next Q80 problem

CM24 work is complete. The next task is the **generic characteristic-zero lift from orbit 1222 onward**: recover generic horizontal sections, lift resolved quotient data, track fields of definition and minimal models, and verify the final seventeen sections on the intended characteristic-zero specialization before any rank claim.
