# Research update — 2026-08-22

## H3 q8 repaired exactly over QQ

The H3 route is again exact through its first two neighbours:

```text
H3 E7+E8/MW2
 --q6--> E8+E6/MW3
 --q8--> D13/MW4.
```

The earlier q8 pause was productive: two independent compiler/marking bugs were identified and repaired.

### Binary-quartic 2-cover factor

The q6 child's binary-quartic covariant map is the degree-two covering map to the Jacobian. The old q8 code treated differences of covariant images as primitive MW differences and doubled them again.

For

```text
Pmap=phi(E7_7)-phi(old_O),
Qmap=phi(E7_7)-phi(affine_E7),
```

the actual geometric height data are

```text
height(Pmap)=32/3
height(Qmap)=32/3
<Pmap,Qmap>=4/3
height(Pmap+Qmap)=24
height(2Pmap+2Qmap)=96.
```

Thus the corrected q8 marked section is

```text
S=Pmap+Qmap,
MW=(-2,-2,0),
S.O=10,
height=24,
collision degree=10.
```

The withdrawn section is exactly `2*S`, explaining its collision degree `46` and height `96`.

### Missing `Dx` in q normalization

For `x(S)=Nx/Dx`, `y(S)=Ny/Dy`, `p=-y(S)/x(S)`, `q=(m-p)/h`, the correct pole cancellation is

```text
R*h*Dy == Ny*Dx mod Nx.
```

The old formula omitted `Dx`, leaving a hidden degree-24 vertical pole. Before the repair the recovered modular pencil had generic branch degree about `100`; after the repair every nonsingular finite level has branch degree `4` at both `p=43` and `p=59`.

### Exact QQ certificate

The corrected q8 global Riemann--Roch system has

```text
h degree = 10
Nx degree = 24
Ny degree = 36
deg(s) <= 6
deg(t) <= 5
ambient = 13
condition rank = 11
kernel = 2.
```

The exact pencil eliminates to a quartic. Its Jacobian has

```text
I9* + 9 I1,
infinity smooth,
root rank 13,
root determinant 4,
Euler number 24,
MW rank 4 (rho=19).
```

Therefore the q8 child is exactly `D13/MW4`.

Reproduce with

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_marking.sage
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --representative component-nef \
  --output artifacts/local/elkies-k3/q8-target-component-nef-audit.json
cmp artifacts/local/elkies-k3/q8-target-component-nef-audit.json \
  elkies-k3/data/fibrations/h3_q8_component_nef_physical_root_target.json
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_corrected2cover_qq.sage
```

The detailed trust/supersession ledger is [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md).

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

The full CM24 equation corridor is exact through the final q6. The specialization sequence after orbit 1222 is

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

The final equation certificate is

```text
data/fibrations/kumar_q80_final_q6_cm24_equation_gf73.txt
```

and the complete stage ledger is

```text
data/fibrations/kumar_q80_cm24_equation_progress.tsv.
```

## Reusable compiler/math results

1. Binary-quartic covariant maps are 2-covering maps; check the MW multiplier before converting point differences to lattice coordinates.
2. Clear the complete rational expression before CRT normalization; section-coordinate denominators can change the residue.
3. Specialize the actual divisor before equation search.
4. Connected ADE corrections are quotient-line/module conditions, not independent component rows.
5. Component labels must be discriminant-group based when root bases are non-Cartan.
6. A3 `(-2,-1,-1)` is the middle-component double-vanishing module.
7. D4 `(-1,0,-1,-1)` is the ramified-chart outer-complement row with residue `0`.
8. The final A5 `(-1,0,-1,-1,0)` resolves to the `+/-4` quotient pair for the two horizontal signs.

## Current execution order

1. **H3 primary:** continue from the exact `D13/MW4` child toward the rootless/high-rank target, keeping the 2-cover and full-residue regressions active.
2. **Q80 secondary:** use the completed CM24 corridor to attack the generic characteristic-zero lift from orbit 1222 onward.
3. Do not revive the historical H3 degree-46, `true1600`, or hand-built `corrected1278` q8 pipelines as canonical constructions.
