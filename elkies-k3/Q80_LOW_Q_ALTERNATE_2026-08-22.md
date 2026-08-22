# Q80 low-q alternate corridor — 2026-08-22

## Status

This note records the current state of the **secondary/fallback Q80 low-q route**. The route is complete at the generic lattice level and reaches a new rootless `MW17` frame. It remains secondary to the corrected H3 source-polarization route.

The main structural result is that **every retained new neighbour from `D7+D5/MW5` through rootless has chamber-reduced old-fibre degree two**, even though the lattice shells alternate between q4 and q6.

The live Q80 task is no longer lattice discovery. It is equation-level algebraization of the certified route. At CM24, explicit equations are now pinned through **q4_1938**; the next live gate is **q4_6855**.

For the detailed current equation ledger, use [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md). Machine-readable CM24 stage status is in [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv).

## 1. Complete generic low-q route

The common Q80 prefix is

```text
E6+D5+A3/MW3
 --q4 (2,2)--> D9+A4/MW4
 --q4 (2,2)--> D7+D5/MW5.
```

The alternate continuation is

```text
D7+D5/MW5
 --q6 (2,3)--> D7+D4/MW6                 [escape]
 --q4 (2,2)--> A6+A4/MW7                 [orbit 424]
 --q4 (2,2)--> A6+A3/MW8                 [orbit 1222]
 --q6 (2,3)--> A4+A2+A1/MW10             [7774]
 --q4 (2,2)--> A3+A2/MW12                [1938]
 --q4 (2,2)--> 4A1/MW13                  [6855]
 --q4 (2,2)--> A1/MW16                   [candidate 1]
 --q6 (2,3)--> rootless/MW17.
```

The exact vectors are pinned in

```text
data/fibrations/kumar_q80_lowq_alternate_prefix.tsv
data/fibrations/kumar_q80_new_lowq_rootless_path.tsv
data/fibrations/kumar_q80_new_lowq_rootless_final_q6.txt
```

and the full chamber geometry is in

```text
data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv.
```

The new endpoint is genuinely non-canonical: the retained `4A1/MW13` and `A1/MW16` frames are not integrally isometric to their canonical-route counterparts.

## 2. Generic equation geometry

The retained new moves have:

| step | move | child | D.F | P.O | MW height | twist | vertical support |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | escape q6 | `D7+D4/MW6` | 2 | 2 | `8` | 1 | 1 fibre / 5 components |
| 2 | orbit424 q4 | `A6+A4/MW7` | 2 | 1 | `13/4` | 1 | 1 fibre / 3 components |
| 3 | orbit1222 q4 | `A6+A3/MW8` | 2 | 1 | `122/35` | 1 | 1 fibre / 2 components |
| 4 | q6_7774 | `A4+A2+A1/MW10` | 2 | 3 | `219/28` | 0 | 1 fibre / 2 components |
| 5 | q4_1938 | `A3+A2/MW12` | 2 | 1 | `47/10` | 1 | 1 fibre / 2 components |
| 6 | q4_6855 | `4A1/MW13` | 2 | 2 | `19/3` | 0 | none |
| 7 | q4 candidate 1 | `A1/MW16` | 2 | 2 | `6` | 0 | none |
| 8 | final q6 | rootless/MW17 | 2 | 4 | `23/2` | -1 | none |

The effective-section convention is essential: chamber-reduce first, enumerate exact shortest lifts in the root coset, choose the unique section nonnegative on the old chamber, and only then compute the vertical decomposition. Arbitrary shortest root-coset representatives gave incorrect earlier L1 values.

## 3. Reusable compiler results

Q80 work exposed two reusable implementation requirements.

First, exact resolved quotient/module intersections must be coefficient-field generic. The compatibility layer

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage
```

removes the earlier hard `QQ` assumption while preserving historical behavior.

Second, a connected vertical ADE divisor should be compiled as a **single resolved quotient-line condition**, not as independent evaluation rows for every listed exceptional component. q6_7774 supplied the decisive counterexample: componentwise A6 evaluation overconstrained a four-dimensional ambient to dimension one, while the connected quotient-line scan recovered the correct two-dimensional pencil and exact lattice target.

A related component-group lesson from q4_1938 is that `root_component_data()` can return arbitrary integral bases of root sublattices. Component labels must therefore be derived from the Smith discriminant group and a minimal endpoint-weight class, not by assuming a Cartan/simple-root basis.

## 4. q6 escape and orbit 424

The q6 escape is the first productive low-q representative in the old q12 MW coset:

```text
q=6, (a,b)=(2,3)
v=(-5,-3,6,6,-8,-4,2,4,-1,8,-16,-1,0,3,5,-2,-2)
child=D7+D4/MW6.
```

At CM24 it has an exact characteristic-zero binary-quartic model over `QQ(sqrt(-6))`, specializing to `D8+D6+2A1/MW2`.

Orbit 424 is

```text
q=4, (a,b)=(2,2)
v=(32,48,-21,28,8,-52,-34,0,18,5,-23,43,9,-18,16,-6,-6)
child=A6+A4/MW7.
```

At CM24 its marked section becomes rational 2-torsion and the child is exactly `2A7/MW4`, with fibres `I8+I8+8I1`. The characteristic-zero equations for both stages are consolidated in [`scripts/verify_q80_lowq_cm24_equations.sage`](scripts/verify_q80_lowq_cm24_equations.sage).

## 5. Orbit 1222: exact CM24 module and splice

Orbit 1222 is

```text
q=4, (a,b)=(2,2)
v=(10,53,-192,-114,29,-256,-170,-12,-14,74,-32,-14,-6,-26,-58,84,-28)
child=A6+A3/MW8.
```

At CM24 the child is `2A6+3A1/MW3`.

The successful Riemann--Roch construction uses the canonical seven-dimensional raw marked-chord ambient for `O+P+2F`. The smooth `P.O=1` collision contributes rank 2, the infinity I8 restriction rank 2, and the finite A7 vertical gate rank 1, giving the required `7 -> 2` pencil.

This child is **exactly the same GF(73) elliptic fibration as the old pinned third-q12 CM24 child**, not merely the same root signature. A simple pinned base transport is

```text
V = 1/(52T+62),
T = (V+20)/(18V),
```

and the two exact splice maps differ by the involution

```text
V -> (V+28)/(23V+72).
```

## 6. q6_7774: specialization-first correction

Generic q6_7774 has `P.O=3`, height `219/28`, and a two-component vertical correction. Its actual CM24 specialization changes sharply:

```text
P.O = 0
height = 8/7
horizontal = old pinned P3
fiber twist = 2
vertical = one A1 plus one connected A6.
```

Exact rank-18 lattice transport forces

```text
special child = A5+2A4+2A1/MW3
root data = (15,74,600).
```

The connected-A6 quotient-line scan checks `876` exact cases and finds exactly four symmetry-related target hits. A pinned representative is

```text
sign=+1
I7=6
I2=5
c7=58
c2=62
m=T(V-6)(V-5)-(4V+64).
```

The pinned Weierstrass model is recorded in [`data/fibrations/kumar_q80_q6_7774_cm24_weierstrass_gf73.txt`](data/fibrations/kumar_q80_q6_7774_cm24_weierstrass_gf73.txt). Its reducible fibres are `2I2+2I5+I6`.

## 7. q4_1938: exact CM24 equation

The actual q4_1938 specialization is

```text
D.F = 2
D.O = 0
P.O = 1
height = 12/5
fiber twist = 1
vertical = one A4 with coefficients (-1,-1,-1,0).
```

The specialized child is forced to be

```text
2A4+2A3+A1/MW3
root data = (15,66,800).
```

The parent MW basis has profiles

```text
P1=(0,0,3,2,4)
P2=(0,1,1,1,3)
P3=(1,1,3,0,2),
```

and the q4 horizontal is

```text
H=-P1+P2+2P3
coordinates=(-1,1,2)
profile=(0,1,4,4,3)
height=12/5
P.O=1.
```

Node-constrained equation-side reconstruction yields one `+/-` pair. A pinned sign has

```text
X_H=(60+40T+34T^2+58T^3+48T^4+30T^5+3T^6)/(19+52T+T^2),
Y_H=(24+15T+57T^3+28T^4+65T^5+15T^6)/(56+57T+5T^2+T^3).
```

The connected-A4 scan checks

```text
2 signs x 2 I5 fibres x 73 residues = 292 cases.
```

Exactly four cases reduce to quartics, and all four hit `(15,66,800)/MW3`. The pinned representative is

```text
sign=+1
I5=14
residue=34
eval_row=(1,14,50,43,46)
kernel=((1,0,51,40,23),(0,1,36,27,48)).
```

The equation certificate is pinned in [`data/fibrations/kumar_q80_q4_1938_cm24_equation_gf73.txt`](data/fibrations/kumar_q80_q4_1938_cm24_equation_gf73.txt).

## 8. Current live frontier: q4_6855

The generic q4_6855 move is

```text
A3+A2/MW12 --q4--> 4A1/MW13
D.F=2
P.O=2
height=19/3
fiber twist=0
vertical correction=0.
```

Do **not** assume those generic data survive CM24. The next exact task is:

1. transport the actual q4_6855 divisor into the full rank-18 CM24 q4_1938 child;
2. chamber-reduce it and determine the true specialized horizontal/vertical decomposition;
3. predict the special child root data;
4. reconstruct the horizontal on the pinned q4_1938 equation;
5. compile any connected vertical ADE correction using quotient-line/module conditions;
6. continue to q4 candidate 1 and final q6 rootless.

## 9. Scope and priority

The Q80 fallback now has both a complete generic all-degree-two rootless lattice route and a CM24 equation corridor through q4_1938. This is a meaningful equation-development path, but it is **not yet a generic characteristic-zero rootless construction**.

The generic characteristic-zero lift from orbit 1222 onward remains open, including field-of-definition and specialization control for the eventual rational sections. Q80 should only replace the corrected H3 route if that lift becomes materially simpler or faster.

Further Q80 work should therefore stay focused on the existing equation corridor. Do not resume broad lattice shell search unless a specific later equation obstruction requires a new neighbor.
