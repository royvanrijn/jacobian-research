# Q80 CM24 equation ledger — 2026-08-22 evening

## Scope

This ledger is the source of truth for the equation-level status of the **secondary/fallback Q80 low-q corridor**. The generic lattice route is already complete to a new rootless `MW17` frame. This file separates:

1. generic lattice/chamber certificates;
2. CM24 lattice-specialization certificates;
3. CM24 equation certificates;
4. the next live equation gate.

The corrected H3 source-polarization route remains primary. Nothing here upgrades Q80 above H3.

## Generic route

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

All eight retained divisors have chamber-reduced old-fibre degree `2`. The machine-readable generic route and chamber geometry are in:

```text
data/fibrations/kumar_q80_new_lowq_rootless_path.tsv
data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv
```

## Equation-level stage ledger

| stage | generic child | CM24 special child | special horizontal | special vertical/module | equation status |
|---|---|---|---|---|---|
| q6 escape | `D7+D4/MW6` | `D8+D6+2A1/MW2` | old q12 section `Q=P1+3P2` | exact E6 jet, binary quartic | **exact over `QQ(sqrt(-6))`** |
| orbit 424 | `A6+A4/MW7` | `2A7/MW4` | rational 2-torsion | exact degree-two chord | **exact over `QQ(sqrt(-6))`** |
| orbit 1222 | `A6+A3/MW8` | `2A6+3A1/MW3` | `P.O=1`, height `25/8` | canonical raw saturated chord module | **exact over GF(73); splice to old q12 child proved** |
| q6_7774 | `A4+A2+A1/MW10` | `A5+2A4+2A1/MW3` | old-model `P3`, `P.O=0`, height `8/7` | one A1 gate plus one connected-A6 quotient line | **exact over GF(73), equation pinned** |
| q4_1938 | `A3+A2/MW12` | `2A4+2A3+A1/MW3` | `-P1+P2+2P3`, `P.O=1`, height `12/5` | smooth `P.O=1` saturation plus one connected-A4 quotient line | **exact over GF(73), equation pinned** |
| q4_6855 | `4A1/MW13` | not yet equation-targeted | generic `P.O=2` | generic vertical correction zero | **next live gate** |
| q4 candidate 1 | `A1/MW16` | not yet equation-targeted | generic `P.O=2` | generic vertical correction zero | pending |
| final q6 | rootless/MW17 | not yet equation-targeted | generic `P.O=4` | generic vertical correction zero | pending |

## Orbit 1222: canonical saturated module

The explicit orbit-424 parent over GF(73) is the `2I8+8I1` model. The successful orbit-1222 construction did **not** come from a guessed three- or five-dimensional `qsat` subspace.

For

```text
L2 = O + P + 2F
```

the correct raw marked-chord ambient is seven-dimensional:

```text
a = A(U)/H^2,  deg(A) <= 4,
b = B(U)/H,    deg(B) <= 1,
f = a + b*m.
```

The exact restrictions are

```text
smooth P.O=1 collision: rank 2   -> 7D to 5D
infinity I8 whole-fibre gate: rank 2
finite A7 vertical gate: rank 1
--------------------------------------
total rank 5                    -> 7D to 2D.
```

For the retained A7 orientation the kernel is

```text
((1,0,2,71,0,70,0),
 (0,1,38,69,0,0,0)).
```

The resulting child is

```text
2I7 + 3I2 + 4I1
root data = (15,90,392)
MW = 3.
```

### Exact splice to the old third-q12 CM24 child

The new orbit-1222 child is the same GF(73) elliptic fibration as the old pinned third-q12 CM24 child. Two exact PGL2/Weierstrass identifications occur:

```text
phi1(T) = (T+22)/(13T+28),  d=3=21^2,
phi2(T) = 1/(52T+62),       d=61=34^2.
```

Their relative base automorphism is

```text
psi(V) = (V+28)/(23V+72),
```

an involution in `PGL2(F_73)`. The simpler pinned transport is

```text
V = 1/(52T+62),
T = (V+20)/(18V).
```

This closes the earlier open question asking whether the two `2A6+3A1` CM24 constructions were actually the same fibration.

## q6_7774: specialization changed the horizontal section

The generic q6_7774 geometry is

```text
D.F = 2
D.O = 1
P.O = 3
MW height = 219/28
vertical support = 2 simple components
fiber twist = 0.
```

The **actual CM24 specialization** is instead

```text
D.F = 2
D.O = 0
P.O = 0
MW height = 8/7
fiber twist = 2
nearest shortest lifts = 441.
```

The horizontal section is the old pinned MW basis section `P3`. The vertical correction is

```text
one A1:  (-1)
one A6:  (-1,0,0,-1,-1,-1)
other fibres: zero.
```

The specialized child is forced by the rank-18 CM24 lattice to be

```text
A5 + 2A4 + 2A1
root data = (15,74,600)
MW = 3.
```

The lattice transport is regression-checked one step earlier by reproducing the orbit-1222 `(15,90,392)/MW3` child exactly.

### Superseded q6 diagnostics

A maximal-height `P.O=3` MW census produced four `+/-` pairs, but these are **not** the specialization of the generic q6_7774 horizontal: the actual specialized section is `P3`, with `P.O=0` and height `8/7`.

Likewise, imposing independent vanishing on every listed A6 exceptional component overconstrains the pencil. It gives rank `2` from the A6 block; together with the A1 gate this leaves only one section in a four-dimensional ambient, contradicting `h0(D)=2`.

### Connected ADE quotient-line rule

The reusable correction is structural: a **connected vertical ADE divisor contributes one resolved quotient-line condition**, not one independent condition per exceptional component.

For q6_7774 the componentwise A6 rows span two directions of the form

```text
(1,r,r^2,c1),
(1,r,r^2,c2).
```

Scanning the 73 quotient lines inside that span, together with the exact intrinsic I2/A1 row, gives `876` cases. Only `24` square classes reduce to quartics, and exactly `4` symmetry-related cases hit `(15,74,600)/MW3`.

A pinned representative is

```text
sign = +1
I7 fibre = 6
I2 fibre = 5
A6 quotient residue c7 = 58
A1 residue c2 = 62
a(V) = 4V + 64
m = T(V-6)(V-5) - (4V+64).
```

The pinned short Weierstrass model is

```text
A(T) =
46*T^8 + 5*T^7 + 16*T^6 + 44*T^5 + 6*T^4
+ 13*T^3 + T^2 + T,

B(T) =
54*T^12 + 58*T^11 + 48*T^10 + 16*T^9 + 42*T^8
+ 67*T^7 + 25*T^6 + 19*T^5 + 27*T^4 + 45*T^3
+ 61*T^2 + 44*T + 49.
```

Its reducible fibres are

```text
I2 at T=46,47
I5 at T=14,42
I6 at infinity,
```

giving exactly `2A1+2A4+A5`.

## q4_1938: specialized marking and equation

Transporting the actual generic q4_1938 divisor into the pinned q6_7774 CM24 frame gives

```text
D.F = 2
D.O = 0
P.O = 1
MW height = 12/5
shortest norm = 6
fiber twist = 1.
```

The vertical correction is supported on one A4 only:

```text
(-1,-1,-1,0),
L1 = 3,
max coefficient = 1.
```

The specialized child is predicted exactly as

```text
2A4 + 2A3 + A1
root data = (15,66,800)
MW = 3.
```

### Specialized parent MW basis

Component-group labels are computed basis-independently: `root_component_data()` can return arbitrary integral root-lattice bases, so the code uses a Smith discriminant generator and selects the minimal endpoint-weight class rather than assuming a Cartan basis.

The exact optimal MW data are

```text
components = (A1,A1,A4,A4,A5)
height =
  [4/15,   0, 2/15]
  [   0, 2/5, 1/10]
  [2/15,1/10, 7/15]

basis profiles =
  P1: (0,0,3,2,4)
  P2: (0,1,1,1,3)
  P3: (1,1,3,0,2)

P1.O=P2.O=P3.O=0,
P1.P2=P1.P3=P2.P3=0.
```

The q4_1938 horizontal is

```text
H = -P1 + P2 + 2P3
MW coordinates = (-1,1,2)
profile = (0,1,4,4,3)
height = 12/5
P.O = 1.
```

Node-constrained polynomial searches on the pinned q6 equation recover the three basis roles. Exact calibration of the two I5 groups and the I6 group at infinity leaves four abstract basis markings; all give one `+/-` pair for `H`.

A pinned sign has

```text
X_H =
(60 + 40T + 34T^2 + 58T^3 + 48T^4 + 30T^5 + 3T^6)
/ (19 + 52T + T^2),

Y_H =
(24 + 15T + 57T^3 + 28T^4 + 65T^5 + 15T^6)
/ (56 + 57T + 5T^2 + T^3).
```

### Connected A4 quotient-line scan: passed

For

```text
L = O + H + F,
L^2 = 2,
h0(L) = 3,
```

a raw five-dimensional chord ambient is cut by the exact `P.O=1` smooth-collision congruence with rank `2`, giving the expected three-dimensional `H0(L)`. The connected A4 correction then contributes one quotient-line condition.

The exhaustive scan is

```text
2 signs x 2 I5 fibres x 73 residues = 292 cases.
```

Exactly `4` cases reduce to quartics, and all `4` hit the hard lattice target `(15,66,800)/MW3`. Their keys are

```text
(+1, I5=14, residue=34)
(+1, I5=42, residue=59)
(-1, I5=14, residue=39)
(-1, I5=42, residue=14).
```

A pinned representative is

```text
sign = +1
I5 fibre = 14
quotient residue = 34
evaluation row = (1,14,50,43,46)
kernel =
  (1,0,51,40,23)
  (0,1,36,27,48).
```

Its child fibres are

```text
I5 at S=-1 and S=-65
I4 at S=-36 and S=-50
I2 at S=-46
one rational I1
one cubic factor giving 3 I1
infinity smooth,
```

hence exactly

```text
2A4 + 2A3 + A1
root data = (15,66,800)
MW = 3.
```

The pinned quartic is

```text
T^4
+ ((20*S^2 + 22*S + 50)/(S^2 + 51*S + 50))*T^3
+ ((53*S^4 + 63*S^3 + 65*S^2 + 45*S + 55)
   /(S^4 + 29*S^3 + 63*S + 18))*T^2
+ ((19*S^4 + 8*S^3 + 49*S^2 + 12*S + 47)
   /(S^4 + 29*S^3 + 63*S + 18))*T
+ (6*S^4 + 27*S^3 + 44*S^2 + 48*S + 19)
  /(S^4 + 29*S^3 + 63*S + 18),
```

with twist

```text
(6*S^4 + 28*S^3 + 13*S + 35)
/ (S^4 + 71*S^3 + 38*S^2 + 36*S + 32).
```

## Current live gate: q4_6855

The CM24 equation frontier is now **q4_6855**. Its generic geometry is unusually favorable:

```text
generic source = A3+A2/MW12
generic child = 4A1/MW13
D.F = 2
P.O = 2
MW height = 19/3
fiber twist = 0
vertical correction = 0.
```

Do not assume those generic horizontal data survive CM24. The next step is the same specialization-first workflow that succeeded for q6_7774 and q4_1938:

1. transport the actual q4_6855 divisor through the full rank-18 CM24 lattice;
2. chamber-reduce it and identify the true specialized horizontal section and any new vertical correction;
3. predict the special child root data before equation search;
4. reconstruct the horizontal on the pinned q4_1938 equation;
5. compile any connected ADE correction as quotient-line/module conditions.

After q4_6855, continue to q4 candidate 1 and the final rootless q6.

The generic characteristic-zero algebraization from orbit 1222 onward remains open. These CM24 modular certificates are a development and verification scaffold, not a substitute for that final lift.
