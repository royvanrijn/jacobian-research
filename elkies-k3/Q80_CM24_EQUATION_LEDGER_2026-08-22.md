# Q80 CM24 equation ledger — 2026-08-22 evening

## Scope

This ledger records the equation-level status of the **secondary/fallback Q80 low-q corridor**. The generic lattice route itself is already complete to a new rootless `MW17` frame. The purpose of this ledger is to separate:

1. generic lattice/chamber certificates;
2. CM24 lattice specialization certificates;
3. CM24 equation certificates;
4. work that is still pending.

The corrected H3 source-polarization route remains primary. Nothing in this file upgrades Q80 above H3.

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

| stage | generic child | CM24 special child | special horizontal | special vertical / module | equation status |
|---|---|---|---|---|---|
| q6 escape | `D7+D4/MW6` | `D8+D6+2A1/MW2` | old q12 section `Q=P1+3P2` | exact E6 jet, binary quartic | **exact over `QQ(sqrt(-6))`** |
| orbit 424 | `A6+A4/MW7` | `2A7/MW4` | rational 2-torsion | exact degree-two chord | **exact over `QQ(sqrt(-6))`** |
| orbit 1222 | `A6+A3/MW8` | `2A6+3A1/MW3` | `P.O=1`, height `25/8` | canonical raw saturated chord module | **exact over GF(73); splice to old q12 child proved** |
| q6_7774 | `A4+A2+A1/MW10` | `A5+2A4+2A1/MW3` | old-model `P3`, `P.O=0`, height `8/7` | one A1 gate plus one connected-A6 quotient line | **exact over GF(73), equation pinned** |
| q4_1938 | `A3+A2/MW12` | `2A4+2A3+A1/MW3` | `-P1+P2+2P3`, `P.O=1`, height `12/5` | one connected A4 with coefficients `(-1,-1,-1,0)` | horizontal pinned; **final A4 quotient-line scan pending** |
| q4_6855 | `4A1/MW13` | not yet equation-targeted | generic `P.O=2` | generic vertical correction zero | pending |
| q4 candidate 1 | `A1/MW16` | not yet equation-targeted | generic `P.O=2` | generic vertical correction zero | pending |
| final q6 | rootless/MW17 | not yet equation-targeted | generic `P.O=4` | generic vertical correction zero | pending |

## Orbit 1222: canonical saturated module

The explicit orbit-424 parent over GF(73) is the `2I8+8I1` model. The successful orbit-1222 construction did **not** come from a guessed three- or five-dimensional `qsat` subspace.

The correct raw marked-chord ambient for

```text
L2 = O + P + 2F
```

is seven-dimensional:

```text
a = A(U)/H^2,  deg(A) <= 4,
b = B(U)/H,    deg(B) <= 1,
f = a + b*m.
```

The exact restrictions are:

```text
smooth P.O=1 collision: rank 2   -> 7D to 5D
infinity I8 whole-fibre gate: rank 2
finite A7 vertical gate: rank 1
--------------------------------------
total: rank 5                  -> 7D to 2D.
```

For the retained A7 orientation the kernel is

```text
((1,0,2,71,0,70,0),
 (0,1,38,69,0,0,0)).
```

The resulting quartic has the expected specialized child

```text
2I7 + 3I2 + 4I1
root data = (15,90,392)
MW = 3.
```

### Splice to the old third-q12 CM24 child

The new orbit-1222 child is exactly the same GF(73) elliptic fibration as the old pinned third-q12 CM24 child. Two exact PGL2/Weierstrass identifications occur:

```text
phi1(T) = (T+22)/(13T+28),  d=3=21^2,
phi2(T) = 1/(52T+62),       d=61=34^2.
```

Their relative base automorphism on the old model is

```text
psi(V) = (V+28)/(23V+72),
```

and is an involution in `PGL2(F_73)`. The simpler pinned transport is

```text
V = 1/(52T+62),
T = (V+20)/(18V).
```

This closes the old open question asking whether the two `2A6+3A1` CM24 constructions were actually the same fibration.

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

The **actual CM24 specialization** is very different. Exact rank-one NS transport and chamber reduction give

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

This lattice prediction is independently regression-checked by replaying the preceding orbit-1222 specialization, which reproduces `(15,90,392)/MW3` exactly.

### Connected ADE quotient-line lesson

A first resolved attempt imposed an independent evaluation condition on every listed A6 component. That produced rank `2` from the A6 block, and together with the A1 gate gave rank `3` on a four-dimensional ambient. Every one of 24 sign/fibre/orientation cases therefore left only one section, contradicting `h0(D)=2`.

The correction is structural: a **connected vertical ADE divisor contributes one resolved quotient-line condition**, not one independent condition per exceptional component.

For q6_7774 the componentwise A6 rows span two evaluation directions of the form

```text
(1,r,r^2,c1),
(1,r,r^2,c2).
```

Scanning the 73 possible quotient lines inside that span, together with the exact intrinsic I2/A1 row, gives `876` exact cases. Only `24` square classes reduce to quartics, and exactly `4` symmetry-related cases hit the required `(15,74,600)/MW3` target.

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

with the remaining discriminant factors simple, giving exactly `2A1+2A4+A5`.

## q4_1938: exact specialized marking

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

A basis-independent discriminant-group labeling was required here: `root_component_data()` supplies arbitrary integral bases of the root sublattices, so component labels must not assume a Cartan/simple-root basis. Using Smith generators and selecting the minimal endpoint-weight class gives the exact optimal MW data

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

### Equation-side horizontal reconstruction

Node-constrained polynomial searches on the pinned q6_7774 equation recover the three `P.O=0` basis roles. Exact calibration of the two I5 component groups and the I6 group at infinity leaves four abstract basis markings, all related by sign/fibre orientation conventions. Exact height polarization with `Pi-Pj` verifies the target Gram. All four markings produce the same single `+/-` pair for `H`.

A pinned sign has

```text
X_H =
(60 + 40T + 34T^2 + 58T^3 + 48T^4 + 30T^5 + 3T^6)
/ (19 + 52T + T^2),

Y_H =
(24 + 15T + 57T^3 + 28T^4 + 65T^5 + 15T^6)
/ (56 + 57T + 5T^2 + T^3).
```

The other candidate is exactly `-H`.

## Current live gate: q4_1938 connected A4 quotient line

For the pinned horizontal `H`, the natural ambient is

```text
L = O + H + F,
L^2 = 2,
h0(L) = 3.
```

A raw five-dimensional chord ambient

```text
A(T)/h^2, deg A<=3,
B/h * m
```

is cut by the exact `P.O=1` smooth-collision congruence with rank `2`, giving the expected three-dimensional `H0(L)`. The connected A4 correction must then contribute **one quotient-line condition**, leaving the required two-dimensional pencil.

The current exact scan is

```text
2 signs x 2 I5 fibres x 73 quotient residues = 292 cases,
```

with the hard target

```text
root data = (15,66,800),
MW = 3.
```

**Status: pending.** Do not record a q4_1938 equation as certified until this scan returns an exact target hit.

## Open after q4_1938

Once q4_1938 is pinned, continue the same specialization-first discipline:

1. transport q4_6855 through the full CM24 lattice before guessing an equation marking;
2. identify the actual specialized horizontal section and vertical divisor;
3. compile connected ADE corrections as quotient-line/module conditions rather than independent component evaluations;
4. continue to the `4A1`, `A1`, and final rootless stages;
5. only after the modular corridor is stable, decide whether a generic characteristic-zero lift is competitive with the corrected H3 route.

The generic characteristic-zero algebraization from orbit 1222 onward remains open. The CM24 modular certificates are a development and verification scaffold, not a substitute for that final lift.
