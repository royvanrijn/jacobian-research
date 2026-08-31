# Fixed-corridor reverse equation lift from the q12 endpoint

<!-- status-consumer: EC-K3-H3-FIXED-REVERSE-A1-2A1-QQ 26bce707a77972c4 -->
<!-- status-consumer: EC-K3-H3-FIXED-REVERSE-4A1-QQ 8f3863b630d27e16 -->

Date: 2026-08-26

Status: exact characteristic-zero construction for the last four arrows only;
the next q4/orbit52 section data are exact over `QQ`, but its Riemann--Roch
pencil and `5A1` Jacobian remain open.

## Result

The exact q12/orbit5867 rootless Jacobian and its determinant-one section-basis
identification with pinned R17 provide an equation source for the historical
fixed corridor.  Inverting the certified fixed transitions gives exact
equations and markings for

```text
2A1/MW15 --q4 orbit981--> A1/MW16 --q6 orbit2247--> rootless/MW17.
```

The same reverse construction now also proves

```text
3A1/MW14 --q4 orbit498--> 2A1/MW15.
```

The construction is executed from right to left.  It does not alter the
selected lattice transitions or claim that the earlier fixed-corridor arrows
are equation-explicit.

## Rootless to A1

In pinned rootless coordinates the reverse A1 fibre is `(4,2,w)` with
`w^2=16`.  Therefore

```text
F_A1 = O + P - 2F,   P.O=6.
```

The exact endpoint-basis word for `P` is

```text
(0,-1,-2,-1,1,0,1,-1,-2,-2,1,-2,1,0,-1,1,1).
```

Exact group law gives `(deg X,deg Y,deg Z)=(16,24,6)`.  The smooth chord
module has ambient dimension 14, condition rank 12, and `h0=2`.  Its binary
quartic has degree four.  The minimal Jacobian has degrees `(8,12,24)`, fibres
`I2+22I1`, smooth infinity, Euler number 24, root rank one, and ADE `A1`.

The prescribed fixed A1 zero is another exact endpoint word.  Its restriction
to the quartic has degree one and the pointed generalized model satisfies

```text
81*A_pointed=A_child,   729*B_pointed=B_child.
```

The transported effective A1 root contracts to the unique I2 support.

## A1 to 2A1

The reverse orbit981 fibre again has the form

```text
F_2A1 = O + P - 2F,   P.O=6.
```

Among all 2,622 signed norm-four pinned-R17 vectors, exactly 21 define
rootless sections meeting the new A1 fibre once.  Their A1 MW tails have rank
15.  An exact Smith solve expresses the orbit981 horizontal as the seven-term
word with coefficients

```text
(1,-1,1,-1,1,1,-2)
```

on selected degree-one pointed images.  Exact group law again gives degrees
`(16,24,6)`.  The smooth chord calculation is `14 -> rank 12 -> h0 2`; its
minimal Jacobian has degrees `(8,12,24)`, fibres `2I2+20I1`, smooth infinity,
Euler number 24, root rank two, and ADE `2A1`.

The prescribed 2A1 zero is the nonidentity component of the old A1 fibre.
An exhaustive exceptional-conic calculation modulo 131 selects the positive
exact quartic ordinate.  The characteristic-zero pointed invariants satisfy
the same `81/729` identities.  Two exact transported root sections contract
to the two distinct rational I2 supports, binding both A1 roots to the
equation.

## 2A1 to 3A1

For the reverse orbit498 fibre, exact short-vector enumeration in twice the A1
height lattice finds 78 sections that have degree one for the 2A1 pencil.
Their 2A1 Mordell--Weil tails have full rank 15.  Integral Smith words, two
successive pointed-quartic transports, and exact group law construct

```text
F_3A1 = O + P - 2F,   P.O=6,
```

with `(deg X,deg Y,deg Z)=(16,24,6)`.  The same calculation constructs two
effective horizontal 3A1 root curves and determinant-one bases on both sides
of the transition.

The smooth chord system is again `14 -> rank 12 -> h0 2`.  Its degree-four
binary quartic gives a minimal Jacobian of degrees `(8,12,24)` with fibres
`3I2+18I1`, smooth infinity, Euler number 24, root lattice `3A1`, and MW rank
14 using the independently proved geometric Picard rank 19.  The old
nonidentity A1 component selects quartic sign `-1` in the exhaustive mod-131
exceptional-conic gate, and the exact pointed invariants satisfy

```text
81*A_pointed=A_child,   729*B_pointed=B_child.
```

The two horizontal roots contract to two distinct rational I2 supports; the
third marked 3A1 root is the remaining vertical component.  Thus both the
prescribed zero and all three root components are bound to the equation.

## Method and cost

The calculations use:

- integral pinned-marking composition;
- exact elliptic group law;
- degree-one pointed-quartic maps;
- a 12-by-14 integer kernel after clearing one global denominator;
- bivariate gcd square-part removal;
- squarefree univariate discriminant gates;
- one exhaustive finite-field exceptional-conic sign test.

No Groebner basis, surface elimination, or nonlinear characteristic-zero
section solve is used.  The first two replays take seconds; the 3A1 horizontal
takes about two minutes.  The largest reported rational coefficient sizes are
79,273 bits for the A1 Jacobian, 217,193 bits for the 2A1 Jacobian, and 568,869
bits for the 3A1 Jacobian.  The corrected 4A1 Jacobian has maximum reported
coefficient size 2,512,351 bits.

## Replay

The commands and artifact paths are listed in
[`scripts/README.md`](scripts/README.md#fixed-corridor-reverse-lift-from-the-q12o5867-endpoint).
The canonical local outputs are:

- `artifacts/local/elkies-k3/fixed-final-a1-horizontal-from-q12-endpoint-qq.json`;
- `artifacts/local/elkies-k3/fixed-final-a1-reverse-rr-qq.json`;
- `artifacts/local/elkies-k3/fixed-final-a1-reverse-pointing-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-2a1-horizontal-from-a1-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-2a1-rr-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-2a1-pointing-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-3a1-horizontal-from-2a1-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-3a1-rr-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-3a1-pointing-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-4a1-horizontal-from-3a1-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-4a1-physical-nef-audit.json`;
- `artifacts/local/elkies-k3/fixed-reverse-4a1-rr-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-4a1-pointing-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-4a1-compact-crossratio-qq.json`;
- `artifacts/local/elkies-k3/fixed-reverse-5a1-compact-sections-crt-qq.json`.

## Boundary

This proves the last four equation arrows and their markings.  For q4/orbit114
the stored abstract fibre was not in the physical 3A1 chamber: its six old-I2
component degrees were `(6,4,-1;-4,-2,3)`.  Four exact affine-Weyl reflections
give the q24 representative

```text
(12,2,12,-5,0,-1,2,2,-2,-3,3,-3,-2,3,-3,-3,3,-3,0).
```

All 77 constructed sections are nonnegative, the exact all-section minimum is
zero, and the complete finite horizontal-wall list is empty.  The corrected
identity is `D=O+P-C1-9F`.  Its compact chord calculation has dimensions
`45 -> 3 -> 2`; the last condition is the exact opposite-point restriction on
the first old nonidentity I2 component.  Exact square stripping gives a binary
quartic, whose minimal Jacobian has `4I2+16I1`.  The third old nonidentity I2
component points the zero, two exact horizontal roots contract to distinct
child I2 supports, and two roots remain vertical.  The reflected full basis
has determinant `-1` and sends the abstract zero exactly to the effective
class `-e4`.

The rejected 57-to-15-to-2 genus-two quotient remains a valid negative
experiment: it compiled a different, nonphysical marked divisor.

For the next q4/orbit52 edge, six exact affine-Weyl reflections give the
physical q18 fibre

```text
(9,2,-3,-6,4,0,2,-3,-2,0,0,0,1,-1,3,0,1,0,0).
```

Its old-I2 component degrees are `(1,2,0,1;1,0,2,1)`, the exact all-section
minimum is zero, the complete finite horizontal-wall list is empty, and the
full `5A1` marking is unimodular.  The prescribed zero is exactly the fourth
old nonidentity component.  A cancellation-preserving fibrewise Abel-word
calculation over `GF(167)` constructs the required old `P.O=15` section and
the three horizontal root sections (`P.O=1,1,3`).  All four interpolate,
satisfy the current `4A1` Weierstrass equation identically, pass holdouts, and
have full-rank coefficient charts after their exact I2-node incidences are
added.

Blind Hensel lifting was then replaced by a compact exact model.  The
cross-ratio sending three rational old `I2` supports to `0`, `1`, and infinity
sends the fourth to `923/3815`; the induced Weierstrass change is an exact
`QQ` isomorphism, not a twist.  It lowers the maximum `A,B` coefficient size
from about 2.5 million bits to 215 bits while preserving `4I2+16I1`.
Fibrewise Abel interpolation at 110 good 26-bit primes, followed by CRT in
this compact normalization, reconstructs all four sections exactly over
`QQ(t)`.  Their compact `(x,y)` numerator/denominator degrees are

```text
P.O=15: (34/30, 50/45)
P.O= 1: ( 6/ 2,  9/ 3)
P.O= 1: ( 6/ 2,  8/ 3)
P.O= 3: (10/ 6, 14/ 9).
```

Literal substitution, replay at every CRT prime, and the withheld prime 167
all pass.  Maximum reconstructed rational size is 816 bits.  Thus the q52
horizontal and three horizontal roots are now exact over characteristic zero;
the remaining completion gate is the exact `D=O+P-C2-6F` two-plane, its
binary quartic/minimal `5A1` Jacobian, and the prescribed-zero/component
pointing.  No Groebner calculation or further large-normalization Hensel lift
is indicated.
