# Exact Picard-20 q=8 chord fibration (2026-08-21)

## Status and scope

The q=8 neighbor of the rational discriminant-43 K3 has now been executed
geometrically and exactly.  Its fiber configuration is

```text
I3 + I4 + I5 + I8 + 4 I1,
ADE = A2 + A3 + A4 + A7,
MW rank = 2.
```

This is a Picard-rank-20 Noether--Lefschetz boundary construction.  Its fiber
class contains the extra section `S`; it is not a q=8 neighbor in the generic
determinant-948 rank-19 lattice.  The generic q=8 hit from the `A13+A1`
continuation is separately known to return to the old `A10` frame.  Neither
fact is an obstruction to using the present surface as an exact CM boundary
anchor, but the two q=8 phenomena must not be conflated.

## Chord and unique degree-one pencil

On the explicit `IV*+I0*+2I3+I2` model, let

```text
z = (y-S.Y)/(x-S.X).
```

The q=8 divisor has old-fiber degree two and horizontal part `3O-S`, so this
is the correct chord coordinate.  The raw line discriminant is a genus-two
quintic in the old base coordinate `t`.

The complete two-plane scan in

```text
<1,t,z,tz>
```

modulo constant base `PGL2` contains 293,090 planes over `GF(23)`.  Exactly one
has the target root system.  Its reduced rows are

```text
((1,0,11,0),(0,1,5,0)),
W = (t+5z)/(1+11z) mod 23.
```

The complete `GF(11)` scan also has exactly one target, with coefficients
`(6,10)`.  Both reductions reconstruct from the rational `A1` branch
`(t,z)=(49/25,-17/3)` to

```text
W = (t+(147/425)z)/(1+(3/17)z).
```

After substitution, the completed quintic contains the exact square factor
`(t-49/25)^2`.  Its odd part is cubic.  Accounting for the constant factor in
the raw line discriminant gives a pointed cubic `v^2=2*C_W(t)`, rather than
the untwisted Jacobian.

## Compact integral model

Apply the base change and Weierstrass scale

```text
U = -(14/3)*(W+21/50)/(W-49/25),
W = (147/25)*(U-1)/(3U+14),
(x,y) -> ((6250/441)^2*(3U+14)^4*x,
          (6250/441)^3*(3U+14)^6*y).
```

Set

```text
q = 9U^2 + 84U + 49,
L = 1815156 U(U-1)^2.
```

The surface becomes the compact integral equation

```text
y^2 = x^3 + (q^2 + 5292 U(U+7)) x^2 + 2Lq x + L^2.
```

Its reducible fibers are

```text
U=28/3: I3,
U=1:    I4,
U=0:    I5,
U=infinity: I8.
```

The remaining discriminant factor is the irreducible quartic

```text
243U^4 + 8316U^3 + 156114U^2 - 2210292U + 45619,
```

giving four `I1` fibers.

## Saturated Mordell--Weil basis

Two particularly small sections are

```text
P = (0, -L),

Q.x = -115248(U-1),
Q.y = -777924(U-28/3)(U-49/9)(U-1).
```

Their sum is

```text
(P+Q).x = -12348U(U-1),
(P+Q).y = -111132U(U-1)(U-28/3)(U+7/3).
```

In the fiber order `(I3,I4,I5,I8)`, exact singular-reduction orders and the
local `A_n` resolution depths give signed component profiles

```text
P = (0,2,1,3),
Q = (1,1,0,3).
```

All three intersections `P.O`, `Q.O`, and `P.Q` vanish on the resolved K3.
Shioda's formula gives

```text
height(P,Q) = [13/40  -3/8]
              [-3/8  17/24],

det(height) = 43/480.
```

This is unimodularly equivalent to the independently recovered saturated
height Gram

```text
(1/120) [34  6]
        [ 6 39].
```

Indeed `P=e2` and `Q=-e1-e2`.  Thus the two displayed rational sections are a
full saturated MW basis, not merely an independent subgroup.

## CM-43 Kumar identification

The Gross vector

```text
beta = -11*i-j-5*k,  norm(beta)=43
```

transports to `(169,167,-128)` in the determinant-948 transcendental lattice.
It has norm `-40764`, divisibility `948`, and perpendicular Gram

```text
[22 1]
[ 1 2],
```

of determinant 43.  The primitive closure of the `H2` Kumar `E7+E8/MW2`
frame retains the `E7+E8` roots and has MW rank three with reduced height Gram

```text
[ 5/2 -1/2  -1]
[-1/2  5/2   0]
[  -1    0   4].
```

Hence the explicit q=8 surface is the `Delta=-43` CM point on the recovered
`X_0^6(79)` Gross/K3 model.  The clean reverse-construction coordinates are
therefore the two height-`5/2` directions in this Kumar frame, not a generic
q=8 deformation of the semistable pencil.

The same CM point is now located directly in the Humbert-8 plane:

```text
r = -1225/722,
s = -93312/442225,
z = +/-(11664/6859)*sqrt(-43).
```

Substitution gives the exact `E7+E8` Kumar equation.  Two rational polynomial
sections are verified coefficient by coefficient and have height Gram
`[[5/2,-1/2],[-1/2,5/2]]`.  Direct point counts at fourteen good primes match
the weight-three CM-43 coefficients, including zero at the tested inert
primes.  This is a strong independent identification fingerprint; it is not
being substituted for the still-missing explicit birational inverse neighbor
between the semistable and Kumar presentations.

## Direct comparison with the Kumar q=8 and q=60 neighbors

The semistable q=8 fibration is not a one-step q=8 neighbor of the CM-43
Kumar `E7+E8/MW3` closure.  Exact `qfminim` gives 1,421,331,656 signed
norm-16 vectors in that frame.  Decomposing them into the saturated MW
projection and dominant `E7/E8` Dynkin labels reduces the complete shell to
303 Weyl orbits.  Of these, 292 give primitive `(a,b)=(2,4)` neighbors.  None
has root invariants `(rank,count,det)=(16,94,480)`, so none is integrally
isometric to the pinned `A7+A4+A3+A2` frame.

This is stronger than a bounded raw-vector search: it is a complete orbit
classification.  It shows that the repeated q=8 is a structural boundary
phenomenon, but not the missing direct Kumar-to-semistable transition.

The same complete orbit table nevertheless gives a useful exact
factorization of the q=60 CM automorphism.  After imposing the genuine
`(P3,Q79)` marking, the best primitive q=8 child already has root frame
`E7+E8` and MW rank three.  In that child the q=60 fiber has `(a,b)=(3,3)`,
so the return to `E7+E8/MW3` is a q=9 neighbor:

```text
E7+E8/MW3 --q=8--> E7+E8/MW3 --q=9--> E7+E8/MW3.
```

The first witness and full second-fiber coordinates are pinned by
[`scripts/verify_kumar_cm43_q8_q9_factor.sage`](scripts/verify_kumar_cm43_q8_q9_factor.sage).
This intermediate is not the semistable `A7+A4+A3+A2` q=8 fibration.

The optimal q=8 fiber has marked horizontal projection `P1-2*P2`, of height
`29/2`; it is not the level-79 direction and is NL-special.  Exact group law
on the CM-43 Kumar equation gives `R.O=6` and
`x(R)=N16/h6^2,y(R)=N24/h6^3`, providing a much smaller equation ansatz than
the pole-58 level-79 section.

The NL qualification is exact and rules out deforming this q=8 fiber along
the generic rank-19 family.  For the marked height Gram, the primitive vector
orthogonal to `P3=(0,0,1)` and `Q79=(4,-5,1)` is `W=(116,92,29)`, of height
`40764`, and

```text
P1-2*P2 = -1/4*P3 + 27/79*Q79 - 1/316*W.
```

The residual height is `129/316=9*(43/948)`.  In the full glue-211 closure
the q=8 witness has CM-summand coordinate `1/316`, whereas q=60 has
coordinate zero.  Therefore q=8 is a low-pole CM-43 boundary factorization;
q=60 is integral in the generic rank-19 lattice, but the complete chamber
calculation below shows that its CM-43 specialization is not a movable pencil.
See
[`scripts/verify_cm43_q8_generic_membership.sage`](scripts/verify_cm43_q8_generic_membership.sage).

There are also exactly ten q=8 orbits with the Humbert root frame `D9+E7`.
For the cheapest root-equivalent orbit the same q=60 fiber needs q=289.  The
actual H2 Humbert-8 two-neighbor is more sharply distinguished by horizontal
projection `+P3`; for that marked class the direct second lattice move has
q=2,466,464.  Thus the explicit `D9+E7` equation is an excellent known
birational chart, but it is not the arithmetic q=8 step in the shortest q=60
factorization.

The q=60 `(a,b)=(5,12)` neighbor behaves differently.  At CM-43 it has
`E7+E8` roots, MW rank three, and is integrally isometric to the original
Kumar closure.  In the Kumar CM Neron--Severi basis its exact fiber class is

```text
(5,12,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,0).
```

The generic q=60 basis specializes with pole profile `(0,1,1)` at this point,
but pole count alone misses fixed CM sections.  In the explicit divisor basis,
the O/E7 reduction gives `D1=Q79+4O-43F`; then `4(P1-P2)` and `P3-P2` are
fixed, leaving exactly `F`.  The q=8 class behaves the same way: after its
O/E7 reduction, `-P2` and `P1-P2` are fixed and again leave `F`.  Thus

```text
D60_raw = F + initial_fixed + 4*(P1-P2) + (P3-P2),
D8_raw  = F + initial_fixed + (-P2) + (P1-P2).
```

Both claimed CM pencils therefore collapse to the old fibration.  This is
certified by
[`scripts/verify_cm43_marked_divisor_transport.sage`](scripts/verify_cm43_marked_divisor_transport.sage).
CM-43 remains a marking and chamber boundary, not an equation seed for q=60.
The construction pivot is to the noncollapsed CM-24 choices: q60 `(4,15)`
retains two generic MW directions, while q80 `(4,20)` retains all three and
has two certified formal branches.

The equation-level marking is now fixed as well.  In the explicit CM basis
`(P1,P2,P3)`, with `P3` the height-four section obtained from the second
point of the `D9+E7` quartic, the level-79 direction is

```text
Q79 = 4*P1 - 5*P2 + P3.
```

It has height `237/2`, is orthogonal to `P3`, and exact group law gives pole
58.  Auditing all eight normalized CM closure signs shows that glues `211`
and `737` preserve the ordered marking `(P3,Q79)`.  The older glue `53`
instead produced horizontal coordinates `(5,-4,2)` in an unmarked reduced
basis.  Consequently the q=60 divisor must use horizontal coordinates
`(4,-5,1)` on the explicit Kumar equation.

## Reproduction

Run the exact characteristic-zero construction and section-lattice verifier:

```bash
sage elkies-k3/scripts/derive_picard20_q8_s_chord.sage
```

Replay the two complete finite Grassmann scans:

```bash
sage elkies-k3/scripts/search_picard20_q8_s_grassmann_gf23.sage --p 11

sage elkies-k3/scripts/search_picard20_q8_s_grassmann_gf23.sage \
  --p 23 --start 0 --limit 73273
sage elkies-k3/scripts/search_picard20_q8_s_grassmann_gf23.sage \
  --p 23 --start 73273 --limit 73273
sage elkies-k3/scripts/search_picard20_q8_s_grassmann_gf23.sage \
  --p 23 --start 146546 --limit 73272
sage elkies-k3/scripts/search_picard20_q8_s_grassmann_gf23.sage \
  --p 23 --start 219818 --limit 73272
```

Replay the CM transport and Kumar primitive closure:

```bash
sage elkies-k3/scripts/transport_cm_delta3_to_k3.sage \
  --target 43 \
  --out artifacts/local/elkies-k3/cm-delta43-k3-vector.txt

sage elkies-k3/scripts/verify_cm43_humbert8_anchor.sage

sage elkies-k3/scripts/classify_kumar_cm_frame_extensions.sage \
  --search-max-q 60

sage elkies-k3/scripts/classify_kumar_cm43_q8_orbits.sage
```

The Grassmann searches are exhaustive only in the declared degree-one chart.
The characteristic-zero script is the exact certificate for the recovered
surface; it does not produce the generic rank-17 rootless fibration or a new
rank-greater-than-30 elliptic curve by itself.
