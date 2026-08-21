# Exact q80 to rootless MW17 path

## Status

This note records a bounded discovery followed by an exact lattice replay.
The shell searches were capped and do not prove path minimality. The retained
six-step path, however, is certified with exact integral arithmetic by
[`scripts/verify_q80_to_rootless_path.sage`](scripts/verify_q80_to_rootless_path.sage).

The path starts from the **generic determinant-948 q80 frame**
[`data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt`](data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt),
not from the discriminant-43 Picard-20 specialization. Its terminal frame is
rootless and explicitly integrally isometric to
[`data/lattice/rank17_gram.txt`](data/lattice/rank17_gram.txt).

## Path

Coordinates in each row refer to the positive frame Gram produced by the
preceding row. For `NS=U+(-M)`, the primitive isotropic class is
`f=(a,b,v)` with `a*b=q(v)=v*M*v/2`.

| step | q | `(a,b)` | roots | MW |
|---:|---:|---:|---|---:|
| 0 | - | - | `E6+D5+A3` | 3 |
| 1 | 4 | `(2,2)` | `D9+A4` | 4 |
| 2 | 4 | `(2,2)` | `D7+D5` | 5 |
| 3 | 12 | `(3,4)` | `A5+A3+3A1` | 6 |
| 4 | 12 | `(3,4)` | `4A1` | 13 |
| 5 | 4 | `(2,2)` | `A1` | 16 |
| 6 | 6 | `(2,3)` | rootless | 17 |

The exact vectors and root rank/count/determinant gates are stored in
[`data/fibrations/kumar_q80_to_rootless_path.tsv`](data/fibrations/kumar_q80_to_rootless_path.tsv).
The terminal-frame isometry is pinned in
[`data/fibrations/kumar_q80_rootless_frame_isometry.txt`](data/fibrations/kumar_q80_rootless_frame_isometry.txt).

The first two q4 transitions raise MW rank `3 -> 4 -> 5`; the next q12
reaches MW6. The second q12 removes seven further root directions at once,
reaching `4A1/MW13`. A q4 step leaves one `A1`, and q6 removes the final root.

### This is not a source-marking shortcut

Transporting the two old Kumar marked frame directions through this same
chain gives the following exact Mordell--Weil quotient heights:

| step | height-4 direction | level-79 direction |
|---:|---:|---:|
| 0 | `4` | `120` |
| 1 | `3156/5` | `148039` |
| 2 | `34964` | `33078501/4` |
| 3 | `3623636/3` | `1721011699/6` |
| 4 | `1632452508` | `775184468553/2` |
| 5 | `1171490429116` | `278133960303042` |
| 6 | `6419904449054724` | `1524206179867493990` |

These are orthogonal-projection heights in each frame quotient; they do not
claim that the same displayed curve remains a section after every neighbor.
They do prove that the small-q rank-raising path makes the *marked source
directions* rapidly more complicated, not less.  It should therefore be used
only after the genuine family is known, not to recover the missing `H237`
orientation.

Nor is the compact q60 chart cheaply adjacent to q80.  Transporting all five
q60 presentations into the q80 NS basis gives direct neighbor norms

```text
source (2,30), (3,20), (4,15), (5,12), (6,10)
direct q   1600,    400,    400,    784,     1600.
```

Thus the best direct q80-to-q60 move still has `q=400`.  Executing it would
discard the principal advantage of both compact charts.  The calculations are
printed by `verify_q80_to_rootless_path.sage` and
`classify_kumar_cm_frame_extensions.sage --print-q60-in-q80`.

## Composite transport

Each neighbor reconstruction returns an exact unimodular 19-by-19 basis
transport. Their ordered product, followed by the pinned terminal isometry,
expresses the pinned rootless basis `[U,rank17_gram]` in the initial q80
Neron--Severi basis. The verifier checks

```text
T * (U + -q80_frame) * T^t = U + -rank17_gram,
det(T)=1.
```

It prints the complete matrix and the canonical row-serialization digest

```text
SHA-256 = 7116a499931bd096ba47fffc377a28690754bb451af0bdd7403f0e50438bd00d.
```

The complete matrix is pinned machine-readably in
[`data/fibrations/kumar_q80_rootless_target_to_q80_ns_transport.txt`](data/fibrations/kumar_q80_rootless_target_to_q80_ns_transport.txt),
and the verifier requires exact equality with it.

For equation work, the sequential local witnesses are preferable to the
large composite coordinates. The verifier also prints every transported
fiber in the initial q80 basis.

## Geometric scope

The lattice certificate proves six primitive `U` embeddings in the geometric
Neron--Severi lattice, hence the corresponding sequence of geometric elliptic
fibrations. The first two rational functions are explicit in characteristic
zero.  For the later steps, exact CM24 equation gates now exist over
`GF(73)`: both q12 steps have marked models, the productive compensated fifth
q4 model is pinned, and the saturated final q6 module is solved.  These
finite-field certificates do not supply characteristic-zero rational
functions, nor is every required divisor proved to be defined over the
eventual rational specialization field. Equation-level execution must still:

1. lift the third and fourth q12 models and the compensated q4/q6 suffix to
   characteristic zero;
2. compute or verify the associated Riemann--Roch pencils on the generic
   family;
3. track sections, fiber components, and Galois fields through each neighbor;
4. verify that the final seventeen sections specialize over `QQ` before any
   rank claim is made.

Thus this closes the missing **lattice route**, not the full elliptic-curve
construction.

## Why the recurrent q8 is not a shortcut

There is a structural source for the repeated number eight: the explicit
Humbert-8 entrance is the `D9+E7` family

```text
Y^2 = X^3 + T(r+(2r+1)T)X^2 + 2rs T^4(T+1)X + r s^2 T^7,
U = (X+sT^3)/T^4.
```

Its root determinant is eight, and the displayed two-neighbor is the exact
entrance to the Kumar `E7+E8` frame. This explains why `q=8` is a useful
normalization signal, but it does not make every q8 witness a generic MW3
generator. At CM43 the short horizontal class has a nonzero CM-only
component and its full fixed part collapses to the old fiber; the generic
level-79 direction is instead `Q79=4P1-5P2+P3`. The cases below therefore
record the frame and chamber outcome, not just the neighbor norm.

The oriented Humbert-8 cover also has the bounded rational chart

```text
z^2 = 2(16rs^2+32r^2s-40rs-s+16r^3+24r^2+12r+2),
r = (m^2-1)/(16(2n^2-1)),
s = m(16r-1)/(32r)-r+5/4+1/(32r),
z = (16r-1)n.
```

This is the preferred low-degree orientation chart for future Möbius and
marked-quotient comparisons. The exact checks presently classify the
nonloop second-child q8 as `E6+A7/MW4` (specializing to
`E7+A7+2A1/MW2`), not MW3; no specific q8 neighbor transformation producing
a new MW3 frame has been found.

The intrinsic rootless marking makes this distinction sharper. It contains
the primitive class

```text
h=(4,4,-1,0,-3,0,2,-2,1,-2,1,1,0,1,0,0,-2,-2,2),
h^2=-4, div(h)=4,
```

whose orthogonal complement has determinant 237 and Smith invariants
`1^17,237`. Exact Cauchy--Schwarz equality proves that there is no q8
isotropic fiber in `h^perp`: the `(1,8)` factor sum is too large, while
`(2,4)` would force the nonintegral vector `2x/3`. The first fibers in this
particular `H237` complement occur at q9. A bounded affine-CVP enumeration
finds exactly thirteen classes up to fiber sign. Their child root ranks are
8, 9, or 10 (20--24 roots), so none is rootless or the Kumar `E7+E8` frame;
the best nine have MW rank 9 and root data `(8,20,144)` or `(8,22,108)`.
The machine-readable classification is
[`rank17-h8-orthogonal-q9-fibers.json`](../artifacts/generated-results/rank17-h8-orthogonal-q9-fibers.json),
SHA256 `8998d250ae0b92f2e7dc891e61dc94614242508f152eab07970351cffa3503de`.
Thus q8 remains useful marking/normalization data, but it is not the direct
fiber that splits off this distinguished height-four summand.

The two generic determinant-948 `q=8`, `(a,b)=(2,4)` witnesses retained in
the earlier `A13+A1` continuation are

```text
(0,-2,-3,-2,2,2,-4,5,-4,5,2,0,-2,0,0,0,0),
(0,0,-1,3,1,0,0,1,0,0,0,0,0,0,-1,0,0).
```

Exact analysis of both child frames gives roots `A10+A2+2A1`, root
rank/count/determinant `14/120/132`, and MW determinant `79/11`. Their
reduced height Grams are

```text
(1/66)[79,-17, 1; -17,106, 19; 1, 19,259],
(1/66)[79,-17,-1; -17,106,-19;-1,-19,259].
```

Negating the first basis vector in the first case, or the second basis vector
in the second case, gives the original A10 target Gram exactly. These q8
neighbors therefore return to the already studied A10 frame type and do not
improve the reduction path. This is separate from the CM43 marked q8 class:
that specialization-only class is non-nef and its complete fixed-component
subtraction leaves the old fiber `F`. Neither occurrence shortens the live
q80 route. In contrast, the first q4 class below is primitive, fully nef, and
has the genuinely different child root system `D9+A4`.

The two generic q8 presentations are not non-nef artifacts. Exact old-chamber
reduction takes respectively 39 and 40 reflections, but both raw witnesses
land on the same primitive nef class

```text
(2,2,-32,23,22,15,-8,-4,-2,-23,6,-4,11,11,-18,17,-15,15,23).
```

It has `D.F=2`, `D.O=0`, MW coordinates `(25,-17,-23)` in the saturated
quotient used by the verifier, MW norm `79/14`, and root norm `33/14`. A
rank-three close-vector bound checks every possible section wall, while root
primitivity excludes a negative bisection. Thus it is a genuine generic q8
fibration, distinct from the old fiber, but both presentations produce the
same A10 return frame. No equation-level marked/global cover for this A13 q8
pencil has been constructed; it therefore supplies no evidence that its
marking is the desired source component.

The generic q8 frame calculations replay with

```bash
sage elkies-k3/scripts/analyze_mw3_branch.sage \
  --frame artifacts/generated-results/elkies-k3-lower-q-path-search/a13-mw3-child-search-frames/frame-007.txt \
  --name generic-q8-from-a13
sage elkies-k3/scripts/analyze_mw3_branch.sage \
  --frame artifacts/generated-results/elkies-k3-lower-q-path-search/a13-mw3-child-search-frames/frame-008.txt \
  --name generic-q8-from-a13-second
sage elkies-k3/scripts/analyze_a13_q8_neighbors.sage
```

## First geometric gate

The first `q=4` class has now been reduced in a deterministic old-fiber
chamber by
[`scripts/analyze_q80_rootless_first_neighbor.sage`](scripts/analyze_q80_rootless_first_neighbor.sage).
Starting from

```text
(2,2,-1,-1,0,0,0,0,-1,-1,0,0,0,0,0,0,0,0,0),
```

seventeen exact component reflections give

```text
(2,2,-6,-8,-11,-16,-13,-10,-13,-20,-26,-39,-8,-8,20,8,-2,-4,2).
```

It remains of old-fiber degree two, has zero intersection with `O`, and has
nonnegative intersection with all `A3+D5+E6` fiber components.  In the
chosen simple-root chamber its exact trivial-lattice expression is

```text
4F + 2O
 -2R4-3R5-2R6-2R7-R8
 -2R9-3R10-2R11-4R12-3R13-2R14.
```

Its projection to the generic rank-three MW lattice is zero and its root
norm is eight.  If a section has MW height `h` and root correction `c`, then
Cauchy--Schwarz gives

```text
D.S >= h + c - sqrt(8c) >= h - 2.
```

Consequently only the exact MW shell `h<2` needs inspection.  It consists of
the two directions `+/-P1`; their intersections with the reduced class are
respectively `1` and `7`.  This proves nonnegativity against **every
section**, using a rank-three MW enumeration rather than the infeasible
rank-17 frame shell.

This also suffices to rule out every multisection wall.  If an irreducible
`(-2)` curve `C` had `D.C<0`, it would be fixed in the effective class `D`,
so `1 <= C.F <= D.F=2`.  Degree one is a section.  In degree two, `D-C` is
vertical, so `C` has zero MW projection.  The analyzer explicitly checks
that the `A3+D5+E6` root lattice is primitive (equivalently the old fibration
has trivial MW torsion), hence

```text
C = kF+2O+r,  r in A3+D5+E6.
```

Writing the root part of `D` as `rho`, the identities `rho^2=8` and `C^2=-2`
give

```text
r^2=4k-6,
D.C=2k-<rho,r>,
||r-rho||^2=2(D.C+1).
```

If `D.C<0`, integrality and positive definiteness force `D.C=-1` and
`r=rho`; then `8=4k-6`, or `k=7/2`, contradicting integrality. Thus the
first q4 class is fully **nef**.

The zero MW projection and degree two put the first equation-level neighbor
in a compensated subpencil of `L(2O+4F)`.  The local calculation is now
complete.  Its fixed vertical coefficients split as

```text
D5: (2,3,2,2,1) = (1,2,1,2,1) + (1,1,1,0,0),
E6: (2,3,2,4,3,2) = (1,2,2,3,2,1) + (1,1,0,1,1,1).
```

At `T=0`, writing `x=T*xi` gives the first exceptional cubic
`(xi-1)^2(xi+2)`.  The bisection meets the two spinor ends, so it follows the
double root `xi=1`.  At infinity the IV* calculation passes through the
triple point and then meets the two outer E6 components.  These are the four
linear gates `c0=0`, `c1=-a`, `c3=0`, `c4=0` on
`a*x+c0+c1*T+...+c4*T^4`.  Therefore

```text
L(D) = <T^2, x-T>,             U=(x-T)/T^2.
```

[`scripts/derive_q80_first_q4_pencil.sage`](scripts/derive_q80_first_q4_pencil.sage)
checks the local transforms, derives the generic binary quartic, and verifies
the Jacobian discriminant.

### CM24 compact-pencil audit

The generic reduced class transports orthogonally to the extra CM24 `A1`:
the effective extra root has `P1.A1=1` and `D.A1=0`.  It was incorrect to
infer that the child root rank can rise by only one.  The new CM class can
combine with old MW directions, increasing the child root rank by three while
its MW rank drops by two.

[`scripts/search_q80_first_neighbor_rr.sage`](scripts/search_q80_first_neighbor_rr.sage)
tests the smallest compensated coordinates on the exact rational CM24 model.
The elementary chart already contains the component-selected hit

```text
U=(x-T)/T^2.
```

Its exact CM24 signature is

```text
I5* + I6 + 2 I2 + 3 I1   (D9+A5+2A1).
```

Over the unrestricted four-parameter ambient family its discriminant is

```text
(U-d+1)^4 R9,
```

giving `D9+A3`.  On each marked rank-19 formal branch, the degree-eight
collision factor vanishes and upgrades the finite `I4` to `I5`, recovering
the lattice target `D9+A4`, MW rank four.  The new verifier certifies that
collision through order four on both exact characteristic-zero branches.
Thus the first geometric neighbor is explicit; the next construction step is
to transport the marked sections and execute the second q4 pencil.

### Second q4 chamber

[`scripts/analyze_q80_second_neighbor_chamber.sage`](scripts/analyze_q80_second_neighbor_chamber.sage)
reconstructs the first child frame and reduces the second pinned witness
against its `A4+D9` components.  It again has old degree two, zero pole at the
zero section, and zero MW projection.  In deterministic simple-root order its
exact expression is

```text
D2 = 4F+2O
     -(1,2,2,1)_A4
     -(1,4,4,4,2,2,3,4,2)_D9.
```

The reflection sequence has fourteen steps and leaves nonnegative pairings
with `O` and every displayed old-fiber component.  The saturated MW lattice
has exact reduced height Gram

```text
(1/20)[19,5,0,-7; 5,35,0,-5; 0,0,80,40; -7,-5,40,171],
det=237/5.
```

The only nonzero MW vectors of height below two are `+/-e1,+/-e2`; their
closest integral section lifts have `D.S=2,6,3,5`, respectively.  The root
lattice is primitive and torsion is trivial.  The same norm identity as for
the first q4 divisor rules out a negative bisection: it would force the root
part to equal the divisor root part and then give the impossible fiber
coefficient `k=7/2`.  Thus the second class is fully nef.

The local equation is also explicit.  Write the first-child Jacobian as
`X^3+A1(U)X+B1(U)`, set

```text
U0=d-1,
v=U-U0,
x0=-3*B1(U0)/(2*A1(U0)),
x1=-A1'(U0)/(6*x0).
```

The double root of the infinity cubic is identically `alpha=3`, and

```text
W=(X-3*v^3-x1*v-x0)/v^2
```

gives an exact cubic identity after clearing `v^4`, with zero residual over
the unrestricted four-parameter ambient field.  The pinned rank-19 frame has
fibers `I3*+I1*+8I1`, i.e. `D7+D5/MW5`.

At CM24, `U0=-3/2`, `x0=-81/2`, and `x1=27/2`.  A complete bounded
degree-four search finds exactly two compensated coordinates.  The simple
infinity root gives the wrong `D12+A3+2A1` chamber.  The double-root choice
above gives

```text
I3* + IV* + 3 I2 + I1   (D7+E6+3A1), geometric MW=2.
```

The finite special fiber has valuations `(ord A,ord B,ord Delta)=(3,4,8)`,
so it is `IV*`, not `I2*`; this agrees exactly with the independently
transported pinned second frame.  Replay the equation and bounded audit with
[`scripts/derive_q80_second_q4_pencil.sage`](scripts/derive_q80_second_q4_pencil.sage)
and
[`scripts/search_q80_second_neighbor_rr.sage`](scripts/search_q80_second_neighbor_rr.sage).

### Third q12 gate

The simple compensated pattern ends at the next step.  In the second-child
`D7+D5` frame the raw pinned q12 witness

```text
(3,4,-2,1,2,1,0,0,1,0,0,0,0,0,0,0,0,0,0)
```

has `D.F=4,D.O=-1`.  Exact `O` and fiber reflections reduce it to

```text
(3,3,-18,-20,18,0,-24,-22,9,6,5,42,-71,-10,6,16,19,-8,-8),
```

with `D.F=3,D.O=0`, root coordinates
`(-6,-3,-3,-3,-2,-5,-1,-1,-2,-1,-2,-4)`, and nonzero MW projection of
norm eight. The projection is already an integral shortest frame lift,

```text
(-5,-4,6,3,-6,-5,3,3,1,9,-15,-1,0,4,5,-2,-2).
```

It defines a height-eight section `S` with `S.O=2`, and the reduced divisor
has the exact clean decomposition

```text
D3 = S + 2O + 2F + root_correction.
```

`S` meets the identity component of both the `D7` and `D5` fibers (all twelve
simple-component pairings vanish, while both affine-component pairings are
one). Thus the next small equation ansatz is the compensated chord through
`S`, generated horizontally by `(y+y(S))/(x-x(S))`, with the displayed root
correction selecting its vertical local branches.

Sequential specialization to CM24 now supplies an exact boundary marking.
The saturated CM MW basis has height Gram

```text
(1/12)[5,-2;-2,8].
```

I2 component profiles identify `e1=P2` and `e2=-(P1+P2)` up to simultaneous
sign, so the abstract q12 coordinates `(2,-1)` give
`Q_CM=P1+3P2`, not the initially guessed `2P1-P2`. Exact group law makes this
section polynomial:

```text
x(Q_CM) = -8/27 W^4 + 22 W^3 - 243/2 W^2 + 729 W - 492075/8,
y(Q_CM) = sqrt(-6)*(16/243 W^6 - 22/3 W^5 + 333/2 W^4
                 - 2025/4 W^3 + 190269/4 W^2
                 - 177147/16 W + 199290375/32).
```

It has CM height three, `P.O=0`, and meets only the conjugate quadratic I2
pair. The pole drop from the generic `S.O=2` is a CM boundary cancellation,
not a change to the generic divisor class. Enumerating the norm-four integral
lifts and imposing the explicit section incidence plus nonnegative pairing
with every old component leaves a unique effective lift. Its exact CM
decomposition is

```text
D = Q_CM + 2O + 4F + root_correction,
A1 coefficients = (0),(1),(0),
E6 coefficients = (1,2,3,2,1,2),
D7 coefficients = (2,4,3,3,5,6,3).
```

The extra `2F` compensates precisely for the pole drop. Only one selected
root of the quadratic I2 pair is a fixed component; the earlier ineffective
shortest lift with A1 coefficients `(1),(2),(0)` must not be used.
Therefore the generic-fiber RR
space is already the three-dimensional chord space

```text
< 1, X, (Y+y(Q_CM))/(X-x(Q_CM)) >,
```

leaving only the vertical `D7+D5` gates. This exact seed is replayed by
[`scripts/derive_q80_third_q12_pencil.sage`](scripts/derive_q80_third_q12_pencil.sage).
The sequential CM lattice transport, the height-three coordinates `(2,-1)`,
and the I2-profile identification of the explicit section are independently
replayed by
[`scripts/analyze_q80_third_q12_cm24_marking.sage`](scripts/analyze_q80_third_q12_cm24_marking.sage).
The raw chord by itself is not the pencil: eliminating `X` gives a squarefree
degree-eight double cover of the old `W`-line, hence genus three. The vertical
gates must select a nontrivial base-dependent combination of `1,X,z_Q` that
reduces this to the required genus-one pencil.

The target at this boundary is now pinned independently. Sequential transport
through the CM24 lattice gives root system `2A6+3A1`, root rank/count/determinant
`(15,90,392)`, and MW rank three.  The smallest restricted correction
`V=z_Q+k2*W^2+k1*W+k0` is exhaustive: lowering the residual degree forces
`k2=-2*sqrt(-6)/9` and `k1=15*sqrt(-6)/4` or
`21*sqrt(-6)/2`.  Their Jacobians have roots respectively
`D8+E7+A1` and `D5+E6+A3+A1`, so neither is the transported pencil.  The
missing vertical gate must use a nonzero coefficient of `X` in the full
three-dimensional space, not merely translate the chord. A second bounded
gate tests all 2,401 expressions
`z_Q+(a1*W+a0)*X+b2*W^2+b1*W` over `GF(7)`. Every nonzero X coefficient has
squarefree trigonal branch degree 17 or 23; none has degree at most eight.
Thus the coefficient must live in the saturated local module, with poles or
component compensation at the old `D7/E6` fibers, rather than in this naive
polynomial submodule. These finite exclusions are checked by
[`scripts/verify_q80_cm24_third_transport.sage`](scripts/verify_q80_cm24_third_transport.sage)
and
[`scripts/derive_q80_third_q12_pencil.sage`](scripts/derive_q80_third_q12_pencil.sage),
with the bounded modular scan in
[`scripts/search_q80_third_q12_xgate_gf7.sage`](scripts/search_q80_third_q12_xgate_gf7.sage).
Allowing one common denominator at the CM24 `IV*` fiber also fails: all 2,352
nonzero-X tuples of
`z_Q+((a1*W+a0)*X+b2*W^2+b1*W)/(W+3)` have branch-support degree greater than
eight. This bounded negative result is
[`scripts/search_q80_third_q12_local_module_gf7.sage`](scripts/search_q80_third_q12_local_module_gf7.sage).
The live construction must impose the separate exceptional-component
valuations at `IV*` and the infinity `I3*` fiber; a single unresolved-fiber
denominator is not the saturated module.

More precisely, the CM ambient coefficient space has dimension nine:
`a(W)+b*X+c(W)*z_Q`, with five coefficients in `a`, one in `b`, and three in
`c`. The unique effective correction predicts two local gates
`c(r)=a(r)+b*Qx(r)=0` at the selected quadratic-I2 root, one singular-point
gate at `E6`, and four gates at the infinity `D7` resolution.  The latter
must not be replaced by unresolved cusp jets: that false approximation has
rank seven but normalizes to genus four at `Vnew=1` and genus three at
`Vnew=7`.

The resolved local normal form is

```text
Y^2 + U^2 Z + Z^6 = 0.
```

The effective correction is its `Y`-valuation cycle, so the complete local
ideal is `(Y,U^2,ZU,Z^3)` and the quotient basis is `1,Z,Z^2,U`. Reducing the
nine ambient sections in this quotient gives the four correct D7 rows. With
the selected-I2 and E6 rows, the exact `7 x 9` matrix has rank seven and a
two-dimensional kernel, certified by
[`scripts/derive_q80_third_q12_local_gates.sage`](scripts/derive_q80_third_q12_local_gates.sage).

Writing its two kernel vectors as `N0,N1` gives the exact CM24 base

```text
Vnew = N1/N0.
```

Clearing `z_Q`, cross-multiplying, and eliminating `Y` gives a quartic in
`X`; the factor `X-Qx` is the expected artifact of clearing the chord
denominator. Dividing it leaves an exact cubic of degree nine in the old
base `W`. With weight three on `X`, its infinity form is exactly a nonzero
scalar times `(xi-3)^2*(xi+6)`, independent of `Vnew`; the simple branch
`xi=-6` supplies a distinguished rational infinity place for the Jacobian
conversion. The resolved plane model is replayed by
[`scripts/derive_q80_third_q12_cm24_pencil.sage`](scripts/derive_q80_third_q12_cm24_pencil.sage).
At the split prime `73` (`sqrt(-6)=33`, `sqrt(-3)=17`) the residual cubic is
irreducible and its normalization has genus one at both `Vnew=1` and
`Vnew=7`. This decisive acceptance gate is replayed by
[`scripts/verify_q80_third_q12_cm24_genus.sage`](scripts/verify_q80_third_q12_cm24_genus.sage).
At `Vnew=7`, plane-curve Brill--Noether at the canonical simple infinity
branch `xi=-6` (Singular place 6) gives semigroup `0,2,3`; the resulting
pole-2 and pole-3 functions yield
`y^2+12xy+27y=x^3+51x^2+40x+26` over `GF(73)`, certified by
[`scripts/analyze_q80_third_q12_cm24_weierstrass_gf73.sage`](scripts/analyze_q80_third_q12_cm24_weierstrass_gf73.sage).
The local analyzer also applies a fail-closed rename of the upstream
`hnoether.lib` local ideal `a`, which otherwise collides with Singular's
algebraic-extension parameter on nonsplit fibers.

Keeping the same `xi=-6` origin at 56 good values of `Vnew` reconstructs the
generic CM24 Jacobian over `GF(73)(Vnew)`. The first 49 values give a
`49 x 50` interpolation matrix of rank 49 for bidegree `(24,24)`; all seven
withheld values pass. The resulting short equation is

```text
y^2 = x^3 + A(V)x + B(V),
A = 6V^8+16V^7+47V^6+33V^5+58V^4+2V^3+63V^2+17V+23,
B = 33V^12+64V^10+61V^9+45V^8+14V^7+20V^6
    +54V^5+8V^4+50V^3+57V^2+47V+43                 (mod 73).
```

Its discriminant factors as

```text
10(V+37)(V+46)(V+17)^2(V+30)^2(V+68)^2
  (V+20)^7(V+67)^7(V^2+20V+67),
```

so its fibers are exactly `2I7+3I2+4I1`, giving the transported
`2A6+3A1/MW3` CM24 marking by Shioda--Tate. The numerator of `j` is `c4^3`
with `deg(c4)=8`, and `c4^3-1728*Delta=c6^2` with `deg(c6)=12`. This stronger
certificate is
[`scripts/reconstruct_q80_third_q12_jacobian_gf73.sage`](scripts/reconstruct_q80_third_q12_jacobian_gf73.sage).
The conjugate embedding `sqrt(-3)=56` gives a second exact degree-24
`j`-map. Exhausting the fractional-linear maps forced by its marked
`I7/I2/I1` pole divisor finds no `phi` in `PGL2(F_73)` with
`j_56(V)=j_17(phi(V))`; consequently there is no base-Mobius descent at this
prime. The equal-`j` cross-polynomial has factor degrees `1+1+6+8+32`. One
linear root, `V=-37`, lies on the original discriminant; the other, `V=-27`,
is smooth for both maps and is the unique smooth `F_73`-rational fixed-`j`
residue. The complete two embeddings, normalized `c4,c6,A,B`, and
cross-polynomial are recorded in
[`q80-third-q12-jacobian-gf73.json`](../artifacts/generated-results/q80-third-q12-jacobian-gf73.json)
and replayed by
[`scripts/analyze_q80_third_q12_galois_descent_gf73.sage`](scripts/analyze_q80_third_q12_galois_descent_gf73.sage).
This is a modular descent obstruction and specialization hint, not yet an
exact coefficient pair over `Q(sqrt(-3))`.

The transported fourth horizontal section supplies a sharper test at the
unique smooth fixed-`j` residue.  Independently enumerating the thirty
polynomial sections on both split embeddings (after undoing the conjugate
short model's constant nonsquare twist) gives 16 compatible markings and four
unique fourth-horizontal candidates on each side.  At `V=-27 mod 73` the two
smooth fibres differ by the nonsquare twist `34`, but none of the 16 candidate
pairings satisfies the required twisted x-coordinate identity.  Hence the
fixed-`j` residue does not align the fourth horizontal section under any
marking ambiguity.  This exact finite-field negative diagnostic is replayed
by
[`analyze_q80_fourth_section_fixed_j_gf73.sage`](scripts/analyze_q80_fourth_section_fixed_j_gf73.sage);
it closes only this specialization shortcut, not the global vertical-gate
construction.

Thus the CM24 third q12 pencil, its Jacobian, and its fiber marking are exact
in characteristic 73. The next gate is lifting this compact conjugate pair to
the CM coefficient field, testing whether the anti-invariant characteristic-
zero coefficient has a linear factor reducing to `V+27`, and transporting
the fourth neighbor.

A bounded exact comparison over `q=4,...,12` sampled 251 sign-pairs per norm.
Its best sampled q8 witness

```text
(0,0,-1,0,0,0,-2,-1,-2,0,-2,0,0,-1,0,-1,0),  (a,b)=(2,4),
```

has child roots `E6+A7` and MW rank four. A sampled q10 gives
`D5+A5+A3+A1/MW3`, whereas the pinned q12 gives
`A5+A3+3A1/MW6`. Thus q8 is a genuine alternate high-root presentation at
this node, not a rank-growing shortcut in the bounded sample. Exact chamber
reduction takes 28 zero/root reflections and leaves the primitive degree-two
class

```text
(2,2,-4,-7,3,-4,-6,-6,1,0,2,14,-24,-5,4,4,4,-3,-2),
```

with `D.O=0`; it does not collapse to the old fiber. Its root lattice is
primitive (torsion one), and one reduced saturated MW height Gram is

```text
[[19/24,-7/24,0,1/4],[-7/24,43/24,0,-1/4],
 [0,0,4,2],[1/4,-1/4,2,17/2]],
```

of determinant `79/2`.  Hence it is neither the old A13-to-A10 q8 loop nor
the collapsing CM43 marking. The norm-four/six section shell is exhaustive
for a negative section and checks 124,008 signed lifts, with minimum pairing
zero; the remaining bisection case would force the impossible integral
coefficient `3/2`. Thus this reduced class is fully nef. At CM24 it remains
noncollapsed and specializes
to `E7+A7+2A1/MW2`.  It may provide a useful degree-two coordinate
normalization, but it loses a generic MW direction (`5 -> 4`) while the
pinned q12 gains one (`5 -> 6`), so it cannot replace the q12 rank-growing
step on present evidence. This scan is not exhaustive and does not prove q12
optimal.

### Fourth q12 gate

The next pinned q12 witness starts in the generic `A5+A3+3A1/MW6` third
child. Exact zero/fiber reduction takes 16 reflections and gives

```text
D=(4,3,-14,-57,-2,19,-8,10,12,-2,-5,-11,-1,-6,17,5,-12,3,12),
D.F=3, D.O=1.
```

Its MW projection has height 15 and coordinates `(-3,-1,0,0,0,-1)` in the
pinned reduced saturated rank-six basis. The effective horizontal section
has `P.O=6`, and the exact generic decomposition is

```text
D = S + 2O - F + root_correction.
```

This class is genuinely nef. If a fixed irreducible `(-2)`-curve is
`C=(a,m,l)`, then effectiveness of `D-C` gives `1<=m<=3`, while

```text
D.C = (3/(2m))*||l-(m/3)d||^2 - 3/m.
```

Thus a negative wall must lie at frame distance less than two from one of
three rational centers. Exact enumeration in the rank-six MW quotient,
followed by rank-eleven root-CVP lifting, exhausts every such integral class;
no iterator hits its fail-closed cap and no negative wall exists. This is
replayed by
[`scripts/analyze_q80_fourth_q12_chamber.sage`](scripts/analyze_q80_fourth_q12_chamber.sage).

At CM24 the same section drops from `P.O=6` to `P.O=2`. The specialized
source is the exact `2A6+3A1/MW3` model above. A node-constrained search needs
only 219 quartic square tests for the two four-node profiles and 389,017
tests for the remaining two-node profile; it recovers all 30 polynomial
sections in under ten seconds. One saturated basis is

```text
H = (1/14)[[2,1,0],[1,6,-2],[0,-2,16]],
Q = -P1-3P2-P3,
Q.O=2, h(Q)=33/7.
```

For this basis

```text
h = V^2+60V+53,
x(Q) = (38+69V+69V^2+7V^3+45V^4+4V^5+4V^6+3V^7+70V^8)/h^2,
y(Q) = (50+62V+5V^2+24V^3+56V^4+68V^5+27V^6+62V^7
        +46V^8+38V^9+56V^10+40V^11+34V^12)/h^3            (mod 73).
```

The CM divisor decomposition is `D=Q+2O+3F+root_correction`, so its
generic-fiber Riemann--Roch space is exactly

```text
<1, x, (y+y(Q))/(x-x(Q))>.
```

The specialized fourth child has roots `A5+A3+4A2` and MW rank two; the
generic child remains the pinned `4A1/MW13`. The lattice marking is replayed
by
[`scripts/analyze_q80_fourth_q12_cm24_marking.sage`](scripts/analyze_q80_fourth_q12_cm24_marking.sage),
the section equations and group-law marking by
[`scripts/search_q80_third_child_polynomial_sections_gf73.sage`](scripts/search_q80_third_child_polynomial_sections_gf73.sage),
and the exact coefficients are in
[`q80-fourth-q12-cm24-gf73.json`](../artifacts/generated-results/q80-fourth-q12-cm24-gf73.json).

The equation-level vertical gate is now exact. Smith saturation of
`{1,x,h z_Q}` has invariant-factor degrees `(0,0,6)`; shifted Popov reduction
gives weights `(0,2,0)` and a ten-dimensional global ambient space. For the
effective marking `candidate=1`, selected `I2=5`, the resolved semistable
component matrix has rank eight and a two-dimensional kernel. Removing its
fixed part gives a bidegree `(T,V,X)=(2,14,3)` moving equation. At `T=1` it is
irreducible and its function field has genus one.

The five finite reducible values are `T=14,25,47,58,67`; their geometric
factor-degree patterns are respectively `(1,2),(1,2),(1,1,1),(1,1,1),
(1,2)`, while the infinity cubic is irreducible. Together with the exact
specialized root lattice `A5+A3+4A2`, this gives the semistable multiset
`I6+I4+4I3`. After removing the persistent square factors of the cubic-in-X
discriminant, the remaining branch sextic has local orders `3,4,6,3,3` at
`T=14,25,47,58,67`. Thus these finite fibers are respectively
`I3,I4,I6,I3,I3`; the remaining lattice `A2` is the irreducible infinity
`I3`. Restricting the old polynomial sections gives nine section-or-vertical
hits, including five nonconstant degree-one sections of the new fibration.
These provide rational points for an explicit algebraic lowering of the
trigonal equation.

The compact artifact is
[`q80-fourth-q12-cm24-moving-cubic-gf73.json`](../artifacts/generated-results/q80-fourth-q12-cm24-moving-cubic-gf73.json),
SHA256
`c6560b3db2d1232866e9996fc727924090aa46293c2482885cf9f9dbf4c21c89`.
The stripped branch-discriminant certificate, including the degree-six
branch polynomial as exact coefficient numerator/denominator arrays, its
complete discriminant factorization, and local orders, is
[`q80-fourth-q12-cm24-discriminant-gf73.json`](../artifacts/generated-results/q80-fourth-q12-cm24-discriminant-gf73.json),
SHA256
`6caa3c9bb83a115a1e40689bf23d58dea1dd7ae1c77795f9970b6d24517a7ef0`.
The wrong selected `I2=43` marking does produce a quadratic model, but its
Jacobian has the wrong reducible-fiber signature and is rejected. Full
normalization and backend-dependent `L(3P)` conversions were stopped; no
fourth Jacobian has yet been claimed. The live next gate is a direct
elimination using one of the five degree-one sections, with local-delta
classification as the bounded fallback.

### Fifth q4 readiness

The next pinned q4 class has now been reduced in the generic `4A1/MW13`
fourth-child chamber. One root reflection gives

```text
D=(2,2,-6,3,-4,0,6,2,-2,-5,-3,-5,-6,4,1,1,-1,-2,-1),
D.F=2, D.O=0, MW norm=6.
```

Its child has `A1/MW16`. This is an exact old-component chamber result, not
yet a full nef proof against every high-MW section.

At CM24, transporting the same class into the specialized
`A5+A3+4A2/MW2` fourth child and applying fourteen component reflections
again gives `D.F=2`, `D.O=0`. In a saturated MW basis with

```text
H = [[1/3,1/6],[1/6,5/12]],
```

its horizontal coordinates are `(-1,1)`, height `5/12`, and `P.O=0`.
Applying the fourth-pencil Weyl word to the complete neighbor basis resolves
the earlier raw/geometric marking ambiguity. For each of the five old
polynomial x-classes, exactly one y-sign becomes a degree-one section; their
new MW coordinates are

```text
(1,1), (0,-1), (0,1), (1,-1), (1,1).
```

They span the full rank-two MW lattice, and the fifth horizontal target
`(-1,1)` is the inverse of the section with x-coefficients
`(12,26,13,39,49)`. Thus the fifth q4 horizontal direction is already
present in the compact section package; a fourth global Jacobian conversion
is not required for this readiness gate.

The smallest equation ansatz is now pinned directly on the compact cubic.
For the correctly oriented section

```text
X = 12+26V+13V^2+39V^3+49V^4,
Y = 52+20V+45V^2+26V^3+44V^4+37V^5+9V^6,
T = (-4V+30)/(V-34),
```

the old Weierstrass identity and moving-cubic identity both vanish exactly.
Its local factor is recorded at every finite reducible support and at the
infinity `I3`; at `T=25` it passes through the intersection of the linear and
quadratic factors. If `A` is the normalization of the compact coordinate
ring and `I_R` is this section ideal, the fifth generic-fiber space is the
two-dimensional fractional module

```text
L(O+(-R)) = Hom_A(I_O I_{-R}, A).
```

The geometric fourth zero is an old-degree-18 curve, so replacing `I_O` by
one of the low polynomial sections would be a marking error. The compensated
nonconstant generator and vertical correction remain open. The ansatz is
pinned in
[`q80-fifth-q4-local-module-ansatz-gf73.json`](../artifacts/generated-results/q80-fifth-q4-local-module-ansatz-gf73.json),
SHA256
`7cadf12e4035dc9325f3249158f906762aba85cebdb2e92cb4de93efc2140d15`.

There is now also a separate, equation-friendly fifth q4 presentation. A
bounded norm-eight shell window found sixteen `A1/MW16` candidates whose
specialized horizontal class `(1,0)` is represented by three pairs among the
five degree-one sections. For each fiber the compensated cubic is constructed
first, then the projection is normalized by two further marked sections to
take the values zero and one. Only after this fiberwise marking is fixed are
the coefficients interpolated. The pair `(1,4)` is the unique tested mode
for which all twelve coefficients reconstruct from fourteen samples and pass
four withheld fibers. Its branch squareclass has degree four over the new
base, hence gives an exact genus-one fifth-child gate over `GF(73)`. The
Jacobian discriminant factorization is

```text
(u+17)^2 u^5 (u-1)^6
  (u^5+53u^4+41u^3+70u^2+64u+22),
```

up to the displayed sixth-power denominator. Restoring the non-square unit
discarded by monic factorization gives the required quadratic twist

```text
65u^4+70u^3+61u^2+20u+1
 = 65(u+51)(u^3+68u^2+u+56).
```

The twist clears every denominator and gives polynomial degrees `(8,12,18)`
for `A,B,Delta`. Since `gcd(A,Delta)=1`, the exact semistable signature is
`I6+I6+I5+I2+5I1`: the finite orders are `6,5,2,1^5`, and degree 18 leaves
order six at infinity. Thus the CM24 root rank is 15 and MW rank is 3, with
root data `(15,82,360)` and ADE type `2A5+A4+A1`. A characteristic-zero
lift remains open. The marked-projection artifact is
[`q80-fifth-q4-marked-projection-pair14-gf73.json`](../artifacts/generated-results/q80-fifth-q4-marked-projection-pair14-gf73.json),
SHA256 `e46c9925c6870a6f9185f36994a5aef682382bba7a9bf8d2adc3d897420988fa`.
The separate twisted-Jacobian artifact is
[`q80-fifth-q4-marked-pair14-jacobian-gf73.json`](../artifacts/generated-results/q80-fifth-q4-marked-pair14-jacobian-gf73.json),
SHA256 `78c654d35acccb907a3b019bc309c84d7d7b705d8d6a521e17f3f169fad67ca9`.
The independent local-minimalization audit also shows why the unit is
essential: the monic odd part alone has four additional `I0*` fibers over the
unit divisor, Euler number 48, and root rank 31, so it is not a K3 model.

The pair-14 equation is now matched to an exact **CM24-only lattice class**.
An affine-CVP parameterization of the oriented norm-eight shell for
`(a,b)=(2,2)` and horizontal coordinates `(1,0)` contains 797,472 signed
classes. The bounded default run stops after its first full hit: it processed
193 primitive presentations, 30 of which had the target root count, before
reaching root data `(15,82,360)`. This certifies the displayed class, but not
uniqueness among all 797,472 presentations. Its raw class is already reduced:

```text
D=(2,2,-98,42,47,21,8,66,-179,29,81,219,30,65,48,23,32,53,-89,-8),
D.F=2, D.O=0, old-fiber degree=49, old-zero pairing=17.
```

The child is exactly `2A5+A4+A1/MW3`, and its explicit rank-18 frame has
determinant 24. The search certificate is
[`q80-fifth-q4-cm24-pair14-marking.json`](../artifacts/local/q80-fifth-q4-cm24-pair14-marking.json),
SHA256 `c3949e37638b138184bf9591e127aa19c93d8d927a4b224b865df7c20dcf6cac`;
the child-frame SHA256 is
`9eb0a7508f08b01bf6a226aee9c0fe51834b3a167dc2d1c449c705ee250357d1`.
This closes the equation/lattice marking at CM24, but it does **not** make the
class a divisor on the determinant-948 generic rank-19 lattice.

The generic deformation comparison is also exact within the three retained
productive shell windows. Their 48 distinct `A1/MW16` candidates collapse at
CM24 to only two nef classes: 32 have old-fiber degree 47, and 16 have degree
43. Every one has child root data `(16,66,2048)`, namely
`D4+3A3+3A1/MW2`; none has pair-14 data `(15,82,360)`. The three bounded
records and hashes are:

```text
10000:20000    0cb51a1593af0a12acbb743ac7f03177849d0e828c8295ef8ffb647101b1f943
20000:72500    ed52099c076218d2c826e69a1a512d0590c321883982f543ba502d85ca1f840f
72500:125000   a9bba4a5f62a499f82022833a0d2d86cf03a1e624318b25bbb7de836cd8b6f60
```

These are bounded indexed-window experiments, not a proof about every
possible neighbor norm or marking. A separate full `(a,b)=(1,4)` shell scan
found no generic root-rank-at-most-one child. The safe conclusion is that
pair14 is a certified CM24 boundary fibration, while the known generic
rootless continuation remains attached to the degree-47
`D4+3A3+3A1/MW2` specialization below.

The apparent pair ambiguity is now resolved in the saturated CM24 root
lattice.  The special fourth frame has order-three torsion, so the free MW
coordinate `(1,0)` has three distinct integral section lifts.  Exact
discriminant-component profiles identify the degree-43 class with explicit
section pair `(0,1)`, and the productive degree-47 class with pair `(2,3)`.
Their full vertical decompositions are:

```text
degree 43: O + P_(0,1) + the identity components of two I3 fibers;
degree 47: O + P_(2,3) + 2*Theta0 + Theta1 + Theta2 on one I6 fiber.
```

There are no residual whole fibers.  Consequently the earlier raw pair-23
projection of genus three was not the desired divisor: it omitted the
single-I6 compensation.  The live fifth-neighbor problem is now a local
compensated pair-23 projection, not a degree-47 global Riemann--Roch solve.
The exact lattice certificate is
[`q80-deforming-fifth-q4-vertical-compensation.json`](../artifacts/local/q80-deforming-fifth-q4-vertical-compensation.json),
reproduced by
[`analyze_q80_deforming_fifth_vertical_compensation.sage`](scripts/analyze_q80_deforming_fifth_vertical_compensation.sage).

The compensated pair-23 equation is now exact over `GF(73)`. Two rational
secant values

```text
mu_A=(68*T^2+56*T+50)/(T^2+29*T+37),
mu_B=(34*T^2+55*T+13)/(T^2+29*T+37)
```

give the coherent gauge `s=(U-mu_A)/(mu_B-mu_A)`. Restoring the binary-quartic
factorization unit gives fibers `I0*+3I4+3I2`, root data `(16,66,2048)`, and
CM24 MW rank two, exactly matching the productive degree-47 lattice class.
The standalone verifier is
[`verify_q80_deforming_fifth_pair23_gf73.sage`](scripts/verify_q80_deforming_fifth_pair23_gf73.sage),
and its generated artifact has SHA256
`23fc49bce2618a6d3c5f5e18ded34b4ffbee220be83523ae250bf7774a91db14`.
This is still a finite-field equation certificate; lifting the gauge along the
characteristic-zero branch remains open.

The other marked-section pairs do not repair this mismatch. Pair `(0,1)`
reconstructs a distinct CM24 K3 with root data `(15,78,384)`, namely
`D4+A5+2A3/MW3`; its exact equation artifacts are
[`q80-fifth-q4-marked-projection-pair01-gf73.json`](../artifacts/generated-results/q80-fifth-q4-marked-projection-pair01-gf73.json),
SHA256 `f992b085aa5740f37649a69b83ff9ee6b2eafba620b60932a940b9a3bb3cb0a4`,
and
[`q80-fifth-q4-marked-pair01-jacobian-gf73.json`](../artifacts/generated-results/q80-fifth-q4-marked-pair01-jacobian-gf73.json),
SHA256 `dbabe70bb39c757c3c20a03c0e6178962956549305745dbef5b2f135096807dc`.
It too has now been joined to the CM24 lattice. The first target hit is

```text
D=(2,2,-99,46,46,21,13,65,-182,29,86,216,35,67,49,27,32,52,-95,-8),
D.F=2, D.O=0, old-fiber degree=47, old-zero pairing=16.
```

The bounded marking artifact is
[`q80-fifth-q4-cm24-pair01-marking.json`](../artifacts/local/q80-fifth-q4-cm24-pair01-marking.json),
SHA256 `56a30031c20754d4a490e2d3985dc369199e095b4e024c035cd758cfd2e20946`;
the extracted determinant-24 frame has SHA256
`e0dbf2249464c6cae253b7c53503d11a80b98f6ad3515e3c48b43c1b8c5d5360`.
A 300-vector horizontal/outside-root q6 sample reaches root rank 14/MW4 but
no lower. The **uncompensated** pair `(2,3)` projection has squarefree cover
degree eight and genus three and is rejected; the compensated pair-23 class
above is the generic degree-47 divisor that carries the certified q6-to-MW17
suffix.

Neither CM equation branch currently beats the generic suffix. For pair14,
300 horizontal/outside-root samples at each of q4, q6, q8, q9, and q10 all
remain at root data `(15,82,360)/MW3`; q12 was stopped at the 60-second cap.
Exact p-neighbor beams tested 2,941 frames for `p=5,7,11` through four
generations, 1,200 direct larger-prime neighbors, and 1,600 neighbors of the
best rank-14 frame for `p=5..29`. Their best root rank is 14/MW4. These are
bounded negative searches, not proofs that the determinant-24 surface has no
rootless fibration.

Two nearby constructions are explicitly rejected. A raw line through the
section leaves degree 13 on the singular bidegree-`(14,3)` model.
Interpolating an independently depressed maximal-order coordinate without
first fixing the marked PGL2 gauge gives a squarefree branch degree 12;
postcomposition cannot remove that squareclass, so it is not a coherent
global pencil. The generic and specialized Singular adjoint probes both
reached their strict 30-second caps and were stopped.

The selected alternate lattice child is not close to the *previously pinned*
rootless frame:
exact transport puts that particular endpoint at direct neighbor norm
`1,382,576`, and the closest of the thirteen intrinsic `h^perp` q9 fibers has
norm `19,631,636,396` in the alternate marking. That does not obstruct a
different rootless fibration. A bounded q2--q6 beam found the exact q6 witness

```text
(a,b)=(2,3),
v=(0,-2,4,2,-1,2,1,-1,1,0,1,-1,1,0,0,0,0).
```

Its child has no roots and therefore MW rank 17. The full original-q80 to
alternate-q4 to q6 NS composite has determinant `+1`. This rootless frame is
not integrally isometric to the old pinned `rank17_gram`, which is harmless:
it is a second rootless fibration on the same determinant-948 K3. The exact
frame and composite transport are in
[`q80-alternate-fifth-q6-rootless-transport.json`](../artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json),
SHA256 `48381d91e288b2cefb85b1484d351d659748f801ea57d190453bd2db0a56eaab`.
The frame and transport matrix hashes are respectively
`fa2d994c3cb78283a872918a54dceaa6807ceb9a6ee1d2b5cb51f213b692b721`
and `c859d5e71d1d6158936e4c1460b38f13edc25484cdbb091e90a05831353d5e57`.

The final q6 pencil is also only degree two after chamber reduction. The raw
class has `D.F=3,D.O=-1`; reflecting once in `O` and once in the sole `A1`
component gives

```text
(3,2,-1,-2,4,2,-1,2,1,-1,1,0,1,-1,1,0,0,0,0),
D.F=2, D.O=1, MW norm=23/2.
```

Its pairings with `O`, the simple `A1`, and the affine component are all one.
At CM24 its selected fifth child has roots `D4+3A3+3A1/MW2`; the
horizontal class has MW coordinates `(-1,-3)`, height `7/2`, and `P.O=1`.
Full nefness is now exact. Write the reduced class as `D=(3,2,v)` in
`U+(-M)`, so `v.M.v=12`. A section indexed by `w in M` has

```text
D.S = (w-v/2).M.(w-v/2)-2.
```

Exact closest-vector enumeration gives minimum distance three and therefore
minimum section pairing one. A degree-two `(-2)` class would instead satisfy

```text
D.C = (w-v).M.(w-v)/2-1,   w.M.w = 2 (mod 4).
```

Negativity would force `w=v` by positive definiteness and even integrality,
but `v.M.v=12` is zero modulo four, a contradiction. Since a fixed curve has
old-fiber degree at most `D.F=2`, the component, section, and bisection checks
exhaust all possible walls. Thus `D` is nef and defines a genuine degree-two
pencil. The remaining equation gate is only to reconstruct that function
after the generic fifth vertical class; no further chamber correction is
needed.

The same calculation identifies the marked horizontal section without a
search. The integral class

```text
S=(5,1,-1,-2,4,2,-1,2,1,-1,1,0,1,-1,1,0,0,0,0)
```

has `S^2=-2`, `S.F=1`, `S.O=4`, and meets the nonidentity A1 component. It
satisfies the exact identity

```text
D = O + S - F.
```

Consequently the final generic-fiber space is the two-dimensional relative
space generated by `1` and `(y+y(S))/(x-x(S))`, but the CM24 specialization
has a nontrivial saturated vertical module. The transported effective section
is identity at all three `I2` fibers, meets endpoints of two `I4` fibers, and
meets a nonidentity `I0*` component. Its corrections are therefore
`3/4+3/4+1=5/2`, as required by height `7/2` and `P.O=1`. The earlier filter
selecting two `I2` components and the identity `I0*` component was wrong; the
corrected finite-field section artifact has SHA256
`fcd61f89daab0a68785a006e6b10dc3829b1c30c24243b67a4e1b80c7d6e6e09`.

In the deterministic special-fifth simple roots, the exact vertical class is

```text
D-O-S = 2F-(R1+R5+R8+R11+R13+R16).
```

The three nontrivial local cycles are:

```text
A1: F-R1                    = affine;
A3 on (R5,R6,R8): F-R5-R8  = affine+R6;
D4 on (R11,R12,R13,R16), marks (1,1,2,1):
    F-R11-R13-R16           = affine+R12+R13.
```

All three local transforms require an affine copy, while the global class has
only `2F`. Thus `D-O-S` is not effective and the constant section does not
survive the specialized local module. This explains why raw chord gauges and
single principal-part corrections fail. The same exact checker shows that the
specialized q6 child has roots `A1+2A3+2A4`, root data `(15,66,800)`, and MW rank three; only the
generic child is rootless/MW17. These data are replayed by
[`analyze_q80_final_q6_cm24_marking.sage`](scripts/analyze_q80_final_q6_cm24_marking.sage).

This also determines the smallest correct linear system. Put
`L=O+S+2F`. Exact component checks and a section CVP give minimum section
pairing one; a completed-square argument gives bisection pairing at least two.
Thus `L` is nef, `L^2=6`, and K3 Riemann--Roch gives `h0(L)=5`. The desired
space is the kernel of the three disjoint elementary transforms above:

```text
H0(D) = ker(H0(L) -> k_A1 + k_A3 + k_D4),
dimensions: 2 = 5 - 3.
```

This `3 x 5` saturated local-module calculation is now solved exactly over
`GF(73)`.  If `H=s-27`, clearing the chord denominator gives a two-column
polynomial numerator module whose Smith invariants are

```text
diag(1,H).
```

The primitive minors have gcd `H`, while the saturated minors have gcd one.
After reducing the Smith generator by a polynomial multiple of the constant
section, the compensated generator is

```text
q_sat = (q0+63)/(s-27).
```

Thus the actual ambient basis is

```text
1, s, s^2, q_sat, s*q_sat,
```

not the unsaturated raw-chord basis.  In this basis the three local gates are

```text
I2:   (s,q0)=(72,6),
I4:   (s,q0)=(64,3),
I0*:  a2=0.
```

Their rank-three kernel has the exact rows

```text
(1,0,0,41,48),
(0,1,0, 6,72).
```

The unit-preserving binary-quartic Jacobian has fibers
`2I5+2I4+I2+4I1`, root data `(15,66,800)`, and CM24 MW rank three, exactly
matching the transported lattice child.  The replay command is

```text
sage elkies-k3/scripts/search_q80_final_q6_local_module_gf73.sage \
  --section-index=0 --write-artifact
```

and the generated certificate
[`q80-final-q6-saturated-module-gf73.json`](../artifacts/generated-results/q80-final-q6-saturated-module-gf73.json)
has SHA256
`7d3866855b9995b733193de9c5d5e3ba1cea6aa1292141e06d0d5011f28975e3`.
The artifact also records the explicit `q0(s,R)` and globally minimal
polynomial `A`, `B`, and discriminant coefficients of the selected child.
This closes the final finite-field equation gate.  The remaining downstream
task is characteristic-zero lifting of the pair-23 fifth model and this
compensated q6 pencil.

This downstream q80 chain must not be confused with the newly recovered
source polarization.  The q80/Kumar route uses

```text
2*H2 = diag(8,237),
```

whereas the exact `H21 intersect H92` source branch uses

```text
2*H3 = [[21,6],[6,92]].
```

Both forms have determinant `1896` and Smith invariants `(1,1896)`, but they
are not integrally isometric: `2*H2` has minimum eight, `2*H3` has minimum
21, and PARI `qfisom` returns zero.  Therefore the q80-to-rootless
construction is a valid downstream fibration chain in the determinant-948
genus, not the original `H3` source polarization.  The recurrent `q=8`
signal belongs to the height-four `H2`/Humbert-8 entrance and to special
Noether--Lefschetz return or collapse phenomena.  The Smith factor `s-27`
above is a local saturation of the final q6 chord module; it supplies no
norm-eight polarization and does not identify the `H3` source.

The exact pullback of the generic horizontal section to the compact
fourth frame is

```text
(29,28,-129,115,-159,144,-58,0,72,13,-86,87,13,-28,-13,-15,14,29,0),
```

of old-fiber degree 28 and old-zero pairing one.  Its saturated three-fiber
local transform is the exact `GF(73)` construction above; only
characteristic-zero lifting remains.

### Marked coefficient field

The normalized slope-`5` surface-coefficient branch is rational over
`GF(7)`, but its first marked section is not. Component conditions force

```text
X1 = T + (d-1)T^2.
```

After substituting the certified rational surface parameter, the leading
marked `y` coefficient has square class

```text
v^2 = (t+3)(t+4)(t+6)(t^2+6t+4).
```

The five finite branch points and infinity give six branch points, hence an
exact genus-two quadratic cover. The second marked polynomial section is
defined over a distinct genus-three quadratic cover. Their branch divisors
overlap in degree five, so adjoining both markings gives a genus-six
biquadratic cover; the third quadratic quotient has genus one. Thus the
rational surface parameter does not rationalize either the marked system or
the first-neighbor pencil.

This modular cover is not yet the known characteristic-zero genus-two source
model. Modulo seven, its absolute Igusa invariants are `(4,2,5)`, whereas the
reduction of

```text
u0^2 = 16t0^6-19t0^4+88t0^2-48
```

has invariants `(4,4,3)`. The mismatch rules out a Mobius identification even
over the algebraic closure in characteristic seven. Together with the
genus-six combined marking, this makes descent through the two chosen short
sections a high-genus detour, not the desired source marking. It does **not**
rule out slope 5 as an unmarked q80 surface component: the actual level-79
marking may involve the third pole section `P3`. The fourteen
surface-ideal generators do continue to vanish through order 90, five orders
beyond their fitting jet, so this rejection is about the marked construction
route rather than a simple order-85 interpolation failure.

Two independent symmetry gates agree. The P1 branch sextic has trivial
projective stabilizer over `GF(49)`, while the known source reduction is
bielliptic via `t0 -> -t0`. The candidate's genus-one third quotient has
`j=0`, whereas the two known source elliptic quotients reduce to `j=1,3`.

The third pole section has now been tested directly rather than inferred from
P1/P2.  After accelerating the formal engine, the slope-5 marked branch lifts
through order 230 in seconds.  Over the rational surface parameter `t`, the
P3 pole `lambda` satisfies a unique quadratic equation of `t`-degree 16 with
twelve withheld coefficients.  Its discriminant has factorization

```text
(t+4)^2 (t+6)^4 (t+3)^11 (t^3+3t+5),
```

and hence squarefree model

```text
w^2 = t^4+3t^3+3t^2+1.
```

This is exactly the genus-one `j=0` third quotient above.  Independently, the
P3 numerator `n0` has a quadratic relation of `t`-degree 43 whose degree-86
discriminant has the same squarefree part.  The coordinates `n1,n2,n3,n5,n6`
are linear over `GF(7)(t,lambda)` within degrees 26--31; the remaining high
coordinates exceed the present order-230 support bound.  Thus the selected
P3 marking does not reveal a new genus-two field: it lies on the already
identified third quadratic cover.  Within this selected marking package the
only quadratic covers are therefore P1 (genus two, wrong Igusa invariants),
P2 (genus three), and P3/product (genus one, `j=0`).  A successful `Q79`
descent would require a different CM marking/orientation or a mechanism not
visible in these three selected covers.

The leading slope-`3` (`1/12` in characteristic zero) branch has now been
continued uniquely through order 230. No centered surface relation occurs in
degrees at most five. Over `GF(7)`, degree six initially gave seventeen
relations at order 214, all passing sixteen withheld coefficients. Split-prime
replay corrects the interpretation: at the ordinary primes
`31,73,79,127,151,193,199` and many others the same matrix has rank 195 and
kernel 15; `p=7` is exceptional with kernel 17, while `p=97,103` have kernel
16. Thus two of the `GF(7)` sextics are accidental and the seventeen-generator
ideal does **not** lift wholesale.

The generic fifteen-sextic modular candidate still has affine dimension one,
Hilbert polynomial `48*n-93` (degree 48, arithmetic genus 94), and an
irreducible `(P,D)` eliminant of total degree 32 and bidegree `(15,32)`.
(`p=31` exceptionally drops to total degree 31 and bidegree `(14,31)`.) This
is a stronger cross-prime modular candidate, not yet a characteristic-zero
ideal.  Singular local-delta calculations at two independent ordinary primes
now determine the normalization genus: at both `p=73` and `p=79`, the affine
delta contribution is 223 and the contribution at infinity is 241.  Thus the
irreducible degree-32 plane curve has total delta 464 and geometric genus

```text
(32-1)(32-2)/2 - 464 = 1.
```

The unmarked coefficient branch is therefore modularly an elliptic curve, but
point counts reject the expected quotient identification. At `p=73` the
normalized candidate has `#E=59`, trace `15`, while the two known source
elliptic quotients have traces `-9,7` (quadratic twists `9,-7`). At `p=127`
the candidate has `#E=106`, trace `22`, while the known traces are `8,0`
(twists `-8,0`). Hence it is not isogenous to either known elliptic factor at
two good primes. Conditional on the very strongly supported fifteen-sextic
finite-jet ideal being the true global branch, slope `1/12` is not the desired
`X(6,79)` component.

This retains slope `8/87` as the leading unmarked q80/H2 reconstruction
route, not as the recovered source polarization. Its unmarked surface
normalization is rational, and its level-79 marking may involve
the third pole section `P3`/`Q79`; the failed P1/P2 covers do not rule out that
different quadratic marked cover. A subsequent bounded order-230 calculation
does, however, show that the *selected* `P3` pole and numerator coordinates
lie on the same genus-one `j=0` quadratic field as the third P1/P2 quotient:
its squarefree model is `w^2=t^4+3t^3+3t^2+1`. Thus the current three-section
marking package still does not recover the known genus-two source. This is a
marking/orientation obstruction, not a rejection of the unmarked rational
surface branch; a different CM marking or direct `Q79` descent is required.
Full normalization of the rejected
genus-one candidate was stopped after sixty seconds; the local-delta and
point-count gates already settle the route comparison conditionally.

The direct `Q79` test has now also been carried out, rather than inferred from
the three separate covers.  Exact lattice transport gives

```text
Q79 = -3*G1 - 2*G2 + 4*G3.
```

At the good elliptic-base value `T=5`, the group law was evaluated directly on
the order-230 slope-`5` formal branch.  For all four independent sign classes
of `(G1,G2,G3)` up to simultaneous negation, `x(Q79)` has

```text
no rational Pade reconstruction with numerator,denominator degree <= 100,
no quadratic relation over the surface parameter of parameter degree <= 70.
```

The quadratic test uses at most 213 monomials, fits 214 coefficients, and
requires all sixteen remaining coefficients to vanish.  Thus simply reversing
the signs of the selected three sections cannot recover the desired genus-two
marking.  This is a bounded finite-jet rejection over `GF(7)`, not a global
characteristic-zero proof against the unmarked slope-`5` surface branch: a
different CM basis/marking, rather than a sign change within the selected
triple, remains possible.  The four machine records are
[`1,1,1`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-q79-cover-1_1_1.json),
[`1,1,-1`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-q79-cover-1_1_m1.json),
[`1,-1,1`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-q79-cover-1_m1_1.json), and
[`1,-1,-1`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-q79-cover-1_m1_m1.json).

The second apparent CM24 seed choice does not escape this gate.  Exact
function-field group law relates the finite-pole candidate 3 and the
infinity-pole candidate 5 by

```text
G2'=-G2,   G3'=G3-G2,
-3*G1-2*G2'+4*G3' = -3*G1-2*G2+4*G3.
```

They are two bases for the same marked `Q79` class, not two source branches.
This identity is now asserted by
[`verify_q80_cm24_seed_gf7.sage`](scripts/verify_q80_cm24_seed_gf7.sage).

The exact finite jets and exceptional `GF(7)` bounded sextic basis are
[`q80-cm24-slope-1-12-gf7-jet214.json`](../artifacts/generated-results/q80-cm24-slope-1-12-gf7-jet214.json),
[`q80-cm24-slope-1-12-gf7-jet230.json`](../artifacts/generated-results/q80-cm24-slope-1-12-gf7-jet230.json), and
[`q80-cm24-slope-1-12-gf7-sextic-ideal.json`](../artifacts/generated-results/q80-cm24-slope-1-12-gf7-sextic-ideal.json).
This is strong modular algebraization evidence, not yet a global
characteristic-zero identity.  No further unstructured jet enlargement is
justified before normalization or a global `Q79` gate.

## Reproduction

```bash
sage elkies-k3/scripts/verify_q80_to_rootless_path.sage
sage elkies-k3/scripts/classify_kumar_cm_frame_extensions.sage \
  --print-q60-in-q80 --print-markings-in-q80
sage elkies-k3/scripts/analyze_q80_rootless_first_neighbor.sage
sage elkies-k3/scripts/derive_q80_first_q4_pencil.sage
sage elkies-k3/scripts/analyze_q80_second_neighbor_chamber.sage
sage elkies-k3/scripts/search_q80_second_neighbor_rr.sage
sage elkies-k3/scripts/derive_q80_second_q4_pencil.sage
sage elkies-k3/scripts/derive_q80_third_q12_pencil.sage
sage elkies-k3/scripts/verify_q80_cm24_third_transport.sage
sage elkies-k3/scripts/analyze_q80_third_q12_cm24_marking.sage
sage elkies-k3/scripts/search_q80_third_q12_xgate_gf7.sage
sage elkies-k3/scripts/search_q80_third_q12_local_module_gf7.sage
sage elkies-k3/scripts/derive_q80_third_q12_local_gates.sage
sage elkies-k3/scripts/derive_q80_third_q12_cm24_pencil.sage
sage elkies-k3/scripts/analyze_q80_fourth_q12_chamber.sage
sage elkies-k3/scripts/analyze_q80_fourth_q12_cm24_marking.sage
sage elkies-k3/scripts/search_q80_third_child_polynomial_sections_gf73.sage --match-mw
sage elkies-k3/scripts/derive_q80_fourth_q12_local_gates_gf73.sage \
  --candidate 1 --selected-i2 5 --genus-gate --generic-only \
  --scan-degenerations --write-artifact
sage elkies-k3/scripts/analyze_q80_fourth_q12_moving_cubic_gf73.sage
sage elkies-k3/scripts/analyze_q80_fourth_q12_moving_cubic_gf73.sage --sections
sage elkies-k3/scripts/analyze_q80_fourth_q12_moving_cubic_gf73.sage \
  --discriminant --write-artifact
sage elkies-k3/scripts/analyze_q80_fifth_q4_chamber.sage
sage elkies-k3/scripts/analyze_q80_fifth_q4_cm24_readiness.sage
sage elkies-k3/scripts/build_q80_fifth_q4_local_module_ansatz_gf73.sage
sage elkies-k3/scripts/reconstruct_q80_fifth_q4_projection_gf73.sage \
  --minimum-withheld 8 --write-artifact
sage elkies-k3/scripts/reconstruct_q80_fifth_q4_marked_projection_gf73.sage \
  --mode pair14 --minimum-withheld 4 --write-artifact
sage elkies-k3/scripts/analyze_q80_fifth_q4_marked_jacobian_gf73.sage \
  --mode pair14 --write-artifact
sage elkies-k3/scripts/audit_q80_fifth_q4_pair14_twist_gf73.sage
sage elkies-k3/scripts/search_q80_fifth_q4_cm24_pair14_marking.sage \
  --start-half 0 --stop-half 797472 --a 2 --b 2 \
  --mw-coordinates 1,0 --mw-up-to-sign
sage elkies-k3/scripts/search_q80_fifth_q4_cm24_pair14_marking.sage \
  --start-half 0 --stop-half 797472 --a 2 --b 2 \
  --mw-coordinates 1,0 --mw-up-to-sign --target-root-data 15,78,384 \
  --output artifacts/local/q80-fifth-q4-cm24-pair01-marking.json
sage elkies-k3/scripts/reconstruct_q80_fifth_q4_marked_projection_gf73.sage \
  --mode pair01 --minimum-withheld 8 --write-artifact
sage elkies-k3/scripts/analyze_q80_fifth_q4_marked_jacobian_gf73.sage \
  --mode pair01 --write-artifact
sage elkies-k3/scripts/analyze_q80_alternate_final_q6_chamber.sage
sage elkies-k3/scripts/verify_q80_alternate_final_q6_nef.sage
sage elkies-k3/scripts/analyze_q80_alternate_final_q6_cm24_readiness.sage
sage elkies-k3/scripts/verify_q80_alternate_fifth_q6_rootless.sage \
  --write-artifact
sage elkies-k3/scripts/analyze_rank17_h8_q9_fibers.sage --write-artifact
sage elkies-k3/scripts/analyze_q80_alternate_fifth_rootless_bridge.sage
sage elkies-k3/scripts/search_q80_first_neighbor_rr.sage \
  --max-k 0 --discriminant-pairs-mod 101 --lift-pair-survivors
sage elkies-k3/scripts/analyze_q80_rank19_marked_cover.sage
sage elkies-k3/scripts/extend_q80_rank19_branches_gf7.sage \
  --slope 5 --order 90 --max-degree 0 --relation-summary-only \
  --validate-ideal artifacts/generated-results/q80-cm24-slope-8-87-gf7-ideal.json
sage elkies-k3/scripts/extend_q80_rank19_branches_gf7.sage \
  --slope 3 --order 230 --max-degree 0 \
  --jet-input artifacts/generated-results/q80-cm24-slope-1-12-gf7-jet214.json \
  --validate-ideal artifacts/generated-results/q80-cm24-slope-1-12-gf7-sextic-ideal.json
sage elkies-k3/scripts/analyze_q80_rank19_branch_ideal.sage \
  --input artifacts/generated-results/q80-cm24-slope-1-12-gf7-sextic-ideal.json \
  --plane-pair P D \
  --exclude-parameter artifacts/generated-results/q80-cm24-slope-8-87-gf7-parameter.json
sage elkies-k3/scripts/extend_q80_cm24_branch_modp.sage \
  --prime 73 --slope 1/12 --order 230 --relation-degree 6 \
  --relation-basis-output \
    artifacts/generated-results/q80-cm24-slope-1-12-gf73-sextic-rref.json
sage elkies-k3/scripts/analyze_q80_cm24_split_sextics.sage \
  --input artifacts/generated-results/q80-cm24-slope-1-12-gf73-sextic-rref.json \
  --plane-pair P D --plane-delta
sage elkies-k3/scripts/extend_q80_rank19_branches_gf7.sage \
  --slope 5 --order 230 --max-degree 0 --parameter-max-degree 0 \
  --jet-input artifacts/generated-results/q80-cm24-slope-8-87-gf7-jet230.json \
  --parameter-input artifacts/generated-results/q80-cm24-slope-8-87-gf7-parameter.json \
  --parameter-algebraic-coordinate lam \
  --parameter-algebraic-coordinate-degree 2 \
  --parameter-algebraic-max-parameter-degree 40 \
  --parameter-algebraic-validation 12
for signs in 1,1,1 1,1,-1 1,-1,1 1,-1,-1; do
  slug=${signs//-1/m1}
  slug=${slug//,/_}
  sage elkies-k3/scripts/extend_q80_rank19_branches_gf7.sage \
    --slope 5 --order 230 \
    --jet-input artifacts/generated-results/q80-cm24-slope-8-87-gf7-jet230.json \
    --parameter-input artifacts/generated-results/q80-cm24-slope-8-87-gf7-parameter.json \
    --parameter-max-degree 40 --q79-evaluation 5 --q79-signs "$signs" \
    --q79-degree-bound 100 --q79-algebraic-max-parameter-degree 70 \
    --q79-validation 16 \
    --q79-output \
      "artifacts/generated-results/q80-cm24-slope-8-87-gf7-q79-cover-${slug}.json"
done
```

Expected terminal line:

```text
Q80ROOTLESSPATH|terminal=rootless|MW=17|det=948|composite_transport_det=1|integral_isometry=1|transport_sha256=7116a499931bd096ba47fffc377a28690754bb451af0bdd7403f0e50438bd00d|status=PASS
```
