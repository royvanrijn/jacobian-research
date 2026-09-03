# The noncyclic `4A1/MW13` equation from published R17 (2026-09-04)

## Result

On the pinned Picard-rank-19 K3, the noncyclic local-bridge marking now has an
explicit Jacobian equation over `QQ`.  In the affine base coordinate `s`, put

```text
E_4A1 : Y^2 = X^3 + 27 A0(s) X + 54 B0(s),
```

where coefficients below are in increasing degree order:

```text
A0 = [
 -7968955746765484566833424282663257363523567616,
 -153919846221500919963993931568646348971756961792,
 -129142324346140457923598416750959664469822635008,
 11856854430268329882477628495139189406928786290176,
 -63309720222430457255163477821658255908032564929920,
 177807653884304174205700570072935249832394824567296,
 -254924855689066400418140235174994303557343394408128,
 42594334501111913045987234379105874303509452716688,
 -3004345002478166518709556809722984087690032712321
]

B0 = [
 711380772410965629657989656618818014147094820363975853986255351054336,
 20610407920549344679927468293689618938131281671777825588650900505755648,
 62494889038397708535098517564164594424322483974386763442826405660983296,
 -850251054016803774513459771444966131266092450730274622424405691867627520,
 -633130068562476676281468643994094210232593012156634079870959979110778880,
 -62139672050561740096107909398780717981424203203033139427862957609181171712,
 728850949667724833184808927668698950858311284337190126317585339999402914304,
 -2752420714212142383473939859499819508928329788058186091597704635630810129152,
 4384829038234510245832354585481506773316664029255203923085345040662070777920,
 -3057081513463253665420576681717545201254749831829029359351316198309814398720,
 1708361167923651471240819383983020688934380800141185460463391530705986715776,
 -110743436948332521165163766616311413566971179202301920828338059997087278552,
 5207445156730776117703110077692175346910536357664810150755200591871698881
]
```

The singular fibres are

```text
s = 0, 1, -146234/269481, infinity : I2,
sixteen conjugate finite bases             : I1.
```

Thus the complete configuration is `4 I2 + 16 I1`, its reducible root system
is `4A1`, and its Mordell--Weil group over `QQ(s)` is

```text
MW(E_4A1/QQ(s)) = Z^13.
```

The machine-readable equation, thirteen rational sections, coordinate maps,
frame, and reverse formulas are in
[`elkies-k3-r17-noncyclic-4a1-direct-fibration-v1.json`](../artifacts/generated-results/elkies-k3-r17-noncyclic-4a1-direct-fibration-v1.json).

## The literal marked `U`

Use the ambient basis `(F,F+O)` followed by the short positive R17 frame.  The
new splitting is

```text
D      = [3,2,-1,0,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,0],
D+O'   = [4,3,-2,0,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,0].
```

Its Gram matrix is `U`, the full splitting transport has determinant `1`, and
the local geometric certificate proves that `D` is nef and `O'` is an
effective physical degree-one curve.  Relative to the old splitting, the
cross-pairing matrix is

```text
A = [[2,3],[3,4]].
```

The saturated bridge has Gram `diag(4,8)`, saturation index one,
discriminant group `Z/4 + Z/8`, and order `32`.  It is maximal.  This is the
noncyclic profile that distinguishes the example from the earlier cyclic
controls.

## Degree-two compilation

The trace section in the short R17 frame is

```text
r = [-1,0,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,0].
```

It has height `12` and intersects the old zero four times.  On Elkies's
published model

```text
y^2 = x^3 - 27 S(t) x + (27/4) T(t),
```

the universal nonzero-trace chord compiler writes sections of `O(D)` as

```text
a(t) + b(t) (y+y(r))/(x-x(r)),   deg(a)<=7, deg(b)<=1.
```

The eight resolved congruences modulo the square of the trace denominator
have rank eight in ten unknowns.  Their two-dimensional kernel gives the
pencil `u=L1/L0` and an exact binary quartic `W^2=q(t,u)`.  Pointing that
quartic at the certified curve `O'` gives its Jacobian with the required new
origin.

A rational `PGL2` change sends three double discriminant bases to
`0,1,infinity`; the fourth becomes `-146234/269481`.  Removing the common
weighted rational square leaves twist class `3` and gives the integral model
displayed above.  Its finite discriminant has factor profile

```text
(linear)^2 (linear)^2 (linear)^2 (irreducible degree 16),
```

the degree-16 factor is squarefree, and `gcd(A,Delta)=1`.  At infinity the
orders of `(c4,c6,Delta)` are `(0,0,2)`.  These checks prove `4 I2 + 16 I1`
rather than merely predicting `4A1` from the lattice.

## Frame and arithmetic Mordell--Weil group

The orthogonal positive frame has determinant `948`.  Exact norm-two
enumeration gives eight signed roots, and four displayed simple roots in the
artifact have Gram `2 I4`; their span is primitive.  Integral quadratic-form
isometry identifies this frame with the local-mutation target.  Separate exact
isometry tests reject both stored historical H3 `4A1` frames: the current
suffix frame and the physical-`q8` orbit-376 frame.

Thirteen old R17 sections restrict through the chord construction to explicit
`QQ(s)` points on `E_4A1`.  The artifact stores every source lattice vector,
base restriction, and new `(X,Y)`, and verifies every equation identity.  The
four simple roots followed by these thirteen section vectors form a
determinant-one basis of the frame.  Their Shioda height matrix has rank `13`
and determinant

```text
237/4 = 948/2^4.
```

Because the root span is primitive, the trivial lattice is primitive and the
torsion subgroup is zero.  Shioda--Tate gives geometric rank
`19-2-4=13`; the thirteen rational sections attain that bound.  Hence the
arithmetic conclusion over `QQ(s)` is exactly `Z^13`, not only a geometric
rank statement.

## Target-free reverse hop

The reverse marked-`U` control requests only `root_rank=0` and
`ADE=rootless`; it supplies no target frame Gram.  In its exact one-witness
catalog it selects the published R17 splitting and restores the published
zero.  The equation certificate independently checks the generic inverse:
the pointed-quartic formulas recover `t` and `W`, the residual chord recovers
the literal old `x` and `y`, and `L1/L0` recovers `u` (therefore `s`).  The
recovered coordinates satisfy the published R17 equation identically.

The certified cycle is therefore

```text
published R17 equation
  -> degree-two noncyclic 4A1/MW13 equation
  -> target-free rootless marked-U selection
  -> literal published R17 equation.
```

## Replay and boundary

```bash
sage -python elkies-k3/scripts/compile_r17_noncyclic_4a1_qq.sage
sage -python elkies-k3/scripts/compile_r17_noncyclic_4a1_qq.sage --check
```

The checker source has SHA-256
`1929e7df84b1fbc5fb94f481ef80015516a0559f0f3d1b59ffc6adb0b818184a`;
the generated artifact has SHA-256
`c8f259c90aae27f33a0f39dc544da541d3c7db5af5002c008d2d46871dbbd32c`.

This is the first explicit equation in the repository for this nonhistorical
noncyclic `4A1` frame.  It does not claim external priority, classify its `J1`
surface-automorphism orbit, or search for rank-jumping specializations.

<!-- status-consumer: EC-K3-R17-NONCYCLIC-4A1-DIRECT-EQUATION f657620e07f8f3f0 -->
