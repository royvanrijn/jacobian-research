# A10/MW3 normalized reduction and two-section seed — 2026-08-20

This note records the current exact finite-field frontier on the preferred
semistable reconstruction branch

```text
ADE = A10 + A2 + A1^2
fibers = I11 + I3 + I2 + I2 + 6 I1
MW rank target = 3.
```

All results in this note are computations over `GF(31)`.  They prove
nonemptiness and smoothness of modular charts; they do not yet give a
characteristic-zero K3 family or a high-rank elliptic curve over `Q`.

## Structural P1 reduction

The former fraction-field chain reached 17 variables but produced a 516 KB
checkpoint after only eighteen eliminations.  The cleaner representation
defines

```text
B = Y^2 - X^3 - A*X
```

from the section identity and uses the local node branches directly.

At the `I3` fiber, the low section/discriminant equations reduce to

```text
(a1 + 6*s0*x1)^2 = 12*s0*y1^2.
```

On the nonzero tangent branch, the remaining Weierstrass scaling normalizes

```text
s0 = 3,
a1 + 18*x1 = 6*y1.
```

At infinity the identical tangent calculation gives

```text
sinf = 3*rho^2,
a7 + 6*sinf*x3 = 6*rho*y5.
```

Parametrizing the nonsingular specialization of `P1` at `t=1`, imposing the
node incidence at `t=lambda`, and solving three sparse infinity equations
reduces the P1 locus to

```text
9 variables / 6 equations / expected dimension 3.
```

The deeper `I11` component-2 condition has a simple blow-up interpretation.
In the infinity chart, it is

```text
y5 = 0,
y4 != 0.
```

After `y5=0`, `Delta_18` vanishes identically and `Delta_17/y4` is
affine-linear in `y3`.  Eliminating `y3` leaves the final chart

```text
variables = rho,r1,s1,lambda,x2,x3,y4
equations = 4
expected dimension = 3.
```

The `I2(1)` derivative equation also contains the excluded node factor
`r1^2-3*s1`; the normalized builder divides it out rather than retaining a
boundary component.

Build the chart with

```bash
sage elkies-k3/scripts/build_mw3_a10_p1_normalized.sage \
  --p 31 --stage component2 \
  --export artifacts/local/elkies-k3/mw3-a10-p1/p31-component2.ms
```

The separate `.open.ms` output records the small nonzero chart factors.  It is
essential: clearing the `y3` denominator without this filter creates hundreds
of numerator-only false hits.

## Complete coordinate slices

Fixing `(rho,r1,lambda)` leaves four coordinates.  A vectorized exhaustive
scan tests at most `31^4` points, filters by the sparse equations first, and is
faster and more reliable here than a saturated Groebner basis.

For seed 1,

```text
(rho,r1,lambda) = (4,18,27),
```

the exact open chart has two points.  Both reconstruct to smooth semistable
surfaces with exact valuations `(3,2,2,11)` at
`(0,1,lambda,infinity)` and a squarefree residual degree-six discriminant.

The second point from a later slice,

```text
(rho,r1,s1,lambda,x2,x3,y4) = (19,8,21,23,0,12,1),
```

is the current two-section seed.  Its surface and first section are

```text
A = 4 + 23*t + 18*t^2 + 28*t^3 + 12*t^4 + 18*t^5
      + 23*t^6 + 20*t^7 + 19*t^8,

B = 23 + 24*t + 27*t^2 + 10*t^3 + 6*t^4 + 8*t^5
      + 23*t^6 + 26*t^7 + 26*t^8 + 16*t^9
      + 19*t^10 + 9*t^11 + 15*t^12,

P1:
X1 = 3 + 9*t + 12*t^3 + 29*t^4,
Y1 = 5*t + 21*t^2 + 12*t^3 + t^4.
```

## Canonical P2

For `P2.O=1`, put `z=t-r` and interpolate the finite nonidentity component
conditions in a degree-six numerator `X2`.  Infinity component
`6 ~ -5 mod 11` forces `Y2` to have degree at most four.  Equivalently,

```text
H = X2^3 + A*X2*z^4 + B*z^6
```

has degree at most eight.  The eight high-coefficient equations are a tiny
four-variable system in `(r,q0,q1,q2)`.

On the surface above, one solution completes to an exact square:

```text
r = 27,
X2 = 17 + t + 20*t^2 + t^3 + 27*t^4 + 27*t^5 + 29*t^6,
Y2 = 21*t + 6*t^3 + 4*t^4.
```

Thus

```text
P2 = (X2/(t-27)^2, Y2/(t-27)^3)
```

is an exact second section.  The component data give self-heights

```text
h(P1) = 79/66,
h(P2) = 106/66.
```

Their ratio is `106/79`, not a square in `Q`; hence `P1` and `P2` are
independent in the Mordell--Weil group modulo torsion.  A direct function-field
group-law regression finds no relation for `|m|,|n| <= 16`.

Verify the complete seed with

```bash
sage elkies-k3/scripts/verify_mw3_a10_p1p2_gf31.sage
```

## Complete canonical P3 audit on both P1+P2 surfaces

The P2 scan has two exact P1+P2 surfaces (the two signs of a given P2 are one
X-solution on each).  The target third profile is `(10,2,0,1)` with
`P3.O=1`.  A two-sided meet-in-the-middle recursion exhausts its coefficient
system in about ten seconds per surface, rather than scanning a large symbolic
system.

```text
surface 2: 96,868,800 joins; 0 raw hits; 0 genuine hits
surface 4: 96,868,800 joins; 58 raw hits; 58 open rejections; 0 genuine hits
```

Every surface-4 raw hit has

```text
(t-r)^2 | X_raw,
(t-r)^3 | Y_raw,
```

so the declared pole cancels.  All 58 reduce, up to sign, to the same genuine
polynomial section

```text
R:
X = 3 + 8*t + 3*t^2 + 9*t^3 + 29*t^4,
Y = 2*t + 26*t^2 + 17*t^3 + t^4 + 16*t^5.
```

The stacked finite-quotient certificate in
[`scripts/verify_mw3_a10_p1p2_gf31.sage`](scripts/verify_mw3_a10_p1p2_gf31.sage)
proves that `P1,P2,R` are Z-independent.  This is an exact rank-at-least-three
statement for the modular surface, but `R.O=0`; it is not the requested
canonical `P3.O=1` and does not realize the target height lattice.

Reproduce the exhaustive gate with

```bash
sage -python elkies-k3/scripts/build_mw3_a10_p3_fixed.sage \
  --surface 2 --out artifacts/local/elkies-k3/mw3-a10-p1/p3-surface2.ms
sage -python elkies-k3/scripts/search_mw3_a10_p3_meet.py \
  --input artifacts/local/elkies-k3/mw3-a10-p1/p3-surface2.ms \
  --lambda 27 --nodes 3,22,4 --sinf 17
```

and analogously with `--surface 4`, `--lambda 23`, `--nodes 3,21,10`,
and `--sinf 29`.  The search kernel enforces pole noncancellation and the
identity-component open condition at `t=1`.

This is a bounded finite-field result only.  It closes the two current fixed
surfaces, not the full P1+P2 locus.

## Exact target triple over `GF(23)`

The prime-31 search was extended with persistent P1 reconstruction, an exact
oriented P2 profile gate, exhaustive P3 square-root slices, and a direct
function-field Shioda-pairing verifier.  Across the completed ranges through
seed 13315 it found 147 canonical P2 records and 59 P3 candidates, but no
target Gram matrix.  This is a bounded negative result over `GF(31)`, not an
obstruction in characteristic zero.

Changing the auxiliary prime to 23 immediately found the missing branch.  In
the first 1000 seeded slices there were 636 P1 points, 12 canonical P2 target
records, and seven P3 candidates.  Seed 17, hit 1 gives

```text
A = 19 + 15*t + 6*t^3 + 21*t^4 + t^5 + 17*t^6 + 2*t^7 + 10*t^8,
B = 8 + t + 12*t^2 + 3*t^3 + 4*t^4 + 22*t^5 + 20*t^6
      + 18*t^7 + 13*t^8 + 10*t^9 + 14*t^10 + 5*t^11 + 9*t^12,

P1:
X1 = 3 + 15*t + 7*t^2 + 17*t^3 + 9*t^4,
Y1 = 13*t + 3*t^2 + 13*t^3 + 5*t^4,

P2, z2=t-11:
X2 = 18 + t + 6*t^2 + 3*t^3 + 18*t^4 + 3*t^5 + 9*t^6,
Y2 = 7*t + 22*t^2 + 10*t^3 + 7*t^4,

P3, z3=t-8:
X3 = 8 + 11*t + 17*t^2 + 14*t^3 + t^4 + 9*t^5 + 9*t^6,
Y3 = 5*t + 18*t^2 + 8*t^3 + 18*t^4 + 11*t^5
      + 5*t^6 + 21*t^7 + 4*t^8.
```

Here `Pi=(Xi/zi^2,Yi/zi^3)` for `i=2,3`.  Exact group-law component
orientation and Shioda heights give

```text
profiles = (2,1,0,1), (6,2,1,1), (10,2,0,1),
P1.O,P2.O,P3.O = 0,1,1,

66 * height_gram =
[ 79   17   -1 ]
[ 17  106   19 ]
[ -1   19  259 ],

det(height_gram) = 79/11.
```

Thus the intended reduced Mordell--Weil lattice is now realized exactly on a
finite-field K3.  This does not yet construct the family or a rank-17 curve
over `Q`.

The combined P1/P2/P3 deformation system has 32 coordinates and 44 retained
equations.  At this point its Jacobian has rank 31, so the tangent space is
one-dimensional and smooth.  The tangent has nonzero `rho` component; fixing
`rho=16` gives a transverse full-rank 31-variable system.  Overdetermined
Hensel lifting succeeds through 20 base-23 digits, but rational reconstruction
at that fixed value does not yield a rational point.  The formal lift is a
constructive local certificate, not yet a global `Q`-point.

Reproduce the local deformation and Hensel audit with

```bash
sage elkies-k3/scripts/analyze_mw3_a10_target_lift_gf23.sage \
  --hensel-digits=20
```

## Next frontier

1. Exhaust the 10,626 normalized `(rho,r1,lambda)` slices over `GF(23)` and
   recover all canonical target triples on this one-dimensional locus.
2. Fit low-degree plane projections of the resulting finite-field curve and
   determine its genus/rational parametrization before doing more p-adic work.
3. Lift a characteristic-zero curve model (not merely a single formal point),
   verify the generic sections, and identify rational specializations.
4. Only then run specialization/rank-jump and conductor searches toward the
   rank-21 low-conductor or rank-31 targets.
