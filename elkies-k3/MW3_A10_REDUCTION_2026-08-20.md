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

## Canonical P3 boundary

The target third profile is `(10,2,0,1)` with `P3.O=1`.  On the fixed
two-section surface, the square-root recursion was exhausted exactly:

```text
840 admissible (pole, top-X, top-Y) slices,
31^3 lower-X choices per slice,
0 P3 hits.
```

Run

```bash
sage elkies-k3/scripts/search_mw3_a10_p3_sliced.sage
```

to reproduce the empty fixed-surface search.  This is a bounded finite-field
result only.  It does not show that the P1+P2 locus lacks P3 elsewhere.

## Next frontier

1. Continue complete P1 coordinate slices over `GF(31)`.
2. Reconstruct only exact semistable survivors and apply the four-variable P2
   high-jet/residual-square test.
3. Run the fast P3 square-root recursion only on P1+P2 hits.
4. Once a P1+P2+P3 point is found, verify the full height lattice and the
   reduced Jacobian, then lift the smooth one-dimensional chart at several
   good primes toward characteristic zero.

