# Orbit-96 `A7+D7` Weierstrass equation and arithmetic representation — 2026-09-02

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT96-A7D7-GALOIS ba008502f0e5533f -->

## Result

The genuine orbit-96 pencil is now equation-explicit.  In the old basis

```text
O,F,E6a_1,...,E6a_6,E6b_1,...,E6b_6,A3_1,A3_2,A3_3,P0,P1
```

its fibre is

```text
D96 = 2*P0
    + (2,1,2,2,1,0)_E6a
    + (1,0,0,0,0,0)_E6b
    + A3_2.
```

Put

```text
D0     = 3*k^2-4,
c      = 2*k/D0,
lambda = -(k^2-4)*D0/4,
H      = 1-t^2,

x0 = -H^2,
y0 = c*lambda*H^2,
m  = ((1-3*k^2/4)/k)*(t^2-1).
```

Then the complete resolved Riemann--Roch space is

```text
H0(X,O(D96)) = <1,z>,

z = (t+1)*(y+y0+m*(x-x0))/(x-x0)^2.
```

The corrected elimination gives a binary quartic in `t`, not the cubic in
the superseded audit.  Its Jacobian is

```text
y^2 = x^3 - 27*I6(z)*x - 27*J9(z),
```

where the coefficients below are listed from constant term upward:

```text
I6 = [
  1,
  (12*k^2-16)/k,
  (36*k^4-192*k^2+64)/k^2,
  (32*k^4-704*k^2+768)/k,
  48*k^4-256*k^2+2560,
  768*k^3-3072*k,
  64*k^6-512*k^4+1024*k^2
]

J9 = [
  -2,
  (-36*k^2+48)/k,
  (-216*k^4+864*k^2-384)/k^2,
  (-528*k^6+5568*k^4-6912*k^2+1024)/k^3,
  (-720*k^6+14208*k^4-45312*k^2+18432)/k^2,
  (-864*k^6+8064*k^4-102912*k^2+116736)/k,
  -960*k^6+768*k^4+21504*k^2+184320,
  1152*k^7-10752*k^5+141312*k^3-466944*k,
  18432*k^6-147456*k^4+294912*k^2,
  1024*k^9-12288*k^7+49152*k^5-65536*k^3
].
```

The exact fibres are

```text
z=0          I8       A7,
z=infinity   I3*      D7,
seven roots  I1.
```

Thus the root data are `(14,140,32)`, the Euler sum is `24`, and
Shioda--Tate on the certified Picard-rank-19 K3 gives geometric Mordell--Weil
rank three.  The root lattice is primitive, so torsion is trivial, and the
free height lattice is

```text
diag(3/8,1,3),       regulator 9/8.
```

The exact replay is
[`scripts/compile_e6a1_rho19_orbit96_rr_galois.sage`](scripts/compile_e6a1_rho19_orbit96_rr_galois.sage),
with generated certificate
[`../artifacts/generated-results/elkies-k3-e6a1-rho19-orbit96-rr-galois-v1.json`](../artifacts/generated-results/elkies-k3-e6a1-rho19-orbit96-rr-galois-v1.json).

## Physical `E6` component marking

Use the ordinary-resolution component order

```text
(leaf,+outer,-outer,+inner,-inner,central).
```

The coordinate valuation cycles are

```text
ord(u) = (2,1,1,2,2,3),
ord(x) = (2,2,2,3,3,4),
ord(y) = (3,2,2,4,4,6).
```

The marked section has `(x/u^2,y/u^2)=(-1,c*lambda)` and therefore meets
`+outer`.  This removes the `E6` arm involution.  The map from the old Sage
`E6` root order to the chart order is

```text
(1,2,3,4,5,6)_abstract -> (2,1,4,6,5,3)_chart.
```

For the unregularized tangent coordinate the exact numerator and denominator
orders are

```text
numerator   (3,2,4,4,5,6),
denominator (4,4,4,6,6,8).
```

Hence its pole cycle is `(1,2,0,2,1,2)` in chart order, namely
`(2,1,2,2,1,0)` in the old root order.  Multiplication by `t+1` subtracts the
`u` cycle at `t=-1` and leaves only `(0,1,0,0,0,0)` in chart order, namely
`(1,0,0,0,0,0)` in the old root order.  At infinity the resolved I4
side charts are regular and the middle chart has order `-1`.  These are
exactly the three vertical pieces of `D96`.

The lattice certificate already proves that `D96` is primitive, nef and
isotropic.  Therefore `h0(D96)=2`, so the two displayed sections are the full
pencil rather than a generic-fibre ansatz.

## Why the earlier audit failed

The tangent slope simplifies to a polynomial in `t`, but the old audit kept it
in `Frac(QQ(k,z)[t])`.  That promoted the quadratic in the old `x` coordinate
to a polynomial over a fraction-field coefficient tower.  In that parent,
Sage's zero-argument `quadratic.discriminant()` returned a different
expression from

```text
quadratic[1]^2 - 4*quadratic[2]*quadratic[0].
```

The spurious expression had squarefree degree three and reproduced the old
`2E6+A3` surface.  Coercing `m` back to the polynomial ring before elimination
makes the two discriminant calculations agree.  The square factor is then
`(t+1)^2`, the squarefree residual has degree four, and its Jacobian is the
`A7+D7` model above.

This is an instance of the general equation-level neighbour discipline used
in explicit elliptic K3 constructions; compare Kumar's explicit fibration
methods in [*K3 surfaces associated with curves of genus two*](https://arxiv.org/abs/math/0701669)
and the geometric fibration constructions of Garbagnati--Salgado in
[*Elliptic fibrations on K3 surfaces with a non-symplectic involution*](https://arxiv.org/abs/1806.03097).
Neither paper contains this orbit-96 calculation.

## Fields of definition

Take the old component `E6a_6` as the new zero.  Exact saturation of the child
root lattice gives three generator section classes with height basis
`diag(3/8,1,3)`.  In the old NS basis they are

```text
G1 = [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],

G2 = [0,0,4,2,4,4,2,1,2,0,-1,0,0,0,1,2,0,4,0],

G3 = [0,-1,6,3,6,6,3,1,3,0,-1,0,0,0,0,3,0,6,1].
```

Over `QQ(k)`, the two split IV* fibres, their components, and the old marked
sections are rational.  The only nontrivial constant-field action swaps the
two side components `A3_1` and `A3_3`; their tangent cone is
`eta^2+3*xi^2`.  On the saturated orbit-96 MW basis this conjugation is

```text
diag(1,-1,1).
```

Consequently

```text
G1  is defined over QQ(k)(z),
G2  is defined over QQ(k,sqrt(-3))(z) and is anti-invariant,
G3  is defined over QQ(k)(z).
```

In representation notation,

```text
MW_barQQ tensor QQ = 1 + chi_{-3} + 1,
rank MW_QQ(k)(z)   = 2.
```

Orbit 96 therefore does not supply a third rational direction.  Its quadratic
character is also the same `-3` character as orbit 103, not a new character
that can be combined independently by specialization or base change.  Its
remaining advantage is geometric: two reducible fibres rather than four, and
Weierstrass degrees `(6,9)` rather than `(8,12)` in the displayed charts.

## Reproduction and boundary

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit96_rr_galois.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit96_rr_galois.sage --check
```

Pinned SHA-256 values are

```text
checker   08c2629837f74aa90f9b09df9fe68104f542b226af892d056bff6c08577c01d0
artifact  f72b312561b01cbce6e9fa47f6c96a6ff6e1ec7aca406169f79434995cc6046f
```

The physical divisor, resolved pencil, quartic, Weierstrass equation, fibre
classification, saturated height basis, generator divisor classes, and
constant-field Galois representation are exact.  Polynomial Weierstrass
coordinates for all three generator sections and ranks of individual rational
specializations are not asserted.
