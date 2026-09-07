# Parity obstruction to the all-IV `E8 + A2^3` lift

## Status

This is an exact obstruction.  It rejects the previously proposed Kodaira
lift

```text
E8 + A2^3  ->  II* + IV + IV + IV + II
```

for the Picard-rank-19 target.  It does **not** reject the exact abstract
`E8 + A2^3` root frame or its Mordell--Weil height lattice.

Run:

```text
sage elkies-k3/scripts/verify_e8_a2_all_iv_parity_obstruction.sage
```

## Obstruction

The proposed family was

```text
y^2 = x^3 + [t(t-1)(t-lambda)]^2(t-mu).
```

Its `A` coefficient is identically zero, so its geometric generic fiber has
`j=0`.  Over the algebraic closure it has the order-three automorphism

```text
rho: (x,y) |-> (zeta_3 x,y).
```

On the free Mordell--Weil group, `rho` satisfies
`rho^2+rho+1=0`.  Therefore `MW tensor Q` is a vector space over the quadratic
field `Q(zeta_3)`, and the geometric Mordell--Weil rank is even.

The target lattice has Picard rank `19` and root rank

```text
rank(E8 + A2^3) = 8 + 3*2 = 14.
```

Shioda--Tate therefore requires

```text
rank(MW) = 19 - 2 - 14 = 3,
```

which is odd.  This is impossible in the all-IV `j=0` family.

The same conclusion is visible directly in the discriminant:

```text
Delta = -432 t^4(t-1)^4(t-lambda)^4(t-mu)^2.
```

Together with the orders at infinity, this is generically exactly
`II* + 3 IV + II`; it is not merely a loose ambient family whose generic
member could evade the CM automorphism.

## What went wrong

The lattice calculation determines the reducible **root lattice**, not the
Kodaira symbol of every fiber.  An `A2` root lattice can come from either an
`I3` fiber or an `IV` fiber.  The CM endpoint has an `IV`, but a deformation
away from that endpoint need not preserve the additive Kodaira type.  Choosing
all three `A2` factors to be `IV` silently forced the entire family to remain
isotrivial.

Thus the exact results that remain valid are:

- the embedding of the determinant-948 rank-19 lattice into the
  discriminant-3 CM Neron--Severi lattice;
- the abstract inherited root frame `E8 + A2^3`;
- the exact reduced Mordell--Weil height Gram
  `(1/3)[[8,-1,0],[-1,10,0],[0,0,12]]`.

The invalid promotion is the equation `A=0` and every section search derived
from it.

## Correct backtrack

For a short K3 Weierstrass model with `II*` at infinity, one may take

```text
deg(A) <= 4,   deg(B) <= 7.
```

Each of the three `A2` factors must then be allowed to be `I3` or `IV`.
The all-IV case is excluded, leaving the distributions with zero, one, or two
`IV` fibers.  Their remaining Euler number is supplied by `I1` and/or `II`
fibers.  The correct component must also specialize to the CM endpoint where
the fibers coalesce to `II* + II* + IV`.

The reconstruction should impose the full rank-three Mordell--Weil
height/glue data on this non-isotrivial `II*` family, rather than impose one
arbitrary section on the obstructed `j=0` chart.  Any candidate should then be
cross-checked against the exact rank17-to-E6 neighbor transport in
[`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md).

The first such ambient family has now been derived exactly, with fibers
`II*+2 IV+I3+3 I1`; see
[`E8_A2_MIXED_FAMILY.md`](E8_A2_MIXED_FAMILY.md).
