# A viable mixed `I3/IV` family through the CM endpoint

## Status

This note gives an **exact two-parameter ambient family** with abstract root
frame `E8+A2^3`, nonconstant `j`, and the correct discriminant-3 CM endpoint.
It is a viable replacement for the parity-obstructed all-IV family.  It has
not yet been proved to contain the determinant-948 `X(6,79)` curve; that is
the next lattice-polarization problem.

Run:

```text
sage elkies-k3/scripts/verify_e8_a2_mixed_family.sage
```

## Direct derivation

Choose the admissible Kodaira distribution

```text
II* + IV + IV + I3 + 3 I1.
```

Put the two `IV` fibers at `t=0,1` and write `D=t(t-1)`.  A normalized short
Weierstrass model has

```text
A = a D^2,
B = D^2 C(t),
```

where `C` is monic cubic.  Apart from the forced fourth powers of `D`, the
discriminant is controlled by

```text
H = 4 a^3 D^2 + 27 C^2.
```

To make `t=lambda` an `I3` fiber, require `H` to have a triple root there.
Introduce `s` with `27s^2=-4a^3`.  Then

```text
H = 27(C-sD)(C+sD).
```

On the generic branch only one cubic factor vanishes at `lambda`; because it
is monic, it must equal `(t-lambda)^3`.  The cusp relation is parametrized by

```text
a = -3r^2,   s = 2r^3.
```

This gives the explicit family

```text
D = t(t-1),
A = -3r^2 D^2,
B = D^2((t-lambda)^3 - 2r^3 D),

y^2 = x^3 + A x + B.
```

No Groebner basis or long search is involved.

## Exact discriminant and fibers

The verifier checks

```text
Delta = -432 D^4 (t-lambda)^3
              ((t-lambda)^3 - 4r^3 D).
```

On the open set where

- `r*lambda*(lambda-1) != 0`, and
- the displayed residual cubic is squarefree,

the fibers are exactly:

```text
t=infinity : II*
t=0,1      : IV, IV
t=lambda   : I3
three roots of the residual cubic : I1, I1, I1.
```

The degrees are `deg(A)=4`, `deg(B)=7`, `deg(Delta)=14`, giving orders
`(4,5,10)` at infinity.  The verifier also checks that the generic
`j`-invariant is nonconstant, so the CM rank-parity obstruction does not
apply.

## CM endpoint

At

```text
r=0, lambda=0
```

the family becomes

```text
y^2 = x^3 + t^5(t-1)^2,
```

the affine Utsumi No.1 model with fibers `II*+II*+IV`.  Geometrically, the
`I3`, one `IV`, and the residual `I1` fibers coalesce into the second `II*`.

## Remaining target test

The family has the right endpoint and the right abstract root lattice, but
those facts do not yet identify its desired one-dimensional sublocus.  The
complete rank-three Mordell--Weil height/component glue is now explicit.  Run

```text
sage elkies-k3/scripts/recover_rank3_mw_via_ns_glue.sage
```

In a reduced basis it gives

```text
height Gram = (1/3) * [[8,-1,0],[-1,10,0],[0,0,12]],

P1 = (1,1,0),
P2 = (0,2,0),
P3 = (0,0,0)
```

for the three `A2` factors.  The first two factors are the two components that
split from the old `E8`; they may be swapped and each component generator may
be inverted.  The third is the persistent `A2` from the CM `IV` fiber.

In the mixed family, the natural geometric assignment is therefore the two
broken factors at `t=0` and `t=lambda`, and the persistent factor at `t=1`.
The component-nonzero counts are `(2,1,0)`.  Since every nonzero `A2` label
contributes `2/3`, Shioda's self-height formula shows that all three reduced
generators have `P_i.O=0`; they may all be sought as polynomial sections.
The off-diagonal height formula then gives

```text
P1.P2 = P1.P3 = P2.P3 = 2.
```

The next exact gate is to impose these three section identities and their
pairwise intersections, then prove that the resulting component has
Neron--Severi determinant `948`
and the `X(6,79)` discriminant form.  A promoted candidate must then agree
with the E6/rootless frames through the certified transport in
[`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md).
