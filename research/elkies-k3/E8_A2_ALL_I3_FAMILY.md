# Universal all-`I3` family for the `E8+A2^3` frame

## Status

This note gives an exact four-dimensional ambient family with generic fibers

```text
II* + I3 + I3 + I3 + 5 I1,
```

the generic Kodaira realization of the abstract root frame `E8+A2^3`.  It
contains the discriminant-3 `II*+II*+IV` surface as a boundary point.  The
determinant-948 `X(6,79)` family should be a one-dimensional
Neron--Severi-enhancement locus inside it, but that locus is not yet solved.

Run:

```text
sage elkies-k3/scripts/verify_e8_a2_all_i3_family.sage
```

## Exact parameterization

Normalize the three prospective `I3` positions to `0,1,lambda` and put

```text
P = t(t-1)(t-lambda).
```

Let

```text
q = q0+q1*t+q2*t^2,
a = a0+a1*t,
b = b0+b1*t,
```

and impose the four coefficient equations in

```text
12*q*b - a^2 = c*P.
```

Define

```text
A = -3*q^2 + P*a,
B =  2*q^3 - P*q*a + P^2*b.
```

The verifier checks the polynomial identity

```text
4*A^3+27*B^2
 = P^3*(9*q^2*c + 4*a^3 - 54*q*a*b + 27*P*b^2)
   + 9*P^2*q^2*(12*q*b-a^2-c*P).
```

Consequently, on the declared coefficient locus,

```text
Delta = -16*P^3*R5,
R5 = 9*q^2*c + 4*a^3 - 54*q*a*b + 27*P*b^2.
```

On the open set where `q(0)q(1)q(lambda) != 0`, the three positions are
distinct, and `R5` is squarefree and coprime to `P`, the fibers are exactly
`II*+3 I3+5 I1`.  The degrees are `deg(A)=4`, `deg(B)=7`, and
`deg(Delta)=14`, so the fiber at infinity is `II*`.

## Why the parameterization is complete on the generic chart

At an `I3` fiber, the singular cubic has a double root.  Interpolate those
three node abscissas by the quadratic `q`.  Then `A+3q^2` vanishes at all
three roots of `P`, so it is `P*a` with `a` linear.  The first derivative of
the discriminant forces

```text
B = 2q^3-P*q*a  (mod P^2),
```

so the remaining quotient is `P^2*b` with `b` linear.  The second derivative
condition is exactly `12qb-a^2=0 (mod P)`, equivalently the displayed
constant-`c` identity.  Thus no generic all-`I3` model is lost by this form.

The variables `(lambda,q0,q1,q2,a0,a1,b0,b1,c)` satisfy four coefficient
relations, leaving five parameters.  Quotienting by Weierstrass scaling leaves
the expected four-dimensional `U+E8+A2^3` moduli.

## CM endpoint and a non-isotrivial check

The boundary specialization

```text
lambda=0,
q=0,
a=0,
b=t,
c=0
```

gives

```text
y^2=x^3+t^5(t-1)^2,
```

the Utsumi discriminant-3 model.  The verifier also checks a rational generic
point with `lambda=-1/2`: its residual quintic is squarefree and its
`j`-invariant is nonconstant.

## Target locus

The exact glue computation gives three required polynomial sections with

```text
height Gram = (1/3)*[[8,-1,0],[-1,10,0],[0,0,12]],
profiles    = (1,1,0), (0,2,0), (0,0,0),
P1.P2 = P1.P3 = P2.P3 = 2.
```

Inside this four-dimensional family, imposing those three independent MW
classes should cut out the one-dimensional determinant-948 locus.  This is
the primary reconstruction problem.  The two-parameter mixed family in
[`E8_A2_MIXED_FAMILY.md`](E8_A2_MIXED_FAMILY.md) remains a useful exact slice,
but it is not assumed to contain the target without passing the same lattice
gates.
