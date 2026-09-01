# Q80 third-q12 exact-lift checkpoint

## Status

The fixed `u=-2` construction has passed the horizontal lift gate.

- The exact `p=19` child remains fully certified: resolved genus-one pencil,
  Jacobian and maps in both directions, minimal `I6+I4+3I2+8I1` model, and
  transported `A5+A3+3A1` marking.
- The complete finite-field pipeline also passes at
  `p=61,67,83,89,103,131` with the same discrete profiles.
- Exact Hensel reconstruction proves that the characteristic-zero horizontal
  is defined over a biquadratic field, not one quadratic field.
- The exact horizontal passes a direct test at the previously reserved good
  prime `p=71` for all four sign conjugates.

The active gate is now the complete connected correction over that exact
biquadratic field, followed by its genus-one Jacobian. No characteristic-zero
child equation or Mordell--Weil rank is asserted at this checkpoint.

## Exact closure operands

The six characteristic-zero polynomial-closure equations are even in the
leading coordinate `l`. Replacing `l^2` by `q` leaves five rational unknowns.
At each of the two `p=19` operand branches, the resulting five-by-five
Jacobian is nonsingular. Newton--Hensel doubling and exact rational
reconstruction produce two literal solutions of the original six equations.

Write their reconstructed leading squares as `q1` and `q2`. Exact square
tests prove that `q1`, `q2`, and `q1/q2` are nonsquares over `QQ`. Hence the
two operands generate the degree-four field

```text
K = QQ(a,b),   a^2=q1,   b^2=q2.
```

The certificate includes exact substitution into all six closure equations
and literal reduction to the original `p=19` operand pair:

```text
artifacts/generated-results/
  q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json
```

Replay it with:

```bash
sage elkies-k3/scripts/lift_q80_third_q12_closure_operands_p19_qq.sage \
  --biquadratic-operands
```

The largest reconstructed coefficient has 36,335 bits. These sizes explain
why a direct low-height rational search did not expose the lift.

## Exact horizontal and held-out prime

The two reconstructed operands can be composed exactly on the original
characteristic-zero Weierstrass model. In the basis `(1,a,b,a*b)`, the
result has

```text
x in span(1,a*b),   y in span(a,b).
```

Literal substitution proves the characteristic-zero Weierstrass identity.
This is a degree-four point in general, although its `x` coordinate descends
to the third quadratic subfield `QQ(a*b)`.

The same certificate reduces the exact object directly at the previously
reserved `p=71`. There both `q1` and `q2` are nonsquares while `q1*q2` is a
square. All four sign conjugates pass:

- `P.O=2`;
- height `8` by both height computations;
- the finite `I1*` identity component condition;
- the expected numerator and denominator degree profiles.

Modulo section sign, the four conjugates give two unsigned classes. The
certificate is:

```text
artifacts/generated-results/
  q80-third-q12-um2-biquadratic-horizontal-qq.json
```

Replay it with:

```bash
sage elkies-k3/scripts/certify_q80_third_q12_biquadratic_horizontal_qq.sage
```

This proves the exact horizontal and an independent held-out-prime
realization. It does not yet prove the exact resolved pencil or child.

## Why the earlier CRT route is retired

The exact square classes have the following local characters, listed as
`(q1,q2,q1*q2)`:

| prime | characters | inert operand |
| ---: | :---: | :---: |
| 19 | `(-,+,-)` | `q1` |
| 61 | `(+,-,-)` | `q2` |
| 67 | `(+,-,-)` | `q2` |
| 83 | `(-,+,-)` | `q1` |
| 89 | `(+,-,-)` | `q2` |
| 103 | `(-,+,-)` | `q1` |
| 131 | `(+,-,-)` | `q2` |

Thus the unique local quadratic target alternates between the two independent
global square classes. Local trace, norm, and coefficient discriminant are
generator-free at one prime, but they are not all reductions of a single
quadratic conjugation quotient across this mixed prime set. The old CRT
integers therefore cannot be used as rational reconstruction candidates.

The seven-prime artifact remains useful as a literal residue and branch
diagnostic. It is now schema-versioned and status-labelled accordingly:

```text
artifacts/generated-results/
  q80-third-q12-um2-frobenius-crt-interface.json
```

```bash
python3 elkies-k3/scripts/compile_q80_third_q12_frobenius_crt_interface.py
```

It accumulates 1,947 ordered local slots modulo `7739891239523`, but its claim
boundary explicitly excludes rational reconstruction.

## Active exact-child path

The immutable `p=19` connected compiler uses only field arithmetic after it
loads the horizontal. The exact adapter rebuilds the same sequence over
`K=QQ(a,b)`:

1. Smith saturation with degrees `(0,0,6)`;
2. the seven-dimensional shifted-Popov ambient;
3. the complete connected `D7` ideal;
4. the finite connected `D5` quotient;
5. the rank-five gate and two-dimensional kernel;
6. the moving equation of degrees `(2,9,3)`.

Run it with:

```bash
sage -python elkies-k3/scripts/compile_q80_third_q12_biquadratic_resolved_pencil_qq.py
```

Once that exact pencil is pinned, the remaining ordered gates are:

1. certify generic genus one over `K`;
2. retain explicit birational maps in both directions;
3. minimize the exact Jacobian;
4. factor its discriminant and certify `I6+I4+3I2+8I1`;
5. transport the old components and zero to certify `A5+A3+3A1`;
6. use finite fields only as independent replays, not to infer the
   characteristic-zero Mordell--Weil rank.

The failed first-marking genus-two-cover hypothesis remains closed: its
quadratic field splits at `p=19`, so it cannot supply the mandatory local
quadratic control there.
