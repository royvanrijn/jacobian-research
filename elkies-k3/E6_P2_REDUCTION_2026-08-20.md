# E6 P2 reduction and corrected GF(31) search state — 2026-08-20

This note records the current transition from the first explicit E6/P1 modular reconstruction to the search for a genuinely independent second Mordell-Weil section.

## First reconstructed E6/P1 point: boundary correction

The first fully reconstructed point over `GF(31)` has

```text
r0 = 4, s0 = 18, x1 = 27
a1 = 4, a2 = 16, a3 = 7, a4 = 6
lambda = 24, mu = 18
s1 = 23, sl = 23, sm = 4
b4 = 6, b5 = 9
```

with

```text
A(t) = 3*t^5 + 6*t^4 + 7*t^3 + 16*t^2 + 4*t + 20
B(t) = t^8 + 5*t^7 + 27*t^6 + 9*t^5 + 6*t^4
       + 9*t^3 + 6*t^2 + 21*t + 8

P1:
X1(t) = 16*t^2 + 27*t + 11
Y1(t) = t^4 + 11*t^3 + 28*t^2 + 19*t + 3.
```

The full section identity and all closed fiber equations verify exactly.
However, a later exact discriminant audit found

```text
ord_0(Delta) = 4,
ord_1(Delta) = 5,
ord_24(Delta) = ord_18(Delta) = 2.
```

Thus the fiber at `t=1` is `I5`, not `I4`.  This point lies on a fiber-
enhancement boundary and is not a point of the required open
`IV* + I4 + I4 + I2 + I2 + 4 I1` stratum.  The reconstruction script has
always printed `delta_mult=5`; the earlier prose failed to enforce the exact-
multiplicity open condition.

The 8x8 Jacobian of the reduced fiber system with respect to

```text
a1,a2,a4,lambda,mu,s1,sl,sm
```

has full rank over `GF(31)`, so this point is smooth on the closed equation
scheme with local parameters `(r0,s0,x1)`.  That Jacobian statement does not
put it in the desired exact-fiber open stratum.

## P3 is dependent on this residue branch

A fixed-surface search found the polynomial section

```text
X3(t) = 28*t^2 + 8*t + 18
Y3(t) = t^4 + 29*t^3 + t^2.
```

Direct group-law computation over `GF(31)(t)` gives

```text
P1 + P3 = -P3,
```

hence

```text
P1 = -2*P3.
```

Thus this P3 is not an independent Mordell-Weil generator on this residue branch. This also explains why first- and second-order deformation tests found P3 to lift automatically with P1.

The search therefore moved to the canonical one-denominator section P2, whose component profile is

```text
P2 = (1,1,2,1,1),   P2.O = 1.
```

## Fixed-surface P2 ansatz

Write the pole of P2 as `r` and `z=t-r`. In numerator coordinates,

```text
P2 = (X2/z^2, Y2/z^3)
```

and the section equation is

```text
Y2^2 = X2^3 + A*X2*z^4 + B*z^6.
```

At each finite reducible fiber `a in {0,1,lambda,mu}`, the component condition requires

```text
X2(a) = s_a * (a-r)^2,
Y2(a) = 0.
```

Let `C(t)` be the unique cubic interpolating these four prescribed X-values and put

```text
F(t) = t*(t-1)*(t-lambda)*(t-mu).
```

The most general degree-5 numerator satisfying the four X-conditions was initially written

```text
X2(t) = C(t) + F(t)*(q0 + q1*t).
```

## Structural reduction: q1 = 0

The coefficient `q1` is impossible on this normalized model.

If `q1 != 0`, then `deg X2 = 5`, hence

```text
deg(X2^3) = 15.
```

But

```text
deg(A*X2*z^4) <= 5+5+4 = 14,
deg(B*z^6)     <= 8+6   = 14.
```

Therefore the right-hand side

```text
H = X2^3 + A*X2*z^4 + B*z^6
```

has degree exactly 15. Since a nonzero polynomial square has even degree, `H` cannot equal `Y2^2`.

Hence necessarily

```text
q1 = 0
```

and the P2 numerator reduces to

```text
X2(t) = C(t) + q0*F(t),
```

so in particular

```text
deg X2 <= 4.
```

This is an exact structural reduction, not merely a search optimization.

For a fixed surface, excluding the four reducible-fiber positions for the pole `r`, the search drops from

```text
27 * 31^2 = 25,947
```

candidate `(r,q0,q1)` triples to only

```text
27 * 31 = 837
```

candidate `(r,q0)` pairs.

The polynomial-square test should be implemented by coefficient recursion rather than repeated factorization.

## First P2 check

The boundary surface above was exhaustively tested with the larger
pre-reduction P2 search and produced no canonical P2 section. Together with
`P1=-2*P3`, this confirms that it is not the desired independent rank-3 seed.

## Exhaustive GF(31) core search

The reduced four-parameter core `(a1,a2,a4,s1)` has only `31^4` possibilities. Exhaustive enumeration tested 893,730 points on the intended `s1 != 0` branch, found 1,853 solutions of the two core equations, 10 squarefree-GCD candidates after saturation, and seven candidates whose common lambda/mu factor splits over `GF(31)`:

```text
(a1,a2,a4,s1) = (2,12,25,30), roots [28,7]
                 (4,16, 6,23), roots [24,18]
                 (14,22,30,21), roots [12,4]
                 (15, 4,17,30), roots [22,2]
                 (17,23, 4,19), roots [19,8]
                 (21,27, 0, 4), roots [20,12,8]
                 (23,11, 8, 5), roots [23,10]
```

These are the complete candidates returned by that reduced common-factor
test, rather than random samples.  They still require reconstruction, exact
Kodaira multiplicities, squarefree residual discriminant, and deeper component
tests; the common-factor test alone is not a promotion certificate.

## Corrected seven-core audit

[`scripts/search_e6_mw3_p2_split_cores.sage`](scripts/search_e6_mw3_p2_split_cores.sage)
reconstructs every candidate that lies in the triangular chart and applies all
of those gates.  The exact outcome is:

| core | reconstruction result | reduced P2 result |
|---:|---|---|
| 1 | exact target fiber stratum; double roots `7,28` | 0 square hits |
| 2 | boundary: `ord_1(Delta)=5` | rejected before P2 |
| 3 | denominator zero in this triangular chart | not tested |
| 4 | only one reconstructed double root | rejected |
| 5 | no reconstructed double root | rejected |
| 6 | only one reconstructed double root | rejected |
| 7 | exact target fiber stratum; double roots `10,23` | 2 signed square hits, both wrong component |

On core 7, the square section has finite component profile

```text
(2,2,1,1)
```

at the two `I4` and two `I2` fibers.  The target is `(1,2,1,1)`.  The exact
test is group-theoretic: `2*P2` meets the identity component at both `I4`
fibers, so both `I4` labels are `2`.  The two signs are therefore rejected.
The diagnostic verifier
[`scripts/verify_e6_mw3_p1p2_gf31.sage`](scripts/verify_e6_mw3_p1p2_gf31.sage)
also proves that this wrong-profile section is independent of `P1`; that fact
does not repair its height-lattice mismatch.

The canonical polynomial `P3` ansatz on core 7 is even smaller: parity reduces
it to 961 quadratic-X cases.  The complete search in
[`scripts/search_e6_mw3_p3_fixed.sage`](scripts/search_e6_mw3_p3_fixed.sage)
has zero hits.

Therefore this coordinate slice currently has two exact target-fiber P1
surfaces and **no canonical P1+P2 seed**.

## Superseding full-height gate

Subsequent all-chart scans showed that moving to more coordinate slices while
retaining these gates was still insufficient.  Every candidate must pass, in
this order:

1. exact `I4,I4,I2,I2` multiplicities and squarefree residual discriminant;
2. the reduced 837-case square test `X2=C+q0*F`;
3. the deeper `I4` component test by specializing `2*P2`;
4. the exact P1/P2 height gate `(P1+P2).O=1`;
5. the 961-case canonical `P3` test;
6. the full height lattice and a smooth Jacobian for the combined section
   equations.

A point passing all six gates would be the modular seed needed to cut the
three-dimensional P1 locus toward the one-dimensional rank-3 family.

Complete rational-chart scans over `GF(5),GF(7),GF(11),GF(13),GF(17)` found
69 sections passing gate 3 and none passing gate 4.  Therefore the next search
must encode the missing intersection condition algebraically rather than
continue the old scan over larger primes.  See
[`E6_MW3_ATTACK.md`](E6_MW3_ATTACK.md) for the exact counts.
