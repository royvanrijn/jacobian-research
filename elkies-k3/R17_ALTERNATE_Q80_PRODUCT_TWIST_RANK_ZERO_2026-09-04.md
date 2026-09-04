# One alternate-Q80 product twist has geometric rank zero

## Result

For the rational rank-one `V4` base at shortlist rank `55`, with character
pair

```text
alternate-orbit-0fda0 : alternate-orbit-1037d,
```

the product-character twist of the direct `norm12-orbit-11952`
alternate-Q80 fibration has geometric Mordell--Weil rank zero:

```text
rank E^(d)(QQbar(u)) = 0.                              (1)
```

This is Success B for this product construction.  It is an exact Frobenius
and Picard-number obstruction, not a bounded section search.  In particular,
this `V4` construction retains only its two defining character directions and
does not acquire a third direction from their product character.

The compact certificate is
[`elkies-k3-r17-product-alternate-orbit-0fda0--alternate-orbit-1037d-p131-toric-frobenius-v1.json`](../artifacts/generated-results/elkies-k3-r17-product-alternate-orbit-0fda0--alternate-orbit-1037d-p131-toric-frobenius-v1.json).

## 1. Exact target and good reduction

Write the direct alternate-Q80 model as

```text
E : y^2 = x^3 + A(u)x + B(u),
deg(A)=8, deg(B)=12, deg(Delta)=24,
```

and let `d=q_0fda0*q_1037d`.  The two factors are the exact irreducible
quadratics stored in the rank-55 shortlist row.  At `p=131`, the verifier
checks all of the following in `GF(131)[u]`:

```text
deg(d)=4;                 d is squarefree;
deg(Delta)=24;            Delta is squarefree;
gcd(d,Delta)=1.
```

Thus the product twist has four geometric `I0*` fibres at `d=0`, the original
twenty-four `I1` fibres, arithmetic genus `chi=4`, and trivial-lattice rank

```text
2 + 4*rank(D4) = 18.                                  (2)
```

## 2. Regular toric model

The computation uses the regular quadratic-twist presentation

```text
d(u) y^2 = x^3 + A(u)x + B(u).                        (3)
```

Over the generic fibre it is isomorphic to the short Weierstrass twist by

```text
X=d*x,  Y=d^2*y.
```

This form is important: the short equation
`Y^2=X^3+d^2*A*X+d^3*B` puts the four `D4` resolutions on a degenerate Newton
face, while (3) is nondegenerate for the toric computation.

The Newton polytope of (3) has vertices

```text
(0,0,0), (12,0,0), (0,3,0), (0,0,2), (4,0,2)
```

in exponent order `(u,x,y)`, with half-space presentation

```text
u>=0, x>=0, y>=0, u+4x+4y<=12, 2x+3y<=6.             (4)
```

The exact finite-field coefficient compiler emits 28 nonzero monomials.
The open-source Costa--Harvey--Kedlaya
[`ToricControlledReduction`](https://github.com/edgarcosta/ToricControlledReduction)
backend is pinned at commit
`74cda9e8148cd8e9a3928fc15a558c9a70b67cc1`.  Its exact Jacobian reductions
pass the nondegeneracy gate and return primitive Hodge vector

```text
[3,30,3]
```

and a monic weight-two Frobenius polynomial `P_toric(T)` of degree `36`.

## 3. Full cohomology and exact boundary correction

The toric primitive factor contains eight boundary classes in addition to
the elliptic degree-28 factor.  They are computed independently as a finite
permutation motive.  Put

```text
D = {d=0},
Z = {d=0, x^3+A*x+B=0}.
```

Write `P_D` and `P_Z` for the characteristic polynomials of `131*Frob` on
`H^0(D)` and `H^0(Z)`.  Then the boundary factor is

```text
P_boundary(T)
 = P_Z(T)/P_D(T).                                      (5)
```

All four branch points are rational at `131`.  The cubic factor-degree
patterns over their residue fields are

```text
1+1+1,  1+2,  1+1+1,  3.
```

Consequently (5) is the exact degree-eight polynomial

```text
(T-131)^4 (T^2-131^2) (T^2+131*T+131^2)
 = (T-131)^5 (T+131) (T^2+131*T+131^2).               (6)
```

The certificate verifies exact divisibility and defines

```text
P_E(T) = P_toric(T) / P_boundary(T),
deg(P_E)=28.                                           (7)
```

It also reconstructs the factors omitted from primitive toric cohomology.
The central and three nonidentity outer components of every `I0*` fibre give
`H^0(D)` and `H^0(Z)` respectively, while the fibre and zero section give two
fixed Tate classes.  Hence

```text
P_triv(T) = (T-131)^2 P_D(T) P_Z(T),       deg(P_triv)=18,
P_ambient(T) = (T-131)^2 P_D(T)^2,         deg(P_ambient)=10,
P_H2(T) = P_ambient(T) P_toric(T)
        = P_triv(T) P_E(T),                deg(P_H2)=46.       (8)
```

All four identities and the complete integral coefficient lists are stored
in the certificate.  Thus the degree-28 quotient is taken from the full
`H^2` polynomial with the exact `U+4D4` Frobenius factor removed; its degree
is not inferred only from a conductor count.

As an independent cross-check, the first two power sums from (7) are

```text
s_1=-119,  s_2=18305,
```

exactly the two values previously obtained by direct fibrewise character
sums.  Thus (7) extends the old two-moment audit rather than replacing it
with an unrelated normalization.

## 4. Weil and complete Tate-root tests

The polynomial (7) has functional-equation sign `+1`.  Normalize its roots
to weight zero by

```text
Q(Z) = P_E(131*Z)/131^28.                              (9)
```

The verifier checks the functional equation exactly.  It then writes

```text
Q(Z)/Z^14 = R(Z+Z^(-1))
```

and uses certified real-algebraic root isolation to prove that all fourteen
roots of `R` lie in `[-2,2]`.  This independently checks the expected Weil
absolute value.  Finally it computes

```text
gcd(Q(Z), Phi_m(Z)) = 1
```

for every `m` with `phi(m)<=28`.  The finite enumeration through `m=1568` is
complete because the elementary bound `phi(m)>=sqrt(m/2)` implies
`m<=2*28^2` in this range.  There are no cyclotomic hits.  Hence the
degree-28 elliptic factor contains no eigenvalue of the form `131*zeta` with
`zeta` a root of unity.

Every geometric divisor class specializes to such a Tate eigenvalue.  The
special fibre therefore has

```text
rho <= 18.                                            (10)
```

Together with the explicit trivial lattice (2), this gives `rho=18` at the
special fibre.  Specialization of Neron--Severi groups and Shioda--Tate now
give

```text
rank E^(d)(QQbar(u))
 <= rho - 2 - 4*rank(D4)
 <= 18 - 18 = 0,
```

which proves (1).

## 5. Integral Tate-cohomology consequence

For this same product squareclass, put

```text
A=E(QQbar(u,sqrt(d))),       A-=ker(1+sigma).
```

The standard twist isomorphism identifies `A-` with
`E^(d)(QQbar(u))`.  Equation (1) therefore gives `A-=0`.  The already-proved
height argument for the pulled-back `48I1` surface makes `A` torsion-free, and
the exact sequence

```text
0 -> Gamma_d -> A-/2A- -> Hhat^(-1)(<sigma>,A) -> 0
```

now has zero middle term.  Consequently

```text
A-=0,       Gamma_d=0,       Hhat^(-1)(<sigma>,A)=0.  (11)
```

Thus the nonzero Tate-class loophole is absent for the rank-55 target.  The
earlier height-eight inversion remains a useful regression, but is no longer
the closing argument for this product character.

## 6. Replay

The complete compiler, pinned open-source build, approximately four-minute
Frobenius computation, and independent quotient verifier are run by

```bash
elkies-k3/scripts/run_r17_product_toric_frobenius.sh \
  'alternate-orbit-0fda0:alternate-orbit-1037d' 131
```

The terminal certificate line is

```text
R17TORICFROB|pair=alternate-orbit-0fda0:alternate-orbit-1037d|p=131|degree=28|tate_degree=0|rho_upper=18|mw_upper=0|status=PASS_GEOMETRIC_PRODUCT_TWIST_RANK_ZERO
```

No Magma installation or license is used.  Sage supplies the exact model,
finite-field, factorization, boundary, functional-equation, Weil, and
cyclotomic checks; the pinned open-source LGPL toric controlled-reduction code
supplies the p-adic Frobenius polynomial.  The certificate pins the rational
model and character-source hashes, tool commit, exact configure/build and
replay invocations, raw input/output and executable hashes, and the exporter,
raw-parser, verifier, and runner hashes.
