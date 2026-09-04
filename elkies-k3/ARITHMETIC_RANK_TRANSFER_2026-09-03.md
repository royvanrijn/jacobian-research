# Galois-equivariant Shioda--Tate balance and arithmetic marking gate (2026-09-03)

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: OP-K3-DIFFERENT-NS-ARITHMETIC-MW17 71f43dc9ef3af620 -->

<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->

## Outcome

The classical Shioda--Tate balance now has a Galois-equivariant quotient
formulation in the repository's marking notation.  Its canonical statement
and proof are Theorem A2 and Corollaries
A2.1--A2.2 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).

For a Jacobian fibration over a characteristic-zero field `K`, put

```text
V = NS(X_Kbar) tensor QQ,
U = <F,O+F>,
W = U^perp,
R = geometric fibre-root space,
M = MW(pi_Kbar) tensor QQ.
```

If the fibre and zero section descend to `K`, then the finite Galois image
`Gamma` fixes `U` pointwise and

```text
M = W/R,
rank MW(pi/K(C)) = dim (W/R)^Gamma
                   = dim V^Gamma - 2 - dim R^Gamma.
```

For two `K`-defined fibrations on the same surface, the stronger
representation-ring identity is

```text
[M_2]-[M_1]=[R_1]-[R_2].
```

Thus a fibration hop redistributes Galois representations, not only
dimensions.
Removing an anti-invariant root creates an anti-invariant Mordell--Weil
direction and does not raise arithmetic rank.

Shioda's papers on Mordell--Weil lattices and Galois representations already
establish the Galois action and height compatibility.  The representation-
ring equality above is a tailored quotient corollary, not a new general
arithmetic theorem.  The project contribution is the fail-closed integral
marking schema, verifier, and the exact control computations.

## Rational-source inheritance

The useful construction consequence is source-first.  If one marked source
supplies `rho(X_Kbar)` independent divisor classes over `K`, then Galois acts
trivially on `NS(X_Kbar) tensor QQ`.  Every further marked `U` represented by
`K`-divisors and passing the primitive-nef-effective-zero gates defines a
Jacobian fibration over `K`, with

```text
rank MW(K(C)) = rho(X_Kbar)-2-rank(R).
```

For `rho=19`, a rootless target therefore has arithmetic rank seventeen.
This rank theorem does not require explicit formulas for all seventeen
endpoint sections.  A saturated integral source basis is still required if
the full integral endpoint Mordell--Weil lattice, rather than only its rank,
is claimed.

## Machine-readable arithmetic marking

The schema is
[`data/arithmetic/arithmetic-marking-v1.schema.json`](data/arithmetic/arithmetic-marking-v1.schema.json).
It records:

1. one common geometric Neron--Severi Gram matrix;
2. integral generators of the finite Galois image;
3. the marked hyperbolic plane for every fibration;
4. the embedded geometric fibre-root basis;
5. expected geometric and arithmetic Mordell--Weil ranks;
6. optional component labels, section classes, fields, and orbit sizes;
7. optional rank-transfer edges.

The exact verifier is
[`scripts/certify_arithmetic_rank_transfer.sage`](scripts/certify_arithmetic_rank_transfer.sage).
It checks that every action is integral, unimodular, and Gram-preserving;
closes the finite group; verifies that each `U` is fixed pointwise and each
root space is stable; computes the rational fixed NS, root, and Mordell--Weil
spaces; and checks the representation-ring identity by its trace on every
group element.  Its output includes the induced matrices on the root and
Mordell--Weil quotient bases, a rational basis of the fixed Mordell--Weil
space, and orbit--stabilizer data for every supplied section class.  Missing
descent evidence is typed `UNKNOWN` and never treated as rank zero or as
arithmetic promotion.

Reproduce the pinned controls with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_arithmetic_rank_transfer.sage --check
```

An additional self-contained marking can be checked with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_arithmetic_rank_transfer.sage \
  --marking path/to/arithmetic-marking.json
```

The generated certificate is
[`../artifacts/generated-results/elkies-k3-arithmetic-rank-transfer-controls-v1.json`](../artifacts/generated-results/elkies-k3-arithmetic-rank-transfer-controls-v1.json).

## Exact controls

| control | geometric rank | fixed rank | arithmetic conclusion |
| --- | ---: | ---: | --- |
| H3 rootless R17 over `QQ` | 17 | 17 | arithmetic rank 17 |
| degree-two alternate Q80, `norm12-orbit-11952` over `QQ` | 17 | 17 | arithmetic rank 17; direct equation and saturated sections now compiled separately |
| unordered E6 incidence over `QQ(k)` | 4 | 2 | two exchanged pairs leave two invariant sums |
| E6A1 orbit 103 over `QQ(k)(r)` | 3 | 2 | `2*trivial + chi_-3` |

The H3 control cross-checks the exact Picard-rank-19 endpoint against the
separate published-target certificate and its seventeen exact `QQ(t)`
sections; their unimodular identification with pinned `R17` supplies a
rational rank-19 divisor span together with fibre and zero.  The E6 incidence
control uses the full integral
rank-19 Neron--Severi Gram: conjugation swaps `P,Q` and `R1,R2`, fixes the
root space, and leaves fixed NS rank seventeen.  Formula (A2.2) then gives
`17-2-13=2`.  Orbit 103 independently checks the one-dimensional
anti-invariant field `QQ(sqrt(-3))`.

## Arithmetic promotion of the alternate-Q80 pencil

The gate now consumes the cheapest exact degree-two copy of the alternate-Q80
frame on the published R17 surface.  In the rational divisor basis

```text
(F,O,Q1,...,Q17),
```

where the `Qi` are the determinant-one pinned combinations of the published
sections, its fibre and mate are

```text
D   = (40,-1,-1,-1,3,0,0,1,-1,-1,3,1,1,0,2,-1,1,-2,-2),
O+D = (40, 0,-1,-1,3,0,0,1,-1,-1,3,1,1,0,2,-1,1,-2,-2).
```

The exact bisection certificate supplies an irreducible smooth genus-one
`QQ`-curve in class `D`; `D.O=1`, so the rational old zero is also the new
zero.  The marked plane is primitive, its complement is the rootless
alternate-Q80 frame, and the rational source basis makes the Galois action on
the full rank-19 NS space trivial.  The gate therefore certifies

```text
rank E_alt(QQ(t)) = 17
```

without compiling the new pencil equation or recovering its seventeen
sections.  The canonical application proof is
[`R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md`](R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md).

The verifier also contains an explicitly labelled representation-only `C2`
regression.  Removing an anti-invariant `A1` raises geometric rank by one and
arithmetic rank by zero; removing a fixed `A1` then raises both by one.  This
tests (A2.3) on every group element and is not asserted to be a new geometric
surface.

## Application to the completed-core NS0024 route

The arithmetic promotion gate has been applied to the completed-core path

```text
D5+E8/MW4 -> 3A1+A2/MW12 -> 3A1+A2/MW12 -> rootless/MW17.
```

Its geometric ranks `4,12,12,17` remain exact.  Its arithmetic ranks are all
currently unknown.  The gate returns

```text
FAIL_CLOSED_GEOMETRIC_ONLY.
```

This is forced by the existing proof boundary: the primes `17,13,7` describe
Kneser neighbours of completed positive frames, not elliptic-neighbour
degrees or a marked sequence of hyperbolic planes on one field-defined K3.
There is currently no characteristic-zero equation for the `D5+E8/MW4`
source, no rational rank-19 divisor span, no Galois component action, and no
field-defined target `U` certificate.

The shortest arithmetic promotion route is therefore:

1. construct one characteristic-zero NS0024 source surface over `K`;
2. exhibit nineteen independent `K`-divisor classes, preferably the source
   `U`, thirteen split `D5+E8` roots, and four rational sections;
3. construct an explicit integral target rootless `U` in that rational
   divisor basis and certify its primitive, nef, and effective-zero gates;
4. apply Corollary A2.1 to obtain arithmetic rank seventeen;
5. compile an endpoint equation and compute the integral section lattice only
   if an explicit family or full-lattice theorem is required.

The existing modular NS0024 searches remain feasibility and Frobenius-ranking
evidence.  They are not a substitute for steps 1--3.

This completed-core path is distinct from the certified thirteen-edge
degree-two lattice corridor beginning at `A3+A4+A6/MW4`. Neither route can
carry the full rational rank-19 marking required over `QQ`, by the later
[`NS0024 rational-marking obstruction`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md).
The live milestone has moved to the
[`different-NS arithmetic MW17 objective`](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md).

<!-- status-consumer: EC-K3-ARITHMETIC-RANK-TRANSFER 3031dd2365a29cd5 -->
