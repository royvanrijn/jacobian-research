# Picard-rank-20 extension of the reconstructed K3 (2026-08-21)

## Status

The reconstructed rational K3 has an exact third independent section.  The
resulting saturated Neron--Severi lattice has rank 20, discriminant `-43`,
and cyclic discriminant group `Z/43`.  Thus the rational surface is a
singular K3 of discriminant 43, and the original
`IV*+I0*+2I3+I2` fibration has Mordell--Weil rank three rather than two.

This corrects the lattice target used by the earlier rank-19 search.  A
rootless fibration on this surface would have MW rank 18, not 17.  The active
positive frame is
[`data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt`](data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt).

## The third section

On the rational model in
[`MW2_FIBRATION_PATH_2026-08-21.md`](MW2_FIBRATION_PATH_2026-08-21.md), the
new section is

```text
S.X = -(625/441)t^2 + (3550/1323)t,

S.Y = (390625/18522)t^4 -(359375/9261)t^3
      +(1875/98)t^2.
```

Literal substitution verifies `S.Y^2=S.X^3+A*S.X+B`.  Resolution at the
reducible fibers gives, in the ordering
`(E6,D4,A2_P2,A2_P1,A1; S.O)`,

```text
S = (2,d2,0,0,0;0),
```

where `d2` is either of the two D4 triality classes different from the class
`d1` met by `P2`.  More concretely:

- at `IV*` infinity, `S` has the same exceptional leading sign as `P2` and
  the opposite sign from `P1`;
- at `I0*`, `S` and `P2` reach distinct first exceptional centers;
- `S` avoids the nodes of both `I3` fibers and the `I2` fiber.

The exact section intersections are

```text
P1.S = 0,   P2.S = 2,   S.O = 0.
```

The raw `P2`/`S` coordinate gcd also has a factor `t`, but that contact is at
the unresolved `I0*` singular point.  The two blowup paths separate there,
so only the smooth roots `24/25` and `32/25` contribute on the resolved K3.

## Height lattice and discriminant

In the basis `(P1,P2,S)`, Shioda's formula gives

```text
height Gram = (1/6) * [ 9  2  8]
                        [ 2 18 -5]
                        [ 8 -5 10],

det(height Gram) = 43/216.
```

Since the reducible-fiber root determinant is `216`, the full NS determinant
is `-43`.  The explicit twenty-class Gram has signature `(1,19)` and Smith
diagonal

```text
1,1,...,1,43.
```

It is automatically saturated: a proper finite-index even overlattice would
have index whose square divides the prime `43`, which is impossible.  Rank
20 is maximal for a complex K3, so the Picard rank is exactly 20.

After splitting off the fiber and `F+O` as a hyperbolic plane, the positive
rank-18 frame has determinant `43`.  Its current root invariants are

```text
root rank/count/determinant = 15/110/216,
ADE = E6 + D4 + 2A2 + A1,
MW rank = 18-15 = 3.
```

The saturated height lattice reduces to

```text
(1/6) * [ 2 -1  0]
        [-1  3  1]
        [ 0  1  9].
```

## Optimal-rank target

MW rank zero is impossible on this NS.  If it existed, its rank-18 ADE root
lattice `R` would satisfy

```text
det(R) = 43 * |MW_tors|^2,
```

and hence would have determinant divisible by 43.  But every irreducible ADE
factor of rank at most 18 has determinant supported on primes at most 19.
Therefore the best possible endpoint is MW rank one.  This argument is
checked by [`scripts/verify_mw0_obstruction.py`](scripts/verify_mw0_obstruction.py).

That optimum is now attained by a two-neighbor exact lattice path, ending at
the `A12+A3+A2` frame with free height `43/156`, trivial torsion, `P.O=0`,
and a semistable candidate presentation.  See
[`PICARD20_MW1_OPTIMAL_PATH_2026-08-21.md`](PICARD20_MW1_OPTIMAL_PATH_2026-08-21.md).

## Reproduction

Run

```bash
sage elkies-k3/scripts/verify_picard20_ns_extension.sage

sage elkies-k3/scripts/analyze_mw3_branch.sage \
  --frame elkies-k3/data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt \
  --name picard20-e6-d4-a2a2-a1-mw3

python3 elkies-k3/scripts/verify_mw0_obstruction.py
```

The first command checks the section identities, resolved intersections,
component profile, height Gram, full discriminant and Smith group, integral
`U` splitting, and pinned rank-18 frame.  The current neighbor search is a
bounded discovery computation and is recorded separately from these exact
certificates.

The q=8 neighbor using the extra section `S` has since been executed exactly.
Its compact `I3+I4+I5+I8+4I1` equation, saturated MW-rank-two basis, complete
degree-one finite-field discovery scan, and CM-43 Kumar identification are in
[`PICARD20_Q8_CHORD_2026-08-21.md`](PICARD20_Q8_CHORD_2026-08-21.md).  That
neighbor is Noether--Lefschetz-special and is not claimed to survive in the
generic determinant-948 lattice.
