# Galois groups of cancellation parameter polynomials through degree thirty

This note classifies the arithmetic Galois group of `M_(m,r)` for every pair
with parameter degree `mr<=30`.  It is an exact finite theorem, not a uniform
Galois-group theorem.  Degrees at most six use resolvents; degrees seven
through thirty use exact modular factorization certificates and Jordan's
theorem.

All polynomials in the table are irreducible by the certificates in
[PARAMETER_IRREDUCIBILITY.md](PARAMETER_IRREDUCIBILITY.md), so every displayed
group acts transitively on the parameter roots.

| `(m,r)` | degree `mr` | Galois group in its root action |
|---|---:|---|
| `(1,1)` | 1 | trivial |
| `(1,2)`, `(2,1)` | 2 | `S_2` |
| `(1,3)`, `(3,1)` | 3 | `S_3` |
| `(1,4)`, `(4,1)` | 4 | `S_4` |
| `(2,2)` | 4 | `D_4` of order 8 |
| `(1,5)`, `(5,1)` | 5 | `S_5` |
| `(1,6)`, `(3,2)` | 6 | `S_6` |
| `(2,3)`, `(6,1)` | 6 | transitive `S_5`, equivalently `PGL_2(F_5)`, of order 120 |
| all `(m,r)` with `7<=mr<=30`, except the next rows | `mr` | `S_(mr)` |
| `(4,3)` | 12 | `A_12` |
| `(2,8)`, `(16,1)` | 16 | `A_16` |
| `(17,1)` | 17 | `A_17` |
| `(1,24)`, `(12,2)` | 24 | `A_24` |

The order and transitivity identify every entry.  In degree four, the
transitive order-eight subgroup of `S_4` is dihedral.  In degree six, a
transitive subgroup of order 120 is the exceptional six-point action of
`S_5`, identified with the action of `PGL_2(F_5)` on the six points of the
projective line.  Thus the two degree-six exceptions are not `S_6` despite
being irreducible and having nonsquare discriminant.

The exact script computes the Galois groups over `Q` using resolvent methods,
then checks degree, order, transitivity, containment in the alternating group,
and the cycle types needed to distinguish the exceptional actions:

```bash
.venv/bin/python scripts/verify_parameter_galois_groups.py
```

## Degrees seven through thirty: Frobenius--Jordan certificates

For each pair in this range, the degree-sieve certificate proves
irreducibility, hence transitivity.  Except for a few prime-degree rows, an
unramified modular factorization of type `(n-1,1)` supplies an `(n-1)`-cycle by Dedekind's
theorem.  A transitive group containing such a cycle is 2-transitive and
therefore primitive.  In the remaining rows, transitivity itself implies
primitivity because the degree is prime.

A second unramified factorization has exactly one cycle of prime length
`ell<=n-3`, while none of its other cycle lengths is divisible by `ell`.
Raising that permutation to the least common multiple of the other lengths
kills every other cycle and leaves an `ell`-cycle.  Jordan's theorem now says
that the primitive Galois group contains `A_n`.  Finally, the discriminant is
a rational square exactly when the group lies in `A_n`.  Thus the group is
`A_n` in the square case and `S_n` otherwise.  The square decision is uniform,
not computational: the closed formula and parity criterion are proved in
[PARAMETER_DISCRIMINANT.md](PARAMETER_DISCRIMINANT.md).

The square-discriminant cases in degrees seven through thirty are

\[
 (m,r)=(4,3),(2,8),(16,1),(17,1),(1,24),(12,2).
\]

For `(4,3)`, the factorization types are `(11,1)` modulo `3` and `(7,5)`
modulo `11`.  For `(2,8)` they are `(15,1)` modulo `41` and `(11,5)` modulo
`2`; for `(16,1)` they are `(15,1)` modulo `37` and `(8,4,3,1)` modulo `5`.
The second type in each case isolates respectively a 5-, 5-, and 3-cycle.
Their discriminants are exact rational squares, so the groups are `A_12`,
`A_16`, and `A_16`.  For `(17,1)`, factorization type `(7,6,4)` modulo `7`
isolates a 7-cycle; prime-degree primitivity and its square discriminant give
`A_17`.  The same certificates give `A_24` for `(1,24)` and `(12,2)`.  All
other pairs in the range have nonsquare discriminant and group `S_n`.  The
exact regression checks every certificate:

```bash
.venv/bin/python scripts/verify_parameter_galois_jordan.py
```

The table shows that parameter irreducibility alone does not force full
symmetric Galois group.  A uniform theorem must account for the dihedral
`(2,2)`, the exceptional degree-six families, and the six alternating cases
before extrapolating from the symmetric cases.
