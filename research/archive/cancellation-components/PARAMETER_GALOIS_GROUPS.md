# Galois groups through degree thirty and an odd square-family series

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

## Higher odd square-family certificates through degree 1057

The degree-thirty script contains exact certificates for the next two members
of the odd square-discriminant family after `(17,1)`:

\[
 \operatorname{Gal}(M_{49,1})=A_{49},\qquad
 \operatorname{Gal}(M_{97,1})=A_{97}.
\]

For degree 49, the modular types

\[
 (22,13,9,3,2),\ (22,10,8,6,3),\ (24,15,10),\ (34,12,3)
\]

at primes `7,11,13,19` have subset-sum intersection `{0,49}`, proving
irreducibility.  The first type isolates a 13-cycle.  A transitive
imprimitive degree-49 group would embed in `S_7 wr S_7`, whose order is not
divisible by 13, so the group is primitive.  Jordan's theorem gives
containment of `A_49`, and the square discriminant gives equality.

For degree 97, the types `(89,6,2)`, `(91,5,1)`, and `(83,12,2)` at primes
`5,19,23` prove irreducibility.  Prime degree gives primitivity, the first
type isolates an 89-cycle, and Jordan's theorem plus the square discriminant
gives `A_97`.

A separate slow exact regression continues with every odd parameter
`m=2a^2-1` for `a=9,11,...,23`:

| `m` | auxiliary primes in the degree sieve | isolated prime cycle | primitivity |
|---:|---|---:|---|
| 161 | `11,13` | 131 | the 131-cycle excludes block sizes 7 and 23 |
| 241 | `5,13` | 149 | prime degree |
| 337 | `5,17` | 109 | prime degree |
| 449 | `19,47` | 127 | prime degree |
| 577 | `5,7,19` | 53 | prime degree |
| 721 | `5,7,11` | 47 | the type modulo 5 excludes block sizes 7 and 103 |
| 881 | `5,17` | 139 | prime degree |
| 1057 | `7,17` | 41 | the type modulo 7 excludes block sizes 7 and 151 |

In every row the subset-sum intersection of the displayed unramified
factorization types is `{0,m}`, proving irreducibility.  The isolated cycle
and the stated primitivity certificate put `A_m` in the Galois group by
Jordan's theorem, while the uniform square-discriminant formula puts the
group inside `A_m`.  Consequently

\[
 \boxed{\operatorname{Gal}(\mathcal M_{2a^2-1,1})
 =A_{2a^2-1}\quad(3\le a\le23,\ a\text{ odd}).}
\]

The complete factor-degree tuples and exact wreath-product compatibility
checks are replayed by

```bash
.venv/bin/python scripts/verify_odd_square_galois_series.py
```

This verifier is intentionally separate from the standard degree-thirty
regression because exact modular factorization in degrees 721 through 1057
is substantially slower.

The table shows that parameter irreducibility alone does not force full
symmetric Galois group.  A uniform theorem must account for the dihedral
`(2,2)`, the exceptional degree-six families, and the alternating square
locus before extrapolating from the symmetric cases.

## Fixed-row growing ramification

There is now a uniform asymptotic statement beyond the finite table.  For
each fixed `r`, all but

\[
 O_r\!\left(\frac{X\log\log X}{\log X}\right)
\]

parameters `m<=X` give an irreducible degree-`mr` polynomial whose Galois
group has order divisible by a prime `p>(log X)^10`.  This follows by
retaining the root-of-unity Newton clusters in the proof of the
Borisov--Filaseta--Lam--Trifonov fixed-derivative theorem.  The exact
translation is in
[the fixed-row ramification note](../../cancellation/FIXED_R_NEWTON_RAMIFICATION.md).

This result deliberately stops short of calling the resulting order-`p`
element a Jordan cycle.  Several clusters can contribute simultaneous
`p`-cycles.  Isolating one cluster and proving primitivity are the two missing
uniform steps toward density-one `A_(mr)` or `S_(mr)`.
