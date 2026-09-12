# Toric Cox--Keller maps with permutation reductions

The ordered-factor Cox maps do not give permutation reductions because their
symmetric fibers have either many rational lifts or none.  This note gives a
toric control family showing that the permutation objective itself is
compatible with constant-Jacobian Cox suspensions.

The construction is elementary, but it isolates the missing ingredient:
exceptional arithmetic monodromy rather than another determinant identity.

## 1. Separated toric construction

Work over a characteristic-zero field `k`.  Fix positive integers
`n_1,...,n_r`, all at least two, and put

\[
 X_r=(\mathbb G_m)^r\times\mathbb A^r
\]

with coordinates `(x_1,...,x_r,z_1,...,z_r)`.  Define

\[
 \boxed{
 F_{\mathbf n}(x_\bullet,z_\bullet)
 =
 \left(
 x_1^{n_1},\ldots,x_r^{n_r},
 \frac{z_1}{n_1x_1^{n_1-1}},\ldots,
 \frac{z_r}{n_rx_r^{n_r-1}}
 \right).
 }                                                    \tag{1}
\]

All displayed quotients are regular because the `x_i` are units.

### Theorem 1.1

The morphism (1) is finite etale of geometric degree

\[
 \prod_{i=1}^r n_i
\]

and has Jacobian determinant one.

Indeed, its derivative is block triangular.  The horizontal diagonal block
has determinant

\[
 \prod_i n_ix_i^{n_i-1},
\]

while the vertical diagonal block has the reciprocal determinant.  Once the
`x_i` are chosen in a fiber, every `z_i` is reconstructed uniquely, so the
geometric degree is the degree of the power map on the torus.

The unit rank is

\[
 \operatorname{rank}\mathcal O(X_r)^*/k^*=r.
\]

Thus (1) is a separated-character Cox ledger with one primitive coordinate
per independent toric boundary direction.

## 2. Compressed construction

As in the linear-factor ledger, determinant cancellation alone needs only
one new coordinate.  Put

\[
 J_{\mathbf n}(x)=\prod_i n_ix_i^{n_i-1}
\]

and define

\[
 F_{\mathbf n}^{\mathrm{compressed}}:
 (\mathbb G_m)^r\times\mathbb A^1
 \longrightarrow
 (\mathbb G_m)^r\times\mathbb A^1,
\qquad
 (x_\bullet,z)\longmapsto
 \left(x_1^{n_1},\ldots,x_r^{n_r},
 \frac z{J_{\mathbf n}(x)}\right).                   \tag{2}
\]

It also has determinant one and degree `prod_i n_i`.  For `r>1`, (2) is a
second explicit disproof of the unrestricted unit-rank suspension bound.

## 3. Exact finite-field criterion

Let `F_q` have characteristic prime to every `n_i`.  The vertical
coordinates in (1) and (2) are nonzero scalar maps once `x` is fixed.
Therefore either map is a permutation of its `F_q`-points exactly when every
power map on `F_q^*` is a permutation:

\[
 \boxed{
 F_{\mathbf n}\text{ permutes }X_r(\mathbb F_q)
 \quad\Longleftrightarrow\quad
 \gcd(n_i,q-1)=1\quad\text{for every }i.
 }                                                    \tag{3}
\]

This produces nontrivial characteristic-zero covers whose reductions are
permutations.

If all `n_i` are odd and `N` is the product of their distinct prime
divisors, every prime

\[
 q\equiv2\pmod N
\]

satisfies (3).  Since `gcd(2,N)=1`, Dirichlet's theorem gives infinitely many
such prime fields.

For example, with `(n_1,n_2)=(3,5)`, the four-dimensional separated map

\[
 (x_1,x_2,z_1,z_2)\longmapsto
 \left(
 x_1^3,x_2^5,
 \frac{z_1}{3x_1^2},
 \frac{z_2}{5x_2^4}
 \right)
\]

has degree fifteen, determinant one, and is a permutation over `F_17`,
because `gcd(3,16)=gcd(5,16)=1`.

## 4. Relation to the factorization programme

This family clears one logical gap:

> constant-Jacobian Cox maps of geometric degree greater than one can have
> permutation reductions at infinitely many primes.

It does not solve the stronger factorization problem.  The toric power cover
has abelian monodromy and no analogue of the three ordered collision
divisors.  The oriented cubic chart has the desired dicritical geometry but
cyclic degree-three fiber behavior which is not exceptional in the required
point-permutation sense.

A construction combining both features needs a non-symmetric exceptional
cover whose boundary characters admit an oriented Cox ledger.
The repository's Davenport--Sunada pair supplies the exceptional group
theory, but its
[three-column boundary obstruction](DAVENPORT_COX_BOUNDARY_OBSTRUCTION.md)
shows why the present Cox suspension cannot yet be filled to affine space.

The
[Danielewski multi-dicritical family](DANIELEWSKI_MULTI_DICRITICAL_FAMILY.md)
now combines dicritical geometry and bijective reductions on a smooth Cox
target by placing the missing divisors in nonsplit Frobenius orbits.  Its
remaining gap is geometric rather than arithmetic: the target is not
affine space.

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/verify_toric_cox_permutations.py
```

The checker verifies the symbolic determinant in arbitrary displayed rank,
the compressed identity, exact point permutations for several fields, and
failure when the gcd criterion is violated.
