# Symbolic discriminant factors do not split the remaining cofactors

The five unfinished panel rows admit a cheap alternative to general integer
factorization: factor the family discriminant over Q[t] first and specialize
its factors. The bounded test is complete and **produces no proper split
of any of the five retained composite cofactors**.

The exact family factorizations are:

| Compact family | Discriminant divisor at finite t |
|---|---|
| R17 074d9 | one irreducible degree-24 factor, multiplicity one |
| R17 11952 | one irreducible degree-24 factor, multiplicity one |
| MW16 a1-fibration-01 | (t+2)² times one irreducible degree-22 factor |

Each statement includes a nonzero rational constant omitted from the table.
The constants, exact factors and finite irreducibility certificates are
retained in the [arithmetic output](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_symbolic_discriminant_v1.json).
Their degrees sum to 24, the full discriminant degree in the compact
Weierstrass models with degrees A≤8 and B≤12. There is no discriminant
zero at infinity in these models.

## Exact experiment

The [protocol](FRESH_SYMBOLIC_DISCRIMINANT_PROTOCOL.json) freezes the same
three families and five incomplete fibres as the preceding
[retained-factor audit](HISTORIC_RETAINED_FACTOR_BOUNDARY.md). It projects
only A(t), B(t), the five parameters and equations, and their unresolved
cofactors. The [input](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_symbolic_discriminant_inputs_v1.json)
contains no sections, exceptional coordinates, rank labels or CT matrices.

For each family the worker computes

\[
 D(t)=-16\bigl(4A(t)^3+27B(t)^2\bigr)
\]

and its factorization over Q[t]. A separate modular certificate proves
irreducibility, rather than merely trusting the factor command. Good
reductions have the following irreducible factor degrees:

| Rational factor | Prime / finite factor degrees |
|---|---|
| 074d9 degree 24 | 131: 8+16; 137: 4+6+14 |
| 11952 degree 24 | 131: 2+3+19; 137: 3+4+6+11; 151: 2+22 |
| MW16 degree 22 | 131: 4+5+13; 137: 8+14 |

Any proper rational factor would have a degree appearing as a proper subset
sum in every good reduction. The displayed intersections are empty.
The linear factor needs no modular proof.

At t=n/d, each factor q(t) is homogenized to d^deg(q) q(n/d), with rational
denominators handled exactly. Its numerator is intersected by gcd with
the already retained composite cofactor N. Exact rational scalings verify
that the specialized compact equations and the frozen panel equations are
the same curves. The independent replay also checks the constant factors,
which are coprime to each N; no omitted constant supplies a proper split.

| Frozen case | Factor receiving the entire unresolved cofactor | Proper splits |
|---|---:|---:|
| new-40, 2818/1535 | 074d9 degree 24 | 0 |
| new-72, 2012/211 | 11952 degree 24 | 0 |
| new-186, 4286/1881 | 11952 degree 24 | 0 |
| new-90, −1867/270 | MW16 degree 22 | 0 |
| MW16 control, −3187/3697 | MW16 degree 22 | 0 |

An irreducible polynomial can have highly composite integer values. This
test does not prove those values hard to factor or prime, nor rule out
other algebraic methods for factoring them. It closes this specific
symbolic-factor route without spending an integer-factorization budget.

## The MW16 linear factor is a fixed singular fibre

At t=−2, the degree-22 factor is nonzero and c4=−48A is nonzero, while
the discriminant has order two. Thus the elliptic surface has multiplicative
type I2 there. The verifier also checks the singular cubic identity

\[
 x^3+A(-2)x+B(-2)=(x-r)^2(x+2r),\qquad
 r=-\frac{3B(-2)}{2A(-2)}.
\]

For the two retained MW16 parameters, n+2d is respectively **−1327** and
**4207**, both nonzero. Their unresolved cofactors are coprime to these
linear values. The high fibre's recorded jump therefore does not occur
at this characteristic-zero degeneration, and the fixed linear factor
does not resolve its missing arithmetic support.

This statement concerns the **elliptic discriminant**, not the discriminant
or resultant of an auxiliary 2-cover, relation curve or simultaneous-lift
carrier. It does not exclude such auxiliary objects becoming reducible,
changing genus or gaining rational points at a smooth elliptic fibre.
In particular, this is not a disproof of the user's proposed simultaneous-
solubility mechanism.

## Implication for the rank-jump work

The panel's local-boundary coverage remains eleven of sixteen. The fresh
MW16 +11 comparison is still incomplete. The new family factor data can
be reused in local arithmetic, but the factorization pattern is shared
by every smooth member of each family and does not distinguish the five
specializations.

No new incidence, solubility or point-visibility feature is established.
The result is a failed computational route to complete incidence data.
The remaining mechanism question still requires independently constructed
additional classes and an explanation of their rationality. Reusing the
known exceptional classes to fill the missing CT matrices would not meet
the masked comparison requirement.

## Verification and scope

All three workers finish within their 30-second caps. The
[portable verification](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_symbolic_discriminant_verification_v1.json)
uses ordinary Python rational arithmetic for polynomial multiplication,
specialization and scaling, and independent finite-polynomial Rabin tests
for every modular irreducible factor. It also verifies all gcd outcomes,
constant-factor coprimality and the MW16 singular-fibre identities. It has
no Sage dependency. An initial verifier string-to-rational conversion error
was corrected before its certificate was produced.

```sh
timeout 30 python3 elliptic-curves/rank-jump/verify_fresh_symbolic_discriminant.py check
```

The scripts read the existing compact atlas coefficients without editing
them. There is no parameter sweep, integer factoring campaign, point search
or change to active-search protocols, outputs or mathematical status.
