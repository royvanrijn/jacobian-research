# Two Mestre directions on one cover, with a genus obstruction

There are now two explicitly constructed independent new sections on **one
quadratic cover** of the published R17 parent. Their height matrix is
`diag(24,24)`, so they give a subgroup of rank at least 19. The covering base
has **genus 21**. It therefore fails the infinite-rational-base requirement
of the [active construction goal](CORRELATED_QUADRATIC_GAINS_2026-09-12.md).

Changing the auxiliary value to an arbitrary rational function does not
repair this construction: the same Mestre identity has covering genus at
least **9** on this parent. This is a proved obstruction for this identity,
not for arbitrary correlated quadratic covers.

## Construction from the generic equation

Let the published equation be `E: y^2=x^3+A(t)x+B(t)`. Its exact degree-8
and degree-12 coefficients are retained in the
[generic-only packet](../artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1/input.json).
No exceptional point, specialization parameter, or specialized rank is used.

For `v=u^2`, write `k=v^2+v+1`. Mestre's identity gives

```
x1=-(B/A)*k/(v+1),       x2=x1/v,
f(x1)=u^6*f(x2),         f(x)=x^3+A*x+B,
D=-A*B*(v+1)*(B^2*k^3+A^3*v^2*(v+1)^2).
```

On `s^2=D` the corresponding ordinates are

```
y1=s/(A^2*(v+1)^2),      y2=y1/u^3.
```

The identity follows by subtracting `v^3*f(x2)` from `f(x1)`; the cubic
terms cancel and the remaining equation is linear in `x1`. The underlying
construction is established, not a new twist method. See
[Rubin–Silverberg, Theorem 3.7](https://www.maths.tcd.ie/EMIS/journals/EM/expmath/volumes/10/10.4/Rubin.pdf),
which attributes it to Mestre. Its theorem concerns a constant elliptic
curve with a moving auxiliary variable. Independence after substituting
the present moving coefficients requires its own proof below.

Select `u=2`, the first positive integer after the degenerate value `1`.
This gives the explicit cover and sections

```
D=-5*A*B*(9261*B^2+400*A^3),
P=(-21*B/(5*A),   s/(25*A^2)),
Q=(-21*B/(20*A),  s/(200*A^2)).
```

Coefficientwise rational identities verify that both points lie on `E`.
The branch polynomial has degree 44 and is squarefree. It is coprime to
the original degree-24 discriminant `4*A^3+27*B^2`; that discriminant is also
squarefree. All these properties have a simultaneous finite-polynomial
witness at prime 131, with no degree drop. Thus the cover is smooth of
genus `(44-2)/2=21`, has no branch point at infinity, and branches only
over smooth original fibres.

## Exact independence modulo the inherited subgroup

Use the Shioda height convention from the
[rank and lift theorem layer](RANK_MUTATION_AND_LIFT_THEOREMS.md).
The quadratic pullback has `chi=4`. Since the branch points avoid the
original singular fibres, its fibres remain irreducible.

At every simple root of `A`, the cover is ramified. The coordinates of both
`P` and `Q` have poles of orders 2 and 3 in a uniformizer of the covering
base, giving one intersection with zero at each of the eight roots. They
have no other finite poles. At infinity their original `x`-coordinates
grow at most as `t^4`, so the standard K3 chart shows no zero intersection
there. Consequently

```
P.O=Q.O=8,
height(P)=height(Q)=2*chi+2*(P.O)=24.
```

The ordinary group law gives a second exact identity:

```
x(P+Q)=-(7/15)*B/A-(20/81)*A^2/B.
```

Its denominator has one simple factor at every root of `A` and `B` and
no cancellation. All twenty points are ramification points of the cover,
so `(P+Q).O=20`. There is again no pole beyond the allowed degree at
infinity. Thus `height(P+Q)=48` and

```
<P,Q>=(48-24-24)/2=0.
```

The resulting positive-definite matrix proves independence. Moreover the
quadratic involution sends each point to its negative. Each is orthogonal
to every section over `Q(t)`, so their independence is also independence
modulo the full inherited subgroup. The parent's certified 17 independent
directions therefore give the stated rank-at-least-19 subgroup. No assertion
of exact total rank is made.

## Why every rational auxiliary function still has large genus

The polynomial `A` is irreducible over Q, as witnessed by its irreducible
degree-8 reduction at 29. It also has the simple root `51` modulo 59.
Hensel's lemma embeds its root field in `Q_59`. Since `59` is congruent to
`11 mod 12`, this field contains neither a square root of `-1` nor a primitive
cube root of unity. Neither can the root field of `A`.

Likewise `B` has an irreducible degree-12 reduction at 59 and the simple
root `93` modulo 107. Its root field embeds in `Q_107`, which contains no
square root of `-1` because `107` is `3 mod 4`.

Let `u` be any nonzero rational function in `Q(t)`. At either irreducible
divisor write `n=ord(u)`. The displayed formula for `D` gives the exact
valuation table

| Divisor | `n<0` | `n=0` | `n>0` |
|---|---:|---:|---:|
| `A=0` | `1+14*n` | `1` | `1` |
| `B=0` | `3+14*n` | `1` | `3` |

Here is the cancellation check. At `A=0` and `n=0`, the residue of `u^2+1`
cannot vanish, and neither can that of `u^4+u^2+1`: the latter would supply
a primitive cube root of unity through `u^2`. Thus the term `B^2*k^3` is
a unit. When `n<0`, its valuation `12*n` is strictly smaller than the other
term's valuation `3+8*n`; when `n>0`, it is again a unit. This gives the
first row including the prefactor.

At `B=0` and `n=0`, the residue of `u^2+1` cannot vanish, so the term
`A^3*v^2*(v+1)^2` is a unit. For `n>0`, the two bracket terms have valuations
`2` and `4*n`, with the former smaller. For `n<0`, they have valuations
`2+12*n` and `8*n`, again with the former strictly smaller. This gives the
second row. The degenerate constant `u=0` has `D=-A*B^3` and the same odd
valuation conclusion, although the second point formula is then undefined.

Every entry in the table is odd. Removing rational square factors therefore
leaves all eight geometric roots of `A` and all twelve roots of `B` in the
branch divisor. Any additional branch points only increase the genus. Hence

```
degree(branch divisor)>=20,
genus(Q(t,sqrt(D)))>=(20-2)/2=9.
```

This covers arbitrary degrees, heights, zeros and poles of `u(t)`. It is
not an empty bounded search. Enlarging the auxiliary-function box or
rewriting the original parameter by a birational change cannot evade it.

By the finiteness theorem for rational points on curves of genus greater
than one, these bases cannot have infinitely many rational points. See
[Faltings, Inventiones 73 (1983)](https://link.springer.com/article/10.1007/BF01388432).
The required positive endpoint therefore remains unachieved even though
the two new directions themselves are explicit and independent.

## Certificates and scope

The [constructor result](../artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1/result.json)
retains the literal points and branch coefficients. The
[branch-gate certificate](../artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1/branch-gate.json)
records the four field primes and the simultaneous smoothness prime.
The [independent replay](../artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1/independent-replay.json)
uses only Python fractions and integer arithmetic. It verifies the point
and addition identities, pole bounds, squarefreeness and coprimality,
Rabin irreducibility tests, and simple local roots. It derives the heights
from the verified pole divisors instead of trusting the producer's matrix.
The written local-height and valuation arguments remain part of the proof;
this is not formal verification.

```
python3 research/elkies-k3/scripts/verify_r17_mestre_shared_twist.py
python3 -m unittest discover -s research/tests -p 'test_r17_mestre_shared_twist.py' -q
```

Discovery used the retained generic packet and no parameter search. The
constructor had a 30 CPU-second cap; the finite-polynomial gate had a
20-second cap. Both used 4 GiB address-space limits and one BLAS thread.
The [execution receipt](../artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1/execution.json)
records component timings and the six failure/control tests.

The exclusion applies only to this Mestre identity on the published R17
parent. Other carrier identities, higher-pole constructions and other
arithmetic parents remain open. No rational specialization rank improvement
or infinite family of rational specializations has been proved here.
