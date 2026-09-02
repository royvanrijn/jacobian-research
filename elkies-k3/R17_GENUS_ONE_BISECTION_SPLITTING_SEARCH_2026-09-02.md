# High-throughput R17 genus-one-bisection splitting search (2026-09-02)

<!-- status-consumer: EC-K3-R17-GENUS1-HIGH-THROUGHPUT-SPLITTING cad3d98ce58c89e7 -->

## Result

The new search removes the rank-28 target fitting from the discovery loop.  It
uses 100 structurally diverse traces selected inside the 1,000
equation-cheapest norm-eight classes and all 43 norm-twelve deep classes.  All
143 compiled branch polynomials are irreducible squarefree quartics over
`QQ`, coprime to their trace denominators and to the degree-24 surface
discriminant.

The production run evaluated 5,474,328 primitive rational parameters: the
complete box

```text
|a| <= 2000,  1 <= b <= 2000,  gcd(a,b)=1
```

together with one million deterministic large-coordinate proposals.  A
separate population contains 900 parameters obtained from multiples on 100
pointed quartic Jacobians.  Twelve Legendre-symbol tables in three disjoint
prime blocks left 77,704 distinct extreme modular collisions.  Every one was
tested exactly, requiring 160,395 integer-square tests.

There is one exact simultaneous split:

```text
t = 1/25
norm-eight class  0x0f6b1
norm-twelve class 0x103b2
```

The norm-eight cover is the neutral member seeded by the sixteenth generic
R17 section at this parameter.  The norm-twelve point is

```text
x = 36075981547811164726251 / 244140625
y = -5898338731136062956741857359589376 / 3814697265625.
```

Exact finite reduction gives baseline rank 17 and combined rank 18 both in a
product of `E(F_p)/2E(F_p)` quotients and independently in a product of
`E(F_p)/3E(F_p)` quotients.  Thus the norm-twelve point is outside the
specialized generic MW17 subgroup.  This proves one quotient direction at
this fibre, hence a rank-at-least-18 specialization; it does not give a new
record rank.

## Trace geometry

For a norm-eight trace

```text
tau = (Nx/h^2, Ny/h^3),  deg(h)=2,
```

regular chord slopes form the pencil

```text
M = M0 + lambda*h^2.
```

The search initializes one member without any exceptional point.  It chooses
a distinct small rational fibre and a generic R17 section on that fibre,
solves the one linear incidence equation for `lambda`, and rejects any member
that is reducible, singular, or badly branched.  The resulting rational point
only initializes the pointed binary quartic.  The actual auxiliary parameter
population is then obtained from group-law multiples on its Jacobian using
the exact inverse pointed-quartic map.  Seed parameters are forced distinct,
so repeated generic initialization cannot create an artificial simultaneous
collision.

For every norm-twelve deep trace the exact closest representative has

```text
deg(h)=4,  deg(M0)<8,
q=(M0^4-6*M0^2*Nx-8*M0*Ny-3*Nx^2-4*A*h^4)/h^6,
deg(q)=4.
```

There is no free `lambda` at this degree: `M0` is the unique regular member.
Complete exact construction verifies all 43 quartics.  The bounded
small-height point pass did not point these curves in advance; the hit above
was found by the common parameter sieve.

## Sieve architecture

The old mixed-trace scanner stored exactly 128 cover bits.  The replacement
uses `ceil(number_of_covers/64)` words throughout, so increasing
`--norm8-count` from 100 through 1,000 requires no source change.  Each prime
table stores, for every `(a:b)` modulo `p`, the simultaneous quadratic-residue
mask for every branch quartic.  A prime at which a clearing denominator
vanishes gives that curve no filtering information; prime selection minimizes
this blindness while balancing the `p^2` table cost.

For each parameter the scanner maintains separate intersections in at least
three disjoint prime blocks.  Its ranking key is, lexicographically,

```text
minimum blockwise F2-rank of surviving trace masks,
minimum blockwise surviving-cover count,
sum of blockwise ranks,
sum of blockwise counts,
all-block survivor count.
```

This distinguishes many locally correlated covers from covers spanning many
classes in `R17/2R17`.  It is still only a local diversity heuristic.  For an
exact hit, the lifted formulas construct every point on the specialized
Weierstrass equation, and the common finite-quotient implementation computes
the image of those points modulo the specialized generic MW17 subgroup.

## Literature boundary

Garbagnati--Salgado explain why special multisections are geometrically tied
to rank jumps, but do not supply this simultaneous arithmetic sieve or its
diversity score:

- A. Garbagnati and C. Salgado,
  [*Rank jumps and Multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159).

The binary-quartic/Jacobian step is the standard degree-two genus-one covering
construction.  Magma's genus-one-model handbook gives an independent
software-level description of the same invariant/covariant map:

- [Magma handbook: genus-one models as coverings](https://magma.maths.usyd.edu.au/magma/handbook/text/1594).

These references justify the geometric and covering framework.  They do not
turn a bounded Legendre scan into a completeness theorem.

## Reproduction

Run the stored production profile with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_elkies_2026_genus_one_bisection_splitting.sage \
  --norm8-count 100 \
  --equation-pool-size 1000 \
  --output artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json \
  --local-directory artifacts/local/elkies-k3/r17-genus-one-bisection-splitting/production-v1
```

The generated artifact has SHA-256
`4745f53993675286173298c9444da023825fd4429e74ede582a1d8c14979d07e`.
The C++ scanner is compiled by the Sage driver.  Passing `--exact-limit 0`
(the default) exact-tests every modular extreme; a positive limit is an
explicit truncated development run and is recorded as such.

## Boundary

The compact box is exhaustive, but the million large-coordinate proposals,
the trace selection, and the Jacobian multiple ranges are bounded.  The
Legendre ranks are search rankings, not Selmer bounds or Mordell--Weil ranks.
Positive finite-quotient escape proves non-membership in the generic subgroup;
failure to escape the chosen quotient products does not prove dependence.
Nothing here supplies the fifteen quotient directions required for rank 32.
