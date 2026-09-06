# Fixed collision primes control lifting a rational product point

For both successful quartets, the remaining native squareclass condition
can be stated using a fixed finite set of primes, determined by the
cover equations before exceptional point coordinates are supplied.
The real condition is automatic on their product curves.

The full carrier is the smooth projective normalization of
`C4: u_i^2=f_i(t)` for four native quadratics. It has genus 17.
Its quotient by the eight even sign changes is the genus-three curve
H4 below, and C4→H4 is unramified of degree eight. H4 has minimal genus
among the native unramified deck quotients; minimality among arbitrary
maps is UNKNOWN. The [carrier and isogeny certificate](FULL_QUARTETS_HAVE_A_GENUS_THREE_DESCENT_TARGET.md)
establishes those assertions. The present calculation makes the remaining
affine lift condition explicit without choosing a known rational origin.

Let

\[
H_4:\ y^2=\prod_{i=1}^4 f_i(t),\qquad
R=\prod_{i<j}\operatorname{Res}(f_i,f_j),\qquad
S=\{p:p\mid R\}.
\]

For either successful system and any finite non-branch rational point
(t,y) on H4,

\[
\boxed{\text{all four native lifts are rational}
\ \Longleftrightarrow\
v_p(f_i(t))\equiv0\pmod2
\quad(p\in S,\ i=1,2,3).} \tag{1}
\]

The complete sets S contain **18 primes for the +7 quartet** and
**23 primes for the +8 quartet**. These are primes of collision between
different native branch factors, not the bad primes of a particular
specialized elliptic curve.

This is an exact **solubility** criterion conditional on a rational product
point. It does not construct that point, certify its original-curve
quotient rank, or produce a rank score. The fixed labels remain
retrospectively selected.

## Why every other valuation is even

Write t=T/Z with coprime integers T,Z, and homogenize the quadratics:

\[
Q_i(T,Z)=a_iT^2+b_iTZ+c_iZ^2.
\]

Because (t,y) is a rational product point,

\[
\prod_i Q_i(T,Z)=(Z^4y)^2.
\]

If p is outside S, no two binary forms Qi,Qj have a common projective
zero modulo p. At a primitive pair (T,Z), therefore at most one Qi(T,Z)
is divisible by p. The product is a square, so that one valuation must
be even. All four valuations are even at every p outside S.

The independent verifier gives integer-coefficient Bezout identities
for each pair on both parameter charts. Each identity has its exact
pair resultant as the constant right-hand side. They certify the
projective noncollision assertion directly, including primes dividing Z.
The resultant is also recomputed by a four-by-four Sylvester determinant.

At a prime in S the product relation still forces the sum of the four
valuation parities to vanish. Thus zero parity for the first three
forces zero parity for the fourth. A positive rational number whose
valuation is even at every prime is a rational square. This proves (1)
once positivity is established.

The obstruction is intrinsically joint: if one native value has odd
valuation at a rational product point, another must have odd valuation
there too. Those two factors vanish at the same reduced parameter, so
their resultant is divisible by that prime. A lone factor's internal
branch collision cannot by itself support this odd-valuation defect.

## The real condition vanishes for both positive systems

The exact order of the real branch roots determines the sign of every
quadratic on every complementary interval. No rational parameter is
sampled. An independent Sturm computation supplies rational isolating
intervals and verifies the root order using sign variations.

For the +7 system use labels A=01333, B=0b2d0, C=13109, D=19e45.
A and D are always positive. The real branch order is

\[
C_-<B_-<B_+<C_+.
\]

The five open-cell sign patterns, in A,B,C,D order, are

```text
++−+   ++++   +−++   ++++   ++−+
```

The negative regions of B and C do not overlap. Whenever the product is
positive, all four factors are positive. At each branch root the other
three factors are positive as well.

For the +8 quartet, in the order 0911e,0a037,1795d,18f5d, the first
quadratic is always positive. Calling the other three B,C,D, the root
order is

\[
B_-<D_-<D_+<C_-<C_+<B_+.
\]

The sign patterns are

```text
+−++   ++++   +++−   ++++   ++−+   ++++   +−++
```

Again, the three negative regions are disjoint. The product has a
negative leading coefficient for both positive systems, so neither H4
has real points above infinity.

Consequently the proper real map C4(R)→H4(R) is surjective in both
positive cases. This is stronger than merely checking that C4 has some
real point, but it is not sufficient for rational solubility.

## Paired results

| System | Pairwise resultants factored | Collision primes | Every real product point lifts? | Retained rational lift |
|---|---:|---:|---|---|
| +7 quartet, 08234-003 | 6 | 18 | Yes | All collision valuations even |
| +8 quartet, 08234-009 | 6 | 23 | Yes | All collision valuations even |
| Obstructed ABCD quartet | 6 | 25 | No | No global lift; previously excluded at 23 |

The control does have real native lifts, including over infinity.
Surjectivity fails because some real product components have two
negative factors: A and B together, or C and D together. Its full
criterion therefore retains the positivity test as well as valuations.

All eighteen resultants were completely factored within the declared
one-worker, sixty-second-per-quartet bounds. Their prime factors were
verified with proof enabled. The complete prime lists and factorizations
are retained in the certificates. Some collision primes are large:
the sets are not merely a list of small primes suitable for cheap
residue-field enumeration. Parity of the valuation of an already given
rational value can nevertheless be checked by repeated exact division,
without enumerating residues modulo a large prime.

## Why this does not contradict the earlier finite-local-test limitation

On the rank-three elliptic pair carrier, a further radicand h(t) can
have an odd valuation at a new prime: there is no product relation
forcing a second factor to vanish there. A fixed finite collection of
local tests cannot characterize all rational lifts on that elliptic base.

On H4, rationality of y supplies the missing product-square relation.
That relation confines every remaining valuation defect to S. This is
why moving to the higher-genus product curve turns the root-lifting
condition into the fixed finite test (1).

The word **rational** is essential. The control's smooth Q23 product
point has four unit values modulo 23, with the C and D units nonsquare.
Even valuation at 23 alone does not make those units square. Equation
(1) concerns rational values, checks all collision primes, and includes
the real signs when needed. Passing a local parity check on a merely
local product point does not invoke the global square criterion.

Nor is the finite valuation test a claim that the product curve satisfies
a local-to-global principle. Finding H4(Q) with the prescribed parity
pattern remains the unresolved global problem. Its Jacobian/Selmer
classes still must meet the embedded curve.

## What the calculation makes available

The arithmetic target for either positive quartet is now concrete:

\[
\text{a rational point on the fixed genus-three product curve}
\ +\ \text{zero valuation defects at the fixed collision primes}
\ \Longrightarrow\ \text{four rational native lifts}.
\]

For the retained specializations the original-curve independence
certificates then give three quotient directions over the rank-17
generic subgroup. Four retained directions on the +7 fibre and five on
the +8 fibre remain outside these quartets.

The real surjectivity and collision support are equation-defined
**solubility features**. They do not demonstrate new Mordell–Weil
incidence, and they do not measure point visibility. The small paired
sample does not validate either as a rank predictor. The useful advance
is a fixed arithmetic specification for the class that a future global
rational-point computation must reach.

## Replay

The [protocol](COLLISION_PRIME_LIFT_PROTOCOL.json) fixes the three existing
systems and contains no rational-parameter or point search. The outputs
are `rank_jump_collision_prime_lift_v1.json` and
`rank_jump_collision_prime_lift_verification_v1.json` under
`artifacts/generated-results/elliptic-curves/`.

```sh
sage -python elliptic-curves/rank-jump/verify_collision_prime_lift.py check
```

The replay checks all prime products, primality, both-chart Bezout
identities, Sturm isolations, sign cells and the two retained parity
certificates. No active search file or external service was changed.
