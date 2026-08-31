# Complete paired-cover geometry after bisection injectivity

## Status

This note records two exact results obtained from the complete 39,120-record
equation-level bisection batch.

1. Every individual quadratic base is a rational conic over `QQ`.  Thus all
   39,120 bisections give parameterizable `P1` base changes of generic
   Mordell--Weil rank at least 18.
2. Every two distinct records have disjoint geometric branch divisors.  Their
   connected biquadratic fibre product has genus one, and the two new sections
   have exact height matrix

   ```text
   [24  0]
   [ 0 24].
   ```

   Consequently all

   ```text
   binomial(39120,2) = 765167640
   ```

   distinct pairs define genus-one base changes of generic Mordell--Weil rank
   at least 19.

This is not a quadratic squareclass collision.  The complete injectivity
theorem in [`BISECTION_COLLISION_SEARCH.md`](BISECTION_COLLISION_SEARCH.md)
still says that no two bisections split over the same quadratic field.  The
rank-two construction instead uses the compositum of two distinct quadratic
fields and its two different nontrivial Galois characters.

## Exact complete classification

Write the branch equation of record `i` as

```text
u_i^2 = q_i(t),             deg(q_i)=2.
```

The exact scan performs the following checks for all records:

- exact Hasse--Minkowski solubility of the projective conic
  `U^2=q_i(T/Z)Z^2` over `QQ`;
- irreducibility of the primitive quadratic defining the geometric branch
  divisor;
- uniqueness of all 39,120 primitive branch quadratics.

The result is:

```text
Q-rational conics                 39120
irreducible branch quadratics     39120
distinct geometric branch sets   39120
anisotropic conics                    0
```

The rational constant squareclass must be retained in the conic-solubility
test.  Removing it is valid for geometric branch comparison but can change a
conic into its anisotropic constant twist.  The verifier deliberately uses
the original `q_i` for Hasse--Minkowski and a primitive polynomial only for
the branch key.

For distinct `i,j`, the two squareclasses are independent and their branch
sets are disjoint.  The compositum is therefore a connected `V4` cover with
four inertia points of order two.  Riemann--Hurwitz gives

```text
2g-2 = 4*(-2) + 4*(4-2) = 0,
```

so `g=1`.  There is no shared-branch genus-zero paired base anywhere in this
catalogue.

On the individual quadratic cover, the previously certified new
anti-invariant section has height 12.  Pulling it through the other quadratic
extension doubles its height to 24.  The two sections transform through
different nontrivial characters of `V4`, so Galois invariance forces their
cross-height to vanish.  This proves the displayed rank-two matrix.  The
invariant rank-17 lattice lies in the trivial character and is independent of
both, proving generic rank at least 19.

Replay the complete classification with:

```bash
/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/analyze_elkies_2026_bisection_pair_covers.sage
```

The result is
[`../artifacts/generated-results/elkies-2026-bisection-pair-cover-geometry-full.json`](../artifacts/generated-results/elkies-2026-bisection-pair-cover-geometry-full.json),
with SHA-256

```text
941b2072d19347d6092ad8eba5ffe2ed7b08c433426360237313241a730172c6
```

The compact output retains the first 100 equation-cheapest exact conic points
and a hash of the complete 39,120-point solver ledger.  The replay itself, not
the sample, is the complete solubility check.

## Immediate rational genus-one bases

There are 104 covers whose actual quadratic leading coefficient is a rational
square and 21 whose actual constant coefficient is a rational square; one
cover has both properties.  Hence inclusion--exclusion gives

```text
binomial(104,2) + binomial(21,2) - binomial(1,2) = 5566
```

paired bases with an immediate rational point over `t=infinity` or `t=0`.
Each is therefore an elliptic curve over `QQ`, not merely a genus-one torsor.
Their Mordell--Weil ranks can be screened through the third quotient

```text
w^2 = q_i(t) q_j(t).
```

The paired curve maps to this quartic by `w=u_i u_j`.  Since the two branch
sets are disjoint, the product involution is fixed-point-free.  After choosing
a rational origin this degree-two map is a 2-isogeny, so the paired curve and
the quartic Jacobian have the same rational rank.

## First new low-complexity rank-19 base

The two equation-cheapest covers in the common-infinity subfamily have orbit
masks `0x06e61` and `0x13f69`, global priority ranks 729 and 893, and equations

```text
u^2 = 41627760409 + 15206854416*t + 2278725696*t^2,
v^2 = 126480025   +   108563070*t +   21650409*t^2.
```

Their common rational point at infinity is

```text
u/t = 47736,       v/t = 4653.
```

The third quotient has quartic coefficients, low to high,

```text
5265080177224330225,
6442600793932536030,
2840374122051096801,
576620074955502864,
49335343317209664.
```

The classical binary-quartic invariants and exact minimization give its
Jacobian

```text
y^2 + x*y + y = x^3 + x^2
                - 10283115414666818054869518594495*x
                + 237013667219266831004461300259717721649810997157.
```

Three displayed rational points are certified independent by exact finite
quotients at `53`, `59`, and `61`; good reduction at `31`, where the group has
order 40, supplies the 3-torsion exclusion needed for infinite descent.  Thus
the Jacobian and the paired base have rational rank at least 3.  In particular
the base has infinitely many rational points and the map to `P1_t` assumes
infinitely many rational values.

The exact replay is:

```bash
/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/verify_elkies_2026_rank19_cheapest_infinity_pair.sage
```

It writes
[`../artifacts/generated-results/elkies-2026-rank19-cheapest-infinity-pair.json`](../artifacts/generated-results/elkies-2026-rank19-cheapest-infinity-pair.json),
with SHA-256

```text
7987ce9e9a40b2693d00abe22689b8e1db9613173f681311d1e073e2c542796b
```

The certificate proves base rank at least 3, not an upper bound.  It proves
generic surface rank at least 19 over the paired base; it does not claim an
exceptional specialization rank.

## Search-policy correction

The new pair is not one of the 8,895,801 norm-four disjointness edges.  The
published rank-19 pair is outside that graph as well.  Distinct Galois
characters, rather than norm-four disjointness of the two bisection classes,
are what prove independence on a biquadratic cover.  Therefore:

- retain the disjoint graph as a geometric or equation-cost priority score;
- do not use it as a hard filter for paired rank-19 arithmetic;
- begin with the 5,566 immediate-point pairs, then widen to locally soluble
  genus-one torsors among all 765,167,640 pairs;
- rank the base Jacobians before specialization scoring, since positive base
  rank supplies infinitely many rational parameters.

## Complete 5,566-base arithmetic census

The planned immediate-point computation is complete.  The exact catalogue
contains a global minimal Jacobian, conductor, root number, quartic, and
mask-attached cover equations for every one of the 5,566 bases:

```text
global root number +1   2823
global root number -1   2743
```

Its SHA-256 is

```text
d27233ba066f311e97b5d9400838417f7aaa687d5a1b0eed11726b41538ea831
```

A complete effort-one PARI point pass followed by exact finite-quotient
independence certificates gives the following lower-bound ledger:

```text
certified lower bound 0   1786
certified lower bound 1   1869
certified lower bound 2   1197
certified lower bound 3    454
certified lower bound 4    182
certified lower bound 5     56
certified lower bound 6     20
certified lower bound 9      2
```

Here lower bound zero means only that the bounded PARI pass returned no point;
it is not a rank-zero theorem.  Every positive number is unconditional because
the displayed points pass exact finite-quotient infinite descent.  The complete
ledger has SHA-256

```text
c28094ba4570b52c784306edef3e90b46e3f4ef5cf3964c157c7a19bb5c0b43b
```

## Rank-at-least-nine base

The arithmetic-complexity-rank-114 leader uses masks `0x0a47e` and `0x0a865`:

```text
u^2 = 4225*t^2 + 38636*t + 289444,
v^2 = 1346816601*t^2 + 7403338254*t + 10921221529.
```

It has the rational point at infinity `u/t=65`, `v/t=36699`.  Its third
quotient is 2-isogenous to the paired base and has global minimal model

```text
y^2 + x*y = x^3
            - 70087047578007713577216*x
            + 3865770423647395544516350651140096.
```

Nine displayed points are independent modulo 3 by the combined exact
finite-quotient matrix at `17,19,53,71,101,107,127,137`; good reduction of
order 28 at 23 excludes rational 3-torsion.  Thus the paired base has rank at
least 9, infinitely many rational points and infinitely many rational
`t`-values.  The surface over its function field retains the exact
`diag(24,24)` new-section height matrix and generic rank at least 19.

This infinitude is now operational rather than only abstract.  Put

```text
u = 65*t + r,       t = (289444-r^2)/(130*r-38636).
```

After substituting in the second conic and setting
`z=(130*r-38636)*v`, the paired base is the pointed quartic

```text
z^2 = 1346816601*r^4 - 962433973020*r^3
      - 309051947898044*r^2 + 168863136988245440*r
      + 46344697121074403584.
```

The reciprocal chart `s=1/r`, `y=z/r^2` has the rational origin
`(s,y)=(0,36699)`.  The exact pointed-quartic transformation gives the paired
base itself, not the binary-quartic third quotient, with global minimal model

```text
y^2 + x*y = x^3
            - 60729194722297004073216*x
            + 5758259762216167074332597509226496.
```

It is the degree-two isogenous companion of the displayed third-quotient
model.  The verifier records the isogeny's rational `x`- and `y`-maps, the minimal-to-pointed
Weierstrass change of variables, the inverse pointed-quartic formula, and nine
exact rational `(t,u,v)` points obtained from the certified independent
third-quotient points.  For example, the shortest stored parameter is

```text
t = -8119579772667928420880724897627
    /1535857632011210964288249442454.
```

The nine transported points remain independent because a nonzero isogeny is
injective after tensoring the Mordell--Weil groups with `QQ`.  Choosing either
sign of an image only shortens its displayed `t`-coordinate and does not alter
that argument.  These values are seeds for specialization searches; no
specialized surface rank is asserted here.

Replay the promoted certificate with:

```bash
/tmp/jacobian-sage-bin/sage -python \
  elkies-k3/scripts/verify_elkies_2026_rank19_rank9_base.sage
```

The output SHA-256 is

```text
a298cd16316bed5584dd3729cbbf63ad9b5e6fbbb2259a70f30ddb233ae7b395
```

The second rank-at-least-nine base has masks `71804:81769`.  Neither leader
lies in the norm-four disjoint graph.  The census and promoted verifier are
independent of the conductor-first and residual-2-Selmer programmes under
`elliptic-curves/`.

<!-- status-consumer: EC-K3-BISECT-BIQUADRATIC-R19 707bffd8b85f8f3e -->
