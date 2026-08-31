# Quadratic-twist rank census for the complete bisection batch

## Status

The complete bisection injectivity theorem distinguishes all 39,120 branch
squareclasses.  It does **not** determine the Mordell--Weil rank of any
quadratic twist.  This note records the exact multiquadratic rank theorem and
a first complete bounded Frobenius census of

- all 39,120 singleton twists `E^{q_i}/QQ(t)`; and
- all 5,566 product twists `E^{q_i*q_j}/QQ(t)` attached to paired bases with
  an immediate rational point.

The census is a heuristic ranking.  It constructs no new section and proves
no new twist rank.

## Exact character theorem

The reusable statement and proof are Theorem F4 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
For independent squareclasses `q_1,...,q_k`,

```text
rank E(QQ(t)(sqrt(q_1),...,sqrt(q_k)))
  = sum over S subset {1,...,k} rank E^{q_S}(QQ(t)).
```

If the quadratic branch divisors are pairwise disjoint, the multiquadratic
base is geometrically connected and has

```text
g = 1 + 2^(k-1)*(k-2).
```

For the published rootless R17 surface, the trivial character contributes
rank 17 and the `k` known singleton characters contribute one direction each.
Their height block after the full base change is

```text
12*2^(k-1) * I_k.
```

Thus the unconditional lower bound is `17+k`.  Every nonempty product twist
contributes its full rank as an additional orthogonal character.  In
particular,

```text
rank E^{q_i} >= 2       => a rational P1 base of generic rank at least 19,
rank E^{q_i*q_j} >= 1   => the paired genus-one base has generic rank at least 20.
```

## Cheap Frobenius census

At a usable prime and finite parameter,

```text
a_p(E^q_t) = chi_p(q(t))*a_p(E_t).
```

The script
[`scripts/screen_elkies_2026_quadratic_twist_ranks.py`](scripts/screen_elkies_2026_quadratic_twist_ranks.py)
precomputes `a_p(E_t)` once on the compact published R17 model and evaluates
every twist by character dot products.  It preserves the rational constant
squareclass of every `q`; a primitive geometric branch polynomial is not a
valid substitute when its discarded scalar is nonsquare.  Singular fibres
have trace zero, and primes at which the twist squareclass or branch divisor
has bad reduction are skipped for that candidate.

For each prime,

```text
A_p(E^q) = (1/p) * sum(t in F_p) a_p(E^q_t).
```

Each disjoint block is scored by

```text
sum_p(-A_p(E^q)*log(p)) / sum_p(log(p)),
```

and candidates are ranked by the weakest of six blocks.  This is the
finite-block form of a Nagao-style search score.  Its value is not a rank
bound.

Replay the complete census with:

```bash
.venv/bin/python \
  elkies-k3/scripts/screen_elkies_2026_quadratic_twist_ranks.py
```

The output is
[`../artifacts/generated-results/elkies-2026-quadratic-twist-rank-census.json`](../artifacts/generated-results/elkies-2026-quadratic-twist-rank-census.json),
with SHA-256

```text
731966cd1dd3d55e4f9e72bec36c0ca0a9540e321d5adcd52f701b00e90c9d81
```

It scores 44,686 twists over 48 primes in six disjoint eight-prime blocks.
All 48 primes are between 211 and 491.  The complete singleton and product
score ledgers are committed by SHA-256 inside the compact output.

## Singleton result: no rank-two signal yet

Every singleton twist already has a known non-torsion direction, so the
39,120 singleton population is an internal rank-at-least-one calibration.
Its block-score distribution is:

```text
                              weakest block       mean of six blocks
median                            0.264                  0.620
90th percentile                  0.475                  0.765
99th percentile                  0.633                  0.886
maximum                           0.864                  1.067
```

No singleton has weakest-block score at least one.  The leaders look like a
particularly stable version of the already forced rank-one signal, not a
separate rank-two cluster.

The best score/equation-cost compromise is mask `18075 = 0x0469b`:

```text
census rank          2
priority rank        448
equation rank        302
six block scores     1.082, 0.853, 1.013, 0.860, 1.008, 1.204
weakest / mean       0.853 / 1.003
```

It is the first singleton to use for a bounded second-section search, but its
score does not itself suggest rank at least two.

## Product result: two weak but uniformly positive leaders

The product-twist population is centred near zero, as expected when most
product characters have rank zero:

```text
                              weakest block       mean of six blocks
median                           -0.525                  0.003
90th percentile                 -0.196                  0.221
99th percentile                  0.043                  0.396
maximum                           0.287                  0.713
```

Only 95 of 5,566 products have positive weakest-block score, only two exceed
`0.25`, and none exceed `0.50`.  Splitting the primes into the first and last
three blocks produced only one common top-100 candidate.  That survivor is
the full six-block leader:

| rank | masks | six block scores | weakest | mean | paired-base rank lower bound |
|---:|---|---|---:|---:|---:|
| 1 | `27431:92937` | `0.401, 0.575, 0.645, 0.539, 0.710, 0.287` | 0.287 | 0.526 | 5 |
| 2 | `63466:79888` | `0.267, 0.390, 0.459, 0.606, 0.597, 0.287` | 0.267 | 0.434 | 4 |

The leading product is

```text
q_27431 = 501490089 + 228711132*t + 34833604*t^2,
q_92937 = 1055919841 + 328702636*t + 14122564*t^2.
```

Its paired base has a rational point at infinity, arithmetic-complexity rank
680, global root number `-1`, and an independently certified rational rank
lower bound 5.  These base facts do not prove anything about the product
twist, but they make the candidate operationally useful if an extra character
section is found.

The two paired bases already known to have rational rank at least nine do not
lead this generic twist census:

```text
42110:43109    census rank 270    weakest -0.106    mean 0.313
71804:81769    census rank 303    weakest -0.120    mean 0.191
```

Base rank and product-twist rank are different invariants.  The rank-nine
bases remain excellent sources of rational specialization parameters, but
the present data do not prioritize them for a generic rank-20 proof.

## Recommended exact follow-up

1. Start with product `27431:92937`.  Search the twist
   `E^{q_27431*q_92937}` first for polynomial sections in the natural
   `chi=4` bounds, then widen to controlled denominators.  Any non-torsion
   section proves generic rank at least 20 on its genus-one paired base by
   Theorem F4.
2. Retain `63466:79888` as the independent second product candidate.  Its six
   blocks are less strong but still uniformly positive.
3. For a rational rank-19 `P1` base, use singleton mask `18075` as the first
   second-section target because it combines the second-best robust score
   with equation rank 302.
4. Before a large section solve, extend the leading candidates to new disjoint
   prime blocks above 491.  The first 24-prime product leaderboard was highly
   unstable; only the six-block survivor should receive expensive work.

An empty section ansatz, a low score, or the absence of a held-out survivor is
only a bounded negative experiment.  Exact rank growth begins only with a new
rational twist section and an independence certificate.
