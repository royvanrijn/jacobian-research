# Quadratic-twist rank census for the complete bisection batch

## Status

<!-- status-consumer: EC-K3-BISECT-MULTIQUADRATIC-CHARACTERS dc58103d8d2494cf -->

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

## Product result in the discovery window

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
`0.25`, and none exceed `0.50`.  These were weak discovery-window signals,
not durable leaders:

| rank | masks | six block scores | weakest | mean | paired-base rank lower bound |
|---:|---|---|---:|---:|---:|
| 1 | `27431:92937` | `0.401, 0.575, 0.645, 0.539, 0.710, 0.287` | 0.287 | 0.526 | 5 |
| 2 | `63466:79888` | `0.267, 0.390, 0.459, 0.606, 0.597, 0.287` | 0.267 | 0.434 | 4 |

The two paired bases already known to have rational rank at least nine do not
lead this generic twist census:

```text
42110:43109    census rank 270    weakest -0.106    mean 0.313
71804:81769    census rank 303    weakest -0.120    mean 0.191
```

Base rank and product-twist rank are different invariants.  The rank-nine
bases remain excellent sources of rational specialization parameters, but
the discovery data alone do not prioritize them for a generic rank-20 proof.

## Strict targeted gate for product 27431:92937

The requested theorem bet uses a target-only gate: no population ranking or
replacement product is allowed to authorize a section solve.  The product
`27431:92937` was scored on two fresh, mutually disjoint 48-prime windows,
`499--821` and `823--1151`.  Together they comprise twelve disjoint eight-prime
blocks, all above 491 and extending through 1151.  The target fails both
holdouts:

| window | weakest block | mean block score |
|---|---:|---:|
| `211--491` discovery | `0.287` | `0.526` |
| `499--821` holdout | `-0.494` | `0.017` |
| `823--1151` holdout | `-0.077` | `0.150` |

This triggers the declared immediate stop.  No `chi=4` polynomial-section
solve and no denominator layer was run for `27431:92937`.  The directly
runnable target-only command is recorded in [`../REPRODUCE.md`](../REPRODUCE.md).
This is a bounded negative heuristic, not a rank-zero result.

### Separate exploratory context

Separate local exploratory artifacts exist outside the strict target-only
replay.  They include a complete rescore of all 5,566 products on `499--821`
and a small multi-target replay on `823--1151`.  They are retained for
provenance, but are not part of the theorem gate and do not authorize a
replacement product search.  The complete rescore's top 100 and the discovery
top 100 have only one key in common, `36197:51178`:

```text
window 211--491   rank 17   weakest 0.108   mean 0.337
window 499--821   rank  3   weakest 0.174   mean 0.493
window 823--1151           weakest -0.200  mean 0.062
```

The discovery leader and the exploratory new-window leader also fail
independently:

| product | `211--491` weakest / mean | `499--821` weakest / mean | `823--1151` weakest / mean |
|---|---:|---:|---:|
| `27431:92937` | `0.287 / 0.526` | `-0.494 / 0.017` | `-0.077 / 0.150` |
| `27431:49826` | not discovery top 100 | `0.291 / 0.539` | `-0.168 / 0.070` |
| `71804:81769` | `-0.120 / 0.191` | `0.096 / 0.489` | `-0.067 / 0.174` |

Thus the exploratory population data also have no stable positive signal
across all three windows.  This remains only a negative heuristic and does
not prove that any product twist has rank zero.

The singleton controls behave differently.  Mask `18075` remains positive
on both holdouts:

```text
window 211--491   weakest 0.853   mean 1.003
window 499--821   weakest 0.313   mean 0.740
window 823--1151  weakest 0.193   mean 0.893
```

This is consistent with its forced rank-one direction, but does not isolate a
second direction.  Mask `24868` is likewise stable on both holdouts, again as
a rank-one control.

The holdout artifacts are:

- [`../artifacts/generated-results/elkies-2026-quadratic-twist-targeted-holdout-p499-821.json`](../artifacts/generated-results/elkies-2026-quadratic-twist-targeted-holdout-p499-821.json),
  SHA-256 `03661d5b772bc1d224ea84d2d7436d0eb47e9173285e7733a84fb8133fc8d30d`;
- [`../artifacts/generated-results/elkies-2026-quadratic-twist-product-census-p499-821.json`](../artifacts/generated-results/elkies-2026-quadratic-twist-product-census-p499-821.json),
  SHA-256 `fedc6e7038ca55d87a8de4f5c8d78a2924a006174d35aaeeebd99b5059130097`;
- [`../artifacts/generated-results/elkies-2026-quadratic-twist-targeted-holdout-p823-1151.json`](../artifacts/generated-results/elkies-2026-quadratic-twist-targeted-holdout-p823-1151.json),
  SHA-256 `43a76b955bbf8c631aaf9eda66d5729fde1119f921b70119b4360a3783ab3210`.

Their argv records are stored in `reproducing_command`; prepend
`.venv/bin/python` because the Python script is not executable.  The strict
target-only replay does not require the complete product-census artifact.

## Exact descent of the singleton control

For mask `18075`, the certified bisection record gives a lifted point
`P(t,u)` on `u^2=q(t)`.  The script
[`scripts/derive_elkies_2026_singleton_twist_section.sage`](scripts/derive_elkies_2026_singleton_twist_section.sage)
forms

```text
R=P-sigma(P),
X=q*x(R),
Y=q^2*coefficient_u(y(R)).
```

It verifies exactly that

```text
Y^2 = X^3 + A*q^2*X + B*q^3,
deg X = 6,   deg Y = 9.
```

This is the known twist section in the full `chi=3`, `P.O=0` polynomial
degree box.  Its corresponding anti-invariant point after quadratic base
change has certified height 12.  The exact coordinate artifact is
[`../artifacts/generated-results/elkies-2026-singleton-twist-section-mask-18075.json`](../artifacts/generated-results/elkies-2026-singleton-twist-section-mask-18075.json).
Its SHA-256 is
`3cc4fdd9651ced2e267ea5a36ab9670a9baf4328ffaf3f84cdd763a051c02c3a`.

Replay it with the bundled Sage runtime path appropriate to the machine:

```bash
sage -python \
  elkies-k3/scripts/derive_elkies_2026_singleton_twist_section.sage \
  --mask 18075 --prime 37
```

## Complete modular polynomial-section scheme for mask 18075

The exporter
[`scripts/export_elkies_2026_twist_polynomial_sections_modp.sage`](scripts/export_elkies_2026_twist_polynomial_sections_modp.sage)
uses the global bounds `deg X<=6`, `deg Y<=9`.  At `p=37` it moves the smooth
fibre `t=0` to infinity.  That fibre has no rational 2-torsion.  Fixing each
of its 48 affine points and recursively eliminating the lower `Y`
coefficients produces 48 systems of nine equations in six `X` coefficients.
The two signs of `Y` give identical systems, hence 24 distinct systems.

```bash
sage -python \
  elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage \
  --singleton-mask 18075 --prime 37

.venv/bin/python \
  elkies-k3/scripts/run_elkies_2026_twist_polynomial_sections_msolve.py \
  --export artifacts/local/elkies-k3/twist-polynomial-sections/\
singleton-18075/p37/export.json \
  --threads 4 --jobs 2 --timeout 180
```

All 24 distinct systems completed.  Twenty-three are empty over
`algebraic_closure(F_37)`.  The remaining system has quotient dimension one,
and covers exactly the leading points `(4,9)` and `(4,28)`.  The exact known
section reduces to the second point, solves all nine exported equations, and
has

```text
(x5,x4,x3,x2,x1,x0) = (2,30,18,29,9,21) mod 37.
```

The solver's linear parametrization is the same point.  Consequently the
complete affine `P.O=0` polynomial-section scheme in this characteristic and
chart consists only of `Q` and `-Q`.  The summary is
[`../artifacts/generated-results/elkies-2026-twist-polynomial-sections-singleton-18075-p37-msolve.json`](../artifacts/generated-results/elkies-2026-twist-polynomial-sections-singleton-18075-p37-msolve.json).
Its SHA-256 is
`41d5dfc76fe1aad67336ec7b2038c9d9865028e64e8a9c2442108863babb1c24`.

This is a complete finite-field scheme computation, not a proof that
`rank E^q(QQ(t))=1`.  A characteristic-zero section can have positive
intersection with the zero section, lie outside this degree box, or have bad
coefficient reduction at 37.  Even another characteristic-zero section in
the same box could specialize to a boundary point of a projective
compactification.  None of those possibilities is excluded here.

## Separate product-system complexity experiment

This separate experiment concerns the different product `36197:51178` and is
not stage two of the `27431:92937` gate.  The only exploratory top-100
overlap was exported at `p=19` in the natural
`chi=4` box `deg X<=8`, `deg Y<=12`.  It gives 12 leading-point blocks and six
distinct systems of 12 equations in eight variables.  A bounded 60-second,
eight-thread run on the first distinct system timed out without a solution
classification.  Its process-tree high-water snapshot at timeout was
2,791,816 KiB (about 2.66 GiB).  The reproducible summary is
[`../artifacts/generated-results/elkies-2026-twist-polynomial-sections-product-36197-51178-p19-msolve.json`](../artifacts/generated-results/elkies-2026-twist-polynomial-sections-product-36197-51178-p19-msolve.json).
Its SHA-256 is
`55f7ce7559ccdfbf7ec4efe5fd973f1af9e6294ef1554a243dc93b961e5b6c2d`.

This is only a complexity experiment.  It neither finds nor excludes a
product-twist section.  Together with the failed third prime window, it says
not to spend the current search budget on the unsliced eight-variable system.

## What remains for an exact twist rank

The Frobenius score is justified as a ranking heuristic by the
Rosen--Silverman/Nagao framework and by later work on twist families; see
[Rosen--Silverman, *On the rank of an elliptic surface*](https://doi.org/10.1007/s002220050238)
and [Kim, *The Sato--Tate conjecture and Nagao's conjecture*](https://arxiv.org/abs/1712.02775).
Neither source turns a finite prime window into a rank bound.

There are now three honest next proof routes:

1. Compute a function-field Selmer upper bound for a singleton twist.  The
   installed Sage build has no rank, Selmer, or two-descent implementation for
   elliptic curves over `QQ(t)`, so this needs Magma or a new descent
   implementation.
2. Bound the Picard number of the `chi=3` twist surface after good reduction
   and apply Shioda--Tate.  This requires a surface Frobenius polynomial (or a
   comparably strong cohomological certificate), not just the `n=1` fibral
   averages used here.
3. Continue constructive searches, but widen in a controlled order:
   singleton sections with `P.O>0` first, and product systems only after a new
   independent prime signal or a dimension-reducing ansatz appears.

At present no new rational section has been found.  Therefore the certified
generic-rank records remain 19 on paired genus-one bases and 18 on each
individual rational quadratic cover; no rank-20 paired base or rank-19
rational cover is claimed.
