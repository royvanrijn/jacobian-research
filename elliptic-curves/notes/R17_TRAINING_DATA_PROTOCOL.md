# Deliberate training data for published R17 specializations

Date: 2026-09-02

## Answer to the rank-27/28 rediscovery question

The existing weakest-block Nagao score places the published rank-at-least-27
and rank-at-least-28 parameters at ranks `422873` and `55387` among all
`121589944` primitive parameters of projective height at most `10000`. Their
population fractions are `0.3478%` and `0.0456%`. Thus a top-one-percent
budget mechanically retrieves both.

That is not yet a blinded rediscovery result. The score implementation and the
four positive-control gate entered the repository together, and the
rank-25--28 fibres were explicitly used to accept the score. The exact claim
is therefore:

> the current score strongly ranks the known controls, but no score trained
> without those controls has yet rediscovered them.

The protocol below creates the missing prospective experiment. All four
published parameters are absent from feature generation, model fitting,
feature selection, hyperparameter selection, stopping rules, and budget
selection. Rank 27 and rank 28 are the primary quarantined replay; rank 25 and
rank 26 remain secondary checks. Because these parameters are already public
and known to the investigators, the eventual result should be called a
**mechanically quarantined replay**, not a human-blind discovery. A truly
blind successor requires an external custodian or a future undisclosed fibre.

## Frozen population and cheap features

Run:

```bash
python3 elliptic-curves/scripts/build_r17_training_dataset.py \
  --count 100000 \
  --height 10000 \
  --seed 20260902 \
  --lane-size 1000 \
  --cover-panel-size 128 \
  --summary artifacts/local/elliptic-curves/r17-training-summary.json
```

The command samples finite primitive `(a:b)` uniformly from the box
`|a|<=10000`, `1<=b<=10000`, rejects duplicates, and hard-excludes

```text
-2/377, -308/251, 2456/135, -9529/5471.
```

The generated feature population and selected cohort are raw experimental
data under `artifacts/local/elliptic-curves/`; they are not mathematical
status artifacts. On the recorded local run the population has `100000` rows and the
five lanes have `1000` rows each, with `4922` distinct selected parameters and
`77` parameters belonging to more than one lane. The deterministic split is
`70092 / 15029 / 14879` for train, validation, and internal test.

The feature levels are explicit:

- **Level 0:** normalized parameter, numerator, denominator, and projective
  height.
- **Level 1:** the existing 102-prime Nagao data in three disjoint blocks,
  including every block score, weakest-block score, mean score, and counts of
  good and bad reductions. The implementation reproduces the four stored
  control block scores exactly, but does not put those rows in the dataset.
- **Level 2a:** a partial local conductor-quality proxy. At the declared small
  primes, the short model is locally scaled and the tame `p>=5` contribution
  to `log|Delta_min|-log(N)` is recorded. This is not an exact conductor; 2,
  3, and untested prime factors remain outside the feature.
- **Level 2b:** the exact dimensions of `E(F_p)/2E(F_p)` and
  `E(F_p)/3E(F_p)` at five small primes. This is a finite-group structural
  code, not the quotient of the specialized Mordell--Weil group by R17. Code
  rarity is fitted on the train split only, with add-one smoothing.
- **Level 2c:** modular quadratic-character diversity across a fixed,
  hash-selected panel of 128 published bisection covers and five primes. It is
  a cheap cover-signature feature, not a count of rational split covers.

These boundaries matter. Raw discriminant radicals are not conductors, finite
reduction codes are not Selmer dimensions, and a modular cover signature is
not a new rational point.

## Stratified deeper cohort

The selector retains five independently auditable lanes:

1. top one percent by weakest Nagao block;
2. highest partial conductor-quality proxy;
3. rarest finite quotient code;
4. highest modular cover-character diversity;
5. deterministic random controls.

Each 1,000-row lane has a frozen `700 / 150 / 150` quota across train,
validation, and internal test, so a rare code confined to a held-out split
cannot consume the training cohort.

Overlaps are retained as multiple lane flags, not silently deduplicated away.
The random lane is immutable and is the population-baseline estimator. Every
job must retain the original lane flags and deterministic train/validation/test
split. Development uses only the train rows; the internal test is opened once
after the model and compute budget are frozen.

Before labels are fitted, exact minimal models and `j`-invariants must be used
to group any repeated isomorphism or quadratic-twist class into one split.
The current parameter-hash split is provisional until that duplicate audit is
complete; row-wise separation alone is not a leakage certificate.

This feature builder is the population-data extension of the existing
[`RANK_JUMP_LABORATORY.md`](RANK_JUMP_LABORATORY.md). New certified labels and
the eventual frozen-score replay should be registered in that laboratory,
not reported through a second metric convention. The quotient-height and
degree-two visibility response schema for the four embargoed controls is
already recorded separately in
[`ELKIES_RANK_JUMP_FINGERPRINTS.md`](ELKIES_RANK_JUMP_FINGERPRINTS.md); none of
those response fields enters this feature population.

The full 39,120-cover split census is appropriate only after this cheap pass.
For each selected fibre it can provide exact split-cover counts, exact points
coming from the known bisections, and finite-reduction lower bounds for the
subgroup they generate. It must not be described as a full rank computation.

That census has now been run on all `4922` selected rows:

```bash
python3 elliptic-curves/scripts/label_r17_training_bisections.py \
  --workers 4 --prime-bound 199
```

It performed `192548640` exact square tests and found 349 nonzero split covers
on 251 fibres. Exact specialization and finite-reduction replay certify 237
new directions on 212 fibres. Six rows are censored because the stored generic
seventeen have finite-quotient rank only 16 in the fixed prime ensemble; a
split point on those rows is not promoted to an eighteenth direction. There
were no ramified covers and no repeated exact `j`-invariants among the 4922
rows. The labels and their compact summary are local raw data:

```text
artifacts/local/elliptic-curves/r17-training-bisection-labels.jsonl
artifacts/local/elliptic-curves/r17-training-bisection-labels-summary.json
```

The lane-level positive rates are `10.5%` for conductor quality, `5.4%` for
Nagao, `2.2%` for cover diversity, `2.2%` for random controls, and `1.6%` for
unusual finite quotient code. These are selected-cohort mechanism-visibility
rates, not population rank probabilities.

## Labels and censoring

Store outcomes as separate fields rather than collapsing them prematurely to
one target:

```text
new_independent_directions_lower_bound
finite_quotient_gain_lower_bound
minimum_exceptional_naive_height
cpu_seconds
peak_rss_bytes
completed_search_bound
residual_two_selmer_dimension
```

The first two are promoted only after exact point membership and independence
replay. `minimum_exceptional_naive_height` is right-censored when the bounded
search finds nothing; it is not infinity. A timeout is missing/censored data,
not label zero. `cpu_seconds / new direction` is undefined when no direction
is found, so the model should use a two-part outcome: probability of at least
one exact gain, followed by cost conditional on gain.

Selection itself is informative. Report results lane-by-lane and against the
random controls. Do not fit an ordinary unweighted regression to the union and
call its output a population probability. The useful optimization target is
expected certified quotient gain per CPU-hour at a fixed review budget.

## Residual-Selmer gate

This programme does not relax the existing fail-closed rank-32 gate. Exact
specialization of the generic 17, evaluation of already certified
multisections, and finite-reduction independence tests can generate training
labels. A new two-cover search, `ratpoints`, slope-box search, or comparable
expensive point search remains forbidden unless the same minimal fibre has a
completed residual two-Selmer quotient of dimension at least 15. Rows that do
not pass that gate have a missing deep-search outcome, not a negative outcome.

## Model and evaluation contract

Start with interpretable baselines before a large model:

1. weakest-block Nagao only;
2. regularized logistic/Poisson or boosted-tree models on Levels 0--2;
3. a cost-aware ranking model for certified gain per CPU-hour.

Freeze the feature list, transforms, weights, hyperparameters, tie-breaks, and
review budget by file hash. Then evaluate the quarantined fibres once. The
primary metrics are their population percentile and recall at fixed budgets
(`top 0.01%`, `0.1%`, and `1%`), plus enrichment over the random lane.
Regression error is secondary because the operational task is prioritization.

A successful answer to the question is:

> a model whose complete development record contains none of the four
> embargoed fibres ranks both `2456/135` and `-9529/5471` inside the frozen
> review budget when the quarantine is opened once.

This criterion was frozen before the first labelled replay below. It was not
met by the learned bisection-gain score.

## First frozen ranker and quarantine result

The first ranker deliberately has no post-holdout tuning. It is the
overlap-adjusted standardized difference between positive and negative class
means on the train split, using the frozen Level 0--2 fields. Its serialized
pre-replay hash is
`dda5ff55a2470ebe42cc8255b7e91d5633d2d9d116c26c7870dc6bf813c0d4b4`.
The exact label semantics are separately pinned by
`b7b3dd079fbf0c786aaae012bc217c2d41c43ce4d489014c02b43649dcc83afc`;
per-row CPU timings are intentionally excluded from that semantic hash.
Run:

```bash
python3 elliptic-curves/scripts/train_r17_bisection_ranker.py
```

The mechanically quarantined replay is pinned in
[`r17_bisection_gain_ranker_quarantined_replay_v1.json`](../../artifacts/generated-results/elliptic-curves/r17_bisection_gain_ranker_quarantined_replay_v1.json).
It answers the primary question negatively for this learned target:

| control | learned rank in sampled 100,000 | weakest-block Nagao rank |
| --- | ---: | ---: |
| rank at least 25, `-2/377` | 580 | 49 |
| rank at least 26, `-308/251` | 232 | 509 |
| rank at least 27, `2456/135` | 784 | 344 |
| rank at least 28, `-9529/5471` | 24,210 | 49 |

Thus the learned bisection-gain score retrieves rank 27 but not rank 28 at a
one-percent budget; the untouched Nagao baseline retrieves both. This is an
empirical rank in the sampled population, not the exhaustive height-10,000
rank quoted at the start of this note.

The failure is informative. The strongest learned coefficients reward small
projective height and partial conductor quality, matching the observed
bisection-splitting response. The published rank-28 fibre is high in the box
and only one of its eleven known exceptional directions is visible to the
complete bisection atlas. Bisection visibility is therefore too narrow a
surrogate for extreme total rank. Because the quarantine has now been opened,
no retuned score may be called a fresh replay on these four controls; further
models are post-unblinding development unless evaluated by a new custodian or
future undisclosed fibres.

## Prospective validation of the narrower target

To distinguish target mismatch from ordinary overfitting, a second cohort was
committed before its labels were calculated. The 5,000 rows are the smallest
SHA-256 keys under salt `r17-prospective-bisection-holdout-v1` among the 95,078
population rows outside the original selected cohort. The committed feature
hash is
`acbd9c534aac5e71090120326eb573275c8fbef0d25c19ca8172845eccc2cce1`.
The model and four comparison methods were frozen at the same time.

The complete atlas then performed another 195,600,000 exact square tests.
There are 94 certified-positive rows and 95 certified directions among 4,999
usable rows; one row is censored. On this genuinely unseen bisection target:

| method | ROC AUC | average precision | top-1% positives | top-1% enrichment |
| --- | ---: | ---: | ---: | ---: |
| frozen learned contrast | 0.7535 | 0.0696 | 6 | 6.38x |
| negative log parameter height | 0.6582 | 0.0408 | 3 | 3.19x |
| partial conductor quality | 0.6599 | 0.0304 | 1 | 1.06x |
| weakest-block Nagao | 0.5874 | 0.0264 | 1 | 1.06x |

At a five-percent budget the learned score retrieves 23 of 94 positives and
23 certified directions, an enrichment of 4.89 over random. Thus the model
does generalize to the response on which it was trained. Its failure on rank
28 is stronger evidence that known-bisection visibility is the wrong target
for extreme total rank, rather than evidence that the fit merely memorized
the stratified training cohort.

The deterministic evaluation is pinned in
[`r17_bisection_gain_ranker_prospective_holdout_v1.json`](../../artifacts/generated-results/elliptic-curves/r17_bisection_gain_ranker_prospective_holdout_v1.json).
The cohort, commitment, exact labels, and timing data remain under
`artifacts/local/elliptic-curves/`. Reproduce the sequence with:

```bash
python3 elliptic-curves/scripts/build_r17_prospective_holdout.py
python3 elliptic-curves/scripts/label_r17_training_bisections.py \
  --input artifacts/local/elliptic-curves/r17-prospective-holdout.jsonl \
  --output artifacts/local/elliptic-curves/r17-prospective-holdout-bisection-labels.jsonl \
  --summary artifacts/local/elliptic-curves/r17-prospective-holdout-bisection-labels-summary.json \
  --workers 4 --prime-bound 199
python3 elliptic-curves/scripts/evaluate_r17_bisection_ranker_prospective.py --check
```

This reaches the natural stopping boundary for the current atlas. More
bisection labels now have diminishing value for the total-rank objective. A
materially broader next round requires either equation-level higher-degree
multisections or a completed residual two-Selmer backend. The current
degree-three census is lattice-only, and the available residual-Selmer runs
remain incomplete, so neither can honestly supply those labels yet.

## Literature context

Mestre--Nagao sums are established rank-search heuristics, but their behavior
depends on the prime cutoff. Bujanovic--Kazalicki--Novak document oscillatory
``murmurations`` and examples where a smaller cutoff classifies rank better
than a larger one: [arXiv:2403.17626](https://arxiv.org/abs/2403.17626).
Kazalicki--Vlah use conductor and Frobenius traces as learned rank features and
compare them with Mestre--Nagao baselines:
[arXiv:2207.06699](https://arxiv.org/abs/2207.06699). Their later work learns
conductor-dependent, multi-cutoff Mestre--Nagao weights:
[arXiv:2506.07967](https://arxiv.org/abs/2506.07967). Those results support
using multiple disjoint prime blocks and learned combinations, but not treating
their outputs as rank bounds. Recent work on twist-class redundancy gives a
separate warning that arithmetic train/test splits must remove structurally
related duplicates rather than only row duplicates:
[arXiv:2605.14288](https://arxiv.org/abs/2605.14288).

The exact R17 equation, known fibres, score calibration, bisection boundary,
and residual-Selmer policy remain canonical in
[`../../elkies-k3/ELKIES_2026_R17_PAPER_IMPACT_2026-08-27.md`](../../elkies-k3/ELKIES_2026_R17_PAPER_IMPACT_2026-08-27.md).
