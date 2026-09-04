# Rank-jump laboratory

## Question and unit of comparison

The laboratory asks one deliberately narrow question:

> Using only information available before a full point search, how highly does
> a ranking rule place already-known exceptional specializations among
> parameters of the same family?

The comparison unit is a rational parameter in one fixed, normalized family
coordinate. Curves from unrelated constructions are not pooled. The primary
outcomes are the exact population ranks of known positive controls, recall at
fixed search-budget fractions, and enrichment over a random ordering of those
same known controls. Classification accuracy is not reported.

The machine-readable inventory is
[`../data/rank_jump_laboratory_v1.json`](../data/rank_jump_laboratory_v1.json).
Its checker verifies the SHA-256 of every certificate and ranking artifact
before calculating any metric.

## Labels are asymmetric

A positive enters the active label set only when both of the following are
available:

1. the arithmetic generic rank of the family is proved; and
2. independent specialized points prove a larger rank lower bound.

The difference is then a certified lower bound for the exceptional quotient
rank. A published rank value without a local independence replay, a numerical
rank estimate, or a raw point list is not sufficient.

Controls have weaker semantics. An ordinary parameter, a high-score fibre on
which a bounded search found nothing, or a low numerical-rank fibre is a
*censored retrieval control*. It is not a theorem that the fibre has no
additional rational points. This is why the laboratory measures retrieval of
known positives and does not use accuracy, specificity, ROC AUC, or a learned
binary classifier that treats every search miss as rank zero.

## Feature boundary

Every ranking run declares the latest stage it uses:

- `parameter_only`: height, denominator, congruence class, root shape, and
  other data available from the normalized parameter;
- `local_arithmetic`: good-prime traces, disjoint-block Nagao statistics,
  root number, discriminant/conductor proxies, and other local data computed
  without searching for rational points;
- `cover_incidence`: counts or local-solubility summaries of a predeclared
  family cover, still before a full point search.

Point-search output and rank labels are forbidden as score inputs. If a method
is tuned on the known controls, the run is labelled `development`. A
`held_out` result requires the score and all hyperparameters to have been
frozen before the held-out family, parameter block, or prime band is opened.
Multiple prime ensembles inside one chosen development score test stability;
they do not by themselves turn the run into a held-out experiment.

## Metrics

For a population of size \(N\), a known positive at rank \(r\) has population
fraction \(r/N\). The reciprocal \(N/r\) is reported as a descriptive
search-space reduction factor. For each predeclared budget fraction \(q\), the
laboratory also reports:

- the number and recall of known positives in the first \(\lceil qN\rceil\)
  candidates;
- the expected number of those same known positives under a uniformly random
  ordering; and
- observed hits divided by that expectation.

These are finite-corpus retrieval statistics, not estimates of the prevalence
of unknown rank jumps. Results are also stratified by the certified minimum
exceptional quotient gain.

## Initial inventory

| Family | Exact admitted positives | Ranking corpus | Current state |
| --- | ---: | --- | --- |
| published R17, compact `t` | gains at least 8, 9, 10, 11 | complete `H<=10000` Nagao population plus sampled 100,000-row bisection-gain replay | development benchmark and first held-out replay active |
| Fermigier/Mestre rank 12 | E22 gain at least 10; rank-20 near miss gain at least 8 | complete normalized 60,815,684-parameter historical global box | retrospective development replay active |
| Nagao section 7 rank 12 | `T=5081/47`, gain at least 8 | complete 18,244,819-parameter historical global box | retrospective development replay active |
| Mestre split-infinity `d`-square families | rank-17/rank-19 fibres exist | no admitted quotient label | exact generic-rank or saturated-subgroup audit required |
| E29 | none | none | excluded until lineage and parameter normalization are exact |

The exclusions are intentional. In particular, “generic rank at least 13” is
not enough to turn a rank-at-least-17 fibre into a certified four-direction
jump beyond the full generic Mordell--Weil group.

## First measured result

The current R17 development rule is the minimum standardized Nagao signal
across three pairwise-disjoint prime ensembles. The complete bounded scan has
121,589,944 primitive projective parameters of height at most 10,000. The
certified rank-at-least-25, -26, -27, and -28 controls occur at ranks

```text
54,624; 593,936; 422,873; 55,387.
```

Thus all four occur in the top `1/200` of the population, where a random
ordering expects about `0.02` of these four controls: known-control enrichment
is about 200-fold. Two occur in the top `1/2000`, where enrichment is about
1,000-fold. The worst control is at population fraction
`593936/121589944 = 0.004884746...`; the two strongest placements reduce the
prefix to inspect by factors of about 2,226 and 2,195.

This answers the laboratory question positively for one family and one
development rule. Because these four fibres were calibration anchors, it is
not yet evidence of out-of-family generalization.

## R17 bisection-gain replay: held-out failure

The separate deliberate-data protocol evaluates the complete preexisting
39,120-cover bisection atlas on 4,922 stratified parameters, attaches exact
finite-quotient gain lower bounds, and fits a fixed train-only class-mean
contrast on Level 0--2 features. All four published controls are absent until
the model serialization hash is frozen. In the sampled 100,000-row population,
the learned score places rank 25--28 at positions

```text
580; 232; 784; 24,210.
```

It therefore retrieves rank 27 but misses rank 28 at a one-percent budget.
The unchanged weakest-block Nagao comparator ranks the same controls at
`49, 509, 344, 49` and retrieves all four. This is a useful held-out failure:
known-bisection visibility overweights small parameter height and local
conductor quality, while only one of the rank-28 fibre's eleven known
exceptional directions is visible to that atlas. The exact protocol, labels,
and no-retuning boundary are in
[`R17_TRAINING_DATA_PROTOCOL.md`](R17_TRAINING_DATA_PROTOCOL.md).

A separately committed 5,000-row random holdout confirms that distinction.
The frozen learned score has ROC AUC `0.7535` and 6.38-fold top-one-percent
enrichment for exact bisection gain, compared with AUC `0.5874` and 1.06-fold
enrichment for Nagao. It learned the narrower cover-visibility response, but
that response does not recover the extreme-rank control. Additional labels
from the same bisection atlas therefore have diminishing value for this
laboratory's total-rank objective.

<!-- status-consumer: EC-K3-R17-TRAINING-EXACT-ARITHMETIC-GROUP-GATE 427bf822e774c81e -->

Further reuse of that learned contrast is now mechanically gated by the exact
R17 arithmetic-group audit.  Across all 100,000 development rows there are no
repeated rational `j`-classes; the labelled selection and prospective holdout
have no twist-class overlap, and the quarantined controls have no twist-class
match in the development population.  The laboratory registry refuses the
learned extractor if the gate is absent, stale, or authorizes another fitted
score.  The historical score and opened evaluations are unchanged.

## Fermigier replay: a negative retrieval result

The normalized Fermigier global box contains 60,815,684 primitive parameters.
The predeclared discovery rank/composite and disjoint held rank/composite
orderings place E22 between positions 2,755,127 and 3,556,250; they place the
rank-20 near miss between 3,070,200 and 12,535,549. Every ordering therefore
has zero recall through a 100,000-candidate budget. The full result and exact
claim boundary are in
[`FERMIGIER_RANK_JUMP_REPLAY.md`](FERMIGIER_RANK_JUMP_REPLAY.md).

This is not a failed theorem and does not turn the remaining population into
negative labels. It is evidence that the historical local score is not a
useful escape-from-generic-lattice retriever in this family. The next
informative evaluation is the normalized Nagao population or a genuinely
separate quotient-aware rule; the Fermigier score should not launch candidate
search or be retuned on these two controls.

The matching recovered corpus is described in
[`FERMIGIER_LABELLED_CORPUS.md`](FERMIGIER_LABELLED_CORPUS.md).  It joins
517,922 canonical specializations, five cheap-feature tables, 927
completed legacy search outcomes, 20 structured numerical-rank outcomes, and
the two exact positives without assigning a negative rank label to any other
fibre.  The positives are quarantined from its deterministic development
splits.  This corpus is now the data gate before any replacement Fermigier
score is fitted.

The subsequent no-fit baseline panel is reported in
[`FERMIGIER_BASELINE_EVALUATION.md`](FERMIGIER_BASELINE_EVALUATION.md).  On the
20,000 finalists of a complete 12,195,252-parameter staged scan, its explicit
family-residual composite puts E22 first and the rank-20 fibre 1,485th.  The
two cumulative residual components put them at `3/147` and `3/109`, whereas
the late-prime window components lose the rank-20 fibre to positions 5,631
and 6,035.  This is retrospective development evidence; E22 influenced the
stage widths, and no discarded fibre received the final-cutoff score.

## Nagao section-7 replay: band instability

The normalized historical box has 18,244,819 positive primitive parameters.
Its frozen training score places `T=5081/47` at position 9,041,935, while the
disjoint validation band places it at 755,065. The latter is a substantial
improvement but still misses a one-percent budget; the former is essentially
median. This is a retrospective external-positive replay, not a prospective
holdout, and all other fibres remain censored. The associated quotient
fingerprint has free rank 8, Smith index 2048, tensor dimensions `19,8,8` over
`F_2,F_3,F_5`, and bounded degree-two visibility spanning all eight free
directions. See
[`NAGAO_SECTION7_RANK_JUMP_REPLAY.md`](NAGAO_SECTION7_RANK_JUMP_REPLAY.md).

## Replay and extension

Run the active inventory and metric replay with:

```sh
.venv/bin/python elliptic-curves/scripts/run_rank_jump_laboratory.py
```

To retain the full JSON report locally:

```sh
.venv/bin/python elliptic-curves/scripts/run_rank_jump_laboratory.py \
  --output artifacts/local/elliptic-curves/rank_jump_laboratory_v1_result.json
```

Adding a family requires a canonical parameter coordinate, proof-bearing
positive labels, censored-control provenance, a frozen pre-point-search score,
and exact population ranks for every admitted positive. Large candidate
tables remain local; a compact ranking artifact should retain the population
definition, feature specification, positive ranks, hashes, and claim boundary.
