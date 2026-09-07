# Fermigier corpus baseline evaluation

## Experimental question

Can cheap information available before a full point search put the two known
exceptional fibres near the front of a same-family search?  The evaluator
applies nine fixed historical/formulaic orderings to the
[`Fermigier labelled corpus`](FERMIGIER_LABELLED_CORPUS.md).  It fits no model
and never treats a censored fibre as rank zero.

The use of family-centred Frobenius residuals follows the same practical
motivation as Mestre--Nagao sieving: local point counts are cheap enough to
rank a large specialization population before expensive Mordell--Weil work.
[Fermigier's original E22 article](https://doi.org/10.4064/aa-82-4-359-363)
also describes the discovery method, not just the resulting curve.  A later
[worked high-rank-family search](https://pmc.ncbi.nlm.nih.gov/articles/PMC6604644/)
states the standard Mestre--Nagao score explicitly and uses it as a sieve.
These precedents justify the feature family, but do not make any finite-cutoff
score a theorem about rank.

## Reproducible commands

The active scorer first required one correctness repair: the integral family
has no good fibre at the fixed family-bad prime 5.  It now skips any prime
whose entire projective fibre locus is bad instead of aborting before the
candidate scan.

```sh
g++ -O3 -march=native -fopenmp -std=c++20 \
  -o /tmp/fermigier-score-sweep \
  elliptic-curves/ecsearch/fermigier_score_sweep.cpp

OMP_NUM_THREADS=32 /tmp/fermigier-score-sweep 100000 200 20000 \
  > artifacts/local/elliptic-curves/fermigier_family_residual_score_h100000_b200_top20000.tsv \
  2> artifacts/local/elliptic-curves/fermigier_family_residual_score_h100000_b200_top20000.log

.venv/bin/python elliptic-curves/scripts/build_fermigier_labelled_corpus.py
.venv/bin/python elliptic-curves/scripts/evaluate_fermigier_corpus_baselines.py
```

The first run scans all 12,195,252 primitive positive `u=a/b` with
`a<=100000`, `b<=200`.  It retains successively 8,000,000, 2,000,000,
500,000, 100,000, and 20,000 candidates at prime cutoffs 100, 200, 400,
1,000, and 2,000.  Each score uses only good-reduction traces and the exact
mean trace over the good projective fibres of this family.

## Exact finite-corpus result

| Ordering inside its materialized cohort | Cohort | E22 | Rank-20 fibre |
| --- | ---: | ---: | ---: |
| staged family-residual composite | 20,000 | 1 | 1,485 |
| cumulative residual `S0` | 20,000 | 3 | 147 |
| cumulative residual `S5` | 20,000 | 3 | 109 |
| last-window residual `S0` | 20,000 | 1 | 5,631 |
| last-window residual `S5` | 20,000 | 1 | 6,035 |
| recovered legacy global score | 5,000 | 2 | 423 |
| recovered legacy multibound composite | 5,000 | 4,704 | absent |
| recovered legacy broad ensemble | 6,071 | 82 | absent |
| recovered hot-neighbourhood score | 5,000 | absent | absent |

For the primary staged composite, the operational full-point-search budgets
are `1/12195252` for E22 and `1485/12195252` for the rank-20 fibre.  This is a
dramatic reduction relative to searching every fibre.  It is **not** an exact
rank under the final-cutoff score on the whole population: fibres discarded
at an earlier stage were not rescored to the final cutoff.

The cumulative components are the most consistent part of the panel.  They
also put the three censored fibres with legacy reported rank at least 17 at
median finalist position 5.  That concordance is descriptive only: the old
point searches were adaptively selected using related scores, and their
reported ranks remain uncertified.

## Leakage and interpretation

This is retrospective development evidence.  E22 was known when the staged
cap widths were chosen, and all nine orderings are being inspected together.
The rank-20 fibre was not used by the score formula, but it is not presented
as a prospective holdout because its existence was already known during this
research programme.

The comparison nevertheless answers the laboratory question at the intended
engineering level: some same-family cheap signals move both exact positives
far toward the front, while other plausible scores fail or omit them.  In
particular, adding the last-prime-window term improves E22 but degrades the
second positive by roughly a factor of ten relative to either cumulative
component.  That instability argues against selecting or fitting weights on
these two labels.

## Natural stopping point

There are only two certified positives in this family.  Additional supervised
fitting would mostly encode those two parameters and would have no credible
uncertainty estimate.  The useful next experiment must therefore be
prospective: freeze one cumulative family-centred rule and its stage caps,
apply it to a disjoint parameter block with no known labels, and open outcomes
only after candidate selection.  Until such outcomes exist, further weight
tuning on this corpus has sharply diminishing evidential value.
