# Empirical source ranking from equation-construction attempts

Date: 2026-09-02.

## Outcome

The proposed feedback loop is useful, but the present data support a triage
model rather than an equation-success probability.  The foundry contains
2,134 MW1 candidates.  Of these, 1,342 have primitive root lattices and hence
complete exact pole, height, and component-correction features.  The other 792
have unresolved torsion/glue, so their pole and section lattice are not safe
model inputs.  Only six MW1 sources have documented equation-ansatz attempts,
covering eight exhaustive declared finite-field charts.  Five charts produce
the prescribed fibre configuration, one produces the full marked section, one
has a positive-dimensional marked local locus, and none has a certified
characteristic-zero equation.

Therefore the first implementation does three narrow things:

1. gives every attempted source a stable composite identity and records its
   exact lattice and finite-field features;
2. distinguishes fibre, marked-section, marked-locus, and equation gates;
3. ranks unattempted primitive MW1 rows by their exact-feature similarity to
   the six attempted sources.

The output is an **equation-precursor analogue score**, not the probability
that an equation exists or will be found.  This distinction is mandatory with
zero MW1 equation successes.

The curated input ledger is
[`data/lattice-foundry/source-equation-attempts-v1.json`](data/lattice-foundry/source-equation-attempts-v1.json).
The deterministic builder is
[`scripts/rank_lattice_foundry_empirical_sources.py`](scripts/rank_lattice_foundry_empirical_sources.py),
and its generated output is
[`../artifacts/generated-results/elkies-k3-lattice-foundry-empirical-source-ranking-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-empirical-source-ranking-v1.json).

## What is now recorded

For every MW1 candidate, the generated table records:

- MW rank, support count, exact minimum pole, generator height, and
  discriminant;
- the exact Shioda correction on each root component, its total and maximum,
  and the number of nonzero corrections;
- primitive versus nonprimitive root lattice;
- expected Kodaira discriminant orders and repeated orders;
- semistable all-`A` compatibility;
- the exact finite-field fibre and marked-section rates when an attempt exists;
- minimum residue-field extension degree among successful marked trials;
- whether a positive-dimensional marked locus appeared;
- the first failed gate and characteristic-zero equation status through the
  linked attempt records.

`source_id` is not globally unique: names such as `NS0011-S005` occur in more
than one prescribed-root shard.  The ledger key is therefore
`(source_artifact, source_id)`.  The reduced-Gram hash is retained as an
additional invariant.  Any future database that joins on `source_id` alone
will silently attach attempts to the wrong lattice row.

The per-component corrections are recomputed from the exact Schur projection
in the root-adapted Gram.  For a root block `R` and generator coupling `b`, the
stored correction is `b^T R^(-1)b`; their sum agrees with the difference
between the frame norm and the MW height.  Nonprimitive rows deliberately get
unknown correction and pole fields rather than values inferred from an
unsaturated root sublattice.

## The current empirical signal

The six attempted MW1 sources are strongly selected, not a random sample.
They were tried because the preceding lattice score already made them look
good.  The exact-chart outcomes are:

| source key | root type | pole | exact fibre rate | exact marked rate | marked locus | equation |
|---|---|---:|---:|---:|---|---|
| group-a / `NS0011-S005` | `A2+A6+A8` | 2 | 1/1 | 0/1 | no | open |
| group-c / `NS0007-S025` | `A1+A3+2A6` | 0 | 0/1 | 0/1 | no | open |
| group-a / `NS0034-S008` | `A2+A3+A4+A7` | 0 | 2/2 | 0/2 | no | open |
| group-a / `NS0043-S005` | `A2+A6+A8` | 0 | 1/1 | 0/1 | no | open |
| group-a / `NS0030-S001` | `2A1+A2+2A6` | 0 | 0/1 | 0/1 | no | open |
| group-b / `NS0048-S030` | `A1+A4+A6+D5` | 0 | 1/2 | 1/2 | yes | open |

Two conclusions already survive the small sample.

First, minimum pole zero is not sufficient.  NS0007 and NS0030 fail before a
squarefree fibre model in their exhaustive characteristic-five charts;
NS0034 and NS0043 reach fibre models but not the prescribed marking.  The
failure stage contains substantially more routing information than one binary
failure label.

Second, the only full marked signal is presently the mixed additive NS0048
profile.  Its `I1*` fibre lowers the short-Weierstrass degree bounds and it
produces smooth one-dimensional marked loci at two primes.  This is a real
reason to keep additive-fibre candidates in the next batch.  It is not yet
evidence that additive profiles are generally superior: the conclusion rests
on one attempted additive source.

The MW2 NS0028 pair is retained outside the fit.  Its fibre configurations are
abundant in characteristics five and seven, while the two required marked
sections never coexist on one model.  That is precisely why a row-level
"fibre success" bit cannot stand in for equation feasibility.

## The intentionally simple model

For each exhaustive declared chart, progress is coded ordinally as

```text
no fibre = 0
prescribed fibre = 1/3
complete marked section = 2/3
positive-dimensional marked locus = 1.
```

Candidates are compared on nine pre-attempt quantities: support count, pole,
log height, total and maximum component correction, number of nonzero
corrections, repeated-fibre excess, semistable compatibility, and log
discriminant.  Each difference is divided by its median absolute deviation in
the complete primitive MW1 population.  The score is a shrunk inverse-distance
average of the six observed progress values.  There are no fitted interaction
terms, no hidden classifier, and no claim that the resulting number is a
probability.

This analogue model is preferable to ordinary logistic regression today.
With one marked success, several predictors separate the data and maximum
likelihood coefficients would be unstable or infinite.  Firth's
bias-reducing penalty is the standard response to separation, while rare-event
logistic probabilities also need special care; neither method manufactures
information absent from six selected sources.  See
[Firth, *Bias reduction of maximum likelihood estimates*](https://doi.org/10.1093/biomet/80.1.27)
and
[King--Zeng, *Logistic Regression in Rare Events Data*](https://doi.org/10.1093/oxfordjournals.pan.a004868).

The selection problem is equally important.  Unattempted candidates are not
negative examples, and the attempted rows were not selected randomly from
the positive class.  Consequently the usual selected-completely-at-random
assumption behind basic positive--unlabelled correction is unavailable; see
[Elkan--Noto, *Learning Classifiers from Only Positive and Unlabeled Data*](https://cseweb.ucsd.edu/~elkan/posonly.pdf).

The leave-one-source-out mean absolute error of the present precursor score is
about `0.217`.  More importantly, when NS0048 is held out, the remaining five
all-semistable attempts cannot predict its marked-locus signal.  The generated
artifact records this failed validation explicitly.

## Which ten to try next

The unconstrained top ten are dominated by repeated reduced-Gram rows in
NS0048 and NS0043.  That is a poor experimental allocation: one hidden
arithmetic or normalization obstruction could remove most of the batch at
once.  The recommended list retains the positive NS0048 lead, excludes already
attempted sources with no marked hit, and takes at most one source per NS
class:

| priority | global analogue rank | source key | root type | supports | pole | height | correction sum | score |
|---:|---:|---|---|---:|---:|---:|---:|---:|
| 1 | 1 | group-b / `NS0048-S030` | `A1+A4+A6+D5` | 4 | 0 | `37/14` | `19/14` | 0.428 |
| 2 | 3 | group-b / `NS0043-S002` | `2A2+A8+D4` | 4 | 0 | `7/3` | `5/3` | 0.387 |
| 3 | 30 | group-a / `NS0046-S014` | `A5+A6+D5` | 3 | 1 | `109/28` | `59/28` | 0.312 |
| 4 | 37 | group-a / `NS0001-S001` | `2A2+A3+D9` | 4 | 2 | `79/12` | `17/12` | 0.310 |
| 5 | 39 | group-c / `NS0011-S013` | `A2+A6+A8` | 3 | 2 | `352/63` | `152/63` | 0.305 |
| 6 | 43 | group-a / `NS0042-S005` | `A2+A6+D8` | 3 | 4 | `71/7` | `13/7` | 0.294 |
| 7 | 45 | group-d / `NS0039-S001` | `A1+A10+D5` | 3 | 4 | `233/22` | `31/22` | 0.292 |
| 8 | 46 | group-d / `NS0034-S001` | `2A2+A4+D8` | 4 | 2 | `19/3` | `5/3` | 0.288 |
| 9 | 47 | group-b / `NS0003-S017` | `A2+A8+D6` | 3 | 3 | `47/6` | `13/6` | 0.288 |
| 10 | 48 | group-b / `NS0047-S004` | `A1+A4+A6+D5` | 4 | 1 | `47/10` | `13/10` | 0.285 |

The values in the last column are comparable ranking coordinates only.  In
particular, `0.428` does not mean a 42.8 percent equation probability.

This list teaches a useful qualitative lesson: after feeding back the actual
construction attempts, the immediate frontier moves away from an all-`A`,
support-first policy and toward mixed additive profiles.  The correct response
is a diversified test batch, not a wholesale rewrite of the lattice score.
The fifth row deliberately preserves a semistable control with the same exact
features as the attempted NS0011 source but a different reduced Gram.

## How to obtain a real success model

Every new run should append one attempt record before inspecting the outcome.
The record needs:

- composite source key and reduced-Gram hash;
- ansatz version, normalization, twist, prime, and extension degree;
- exhaustive versus bounded coverage and the exact exposure count;
- fibre-model count, raw-section count, fully marked-section count, and first
  failed gate;
- Jacobian rank/local dimension at every marked point;
- wall time, memory, and termination reason;
- characteristic-zero lift, reconstruction, substitution, and equation status.

The next batch should contain the ten rows above plus matched controls selected
before running: several all-`A` low-pole rows, several additive rows outside
the top score, and random primitive rows from score deciles.  This breaks the
current selection loop and makes failures interpretable.

After at least five independent marked-locus successes, fit two small gate
models rather than one monolith:

```text
P(prescribed fibre | static lattice features)
* P(marked section/locus | prescribed fibre, static features).
```

Use no more than three or four predeclared predictors initially, Firth or a
weakly regularized Bayesian logit, and leave one NS class out at a time.  Only
after multiple characteristic-zero equations succeed should the response be
renamed from "equation precursor" to "equation success".  Until then, a top
ten is an allocation rule for the next computation, not a mathematical
existence claim.

## Replay

Generate and check the deterministic ranking with

```bash
python3 elkies-k3/scripts/rank_lattice_foundry_empirical_sources.py
python3 elkies-k3/scripts/rank_lattice_foundry_empirical_sources.py --check
```

No Sage, Gröbner-basis, or new finite-field search is run by this analysis.
It consumes the pinned lattice ranking, pole audit, and exact/bounded attempt
artifacts already present in the repository.
