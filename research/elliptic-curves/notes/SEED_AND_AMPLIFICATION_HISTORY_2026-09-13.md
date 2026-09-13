# Seed production and amplification in retained search history

**The history analysis is complete. Cheap features help identify first seeds,
but no tested deep-amplifier rule passed the declared transfer gate. Searching
harder on the current stalled rank27/28 curves is retired as the default
rank32 strategy. No new point search or parent construction was launched.**

The best remaining feature signal is smaller displayed model complexity in
older score-selected R17 cohorts. That signal is worth distinguishing from
generic rank or a parent name, but its apparent deep-tail performance weakened
substantially when a fitted rule was transported to a later cohort.
This does not exclude better features, other detectors, or new parents.

## Dataset and scope

The [dataset](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/fibres.csv)
contains **4,482 cohort records in 33 retained cohorts**, with 20 recorded
parent/fibration identifiers. These are not 4,482 distinct elliptic curves or
20 independent surfaces. Repeated equations and blind catalogue rediscoveries
remain flagged. A common j-invariant is used conservatively to exclude possible
training leakage; it is not used to identify rational twists as the same curve.

Included are all 2,080 slots in the completed broad-rank run, the 60-fibre R17
seed/complement panel, 1,221 rows in 28 earlier MW16/R17 full-cohort exports,
922 rows in the last inherited high-rank-foundry ledger, 151 in the last
parent-foundry ledger, and the 48 completed determinant1092 record-scale
attempts. Foundry administrative stops and preparation failures are retained
separately from completed fixed exposures. Inherited ledgers are read once,
not added repeatedly across controller versions.

This is a census of the enumerated full-cohort records, not a claim that every
historical run directory supplies a complete prospective denominator. Early
success-only compact-R17 exports, public record controls, commissioning,
non-dispatched proposals and individually targeted follow-up experiments do
not enter the predictive denominator. The last category has different
selection and budgets; its results remain in the original notes, including
the [completed next-direction benchmark](NEXT_DIRECTION_RANK32_BENCHMARK_2026-09-12.md).
Here `final_rank` always means the endpoint of the named cohort, not the
largest bound ever obtained later on that equation.

The export retains exact parameters and models, initial/final certified
subgroup bounds, source and generic Gram bindings, frozen selection fields,
call counts, acquisition histories, continuation policy, censoring and timing
scope. Eight hundred ninety-two rows also join a full selection record by its
frozen protocol hash. Earlier and later score cutoffs are kept separately:
some historical protocols reuse the same field name for different cutoffs.
All 2,080 broad-run prime vectors are materialized by lookup in the frozen
residue tables; no new finite-field point counts are computed.

The [full generic Grams](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/lattices.json)
are Mordell–Weil data, not replacements for integral Neron–Severi frame lattices.
Parent identifiers and generic rank are descriptive strata, not one-hot
predictors of an unseen family. The broad X948 cohort has generic rank17
throughout and cannot itself test whether MW16 is preferable to MW17.

## Definitions that separate the mechanisms

A **seed** is a certified direction beyond the certified starting subgroup.
A **later gain** requires another originating point call to contribute a
direction. Admissions from one saved cloud are merged, including subsequent
arithmetic reconciliation of that cloud. Counting two points from one call
as two amplification steps would change the question.

A **deep amplifier** has an acquiring call whose directions from strictly
earlier calls already certify rank at least 23. A first cloud jumping directly
from 17 to 24 is not a deep amplification event. Originating-call attribution
does not assert that later-reconciled directions were already available to
the search scheduler at that time.

The mutually exclusive trajectory classes are no seed, one acquiring call
(seed only), multiple acquiring calls (amplifier), and deep amplifier, with
the last taking precedence. Endpoint rank18/19 is a separate bucket: it does
not establish a stall. A final lower bound is never an upper bound.

The conditional later-gain denominator contains seeded fibres with recorded
search after the first acquiring call. The broad run and R17 panel supply
that exposure for every seeded fibre. Other cohorts retain absent exposure
and administrative censoring explicitly. There are 4,338 complete acquisition
trajectories and 144 unclassified preparation records. Of the complete
trajectories, 1,005 have censored or administratively stopped exposure; their
positive certificates are retained but they are not silently pooled into
the completed-policy model fits.

## What the completed cohorts show

| Cohort | First seed | Later acquiring call, given seed | Deep amplifier |
|---|---:|---:|---:|
| Broad run | 839/2,080 | 494/839 | 5/2,080 |
| Earlier compact192 R17 | 156/192 | 118/156 | 36/192 |
| R17 seed/complement panel | 48/60 | 38/48 | 16/60 |
| X1092 record-scale fixed exposure | 0/48 | Undefined | 0/48 |

The broad run's four trajectory classes contain 1,241, 345, 489 and five
rows respectively. Its final rank histogram alone conceals the twelve
multi-direction seed clouds in the one-call class. In the R17 panel, the
rank24 fibre obtained entirely from its seed cloud is likewise excluded from
the deep-amplifier count.

The broad-run presentation table is:

| Presentation | Seeds | Later gains / seeds | Deep amplifiers |
|---|---:|---:|---:|
| X948 074d9 | 144/320 | 94/144 | 0 |
| X948 07ca9 | 169/320 | 104/169 | 1 |
| X948 08234 | 53/320 | 30/53 | 0 |
| X948 08f72 | 148/320 | 81/148 | 1 |
| X948 103b2 | 142/320 | 84/142 | 0 |
| X948 11952 | 135/320 | 74/135 | 3 |
| X1092 original | 48/80 | 27/48 | 0 |
| X1092 class1 | 0/80 | Undefined | 0 |

Thus “X1092 does not seed” would be incorrect even descriptively. Its original
presentation seeds often in this run but does not reach the deep endpoint.
X948 and X1092 use different initial detectors, so this is a comparison of
parent and detector policies, not a causal surface comparison. The small
three-event lead for 11952 is not a validated deep-tail selection rule.

The older cohorts and broad run have different parameter selection and
adaptive exposure. Their rates are observations under those policies, not
interchangeable estimates of an intrinsic family success probability.

## Whole-family validation

The [protocol](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/protocol.json)
fixed ridge logistic regression and a depth-two decision tree, without
hyperparameter search. All fits hold out a complete fibration; preprocessing
uses training rows only. Model inputs are parameter size/sign, displayed
coefficient/discriminant size and recorded cheap score summaries, with the
broad run's existing local-prime count. Post-search arithmetic profiles,
public ranks and timings are not predictor inputs.

The primary native-X948 run contains 1,920 rows. Logistic regression improves
held-out seed Brier score from **0.2481 to 0.1916** and conditional later-gain
score from **0.2426 to 0.2250**; lower is better. Its predicted top quarter
contains 331 seeds in 480 rows, versus 460 in1,440 elsewhere. This is not just
a ranked-versus-control effect: within the ranked arm the corresponding
seed counts are 276/384 versus403/1,152. The conditional-gain tree does not
improve the baseline.

For the five broad deep events, logistic regression puts three in its top
480 and two in the other1,440; the tree puts only one in its top quarter.
That event count does not support a reliable tail rule. AUC with five
positives is especially easy to overinterpret; the exact counts are primary.

There is a positive *within-cohort* signal in older R17 data. Leaving families
out of compact192, the logistic top quarter contains 22 deep amplifiers in 48
rows, versus 14 in144. A separately fitted model on the later R17 panel puts
nine in18 versus seven in42. Smaller model/discriminant coefficients recur
in the fitted rules. Those fits are not independent prospective deployments.

To test reuse of the rule, the same fixed models were then trained on
compact192 and applied to the later R17 panel, always omitting the test
fibration from training as well as possible equation aliases:

| Transported deep-gain model | Top quarter | Remainder | Brier score | Training-prevalence baseline |
|---|---:|---:|---:|---:|
| Logistic | 5/18 | 11/42 | 0.2540 | 0.2041 |
| Depth-two tree | 6/18 | 10/42 | 0.2024 | 0.2041 |

The logistic enrichment collapses; the tree's enrichment is only 1.4-fold
and its log loss is worse than the baseline. The two cohorts differ in
search policy, so this is a failed transfer under the tested policies,
not a theorem that coefficient size is irrelevant. No deep-gain model/family
test survives the reported multiplicity correction. Seed prediction is
more reusable than prediction of repeated or deep amplification.

All 772 model/fold records, including failed and one-label folds, are retained.
Per-study and global Holm adjustments are reported. Concatenated held-out
predictions share training data; no pooled Fisher p-value is presented as
independent replication. The forward cohort check and data-audit corrections
are disclosed in the [audit record](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/audit_corrections.json).
This is retrospective analysis with known historical outcomes, not a
preregistered new search experiment.

## Decision and reproduction

No tested rule earns a fresh deep-amplifier campaign. The analysis stops
without feature expansion, new scoring or more point search. This is not
“nothing predicts anything”: seed incidence is predictable to a useful
extent, and smaller models remain a descriptive visibility signal.
Neither result establishes improved probability of rank32.

The operational default is now candidate generation and parent/fibration
diversity, followed by explicitly bounded seed and amplification stages.
The next proposal must identify a new parent/fibration input or a materially
different candidate population; it must not quietly restart the stalled
27/28 follow-ups. New parent construction still needs its own mathematical
scope and budget. Arithmetic-constructor, carrier and fixed-word routes
remain parked. This analysis launches none of them.

Files:

- [Fibre CSV](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/fibres.csv) and [JSONL](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/fibres.jsonl).
- [Acquiring calls and timing bounds](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/gain_events.csv), [frozen prime vectors](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/prime_features.jsonl), [generic Grams](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/lattices.json).
- [Cohort/family tables](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/cohort_summary.csv), [held-out predictions](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/heldout_predictions.csv), [every fitted rule](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/heldout_folds.json).
- [Input bindings and coverage](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/extraction.json), [summary](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/summary.json), [separate accounting check](../../artifacts/generated-results/elliptic-curves/seed_amplification_history_v1/accounting_audit.json).

The accounting check validates 5,026 merged acquisition events, all 4,338
complete trajectories, 2,080 frozen prime-vector sums and 217 summary groups.
It checks retained packet/replay bindings and data accounting, not a fresh
elliptic arithmetic replay or a new independence proof. Primary timing uses
process-tree CPU and elapsed time separately; a gain without an exact timing
receipt receives its enclosing batch interval. Earlier point-only CPU sums
are labelled as such and never substituted for complete construction cost.

From the repository root, using the already installed NumPy/SciPy environment:

```sh
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python research/elliptic-curves/cas/analyze_seed_amplification_history.py analyze
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python research/elliptic-curves/cas/analyze_seed_amplification_history.py audit
```

`extract` reconstructs the exported dataset from retained files and performs
the frozen-selection joins. It reads discovery evidence but starts no
research calculation. No missing artifact is rebuilt. The original cohort
protocols and certificates remain unchanged.
