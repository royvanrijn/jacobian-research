# Fermigier labelled-corpus protocol

## Purpose

This corpus recovers cheap score tables and later search outcomes for the
same Fermigier--Mestre family before any rank-jump model is fitted.  It is a
within-family retrieval dataset, not a comparison between unrelated record
curves.

The builder currently emits 517,922 distinct specializations.  Exactly two
rows are positive controls:

- Fermigier's E22 fibre, with certified rank at least 22 and exceptional
  quotient rank at least 10 over the generic rank 12;
- the `u=28917/20` near miss, with certified rank at least 20 and exceptional
  quotient rank at least 8.

The remaining 517,920 rows are **censored controls**.  They are not assigned
rank 12, quotient rank zero, or “no jump.”  A bounded point-search miss cannot
establish any of those statements.

## Coordinate and deduplication rule

Historical files use both Fermigier's adapter coordinate `u` and the literal
symmetric-shift coordinate `T=2*u`.  The implemented tuple model depends on
`T^2`.  Every row is therefore keyed by the reduced, nonnegative `T`; the
output also records `u`, `T^2`, and the primitive projective pair for `T`.
This exactly removes coordinate aliases and the known sign symmetry.  It is
not a global rational-isomorphism or quadratic-twist classification.

## Recovered inputs

The population starts with the 487,250-row local hot-neighbourhood table and
adds five pre-point-search feature sources:

- the 5,000-row historical global-score tranche;
- the 5,000-row hot-neighbourhood score tranche;
- the 5,000-row multibound score table;
- the 6,071-row broad ensemble.

The current explicit family-residual scorer adds a fifth table: the top
20,000 operational finalists from all 12,195,252 primitive `u=a/b` with
`1<=a<=100000` and `1<=b<=200`.  Its C++ implementation and run log are
hashed with the table.  The scorer removes the exact mean trace over the good
projective fibre locus at each prime; fixed family-bad primes contribute no
ordinary trace score.

It then attaches structured observations from fifteen archived/active
Fermigier experiments and the Fermigier lane of the conductor-first anchor
pilot.  Four legacy logs add bounded point-search observations.  The summary
records the SHA-256 of every input and both the pre-deduplication record count
and post-deduplication row count for every source.

The separate 487,262-row neighbourhood-provenance table is also hashed and
checked to cover exactly those 487,250 distinct parameters.  It can be joined
by canonical `u` when seed/offset lineage is needed; duplicated provenance
rows do not duplicate corpus fibres.

The local TSVs and logs survived without their historical generators.  Their
hashes make this recovery repeatable in this workspace, but that missing
lineage is an explicit limitation for the four older feature tables.  The new
family-residual table has complete code/run lineage.  The archived global search's complete
60,815,684-fibre box remains a compact population formula and exact ranking
replay; it is not expanded into millions of invented negative labels.

## Output and leakage boundary

Run:

```sh
.venv/bin/python elliptic-curves/scripts/build_fermigier_labelled_corpus.py
```

The command writes deterministic gzip-compressed JSONL to
`artifacts/local/elliptic-curves/fermigier-labelled-corpus-v1.jsonl.gz` and a
hashed summary beside it.  These are raw recovered data and therefore remain
under the ignored local-artifact tree.

Every non-positive row receives a deterministic train, validation, or
internal-test split.  The two positives are quarantined in
`positive_holdout`; their parameters and labels must not be used to invent or
tune a score.  Later evaluation may use them once to measure retrieval rank,
top fraction, recall at fixed budgets, enrichment, and comparisons with the
historical score columns.

Legacy `rank>=` log fields and structured stable numerical ranks are retained
under names ending in `_uncertified` or as numerical observations.  Only the
attached exact certificate changes `label.state` to `certified_positive`.

## Current control coverage

The recovered corpus contains:

- 927 rows with a legacy reported rank field from a completed bounded search;
- 704 rows whose completed legacy search reported zero quartic points;
- 223 rows whose completed legacy search reported at least one quartic point;
- 20 rows with a structured stable numerical height rank, including the
  rank-15 and rank-16 triage fibres;
- 500 logged failed attempts, retained as failures rather than negatives.

These counts are corpus metadata, not mathematical rank claims.  Rebuilding
the corpus is the gate before fitting or selecting any new ranking rule.
