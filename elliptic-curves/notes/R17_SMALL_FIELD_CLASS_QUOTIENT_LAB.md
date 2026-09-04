# R17 small-field prospective class-quotient laboratory

## Current status

Phase 0 is frozen.  The laboratory contains exactly 100 rank-blind ordinary
fibres of the published R17/MW17 family.  No class quotient and no detector
outcome has yet been promoted: every feature cell is null and every point-search
cell is sealed in the cohort artifact.

The commitment is
[`../../artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-cohort-v1.json`](../../artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-cohort-v1.json),
with candidate-list SHA-256

```text
47cb093c0cbf4803e2cb6c176a45579f241ecf574a48c7c742dce2f053c4d95e
```

This is an experiment protocol, not a mathematical-status result.  In
particular, `dim Q_t` remains `UNKNOWN` until every row has an unconditional
BNF certificate and the Phase-1 ledger freezes.

## Frozen question

For the completed-square irreducible cubic field `K_t`, let `S_t` contain all
prime ideals above rational primes dividing `2 Delta_min(E_t)`.  For each of
the seventeen specialized generic sections, the exact localized Kummer
half-ideal gives a class `c(G_i)` modulo the `S_t`-classes.  The pre-search
feature is

\[
 Q_t=
 \frac{\operatorname{Cl}(K_t)}
 {2\operatorname{Cl}(K_t)+
  \langle S_t,c(G_1),\ldots,c(G_{17})\rangle}.
\]

The frozen primary question is whether larger `dim Q_t` predicts a larger
integer gain from the later fixed Stage-A half-lattice detector.  Total
`dim Cl(K_t)/2Cl(K_t)` is retained as the predeclared negative-control
predictor.  The primary statistic is Kendall tau-b; its one-sided randomization
test uses 100,000 deterministic permutations within cubic-signature by
field-discriminant-quartile blocks.  There is no trained threshold.

The endpoint is deliberately “future exactly certified detector-visible
quotient gain,” not true rank jump.  Every counted direction must satisfy the
curve equation and a full finite-reduction independence certificate.  A
bounded miss proves neither rank 17 nor absence of an escape.

## Rank-blind cohort

The universe is every reduced finite `a/b` in the published R17 coordinate
with

```text
b > 0,  gcd(|a|,b)=1,  max(|a|,b) <= 24.
```

There are 719 such parameters.  Structural filtering requires a nonsingular
fibre and irreducible completed-square 2-division cubic and removes duplicate
exact rational `j`-invariants.  The first 100 rows are then taken in increasing
order of the exact explicit cubic-order discriminant

```text
disc(z^3 + 16*A_h(a,b)*z + 64*B_h(a,b)),
```

with projective height, denominator, and numerator as fixed tie breakers.  No
rank, exceptional point, Nagao, local-cylinder, Selmer, or previous search-hit
label is loaded.  The selected order-discriminant sizes range from 286 to 353
bits.  Maximal-order field discriminants are Phase-1 outputs rather than
selection-time assumptions.

“Ordinary” here means nonsingular with irreducible 2-division cubic.  It does
not make the stronger, presently unauditable claim that no one has ever
computed any selected fibre.

## Enforced phase boundary

The execution order is:

1. freeze the rank-blind cohort;
2. compute the complete class group, `S_t` image, generic-MW17 Kummer image,
   units, local data, and `dim Q_t` for all 100 rows;
3. freeze the detector protocol against the whole-file feature hash;
4. run the detector from a redacted manifest which contains parameters but no
   feature values;
5. join the two ledgers once and apply the predeclared analysis.

The feature worker first transports the completed-square root through an exact
`polredabs` field isomorphism.  It then requires both `nfcertify` and
`bnfcertify == 1`.  Class coordinates of every `S_t` prime and every generic
half-ideal are reduced only on even class-group cyclic factors; exact binary
rank gives `dim Q_t`.  A timeout or backend failure leaves `dim Q_t=null` and
keeps all point search sealed.  A GRH-only BNF can never unlock the detector.

The first-row implementation probe reached BNF certification only after a
four-gigabyte PARI stack allocation; its large Zimmert bound made certification
long enough that the probe was manually stopped.  Therefore Phase 1 is a
checkpointed campaign with a one-hour per-row envelope, not a claimed cheap
replay.  This is precisely why the frozen artifact does not pretend that the
quotient values already exist.

## Recovered MW16 expansion boundary

The recovered curve-398 `A1`/MW16 pencil is kept as an explicit expansion
lane, but is not included in v1.  In its current exact `lambda` gauge, even the
`lambda=0` maximal-order reduction exceeded a 30-second feasibility probe.
That is an operational observation, not an exclusion theorem.  The family may
be added only after a rank-blind rational reparameterization demonstrates a
complete unconditional BNF and generic-MW16 class-image computation under the
same resource policy.  Mixing an unproved “small-field” MW16 lane into this
commitment would make the missingness family-dependent from the start.

## Reproduction and execution

Freeze or replay Phase 0:

```bash
sage -python \
  elkies-k3/scripts/build_r17_small_field_class_quotient_cohort.sage --check
```

Run Phase 1 in sixteen deterministic shards.  Each worker may reserve four
gigabytes, so on a 32-GiB host run no more than four shards concurrently:

```bash
sage -python \
  elkies-k3/scripts/run_r17_small_field_class_quotient_features.sage \
  --chunk-index I --chunk-count 16

sage -python \
  elkies-k3/scripts/run_r17_small_field_class_quotient_features.sage \
  --merge --chunk-count 16
```

The merge emits `FROZEN_COMPLETE_UNCONDITIONAL_PRE_SEARCH_FEATURES` only when
all 100 certified rows are present.  Then, and only then:

```bash
sage -python \
  elkies-k3/scripts/freeze_r17_small_field_class_quotient_detector_protocol.sage

sage -python \
  elkies-k3/scripts/run_r17_small_field_class_quotient_detector.sage \
  --chunk-index I --chunk-count 16

sage -python \
  elkies-k3/scripts/run_r17_small_field_class_quotient_detector.sage \
  --merge --chunk-count 16

.venv/bin/python \
  elkies-k3/scripts/analyze_r17_small_field_class_quotient_experiment.py
```

The confirmatory statistic is set to null if even one scheduled Stage-A row is
censored.  Complete-case summaries may still be emitted, but only as
exploratory diagnostics and never by treating a timeout as zero gain.

