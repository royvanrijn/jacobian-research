# Curve302: index, nontrivial basins, and short-vector-only controls

## Question and existing evidence

This is a follow-up to the completed integral-core and short-vector experiments.
The reported enumeration has 1,288,441 sign-canonical primitive directions and
covers the 180 acquisitions in fourteen historical V3 runs. The short-vector
filtration matches the observed common cores at ranks one and three, but not
the whole chain. Static quotient height beat the tested adaptive predictors.
These findings are NOT overwritten or reclassified by this experiment.

The new question is narrower: after conditioning on the preference for short
vectors, is the observed common *integral* core still unusual? A positive
reference-ensemble match supports the simpler sampling explanation. A mismatch
leaves a residual to investigate at chart level. Neither outcome explains why
the original fibre has rank at least 31.

This patch's own tests are synthetic. It does not claim a new Curve302 result.
The input adapter targets the sealed output contracts in the preceding supplied
patch, not uninspected changes in the subsequent optimization commits. Changed
schemas or evidence fail closed; no guessed aliases or positional RREF parsing
are used.

## Run from the repository root

Use the same Python/SymPy environment as the previous exact lattice experiment.
The implementation is tested with SymPy 1.14.0; it does not require Sage or a
point-search backend.

```sh
python3 -m unittest discover -s research/elliptic-curves/tests \
  -p test_curve302_short_core_controls.py -v

python3 research/elliptic-curves/cas/run_curve302_short_core_controls.py run
python3 research/elliptic-curves/cas/run_curve302_short_core_controls.py check
```

Automatic input discovery searches local and committed generated results and
requires exactly one completed short-vector-core bundle. With multiple copies:

```sh
python3 research/elliptic-curves/cas/run_curve302_short_core_controls.py run \
  --source /absolute/path/to/completed-short-vector-core-results
```

That directory must include `REPORT.json`, `enumeration.json`, `filtration.json`,
`ranks-basins.json`, `primitive-directions.tsv`, and the original `plan.json`.
The prior plan identifies its closure-structure input. If that input moved,
pass `--structure /new/path/to/closure-structure-results`; the files must match
its recorded source hashes. Nothing in the original result directories changes.

A new output directory is used:

```text
research/artifacts/local/elliptic-curves/curve302-short-core-controls-v1/
  plan.json, manifest.json
  inputs/                       private copies of all consumed evidence
  code/                         private copies of both new Python files
  index.json                    first-29 Smith/index certificates
  basins.json                   corrected deficit-one summaries
  bands.json                    frozen per-run, per-event rank bands
  panels/00000.json, ...         accepted indices and every censored outcome
  sampling.json                 reference ensemble and complete denominators
  phases/                       bounded-worker logs, receipts, output seals
  SUMMARY.md, REPORT.json
```

`prepare` freezes without starting analysis. `resume` reuses only sealed, matching
stages and starts only untouched stages. `status` reports the stage records.
A failed/unsealed stage is never silently retried; use a new folder after review.
Failed verification attempts retain their logs under `checks/`. If the repository
controller changes later, invoke the frozen `code/run_curve302_short_core_controls.py`
with the original `--folder`; changed live controller code is refused.
The `check` command recomputes all THREE NEW analyses and compares every output,
including every panel and its random draws, byte-for-byte. It does NOT rerun the
previous 1,288,441-direction enumeration or prove completeness afresh.

## 1. Literal first-29 integral generation

Compute the Smith invariants of the original 29 integer row generators. Report
rank, index in saturation, and—when full rank—index in the ambient Z^14. Check the
Smith transformation identity and unimodularity, and cross-check the index with
an independent HNF determinant.

For every coordinate axis, give the smallest positive multiple in the generated
lattice and an explicit integer combination of the ORIGINAL 29 rows producing
it. Thus index one comes with witnesses, not merely a rank statement.

Norm ties are not concealed: compute a separate index after including the full
29th norm shell. The literal 29-row question and the shell-completed question
are separate outputs. The latter avoids allocating a quadratic-size Smith
transformation matrix when a shell is large.

## 2. Basins with automatic successes removed

Read the sealed per-prefix basin counts and independently recheck each reported
core deficit using the exact acquired prefix lattice. Deficit-zero cases are
counted only as already acquired. The primary aggregate contains deficit-one
cases only, and reports both numerator and denominator at:

- the actual next-acquisition norm;
- the full certified enumeration ceiling.

Both the pooled fraction and the equal-case mean are included, as they answer
different questions. A stratum without deficit-one cases has no estimate (`null`),
not a fabricated zero or one.

This reuses the certified source COUNTS; it does not recompute the large basin
census. Their original denominators include rationally independent extensions
which might be nonprimitive. They are saturated-span basin fractions, NOT the
primitive-only sampling law below and NOT V3 success probabilities.

## 3. Frozen conditional short-vector-only ensemble

Default: **128 panels of fourteen trajectories**. The same fourteen seed axes,
per-run number of acquisitions, and shorter rank-29 trajectory are preserved.
For each observed acquisition, use its static-height tied rank interval in the
complete vocabulary. Its WORST tied rank selects a predeclared decade:

```text
1..10, 11..100, 101..1000, 1001..10000, ...
```

Expand either boundary as necessary to include whole equal-norm shells; cap the
last band only at the already-sealed enumeration ceiling. No exact norm is used
as a singleton target and no band is retuned after seeing simulations. Very
small/forced bands and empty admissible sets remain limitations of this null.

At each step sample UNIFORMLY among the directions in the frozen band whose
addition is a primitive rank-one lattice extension. The selection code receives
no target-core features. Actual observed vector identities determine bands only;
they are not inserted into a changing candidate pool or favored during sampling.
Each run and panel has its own SHA-256-derived random stream.

### No implicit saturation

Maintain an exact surjective integer quotient map A with kernel the acquired
prefix S. A candidate v is accepted precisely when

```text
A*v != 0 and gcd(entries of A*v) == 1.
```

An exact unimodular Bezout reduction updates the quotient map. Its kernel is the
LITERAL new lattice S + Z*v. Nonprimitive candidates are rejected; saturation is
not used to create unseen directions. Tests compare this admission rule with
independent Smith calculations, including primitive ambient vectors which would
produce nonprimitive extensions.

### Censoring is part of the result

Uniform rejection sampling has a frozen 4,096-proposal cap per step. After 64
failed proposals, bands of at most 4,096 directions get an exact reservoir scan;
this preserves uniformity on the admissible subset. It distinguishes:

- `CENSORED_EMPTY_ADMISSIBLE_BAND`: exact scan proves no admissible extension;
- `CENSORED_DRAW_CAP`: the bounded large-band search did not produce one.

Neither widens the band, restarts the run, replaces the panel, nor counts as a
negative mathematical result. Other runs continue, allowing earlier prefix
comparisons to remain usable. Completed-only probabilities are NOT reported.
For every Boolean statistic all planned panels remain in the denominator;
unknown panels yield the finite empirical bounds
`[known_true/N, (known_true+unknown)/N]`. These are censoring bounds, not confidence
intervals or formal p-values. Excessive censoring means this particular frozen
null is inconclusive, not a reason to quietly change it.

### What is compared

Reconstruct the 194 observed primitive prefixes and independently verify each
reported common integral core by containment, rational annihilator rank, and
saturation. RREF serialization order is irrelevant.

For simulations measure all-run arrival of L2, L1, L4; common intersection ranks;
and containment/equality of the observed cores. Because every simulated prefix
is primitive, rank and containment establish INTEGRAL core equality. Dimensions
through 12 compare fourteen runs; dimensions 13 and 14 explicitly compare only
the thirteen corresponding survivors. Zero cores and the forced full-rank
endpoint are tagged trivial and excluded from primary comparisons.

The observed named cores and arrival dimensions came from these same histories.
The reference ensemble is therefore retrospective model checking, not an
independent discovery test, a fitted propagation theorem, or a prospective rank
predictor. A matching short-vector null should be allowed to simplify our story.

## Safety and performance

The vocabulary is loaded once per relevant worker into compact int64 STORAGE,
while all arithmetic uses unbounded Python integers/rationals. Coordinates too
large for storage fail explicitly; no wraparound occurs. Norms are streamed,
not duplicated as 1.3 million large Fraction objects. Only shell offsets and
needed observation positions are indexed. Witness norms are rechecked exactly;
other rows retain the preceding sealed enumeration's trust boundary.

Primitive extensions use small quotient maps, not a fresh 14-dimensional Smith
normal form for every proposal. Trials are separately recorded. The default
worker limits are one hour and 4 GiB address space, configured and frozen at
prepare with `--stage-seconds` / `--memory-mib`. The runner locks its output
folder, snapshots source and code, and kills only its own process group on a
supervisor timeout. Existing V3 campaigns are never touched.

No existing source, result, ledger status, or historical certificate is modified.
