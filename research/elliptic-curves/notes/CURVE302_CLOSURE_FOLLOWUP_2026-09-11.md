# Curve302: persistence, filtered subspaces, and held-out next moves

## Implementation status and provenance

**Execution update:** the actual sealed 14-run/180-acquisition bundle has now
completed all three follow-ups and deterministic recomputation. See the
[execution ledger](#actual-data-execution-ledger-2026-09-11) below. The following
implementation-session provenance is retained as history.

This additive runner is based on the source schemas inspected in the preceding
conversation and the user's report of the completed first closure experiment.
**The newly committed result files were not fetched in the implementation session:**
the GitHub connector was disabled, direct Git access failed, and the new files
were not available through the public web fallback. No claim is made to have
replayed or analyzed that new commit here. The delivered Git commit/patch was
prepared locally; publication to the repository's `main` was unavailable.

The runner consumes the actual exported data when run in the repository. It
requires the previous report's hashes to match its three files and the closure
report's landscape hash to match the full original cube. It never substitutes
fixtures for missing data. This is experiment infrastructure, not an update to
`research/MATH_STATUS.json` or a new mathematical-status claim.

## Run from the repository root

Python 3.10+, NumPy, psutil, and a POSIX host are required. The usual Sage Python
can also run this script. It deliberately needs neither a point-search binary
nor new Sage elliptic-curve arithmetic: the source trajectories were already
coordinate-verified by the preceding experiment.

```sh
python3 -m unittest discover -s research/elliptic-curves/tests \
  -p test_curve302_closure_followup.py -v

python3 research/elliptic-curves/cas/run_curve302_closure_followup.py run
```

The default input finder searches generated-results for a directory containing
`closure-laws.json`, `quotient-relations.json`, `trajectories.json`, and
`REPORT.json`. More than one distinct complete bundle is an error, not a reason
to select the newest or most successful result. Identical copies are harmless.
The original `curve302_exceptional_subgroup_landscape_v1.json` is a separate input.

An explicit invocation is available when the exported bundle is elsewhere:

```sh
python3 research/elliptic-curves/cas/run_curve302_closure_followup.py run \
  --results /absolute/path/to/previous-complete-bundle \
  --landscape /absolute/path/to/curve302_exceptional_subgroup_landscape_v1.json
```

The new default output is
`research/artifacts/local/elliptic-curves/curve302-closure-followup-v1/`.
`run` requires a fresh output folder. For a separate preflight, use `prepare`
then **`resume`**, not `prepare` followed by `run`.

```sh
python3 research/elliptic-curves/cas/run_curve302_closure_followup.py status
python3 research/elliptic-curves/cas/run_curve302_closure_followup.py check
```

`check` recomputes all three analyses from the frozen data and source, under the
same resource limits, then compares their complete JSON values. It preserves
check logs and receipts. This is deterministic recomputation, **not an independent
implementation or a new replay of the earlier elliptic-group identities**.

## Experiment 1 — exact threshold persistence

The previous two thresholds already closed the empty seed to the entire cube.
Their zero minimum-generator size was therefore inevitable, not evidence of an
interesting matroid or antimatroid. This runner replaces those snapshots with
complete growth functions for the empty seed and all fourteen singleton seeds.

For each start S, dynamic programming computes the minimum bottleneck needed to
reach every target, retaining concrete path witnesses. Because retained costs
are monotone under inclusion, the union of reachable sets is reachable: execute
one path, then another, whose costs cannot increase in the larger subgroup.
Consequently a direction belongs to the threshold closure exactly when its
minimum arrival threshold is at most T. These arrival values encode the entire
step function at **every distinct retained cost**, without naively rerunning a
closure calculation at each cost. A separate fixed-point algorithm checks both
sides of every growth breakpoint.

Outputs include each seed's reduction in the full-closure threshold, exact
closure-size advantage area up to the empty-seed threshold, subcritical closure
sizes, and arrival witnesses. The empty seed is always included. No closure
axiom or percolation interpretation is assumed. This is a weighted hypergraph
on one chosen set of fourteen axes, not an intrinsic closure on all rational
quotient directions.

## Experiment 2 — strict combinations and subspace intersections

For each actual acquired point, use its original integer quotient word, not just
its individually normalized primitive label. In every growing mod-2 subspace S,
compute a = dim(S intersect V_strict) and b = dim(S)-a. Verify a+b=dim(S).

Every basis vector of S intersect V_strict carries an exact F2 witness expressed
in the initial seed and previously acquired quotient vectors. This detects new
strict combinations even when the newly acquired vector is itself mixed.
Strict here means membership in the specified quotient kernel; a representative
may require adjustment by the generic subgroup. No new local or ideal-class
calculation is inferred from these combinations.

The rank-29 exception's terminal (a,b) is measured, **not hardcoded as (8,4)**.
The output records the first stage filling the local four-dimensional quotient,
the strict/local deficits at every endpoint, pairwise intersection dimensions,
and the common intersection across all runs at each attained dimension.

Pairwise intersections include their excess over max(0,2r-14), the unavoidable
ambient-dimension lower bound. Equality of full fourteen-dimensional endpoints
is explicitly tagged tautological. Equal subspaces mod 2 need not be equal
integral subgroups or equal rational subspaces; no such promotion is made.

## Experiment 3 — actual next acquisitions against frozen vocabularies

The old 0.996423 correlation compared two quantities built from the same height
form. It is informative about a metric approximation, but not independent
validation that geometry predicts V3 acquisitions. This experiment evaluates the
actual next acquired directions using only their current prefix for projection.

Two separate candidate sets are retained:

1. A fixed 18,760-vector set, independent of the acquisition outcomes: primitive
   signed axes; support two/three with coefficients in {-2,-1,1,2}; support four
   with coefficients in {-1,1}. Sign/gcd duplicates are removed.
2. Leave-one-run-out: primitive directions from the other thirteen runs plus
   the fourteen axes. A held-out target is never inserted into its own list.

Candidates already in the admitted mod-2 span are excluded. This reflects the
existing admission condition. A target absent from a candidate vocabulary is
reported as out of vocabulary, keeps its place in the all-event denominator,
receives percentile loss one, and scores no top-k hit. Coverage is reported
separately. This is not a test over the unknown set of all possible points.

Scores are projected residual norm, fractional unlock summed over the remaining
named axes, and static full quotient norm as a nonadaptive baseline. Positive
adaptive improvement means lower run-balanced loss than the static baseline.
QR in whitened coordinates avoids cancellation in repeated Schur subtraction.
Exact rational LDL first checks the supplied rounded-height form is positive
definite. Ill conditioning or numerically vanishing eligible residuals stop as
UNKNOWN, not as a silent candidate deletion.

The unlock score is necessarily tied in residual real dimension one. All ties
are reported as best/worst rank intervals; top-k counts use the conservative
worst rank. This prevents a final-stage degeneracy from becoming a fake top-one
success. The score depends on the named-axis frame, and the vocabulary compares
literal primitive vectors, not all equivalent rational lines modulo the prefix.

Seeded Monte Carlo draws from eligible candidates provide **conditional reference
intervals**, not population p-values. The fourteen starts share a single curve,
form, and much of their point vocabulary. Their acquisitions are not independent
trials. Neither these intervals nor leave-one-run-out removes the retrospective
use of the known M31 quotient metric. A high score also does not show that V3's
actual chart scheduler offered that candidate. A propagation theorem remains a
separate target.

## Operational isolation and failure behavior

No V3 call, point search, remote request, workflow change, or existing campaign
change is made. Only a fresh output directory is written. The controller copies
its own source and the five source files into that directory, hashes them,
and pins the Python/NumPy/psutil environment. Later repository commits do not
mutate these snapshots. A process lock prevents two controllers from using the
same output directory.

Each stage has its own log, start marker, supervisor receipt, and output seal.
Default caps are one thread, 3 GiB process-tree RSS, and 1,800 seconds per stage.
`--phase-seconds` and `--memory-mib` may be set at preparation; later resumes use
the sealed policy. A `STOP` file stops further phases and cancels the current
owned stage process group. No process belonging to another search is signalled.

A failed, censored, or unreceipted stage is preserved and cannot be silently
restarted. Investigate it and use a new folder for a changed experiment. A clean
resume skips sealed stages. Existing results are never overwritten by a fresh
`run`. Status labels seals as unchecked until `check` verifies them.

## Smoke-test evidence and experiment ledger

The implementation-session validation used **synthetic fixtures**, explicitly
labelled `SYNTHETIC_TEST_FIXTURE` in their plans/reports. The 14-run/180-acquisition
shape was chosen to exercise the source contract, not to reproduce real points.

- 37 regression tests: exact binary intersections and combination witnesses;
  minimax versus independent reachability on all small cutoffs; seeded advantage;
  zero-cost edges; exact LDL; projected norm versus rational Schur computation;
  rank ties; out-of-vocabulary accounting; deterministic randomness; corrupt
  inputs/plans/source; changed external originals versus immutable snapshots;
  lock conflict; STOP; timeout; failed-stage sealing; interrupted-run rejection.
- A full-size synthetic pipeline: all 16,384 states, fourteen runs and 180
  acquisitions, the full 18,760-vector policy and 1,000 reference draws.
- Deterministic full recomputation, plus clean resume with unchanged output
  hashes. These tests validate software behavior, not the proposed mechanism.

The actual committed result bundle remains **NOT RUN / NOT INSPECTED in the
implementation session**. This new suite should be recorded once with its
input hashes, outputs and interpretation, then retained rather than relaunched
under a new name. A new threshold or merely rerendering these summaries is not
another scientific experiment.

## Actual-data execution ledger (2026-09-11)

The supplied patch was applied with `git am` as `76b20e01`. All 37 supplied
tests passed. Real-input preparation then rejected `noncanonical source RREF`:
the original closure exporter sorts reduced binary rows numerically, while
the new reader required pivot order. The reader now accepts precisely these
two orderings of the same reduced basis and returns pivot order. It still
rejects unreduced, redundant, zero and noncanonical row sequences. Two added
regressions pass, bringing the specified test suite to **39 passing tests**.
The complete real source bundle also passes validation (14 runs, 180 gains).

The failed preparation is retained at
`artifacts/local/elliptic-curves/curve302-closure-followup-v1/`; no stage ran
there. The successful run uses the fresh sibling
`curve302-closure-followup-v1-rref-compat/`. From the repository root:

```sh
sage -python research/elliptic-curves/cas/run_curve302_closure_followup.py run \
  --results research/artifacts/local/elliptic-curves/curve302-closure-structure-v1 \
  --folder research/artifacts/local/elliptic-curves/curve302-closure-followup-v1-rref-compat
sage -python research/elliptic-curves/cas/run_curve302_closure_followup.py check \
  --folder research/artifacts/local/elliptic-curves/curve302-closure-followup-v1-rref-compat
```

These commands record the completed invocation; do not rerun `run` in that
existing folder. Defaults were retained: full vocabulary, 1,000 reference
draws, seed 3020911, one thread, 1,800 seconds and 3 GiB per stage. No search
was launched. The report status is `COMPLETE_THREE_RETROSPECTIVE_FOLLOWUPS`;
the saved check reports `PASS_DETERMINISTIC_RECOMPUTATION` for all three JSON
outputs, using the frozen source and inputs rather than repeating EC proofs.

| Binding | SHA256 |
| --- | --- |
| Frozen runner | `3e327f5384b2e5548c07ff3e50840313e68d132b86ed05f07142525eb493b43c` |
| Source REPORT | `a3e21d661c29dc30c4334a4cb3bd20ba55ad505d3380eacdbe2584cfd03ec358` |
| Landscape | `1f6f67d6575bb67842f93c8de2b3d7e14671704ad1118cbc67431cbcf861f38b` |
| Persistence output | `4858362cb78b15e50eff4e793cdb8656a36bc1929f0df438c29311cb33eb8792` |
| Strict-filtration output | `e3a73c468bfe91db48561b2834b9c80c0d68b18a7dc20984e62d7116b96c39b0` |
| Next-moves output | `e31fb3598cd715b475f8f206f76982e63f1a6814a6fa0c18beb24c7f78458062` |

The plan retains all individual source-output hashes, software versions and
absolute origins; the source report binds the preceding three outputs.

Persistence encodes all 5,930 distinct edge cutoffs. Seven singletons strictly
lower the full-closure threshold. The empty threshold is 30.94463075 in the
preceding atlas's units (numerator divided by 4,000,000); recovered-strict-02
and recovered-strict-03 attain 25.63565575, a 17.1564% reduction. Unchanged full
thresholds do not imply unchanged arrival curves.

All runs fill the local quotient by quotient dimension 4–7. The rank 29
exception ends at strict/local dimensions (8,4), so its two missing dimensions
are strict. Across all acquisitions there are 128 strict-dimension increments;
125 are supplied by individually mixed vectors, with exact F2 combination
witnesses retained. Thus the earlier count of only three individually strict
acquisitions does not measure strict-subspace growth.

Different intermediate subspaces contain a common core: at quotient dimension 6
all fourteen contain local-02; at 8 they contain local-01 and local-02; at 9 they
also contain local-04. Their common intersection at 12 has dimension 6. These
are statements about mod 2 subspaces, not equality of rational or integral
subgroups. Full-dimension endpoint convergence remains tautological.

Neither adaptive score improves on static quotient norm in either vocabulary:

| Vocabulary | Coverage | Static loss | Residual loss | Unlock loss | Static / residual / unlock top-1 hits |
| --- | --- | --- | --- | --- | --- |
| Fixed 18,760 vectors | 49/180 | 0.728571 | 0.732351 | 0.902795 | 23 / 0 / 0 |
| Leave-one-run-out | 70/180 | 0.632527 | 0.654158 | 0.749235 | 30 / 7 / 0 |

Losses are run-balanced mean percentile losses with absent targets penalized
by one; lower is better. Top-1 counts use conservative tie handling and all 180
events. These bounded literal-vector vocabularies have limited coverage.
This result supports early local completion and shared subspace content, but
does not support the tested adaptive next-move scores over static height.
It establishes no prospective predictor, propagation theorem or rank bound.
