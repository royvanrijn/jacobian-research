# Retained rank-26 source gaps and completed exposure

Two already-retained curves completed **98 new point-search boxes with zero
rank gains**. Their full clouds contain 477 rational points and still certify
rank at least26. This is a completed visibility experiment on the existing
X948 parent, not a parent expansion or a new parameter scan.

Authority: `EC-RETAINED26-SOURCE-GAP-EXPOSURE-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json). The
[report](../../artifacts/generated-results/elliptic-curves/retained26_gap_trial_report_v1.json)
binds source evidence, frozen selection, completed exposure and independent proofs.

## Source audit and fixed selection

InventoryV22 contains18 rank-26 curves. The audit reads1,642 `result.json`
files at depths1,2,3 under the local elliptic-curve artifact directory,
excluding paths containing `portable`, `workspace` or `standalone`.
There are53 exact family-and-parameter matching files. Eleven curves already
have later own26 centre records and receive no new exposure. Historical
RUNNING labels in those records are not evidence of live processes.

The remaining seven have completed discovery records whose centres use only
their first17 generic basis directions. Exact seed-prefix and centre checks
are retained in the
[source audit](../../artifacts/generated-results/elliptic-curves/retained26_source_gap_v1.json).
This is affirmative evidence about those source records and a bounded layout
review. It is not repository-wide absence across aliases, isomorphic models,
deeper directories or archives.

| Stable ID suffix | Family | Parameter | Completed source boxes | New boxes |
|---|---|---|---:|---:|
| 73 | 074d9 | 1175/1901 | 43 | 0 |
| 63 | 074d9 | 2588/1603 | 43 | 0 |
| 42 | 074d9 | 808/2259 | 43 | 0 |
| 74 | 07ca9 | 613/828 | 43 | 49 |
| 49 | 103b2 | 2773/962 | 43 | 0 |
| 50 | 11952 | -1173/127 | 49 | 49 |
| 75 | 11952 | -1376/455 | 49 | 0 |

IDs have prefix `new-20260906-`. Selection was fixed before scoring: decreasing
training score, then coefficient bits and ID; take the first two distinct
fibration labels. The [scores](../../artifacts/generated-results/elliptic-curves/retained26_gap_scores_v1.json)
use all3,510 primes5..32749, for24,570 traces and336 independent direct
character-sum checks. Corrected local good-prime reduction is enabled;
none of these seven models requires a local fourth/sixth-power reduction.
The integral score model differs from the source short model by scale12,
invertible at every scoring prime. Validation primes65537..131071 are neither
computed nor used. IDs74 and50 are selected; distinct fibration labels do
not mean distinct parents. This is not a randomized score-efficacy comparison.

## Equal, completed point exposure

Both map sets precede all point searches. Each uses2,048 distinct deterministic
SHA256 parity masks on the certified26-point seed, restricted to nonzero
coefficients above the first17 directions. The fixed sample domain is
`full11952-specialized-followup-v1`. Canonical heights use384-bit precision,
rounded at scale10^6, followed by unimodular LLL and numerical CVP. The49
largest computed norms supply factor-free Gauss/hyperellred charts.
The finite-mod2 independent26 seed proves these classes lie outside the old
generic17 image in E(Q)/2E(Q). No covering or global CVP optimality is claimed.

Each curve receives49 boxes, height125000 and ten seconds per box, with no
rank stop, refill or following wave. Stage limits are180 seconds for maps,
600 seconds for the point worker and300 seconds for history replay, one
worker and2GiB. All98 boxes complete. Exact replay checks4,096 parity/norm
transports,98 rational maps and point provenance.

| ID suffix | Retained cloud | Initial lower bound | Final lower bound | Gains |
|---|---:|---:|---:|---:|
| 74 | 309 | 26 | 26 | 0 |
| 50 | 168 | 26 | 26 | 0 |

Both full clouds have finite column rank26 modulo2,3,5. The
[standalone verifier](../cas/verify_retained26_gap_rank.sage), copied with the
two mod2 certificates into a fresh directory, independently enumerates complete
finite groups, checks all477 points, and certifies independence and torsion
exclusion without repository imports. Its transcript is
`PASS independent retained26 full clouds [(309, 26), (168, 26)]`.

Total supervised time including seven-curve scoring, maps, workers, histories,
cloud proofs, geometry and independent replay is139.02700246241875 seconds.
This counts constituent cloud jobs once, excluding the enclosing wrapper's
11.661024059052579 seconds to avoid double counting. Layout discovery,
report-generation CPU and total session time are outside that accounting.

There is no inventory change, no exact-rank upper bound, and no proof of
absence of additional rational points. Five source gaps remain unscheduled.
The null result does not automatically authorize another wave.

## Replay

```sh
python3 elliptic-curves/cas/audit_retained26_source_gap.py --check
python3 elliptic-curves/cas/score_retained26_gaps.py check
python3 elliptic-curves/cas/report_retained26_gap_trial.py --check
sage -python elliptic-curves/cas/verify_retained26_gap_rank.sage --input \
  artifacts/generated-results/elliptic-curves/retained26_gap_new_20260906_74_mod2_v1.json \
  artifacts/generated-results/elliptic-curves/retained26_gap_new_20260906_50_mod2_v1.json
```

Source audit check mode replays frozen affirmative bindings; it does not rescan
changing directories. Protocols, maps, results, checkpoints and transcripts
are under `artifacts/local/elliptic-curves/retained26-gap-scores-v1` and
`retained26-gap-trial-v1`. Frozen sources and completed searches are preserved.
