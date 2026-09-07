# Full-score selection on untouched 11952 outer slices

**The fixed parameter scan and score-component replay pass. All 4,096 scalar scores, exact
selection replay and disjoint 64-finalist validation pass. All 3,136 point boxes are complete; all local independent
point proofs pass, with bounds 18 on two fibres and 17 on the remaining 62.
All 194 isolated point-only replay stages pass. No near-record result is claimed.**

The previous outer-annulus trial retained 4,096 short-score candidates per
slice before using the longer projective trace cache. An exact benchmark on
its first negative slice instead scores all 12,690,811 primitive addresses
through prime 32749 before keeping 512. Of those 512 leaders, **366 were absent
from the previous short-score survivors**. Complete signed small-frame tests,
heap ordering and independently compiled cached-reader checks pass. The full
slice takes 20.75 seconds, below the frozen 45-second cost gate. This measures
retention loss for the stated score, not rank density or exceptional-point
sensitivity. The benchmark starts no point search.

## New parameter population

Eight new denominator residues modulo 16384 are chosen by SHA256, one in
each quarter per sign with balanced parity, excluding the earlier eight
outer slices. Every address satisfies

\[131072 < H=\max(|n|,d)\le524288,\qquad \gcd(n,d)=1,\quad d>0.\]

The **76,268,467 primitive addresses** are disjoint from both the previously
complete height-131072 square and the earlier outer trial. All 3,510 primes
through 32749 contribute before heap admission. Each slice retains 4,096,
giving 32,768 candidates. All eight scans, complete signed frame regressions
and independent replay of both cached score components pass. This is a
stratified annulus sample, not complete coverage.

## Frozen new-fibre trial

After the earlier MW16 and 11952 cohorts finish their exact proofs, freeze
1,101 previous equation entries and the pinned 593-equation catalogue as
exclusions. Score 4,096 distinct new equations with fresh scalar traces
through 65521; require agreement with the first cached extension and two
independent character sums per equation. The first eight calls have a fixed
runtime gate, followed by two workers and 16-row checkpoints. Freeze 64
finalists only after exact score replay. Disjoint primes 65537 through 131071
provide validation and cannot change selection.

The point trial uses only the 17 transported generic sections, 49 charts per
fibre, height 125000 and ten seconds per chart. All maps precede all points.
There is no adaptive refill. Exact admission histories, rational geometry,
full point-cloud bounds modulo 2, 3 and 5, final equation comparisons and a
standalone point-proof replay are separate stages. Scores do not prove rank.

Sources: `../cas/newfamily/scan_joint_cache_annulus.cpp`,
`../cas/benchmark_joint_11952_annulus.py`,
`../cas/scan_joint_11952_fresh_annulus.py`,
`../cas/score_joint11952_fresh.py`, and the `joint64` point/proof controllers.
The benchmark certificate is `joint11952_annulus_benchmark_v1.json` under
`artifacts/generated-results/elliptic-curves`. Frozen protocols and raw logs
are under `joint11952-annulus-benchmark-v1`, `joint11952-fresh-annulus-v1`,
`joint11952-fresh-scores-v1` and `joint64-r17-pari-v1` in local artifacts.

The independent homogeneous specialization audit and its replay also pass:
all 32,768 prospective models have distinct j-invariants. Its certificate is
`joint11952_parameter_models_v1.json`. This proves within-pool distinctness,
not universal novelty.

All 64 point attempts finish within their declared boxes. The recorded final
subgroups have sizes 18 on two fibres and 17 on the remaining 62. The two
18-point candidates are at 409137/109156 and 415742/514523. Independent
generic transports, rational geometry, point-cloud and final equation
certificates all pass; there
is no high-rank inventory addition from these search outcomes.

The final results are `joint64_r17_results_v1.json` and
`joint64_experiment_v1.json`. All 64 equations are unmatched against the
frozen 1,101 prior equations and 593 catalogue equations. The 59,067,118-byte
point-only archive is complete; all 194 isolated replay stages pass.
The inventory remains at 195 curves with certified lower bound at least 22.
