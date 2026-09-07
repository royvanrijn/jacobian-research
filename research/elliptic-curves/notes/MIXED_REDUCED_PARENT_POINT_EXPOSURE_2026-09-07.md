# Matched score strata and remaining rank26 exposure gaps

The experiment is **complete**: all441 point boxes and nine exact histories
finished, and all nine clouds have independent certificates. Four sibling
fibres gain six directions in total. The immutable intake, source hashes,
selector output and seed replay are packaged in
[`mixed_reduced_parent_intake_v1`](../../artifacts/generated-results/elliptic-curves/mixed_reduced_parent_intake_v1/manifest.json).
The [terminal report](../../artifacts/generated-results/elliptic-curves/mixed_reduced_parent_exposure_v2.json)
records final lower bounds18/19/19/18/15/16 and26/26/26. No curve reaches a new
near-record bound or enters the high-rank inventory.

| Stratum | Initial bounds | Final bounds | Certified gains | Worker seconds | Geometry, worker, history and proof seconds |
|---|---|---|---:|---:|---:|
| Strong |17,17|18,19|3|498.60|1060.77|
| Moderate |17,17|19,18|3|217.09|404.29|
| Fixed lower sample |15,16|15,16|0|310.12|598.45|
| Retained26 gaps |26,26,26|26,26,26|0|155.73|198.27|

The moderate pair delivers the same three directions with less computation.
This supports retaining a broader score portfolio for the next decision;
two curves per stratum do not establish a calibrated success law. The lower
stratum also has smaller initial section spans. All47361 retained point
witnesses are checked; lower bounds agree modulo2,3,5. None of the six sibling
equations matches the pinned630-entry catalogue or201-curve inventory.

## Mathematical and experimental gate

The [full blind factor-free control](BLIND_FACTOR_FREE_CONTROL_AND_PROSPECTIVE_EXPOSURE_2026-09-07.md)
recovers an independent28th direction from the original27-point subgroup,
without the public exceptional point. The new
[determinant1092 parameter chart](DET1092_REDUCED_PARAMETER_CHART_2026-09-07.md)
has an independently verified small integral equation and seventeen generic
sections. Its first two fixed fibres completed their point exposure, but
were not selected for arithmetic score.

This experiment combines a small height-matched score comparison on that
new parent with the last three retained26 source-centre gaps on X948.
Different X948 fibration labels are not different parent surfaces; the
portfolio contains **two parent surfaces**.

The latest parallel
[class-block and rational-lift proof](../rank-jump/CONSTRUCTED_CLASS_BLOCK_AND_RATIONAL_LIFTS.md)
constructs two independent strict classes from equation/generic inputs, then
proves both covers rationally soluble using known-point witnesses admitted
retrospectively. That supersedes the earlier unresolved solubility statement.
It does not yet supply autonomous point recovery from the constructed covers;
neither its oracle-labelled representatives nor exceptional points enter this
frozen selector or worker.

## Frozen selection

The entire new-parent population has **87 rational addresses**, with
|m|,n≤8 in the reduced chart. The previously searched coordinate1 is omitted,
leaving86 equations. The rule also excludes previous original-chart fibres
and the construction anchor and conservatively deduplicates internal
j-invariants. There is no catalogue-driven selection or replacement.

For every equation, the existing corrected scorer removes exact p⁴/p⁶ short
model scalings before deciding whether reduction is good. The score is
the sum of rounded10¹²(2−a_p)log(p)/(p+1−a_p) over training primes5..32749.
Independent direct character-sum checks include every prime through199,
hence5 and13, and fixed larger checks. Validation primes65537..131071 are
neither computed nor read. Scoring workers consumed7.699501309seconds;
this does not include Python preparation and direct-check overhead.

Choose the most populated64-bit j-numerator-height bin between256 and639,
breaking ties toward the smaller bin. The chosen bin is256..319 with71
equations. Within descending training-score order, select the first two,
positions floor(N/3) and floor(N/3)+1, and two SHA-selected entries from the
bottom third. Those rules were frozen before scores and point searches.

| Stratum | Reduced coordinate | j numerator bits | Certified initial span |
|---|---|---:|---:|
| Strong | −3 | 271 | 17 |
| Strong | 4/3 | 295 | 17 |
| Moderate | 7/2 | 275 | 17 |
| Moderate | −5 | 316 | 17 |
| Lower fixed sample | −1 | 262 | 15 |
| Lower fixed sample | 1/7 | 283 | 16 |

All six reach large original parameter addresses through the certified
rational base matrix. Exact group relations prove the specialized section
losses; no fibre is replaced because its initial span is smaller.

The retained arm uses precisely the still-unexposed IDs73,75,63, in the prior
frozen training-score order. IDs74/50 and42/49 have already completed their
own-subgroup follow-ups. The three remaining inputs have certified26-point
bases and retain those points unchanged. No additional parameter search or
rescoring occurs for this arm.

## Exposure and certification

Every curve receives2048 deterministic parity samples and49 largest computed
norm centres. New-parent centres use the certified specialized section span;
retained26 centres have nonzero parity above the inherited generic17 subgroup.
Numerical heights propose a metric, with exact rounded norm/parity transport;
there is no covering or optimality claim.

All441 factor-free maps were prepared before any point search. Every chart
receives height125000 and ten seconds, with no rank stopping, replacement,
adaptive wave or automatic enlargement. Single-worker stage limits are180
seconds for geometry and1200seconds each for point processing and exact history
replay, with2GiB RSS. The larger processing limits accommodate the many point
witnesses in small arithmetic fibres; they do not enlarge the point boxes.

Seed preparation uses exact finite-independent subsets. For missing generic
directions,384-bit heights propose relations with bounded denominator and
coefficient at most64, accepted only after exact rational group equality.
A standalone Sage replay checks all specialized source sections, point
memberships, score-model transports, original parameter addresses, retained
input equality, and rational-span identities before geometry begins.

Completed clouds must pass mod2 and mod3/mod5 construction and replay,
independent exact geometry, and a standalone finite-group rank check: six
certificate stages per curve. The final report must record completed exposure,
certified gains, elapsed computation and post-search catalogue comparison.
Finite-image ranks are lower bounds, not upper bounds on the curve or cloud.

The initial section spans differ between strata despite matching arithmetic
height. This small comparison is descriptive: report gains from the actual
15/16/17-dimensional starting spans and final delivered bounds separately.
It cannot identify a causal score effect or prove parent superiority.

## Replay and continuation

```
python3 elliptic-curves/cas/package_mixed_reduced_parent_intake.py --check
python3 elliptic-curves/cas/report_mixed_reduced_parent_exposure_v2.py --check
```

The original one-worker driver stopped during certification when parallel
repository maintenance moved its working files under `research/`. At that
instant all441 point boxes, all nine histories, and four whole-cloud proofs
had finished. Original interrupted ledgers and both moved/stale-path process
records are preserved. `complete_mixed_relocated_certificates.py` copied the
unchanged seed/map/result inputs and ran certificate-only replays for the five
unfinished rows. No point search was repeated or rule changed. The continuation
ledger lives in `artifacts/local/elliptic-curves/mixed-reduced-parent-relocation-replay-v1/`.
The terminal report charges31.217197454seconds for the interrupted certificate
driver in addition to completed stages. Total recorded stage time, including
scoring workers and preparation, is2307.008779059seconds. Selector preparation
and direct character-check overhead are not included in that figure.

A separate [exact centre-class criterion](POINTED_QUARTIC_CENTRE_CLASSES_2026-09-07.md)
now explains why newly certified directions can supply different degree-two
presentations for a later fixed-subgroup follow-up. It does not change this
frozen campaign or guarantee additional point discovery.

The [separate focused trial](CURVE302_FOCUSED_POINT_EXPOSURE_2026-09-07.md)
uses302 as a generic17 recovery control and full31 direct target, with the
two certified19 siblings at reduced s=4/3 and7/2. It completes196 boxes,
recovering19 from17 on302 and no stronger public31 or sibling19 bound.
The score comparison above remains unchanged.
