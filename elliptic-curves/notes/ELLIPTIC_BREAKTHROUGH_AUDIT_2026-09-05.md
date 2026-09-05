# Elliptic-curve machinery audit and prospective results

The [latest fixed MW16 follow-on](NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md)
adds four curves with certified lower bounds 22–24, bringing the inventory to
36 distinct curves absent from the refreshed 586-row ICARM snapshot.
Its `3/17` curve has rank at least 22 and a proved exact 76-digit conductor,
third among recorded rank-at-least-22 conductors in that snapshot. Four such
catalogue entries lack conductor data. This is a low-conductor near-record
result; exact rank and new rank-at-least-28/32 targets remain open.
The [index](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v2.json)
contains all 36 equations and point certificates.
The original audit and earlier experiments below retain their own scopes.

The first concrete outcome is [fifteen newly certified prospective curves of
rank at least 22–24](NEW_COMPACT_R17_CURVES_2026-09-05.md), absent from the
pinned 584-equation ICARM snapshot. This does not meet the rank-at-least-28
near-record or rank-at-least-32 record target. Those targets remain open.

## Scope and limits

The starting inventory pins 4,726 tracked files: 765 elliptic-programme files,
2,823 K3 files, 478 generated-result files and 660 archived files, totalling
1,574,420,667 bytes. It identifies 321 relevant mathematical-status entries.
The starting commit is `300ee96`; pre-existing untracked visibility work was
preserved. Other work continued in the shared checkout during this audit.

This is a complete file inventory and programme-level triage, with deeper
examination of the selection, specialization, point, independence, cache and
checkpoint paths used below. It is **not** a claim that every historical line
or large K3 enumeration was independently checked. The prior
[external audit](EXTERNAL_AUDIT_2026-09-04.md),
[structural reassessment](RANK_JUMP_REASSESSMENT_2026-09-05.md), and
[subsequent diagnostics](RANK_JUMP_DIAGNOSTICS_2026-09-05.md) retain their own
evidence and claim boundaries. Their work is not relabelled as this audit.

The raw inventory, bounded-test logs, downloaded source snapshot, exact
recognition ledgers, failed-worker logs and search manifests are retained in
`artifacts/local/elliptic-curves/breakthrough-audit-20260905/` and the sibling
`compact-r17-*` directories. Mathematical status remains in
[`MATH_STATUS.json`](../../MATH_STATUS.json).

## Findings that change the next computation

The later [recorded-point admission audit](RECORDED_POINT_ADMISSION_AUDIT_2026-09-05.md)
finds another concrete gap: prime 257 certifies a 26th direction already
recorded by the blind MW16 control search, whose admission bound was 251.
All 202 retained-cloud proofs replay; only that known control improves.
The same audit identifies avoidable per-point cache storage and validates an
optional quotient-only cache for future workers.
The later [small-conductor follow-up](SMALL_CONDUCTOR_FOLLOWUP_2026-09-05.md)
measures a separate quadratic cost: every state construction revalidates all
old point observations. An optional bounded exact-membership cache preserves
the original serialized states and reduces the measured long-history cost.

**1. The first selection stage can discard the strongest curves.**
On the exact 20,400,078-parameter height-4096 population, the global
short-prime cutoff for the best 1,024 parameters is `61.094343749031`.
The published rank-25/26/27 controls in that box have prefix scores
`60.382592278581`, `60.829596942110`, and `59.192715807466`.
All three are lost before the extension primes can help. The rank-28 control
is outside this height box and is not counted as a lost in-box observation.

This is a scheduling failure, not a mathematical rank exclusion. It was not
fixed by moving the failed cutoff retrospectively. A separate experiment
evaluates all 562 primes from 5 through 4093 before selection. It retains the
published rank-25 control second, and independently rediscovers ICARM curves
414 and 417 with 25 and 26 independent points. The top-64 continuation is
separately frozen. This does not establish a general enrichment theorem.

The [completed full-prime control audit](../../artifacts/generated-results/elliptic-curves/compact_r17_selection_control_audit_v1.json)
shows the remaining failure: the rank-26 parameter `-308/251` is only 667th
at height 4096, and the rank-27 parameter `2456/135` misses even its top 1,024.
Both miss the larger box's top 1,024. The rank-25 control is second/eleventh
in the two boxes, while the rank-28 control is outside the smaller box and
third in the larger. Thus the repair removes an avoidable truncation without
making full-prime score a reliable rank classifier.

**2. Prospective coefficient size and exposure were not matched to the
record calibrations.** The earlier fresh MW16/MW18 fibres were often far
larger than their successful controls. The compact published-R17 trial
produces real gains on fresh fibres using the same pointed search engine:
for example `33/119` gives seven new independent directions beyond MW17.
This demonstrates a productive prospective route, not that arithmetic height
alone causes rank or that every old null was a detector failure.

**3. The production default still selected shallow centres.**
`run_mw_search.py` uses `VoronoiIterator.next_holes`, which enumerates shallow
classes first. The explicit deep-centre machinery is separate. A migration
to that default therefore does not preserve the original 43-deep-class
experiment. The new trials explicitly retain those 43 generic parity labels
and choose their specialized representatives using numerical heights. These
representatives are scheduling data; exact rank proofs do not use their
optimality. No claim of an exact specialized covering radius is made.

**4. Chart checkpoint addressing omitted state that replay required — fixed.**
`PointedQuarticSearch.search` keyed a saved chart by source, curve, generators,
centre and budget, but `verify_record` also required equality of the complete
`MWState` key. Retrying an incomplete chart can add observations without
changing that curve or basis. A later identical chart then hits an old cache
entry that replay rejects as belonging to another state. This occurred in
the `33/119` worker after its first wall-time stop.

The checkpoint key now includes `mw_state_key`. The regression constructs
two states with the same basis and chart but different known-point observation
histories, searches both, and replays each without enumeration. The compact
runner also restores a retained partial state instead of overwriting it with
a restarted initial state. Generation-time source bytes are preserved in the
existing source-snapshot archive; old mathematical witnesses still replay.

**5. Cheap point-signature caching can exhaust the filesystem — addressed
in the new workers.** Millions of per-point/per-prime object, index and lock
files exhausted all 4,086,745 tmpfs inodes. There was still byte capacity;
the resulting `ENOSPC` failures were not point-search nulls. The old cache
was preserved in a Zstandard archive, integrity-tested, and hashed before
the temporary directory was removed. Its archive hash is
`6c1c05b6ae9523e012131b9318af64c5b89f58503ce7428b68dc1027cd85043a`.

The new `MemoryFactStore` retains identical immutable keys and portable
snapshots inside each bounded worker, without one filesystem file per cheap
fact. Tests compare its finite quotient witnesses and `MWState` with the disk
backend, replay across backends, and reject tampered/conflicting snapshots.
Durable per-chart checkpoints remain. Expensive shared number-field facts
still use the disk backend. Point-admission trials now stop at prime 251;
every ambiguous point is retained, and the bound is recorded per chart.
This is a smaller proof-search allowance, not a dependence assertion.

**6. The external comparison data had advanced.** The live snapshot contains
584 curves, including new public lower bounds 30, 29 and 28 beyond the
repository's complete 573-row sweep. An exact equation-only intake tests all
eleven new rows against the six norm-twelve classes and five retained A1/MW16
fibrations. Curve 579 hits `0e80b`, 580 and 581 hit `074d9`, and 583 hits
`08234`, all untwisted over Q. Curve 582 misses all eleven tested families.
These finite-family exclusions neither identify its parent nor exclude other
fibrations. Its public rank metadata is not independently certified by this
equation-only intake. A cheap Sage rational-isogeny check on all ten public
rank-at-least-29 curves returned no prime-degree isogeny; it supplies no new
curve by an untested isogeny shortcut.

**7. Ten pre-existing checker-source hash mismatches were audited and rebound.** They affect
the rank-28 bad-place/local/S-class/residual helpers, the fixed-cubic local
and rank-one/CT helpers, the product-twist rank-zero/regulator checkers, and
the anchored-MW18 checker. Their exact identities and old/current hashes are
in the initial inventory. All ten generation sources were found in the
existing source archive and matched their recorded hashes. Six exact
certificate replays passed. The four rank-28 helper paths have no cheap full
replay mode: their changes and narrow regressions were checked, without
rerunning their historical large local/class-group/descent calculations.
Two implementation-text assertions in the relative-descent tests were stale;
the updated test executes the worker against a controlled arithmetic adapter
and checks unconditional-purpose forwarding and failure propagation.

The [compatibility ledger](../../artifacts/generated-results/elliptic-curves/checker_source_compatibility_audit_v1.json)
and [portable logs/diffs](../../artifacts/generated-results/elliptic-curves/checker_source_compatibility_audit_v1.zip)
separate these verification modes. Active checker-code hashes were updated;
generation-time evidence, all mathematical scopes, and all partial states
were retained. This metadata repair is not a new completed rank-28 descent.
The status renderer and its `--check` now pass. The new curve lower bounds
are recorded separately as `EC-COMPACT-R17-NEW-CURVES-20260905`.

## Programme gaps after this audit

| Programme | What is available | Missing endpoint that matters |
|---|---|---|
| Compact R17/MW17 | Six compact families, 102 exact generic sections and new lower bounds through 25; native known-rank controls | A distinct curve of rank at least 28, then 32; full fresh generic parity censuses for all six presentations and balanced wider searches remain separate unrun experiments |
| A1/MW16 | Five compact fibrations, 80 exact sections and fifteen new curves with lower bounds 22–25 across the completed pilots and next-twelve follow-on | Prospective gains currently reach nine beyond MW16; rank 32 needs sixteen certified quotient directions. The new small-conductor rank22 fibre needs one further direction to beat the recorded rank23 conductor minimum |
| Anchored MW18 | Exact rank-18 specializations and deep-centre diagnostics | Adaptive rather than initial-only recovery, and a point-producing extension that does not merely re-count the same anchor |
| Fixed cubic / 2-Selmer / CT | Exact transported cohomology and insolubility results | Global rational points on surviving covers; a fixed cubic algebra does not preserve the old high-rank group |
| Rational/quartic base changes | Rank-19 carriers and bounded twist calculations | Independent characteristic-zero sections plus the required character/glue proof for a rank-20 construction; the seventeen excluded product twists remain excluded |
| Different-NS foundry | Arithmetic marking obstructions and a remaining unknown queue | A rational non-CM point on the full marked source before equation work; another determinant-948 chart does not meet this endpoint |
| Curve-302 and new-582 parents | Exact recognition misses and point-cloud diagnostics | A constructive parent with transported independent sections; a large displayed jump or another bounded miss is insufficient |
| Older Mestre/Fermigier/Nagao/Kihara/newfamily lanes | Reproducible lower-rank families and some exact ranks | A demonstrated advantage for record-scale construction, not unbudgeted reopening of all old sweeps |
| Record/conductor comparison | A new rank22 curve with a proved global minimal equation and exact76-digit conductor; refreshed586-row comparison | Exact rank remains open. Four catalogue entries at rank at least22 lack conductors, and database absence remains a finite novelty statement |

## Verification and completed experiments

The bounded checks of residual quotient algebra, production proof gates,
fixed-cubic geometry, observability, arbitrary-prime independence, pointed
search and finite-quotient escape all passed. After the checkpoint repair,
all 13 pointed-search tests and all three lazy-search tests passed; both
memory-store tests passed. The 24-point `33/119` full `MWState` transcript
also replayed through all 43 retained charts, including its mixed-generation
history. The compact point certificates have separate Sage-free replays.
These checks do not certify the unrelated whole repository.

The top-64 fullscore continuation fixes height 4096, the retained full-prime
ranking, 43 centres per fresh fibre, point height 100,000, four seconds per
chart, four concurrent workers, 300 seconds and 1.5 GiB per worker, and a
one-hour submission limit. Completed measurements are reused; exact matches
in the pinned public snapshot are skipped. A certified rank-at-least-28 hit
stops further dispatch for independent verification. Unrun or stopped rows
stay censored. The original short-prefix and height-256 populations remain
unchanged.

The two initial rank-24 curves also completed 128 sparse adaptive charts
each, with no certified gain. All 256 chart/admission records replay. Audits
of all retained distinct ambiguous points modulo 3 and 5 through prime 997
also stop at finite column rank 24. This is a bounded negative result, not an
exact-rank claim. The canonical curve note links the portable transcripts and
the additional finite-quotient certificates.

A separate height-16384 expansion completed 326,397,350 parameter scores,
retaining the same 562-prime policy and the top 64 for bounded measurement.
The eight scanner shards each completed inside 180 seconds, with two at a
time. The published rank-28 parameter is now inside the box and ranks third.
Its known equation is skipped prospectively and receives a separate MW17-only
control search. Four workers per cohort remain bounded by the same 300-second
and 1.5-GiB limits; the end of the smaller cohort overlaps the beginning of the
larger one. The expanded box is a disclosed exploratory continuation, not
evidence that enlarging parameter height guarantees higher rank.

The smaller top-64 cohort finished: 53 fresh workers, seven reused completed
measurements and four exact public-equation matches. Eight additional curves
of rank at least 22–24 were independently exported, giving fourteen distinct
prospective curves across the retained experiments.

The rank-28 control recovers rank 26 initially, then rank 28 after only eleven
charts under the broader 301-class rule from the earlier blind ladder. Its
[independent point certificate](../../artifacts/generated-results/elliptic-curves/compact_r17_blind_rank28_recovery_v1.json)
matches ICARM curve 11, so this is calibration. The three new rank-24 curves
receive the same broader follow-up. A separate four-fibre cohort selects the
highest full-prime scores with **any** certified initial gain: `-5284/455`,
`-636/3209`, `-4253/2347`, and `-9795/13063`, starting at ranks 19,18,19,18.
This avoids requiring initial rank 24 when the earlier ladder already showed
large adaptive gains from initial ranks 20/21. Each follow-up fixes at most
301 charts, four seconds per chart, 1,500 seconds and 1.5 GiB per worker.
These are finite follow-ups, not claims that every ambiguity becomes a point.

All seven prospective follow-ups completed all 301 charts without certifying
another direction. Two rank-24 workers first reached the 1,500-second limit;
separate 300-second continuations completed only their remaining frozen charts.
Every chart and admission replayed. The first read-only replay of `33/119`
reached its 120-second checker cap; a separately recorded 300-second replay
passed. All seven final point clouds were checked modulo 3 and 5 through
prime 997, with no higher certified finite column rank.

The larger top-64 cohort completed 58 fresh workers, three reused measurements
and three public matches. Its unique strongest fresh initial measurement was
`-7540/2317`, at rank at least 21. A separately frozen 301-chart follow-up
completed in 1,296 seconds, gained one independent direction, and exported the
fifteenth distinct curve, with rank at least 22. The canonical curve note
links its independently replayed point certificate. All 301 chart and
admission records also replayed.

**Coverage gap measured and tested.** All 2,279 chart records from the 53
fresh height-4096 workers stopped before full-box completion. Their mean
completed-denominator fraction was about 0.715. This is bounded coverage,
not 2,279 exhaustive boxes. A separate tail-completion protocol selected the
four highest full-prime scores with no initial gain in the larger population:
`-919/2282`, `-1121/9525`, `11744/9849`, and `10126/15275`. It verifies the
same quartic, basepoint and coordinate maps, and resumes after each retained
prefix's last completed denominator. All 172 boxes then completed through
height 100,000, still with only the seventeen certified generic directions.
All tail records replayed. The
[manifest](../../artifacts/generated-results/elliptic-curves/compact_r17_tail_replay_v1.json)
and [portable witnesses](../../artifacts/generated-results/elliptic-curves/compact_r17_tail_witnesses_v1.zip)
retain both prefixes and tails. This closes that finite coverage gap for these
four fibres; it is not an exact-rank bound or a global detector guarantee.

The [population manifest](../../artifacts/generated-results/elliptic-curves/compact_r17_population_manifests_v1.json)
and [portable protocols and scan ledgers](../../artifacts/generated-results/elliptic-curves/compact_r17_population_manifests_v1.zip)
preserve the successive, separately declared selection policies. Increasing
the parameter box did not improve the best prospective lower bound. The
remaining bottleneck is selecting and exposing exceptional directions on new
fibres, not a demonstrated shortage of raw parameters.

The [complete initial-measurement certificate](../../artifacts/generated-results/elliptic-curves/compact_r17_initial_measurements_v1.json)
independently checks all 111 fresh workers in the two full-score cohorts,
including every rank-at-least-17 result. Its portable replay is
`python3 elliptic-curves/cas/audit_compact_r17_measurements.py --check FILE`.
The 53 smaller-box workers give measured lower-bound counts
`17:23, 18:9, 19:8, 20:4, 21:1, 22:1, 23:6, 24:1`; the 58 larger-box
workers give `17:45, 18:4, 19:8, 21:1`. All 4,773 initial chart boxes were
incomplete at their recorded stops. The measured mean denominator coverage
was 0.715 and 0.682, respectively. Scores and coverage remain pinned worker
metadata; the independent checker redoes point membership, finite-quotient
independence and generic-section transport, not enumeration completeness.

The [comparison figure](../../artifacts/generated-results/elliptic-curves/compact_r17_initial_measurements_v1.png)
and [PDF](../../artifacts/generated-results/elliptic-curves/compact_r17_initial_measurements_v1.pdf)
plot these initial measurements against score and short-model coefficient
size. Reused measurements and known curves are omitted. The lower bounds
are not an exact-rank distribution, and the score/height pattern does not
establish causation. The separate adaptive results above remain separate.

A cheap [constant-scaling audit](../../artifacts/generated-results/elliptic-curves/r17_constant_scaling_audit_v1.json)
tests whether the large coefficients of the six literal compiled R17
presentations are merely a removable common Weierstrass scale. It strips
weighted-content primes through 997 and remaining exact twelfth powers,
then checks every coefficient identity under `x=u^2*X`, `y=u^3*Y`.
The largest numerator/denominator bit sizes fall from 1584 to 1207 (`074d9`),
1621 to 1289 (`07ca9`), 1603 to 1237 (`08234`), 1887 to 1503 (`08f72`),
1371 to 1141 (`11952`), and 991 to 786 (`103b2`, representing `0e80b`).
This removes some scale but does not establish a compact alternative family.
It is not a minimal-model certificate: base-parameter changes and other
coordinate reductions remain separate possible improvements. None of these
identities constructs a new curve or increases generic rank.

The subsequent [base-coordinate compactification](COMPACT_SIX_R17_ATLAS_2026-09-05.md)
does close the large-family-coefficient preparation gap: all six compiled
models now have coefficients of 141–169 bits and all 102 generic sections
transport exactly. The decisive missing primes were recovered by perfect-power
extraction from a common auxiliary invariant before bounded factorization.
The new small-parameter boxes are different populations, so this does not
retroactively change any old point-search result. A separately frozen
height-1024, four-finalist-per-family pilot tests the practical consequence.

Primary context: [Elkies's explicit R17 paper](https://arxiv.org/html/2608.25406v1),
[Elkies's high-rank search lectures](https://arxiv.org/abs/0709.2908), and the
[ICARM leaderboard](https://elliptic-rank.icarm.cloud/). The broader theoretical
sources and their hypotheses remain indexed in the structural reassessment.
