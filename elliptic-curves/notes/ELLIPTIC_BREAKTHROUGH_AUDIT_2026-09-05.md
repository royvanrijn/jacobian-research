# Elliptic-curve machinery audit and prospective results

The first concrete outcome is [fourteen newly certified prospective curves of
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
| Compact R17/MW17 | Exact generic sections; successful fresh ranks 22–24; calibrated recovery of known 25/26 | A distinct curve with certified rank at least 28, then 32; adaptive recovery on genuinely promising new fibres |
| A1/MW16 | Five distinct fibrations; 55/55 retrospective recovery; fresh bounded nulls | Selection and affordable exposure that transfer prospectively; a witness producing the required 16 quotient directions for rank 32 |
| Anchored MW18 | Exact rank-18 specializations and deep-centre diagnostics | Adaptive rather than initial-only recovery, and a point-producing extension that does not merely re-count the same anchor |
| Fixed cubic / 2-Selmer / CT | Exact transported cohomology and insolubility results | Global rational points on surviving covers; a fixed cubic algebra does not preserve the old high-rank group |
| Rational/quartic base changes | Rank-19 carriers and bounded twist calculations | Independent characteristic-zero sections plus the required character/glue proof for a rank-20 construction; the seventeen excluded product twists remain excluded |
| Different-NS foundry | Arithmetic marking obstructions and a remaining unknown queue | A rational non-CM point on the full marked source before equation work; another determinant-948 chart does not meet this endpoint |
| Curve-302 and new-582 parents | Exact recognition misses and point-cloud diagnostics | A constructive parent with transported independent sections; a large displayed jump or another bounded miss is insufficient |
| Older Mestre/Fermigier/Nagao/Kihara/newfamily lanes | Reproducible lower-rank families and some exact ranks | A demonstrated advantage for record-scale construction, not unbudgeted reopening of all old sweeps |
| Record/conductor comparison | Exact point certificates and refreshed public equations | Exact global-minimal conductor data before claiming a conductor record; database absence remains a bounded novelty statement |

## Verification and current continuation

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

Primary context: [Elkies's explicit R17 paper](https://arxiv.org/html/2608.25406v1),
[Elkies's high-rank search lectures](https://arxiv.org/abs/0709.2908), and the
[ICARM leaderboard](https://elliptic-rank.icarm.cloud/). The broader theoretical
sources and their hypotheses remain indexed in the structural reassessment.
