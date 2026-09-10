# Autonomous high-rank search foundry

The foundry is a persistent Python/Sage/PARI controller for finding new elliptic
curves over Q with certified lower bounds 32, then 31, 30, 29 and 28. It does
not use model calls, Codex automations, or human decisions between batches.
The target is an independent rational point subgroup; no exact-rank upper
bound is required. Production launch and commissioning receipts are linked
below.

## Mathematical and empirical basis

The [sixty-fibre panel](R17_SIXTY_SEED_COMPLEMENT_PANEL_2026-09-09.md) supplies
the strongest prospective evidence: 48 seeds and 38 further cascades in 6,090
point calls, including one M27. This contrasts with the completed
[fixed-fibre norm12 campaign](UNATTENDED_NORM12_SEARCH_2026-09-09.md), whose
6,196 calls on four M27/M28 states supplied no gain. The first, second and
lower-height fresh cohorts reinforce both adaptive amplification and the
importance of reconciling every returned cloud. The independent curve302
seed panel and determinant1092 blinded reconstruction remain calibrations,
not estimates of fresh-fibre success rates.

The mathematical gate for every fresh candidate is exact native specialization
of one of the six compact MW17 families, point membership, a finite-reduction
independence certificate, and rational 2-torsion exclusion. Scoring, rounded
canonical heights and CVP geometry schedule searches; they do not certify
elliptic rank. The constructive search uses pointed quartics, exact rational
maps, retained PARI rational-point witnesses, and enlarged independent
subgroups. A bounded miss never excludes higher rank.

The review also covered `MATH_STATUS.json`, the structural rank-jump
reassessment, breakthrough audit, productive-parent reassessment, actual panel
trajectories and native point packets, the inventory of 291 curves, and its
195-exact/96-unresolved conductor publication cutoff. The newer live conductor
pass belongs to a separate controller. A fresh equation-only ICARM download
contains 687 curves, compared with the previous pinned 630. Its original bytes
and hash are retained in the commissioning input directory and its parsed
contents in every runtime snapshot. It supplies dated comparisons, not a
worldwide novelty guarantee.

## Allocation policy

Every scheduling choice is saved before execution, including the candidate,
heuristic utility, exposure totals, rejected aliases, parent bank, limits and
intake cursor. The defaults are explicit engineering choices informed by the
above experiments, not an optimized or prospectively proved policy.

- Reserve 60% of charged search worker time for fresh fibres, counting their
  first amplification pass as exploration. The other share follows recent
  gains, directions per point call, and distance to ranks 28–32. In-flight jobs
  count toward scheduling exposure. Bounded jobs can temporarily move the
  realized share away from its target.
- Balance accepted fresh candidates across all six families. In eight rounds,
  three draw from the low-height retained score order, one from the higher
  order, two from a deterministic score-independent permutation of the entire
  retained pool, and two from reduced signed rational parameters in increasing
  height. Excluded addresses do not consume a family/height slot. Exhausted
  retained streams fall back to the parameter producer. No score cutoff removes
  an equation from all future consideration.
- Deduplicate by exact rational isomorphism after a j-invariant index. Equal j
  alone does not merge quadratic twists. Exclude the pinned historical rosters,
  current inventory, public snapshot, and earlier verified foundry outputs.
  Existing prospective panel watchlist packets enter as declared historical
  starting subgroups, not fresh discoveries.
- A fresh seed pass tries the complete 43/49 maximum generic classes, stopping
  at the first certified M18 and retaining its entire returned cloud. Each
  seeded fibre then gets at most 100 amplification point calls, shared across
  all subgroup enlargements and reconciliation restarts.
- Parent banks start with disjoint blocks of the exact maximum catalogue, then
  use deterministic samples from the full 17-dimensional parity space. Both
  integer and rational exact CVP solvers check generic minima. Lower-shell
  selections first extend the maximum-class span, preserving directions beyond
  its hyperplane where applicable. Every exact minimum representative is saved;
  the chosen representative is reproducible. Both bounded coordinate policies
  and signed box deduplication are retained.
- A verified gain automatically permits further amplification, including at
  ranks 28–31. Warm continuations can inherit sealed landscapes and chart
  receipts without repeating completed point calls. A stalled batch changes
  parent bank. Cooling occurs after 200 no-gain calls below rank25, 300 at
  ranks25–27, 600 at ranks28–29, and 900 at ranks30–31. New certified directions
  reset staleness using their originating point calls, including cloud gains.
- A minority revival slot becomes available per 40 newly completed fresh
  fibres, at most twice per cooled curve. M17 revivals use new full-space
  generic masks. Cooling and retirement express a resource decision, never a
  rank upper bound. Reaching rank32 seals the result and lets other searches
  continue.

The initial retained population has 6,144 scored addresses. The independent
height-shell producer supplies continued intake after these rows are consumed;
it is not a precomputed finite queue. No new K3 construction, Selmer campaign,
or broad determinant1092 sweep is silently launched.

## Exact evidence and conductor gate

Search replay redoes the rational CVPs, maps, point witnesses and cloud audits.
Every complete returned cloud is reconciled and replayed again before another
landscape starts. If odd-prime cloud evidence exceeds the resolved subgroup,
the certified portion is preserved and the unresolved case is quarantined.

Portable certificates contain the exact equation, native parameter, all
independent points, torsion witness, and finite quotient signatures. The final
certificate retains only rank-increasing prime blocks; **both independent
finite-group implementations replay that smaller witness before publication**.
All original larger cloud certificates remain in the raw history. This removes
redundant verification work, not proof obligations.

Every successful job exports a portable certificate and a manifest binding its
search/replay evidence. Rank28–32 hits and conductor threshold improvements
get separate event files. Public matches with an already reported bound are
labelled reproductions. The source runtime, software versions, CAS executable
hashes, initial data, raw point clouds, exact CVP minima, cursors, failures and
resource receipts are preserved. A live report keeps trajectories and
family/selection-arm summaries; adaptive censoring prevents interpreting these
as a randomized causal comparison.

Conductor work is secondary and receives at most a 2% scheduling allowance.
A cheap exact small-prime scaling and discriminant screen identifies plausible
small arithmetic models. A large discriminant merely lowers priority. Selected
cases get a 30-second factorization attempt and a 120-second independent replay,
using the existing Sage generic Tate / PARI local-reduction agreement and prime
certificates. Rational models receive exact denominator-clearing transport.
Timeouts retain certified divisors/bounds and `UNKNOWN`. Only an exact conductor
can beat a reported rank-threshold benchmark; missing catalogue conductors and
world-record status remain explicit uncertainties.

Automatic exports do not rewrite `MATH_STATUS.json`, the curated inventory,
or historical proof notes. These machine-certified packets are available for
independent verification and subsequent registration, without holding up the
next search. No external submission or message is sent.

## Process and resource contract

The launcher copies source files and required inputs into an isolated runtime
snapshot and hashes every file. Ordinary edits in the checkout do not change
an active search. The Python controller uses a WAL SQLite ledger with FULL
synchronization, durable dispatch intents, exclusive locks and process-start
leases. A restarted controller adopts live workers; interrupted workers replay
their own saved prefixes. A durable point-call intent is written before PARI;
its complete return is journalled before point admission. If interruption falls
between that return and chart publication, the retry reuses the exact saved
return. Journal contents are bound into the exported evidence manifest. Every
retry preserves earlier supervision evidence.
Two crash retries are allowed per job; repeated proof/engineering failures halt
new dispatch after four consecutive failures. A detached guardian restarts a
crashed controller, with bounded backoff and a repeated-failure stop.

Fresh specialization has a separate bounded seed gate: the17 specialized points
must supply17 independent finite columns through prime1000 and a rational
2-torsion exclusion witness through prime200. Missing either witness produces
`UNRESOLVED_GENERIC_SEED`, not a claimed rank17 subgroup. The candidate remains
`UNCERTIFIED_SEED` with unknown rank and no point search, while intake continues.
These input outcomes neither increment nor reset the engineering-failure streak.
Equation mismatches, invalid points, backend failures and failed certificate
replays still follow the engineering halt policy.

Default limits:

| Resource | Limit |
| --- | ---: |
| Concurrent foundry workers | 4 |
| Point-search height | 125,000 |
| One point call / map construction | 10 seconds / 5 seconds |
| One search or replay phase | 1,800 seconds, 3 GiB process-tree RSS |
| One whole job | 10,800 seconds |
| Completed search charts in one job | at most 198 fresh; 100 continuation |
| Additional interrupted-invocation reserve | 2 per search job |
| Renewable UTC-day point budget | 60,000, including failure reservations |
| Renewable UTC-day worker-time budget | 172,800 seconds, including reservations |
| Free-disk / free-inode reserve | 20 GiB / 100,000 |

Daily budgets renew without an AI turn. Resource pressure waits without
claiming a mathematical miss. There is no fixed overall campaign-duration
stop. Graceful stop drains current jobs through proof replay. The detached
process survives the end of this Codex task. WSL or machine shutdown terminates
processes; `resume` restores work from disk. No operating-system boot service
is installed.

The ledger distinguishes completed logical chart calls from physical invocation
intents. Reused complete returns do not start a second PARI call; an interruption
before a return is durable can consume a retry. Failed jobs consume their whole
reservation. Daily worker charges conservatively include recovery time and
unfinished jobs carried over midnight.

## Commands

Run from the repository root:

```sh
python3 research/elliptic-curves/cas/run_high_rank_foundry.py status
python3 research/elliptic-curves/cas/run_high_rank_foundry.py stop
python3 research/elliptic-curves/cas/run_high_rank_foundry.py resume
```

`status` checks process-start tokens. `stop` waits for current bounded jobs and
replays. `resume` retains the source snapshot, ledger and completed results.
Changing the policy requires a separately named prepared campaign; do not edit
its frozen configuration. A new campaign automatically excludes and can import
verified prior foundry outputs.

Portable arithmetic replay, without raw search directories or a point search:

```sh
python3 research/elliptic-curves/cas/verify_high_rank_foundry_certificate.py CERTIFICATE.json
```

The launcher is [run_high_rank_foundry.py](../cas/run_high_rank_foundry.py).
The policy, intake, bounded job, arithmetic adapter and certificate compactor
are separate modules in the same directory. Focused regressions are in
[test_high_rank_foundry.py](../tests/test_high_rank_foundry.py).

## Commissioning and production

The [sealed commissioning packet](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-commissioning-v1/summary.json)
records 33 passing focused regressions and four preserved integration snapshots.
Four initial setup failures exposed an omitted native source dependency; two
more exposed a missing fresh-case directory. Both were corrected before any
point invocation in those failed jobs. All failures and their runtime snapshots
remain available. The subsequent 14 successful jobs replayed 1,376 completed
logical point calls, including cached continuations and repeated cloud restarts.

Six distinct fresh fibres have certified lower bounds 17, 17, 22, 23, 23 and 25.
These are new relative to the frozen repository and 687-curve public inputs;
worldwide novelty is not established. The best fresh commissioning result is
[`07ca9`, parameter `-943/412`, rank at least25](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-integration-v4/job-0000003-f-0c58f3edab4ef28a34de.json).
Its 101-call trajectory was
`17 → 18 → 19 → 20 → 21 → 22 → 23 → 24 → 25`.
One returned cloud supplied three independent directions, including two recovered
by reconciliation. The complete packet passed a separate portable replay.
There is no rank28+ discovery or new conductor-record claim from commissioning.

The recovery controls deliberately killed an owned search worker and, separately,
the detached guardian's controller. The worker resumed with all59 pre-existing
chart receipts unchanged and completed its rank23 proof. The guardian restarted
its controller and adopted both live workers without replacing their processes.
A real retained-chart control obtained18 points with one PARI invocation and
reused the complete return on its second call. An existing conductor11 curve
passed the bounded conductor builder and independent replay. These are engineering
controls, not fresh discoveries.

Production uses the separately frozen `high-rank-foundry-v3` runtime, four
workers, and the limits above. It excludes all completed commissioning curves
from fresh intake and imports their useful certified subgroups as historical
continuation states. Commissioning outcomes do not become labels in the fresh
score-independent streams.

- [Live report](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-v3/REPORT.md)
- [Frozen configuration](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-v3/config.json)
- [Source and input manifest](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-v3/manifest.json)
- [Detached launch receipt](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-v3/launch-receipt.json)

The run has no finite job limit. Python selects the next family, parameter,
continuation, parent bank, cooling/revival decision and eligible conductor job;
Sage/PARI execute the bounded arithmetic. No AI wakeup or scheduler is attached.

### Seed-gate recovery, 2026-09-09

The original production run completed67 jobs and6,706 point calls, certifying30
fresh curves. Its best fresh result was
[`07ca9`, `-1508/909`, rank at least26](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-v1/job-0000067-f-cba789931b68e747133c.json),
reached at call24; no new rank28+ or conductor-record hit was recorded. The run
then reached `t=-1` in four families. Their initial17-point seed witnesses did
not close, and the first controller incorrectly counted those bounded input
misses as four engineering failures. It drained the remaining jobs and stopped
at14:43 Amsterdam time. Its runtime,71 job records and all evidence remain intact.

Version2 separates those seed outcomes as described above. The repair passed44
focused regressions, including all six `t=-1` specializations, a valid M17 seed,
continued scheduling after repeated seed misses, and retained halting on real
backend failures. The
[four original cases also passed a supervised controller replay](../../artifacts/generated-results/elliptic-curves/high-rank-foundry-seed-recovery-controls-v2/verification.json)
with zero point calls and four explicit `UNKNOWN` seed outcomes.

The new snapshot inherits the entire stopped ledger and exact intake cursor.
All48 existing certificate packets replayed during preflight. Curves retain
their identities, lower bounds, trajectories, cooling and revival state; the
charged job history continues to count toward daily budgets. Cached landscapes
bind the old arithmetic source, so their raw evidence stays in version1 and
their next continuation uses a fresh parent bank. Previously failed inputs
remain quarantined with unknown rank. New job numbering continues after71.
The current report includes the inherited totals.

This migration is reproducible with `prepare --inherit-state OLD_FOLDER
--folder NEW_FOLDER`; the parent must be stopped. It creates a new frozen
source/input manifest without modifying the previous one. The standard
`status`, `stop` and `resume` commands above now select version3 by default.

### Unlimited daily budget and automatic continuation, 2026-09-09

At the user's direction, version3 removes the UTC-day point and worker quotas.
Per-job and process-tree limits remain finite, disk reserves remain active, and
physical point returns are still journalled before admission. Work grows in
bounded steps: 100 amplification calls initially, then +25 after each 100 fresh
fibres, capped at 300. At migration the next allowance is 175 calls. The
controller therefore keeps progressing through the queue instead of waiting at
midnight, while resource safety remains explicit.

The detached guardian also restarts a controller that exits cleanly without an
explicit stop or halt state. An explicit `stop` still drains and ends the run.
This makes a completed bounded controller invocation an automatic continuation
point while preserving the same ledger and frozen evidence.
