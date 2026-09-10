# Broad rank search: runnable, finite and token-free between batches

**Operational status, 2026-09-10:** both real backends passed commissioning, and
`artifacts/local/elliptic-curves/broad-rank-v1` is launched with the default eight
presentations and four workers. Its [live report](../../artifacts/local/elliptic-curves/broad-rank-v1/REPORT.json)
records progress; the campaign is not a completed search result. The older
two-parent MW17 campaign was never dispatched and is stopped, with its inputs
and smoke evidence preserved. Mathematical status is unchanged. Use a new
directory when preparing another campaign.
The [historical ledger](HIGH_RANK_CONSTRUCTION_AND_SEARCH_LEDGER_2026-09-10.md)
explains why the productive six-presentation R17 route gets most of this exposure.

## Default experiment

| Group | Presentations | Score-selected per presentation | Controls per presentation | First search stage |
|---|---:|---:|---:|---|
| X948 native R17 | 6 | 256 | 64 | Complete43/49 maximum-class first-seed routine, at most98 calls; then at most100 complement calls after **any** certified gain |
| X1092 original/class1 | 2 | 64 | 16 | Existing generic-parent V3 worker, at most198 calls |

Every presentation scores the **same262144 signed rational addresses**, at indices
131072 through393215. This gives2097152 feature rows and2080 selected slots,
before exact singularity/isomorphism/exclusion attrition. The addresses are unique
within each stream, not uniformly distributed over rationals or heights. Same t
on different fibrations is not a matched elliptic curve.

Default four workers; `--workers` permits1–8. The six native IDs are `074d9`,
`07ca9`, `08234`, `08f72`, `103b2`, `11952`. Comparisons are `x1092-original` and
`x1092-class1`. `--parents` restricts the roster without requiring unused inputs.
These are eight working presentations on **two** K3 surfaces, not eight surfaces
or eight independent geometric frame types.

## Start here: real backend smoke

Run from the repository root, with Linux, Sage, PARI at `/usr/bin/gp`, and NumPy.
The Python tests require NumPy too. Source and input files come from this checkout;
no historical ignored point-search directory is required by the new launcher.

```sh
python3 -m unittest discover -s research/elliptic-curves/tests -p test_broad_rank_search.py -v

SMOKE=research/artifacts/local/elliptic-curves/broad-rank-smoke-v1
python3 research/elliptic-curves/cas/run_broad_rank_search.py prepare \
  --folder "$SMOKE" --workers 1 --parents 074d9 x1092-class1 \
  --window 256 --ranked 1 --controls 1 \
  --comparison-ranked 1 --comparison-controls 1 --rounds 0
python3 research/elliptic-curves/cas/run_broad_rank_search.py run --folder "$SMOKE"
python3 research/elliptic-curves/cas/run_broad_rank_search.py status --folder "$SMOKE"
```

`--sage /absolute/path/to/sage` can be supplied during preparation. The smoke must
include at least one accepted fibre from **each backend**, with completed uncensored
point search and a sealed two-implementation rank replay. A zero-call preflight,
a synthetic fixture, an exit code or an isolated PASS string does not satisfy this.
Inspect `REPORT.json`, `batch-000/result.json`, `packet-verified.json` (native) or
`verified.json` (generic), and their seals. A generic finite-gate UNKNOWN is not
an engineering failure or a rank statement; retain it, but do not call that
backend's end-to-end smoke passed on its account. Run a fresh small fixed smoke
window if the first one contains only aliases or generic-gate misses.

**Validation performed here:**36 Python tests pass. They include a separately
written enumeration of every finite point for small-prime score tables, score/
control invariance, rational scaling, mocked-CAS controller dispatch and resume,
certificate accounting, censorship, STOP and evidence tampering. Python compilation
also passes. Sage/GP were unavailable during the original implementation.

**Local backend commissioning, 2026-09-10:** the exact four-slot smoke above
completed all four accepted fibres. Both native `074d9` fibres have certified
lower bound19 after101 calls each; both class1 fibres retain lower bound17 after
198 calls each. All598 calls have individual attempt and result receipts, with
zero censored exposures or engineering failures. Every final packet passed both
finite certificate implementations; all seven job seals and derived state
bindings verified. The [smoke validation receipt](../../artifacts/local/elliptic-curves/broad-rank-smoke-v1/SMOKE_VALIDATION.json)
retains parameters, hashes and software versions. The prepared full runtime
matches all3,940 shared arithmetic files and the executables tested by the smoke;
additional parent inputs and unrelated new source files are recorded in its
`smoke-runtime-comparison.json`. The built-in detached `launch` command is used.
These are commissioning lower bounds, with worldwide novelty unexamined.

## Prepare and launch the full frozen campaign

Once the smoke passes:

```sh
RUN=research/artifacts/local/elliptic-curves/broad-rank-v1
python3 research/elliptic-curves/cas/run_broad_rank_search.py prepare --folder "$RUN" --workers 4
python3 research/elliptic-curves/cas/run_broad_rank_search.py launch --folder "$RUN"
python3 research/elliptic-curves/cas/run_broad_rank_search.py status --folder "$RUN"
```

`launch` creates a detached Linux process; the search itself makes **no model calls**.
`run` is the equivalent foreground operation. The frozen controller, arithmetic
sources, parent projections and Sage/GP executable hashes are checked before work.
Ordinary checkout edits do not silently change the active runtime. The whole
address universe and selection counts are fixed before scoring; no online
retuning, extra frame realization, class-group factorization or automatic
parameter-window expansion occurs.

Both source comparison parents already have replayed equation/section certificates.
Preparation rejects a changed class1 parent hash. A new legitimate parent requires
an explicit adapter/certificate review, not updating a hash merely to silence a gate.

The current inventory supplies **equation-only** exclusions. Additional frozen
exclusion files can be added at preparation with repeatable `--exclude-equations`;
use `{"models": [[a1,a2,a3,a4,a6], ...]}` or a database-style `curves/ainvs` list.
No exceptional point coordinates or rank labels are copied for candidate selection.
Coverage/freshness is only relative to these snapshots and this queue. Previously
searched low-bound curves absent from the inventory are not automatically excluded;
export their equations explicitly when those local ledgers are available. This
avoids pretending that an inaccessible historical roster was checked.

## Why the first stage is not just24 calls

The inspected R17 panel found seeds late, amplified M18/M19 states substantially,
and recorded several gains at complement call75 or later. Here every successful
native first seed receives up to100 complementary calls, regardless of whether it
already reaches20. The X1092 comparison uses a broader198-call initial exposure.
It is a different detector, so compare parent-plus-detector and per-CPU yield,
not an alleged causal effect of choosing a particular fibration.

After all eligible baseline tasks have received a stage, bounded continuation
rounds give at most100 additional calls per task. Total call caps are398 below20,
512 at20–22,1024 at23–24,2048 at25–26 and8192 at27–31. The native initial seed calls
are included in these totals. Staleness caps are200 below25,300 at25–27,600 at28–29,
and900 at30–31. A new certified gain resets staleness at its originating point call.
These are scheduling choices, not rank upper bounds.

A fixed hash predicate selects approximately1/8 of slots for one **rank17 rescue**,
with a different full-space bank. It is not chosen after inspecting scores or
results. A rescued M18/M19 then uses the same continuation rules. Rank32 ends that
fibre's exposure; other scheduled fibres remain eligible. At most96 continuation
rounds are permitted by default. The job universe is finite even if gains keep
occurring. `--rounds 0` requests baseline stages only.

Scored and control arms use the same policy within each parent. Report baseline
endpoints separately from adaptive endpoints; adaptive depth is outcome-dependent.
Controls are chosen in a deterministic score-independent hash order and include
both signs. Known equation exclusions and cross-presentation aliases are retained
as attrition, not replaced by score-dependent candidate refills.

## Evidence, cost and interruption behaviour

The runner stores compressed JSONL for **all** cheap features, per-prime residue
tables with traces and smoothness flags, fixed selections and exact equations.
Individual a_p vectors can be reconstructed from a parameter and the tables
without repeating point counts. Cutoff scores at97,257,997 are saved when included.
Integer-weighted rational coefficient scaling is exact; scoring performs only
coefficientwise p-minimization. Parameter-specific extra minimization and exact
conductors are not claimed. Singular reductions are omitted score terms, not
mathematical exclusions.

Every backend result is bound to its request, model, parameter, returned packet
and independent replay receipt. All raw search journals/clouds and output files
are hash-sealed. Result integrity is not a substitute for mathematics: the
existing workers perform the actual independent rank replays. The controller
checks prefix preservation, rank monotonicity, gain accounting and budgets before
promoting a packet. A known-point class or generic arithmetic signature never
bypasses those gates. Worldwide novelty remains `NOT_CHECKED`.

Failed or censored exposure retains its logs and cost, not a fabricated rank17.
A replayed lower bound survives censorship with that status explicit. An
unreceipted interrupted dispatch is **not repeated automatically**; a live bounded
driver is adopted on resume, while an orphaned/incomplete one remains UNKNOWN.
This conservative behaviour avoids double-spending point calls without durable
return receipts. Recovering such a raw case is a separate manual decision.

Per-job CPU includes the worker and descendants using an isolated Linux subreaper,
including independent replays and driver work. Scoring, preflight, source
preparation and controller CPU are reported separately. Cost of interrupted jobs
without a seal remains unknown. Per-population-row scoring charge is the parent's
score-stage CPU divided by its frozen population size; first-wave total cost also
includes preparation and preflight. Do not count parent creation again as a new
cost if reusing an already certified parent, but retain its historical cold cost
when comparing complete construction strategies.

Each job is bounded by10800 seconds and3GiB process-tree RSS; individual old worker
phase and point/map limits remain in force. Default free-disk reserve is20GiB.
Four successive engineering failures halt new dispatch; mathematical generic-gate
misses and resource-censored tasks are distinguished. A controller crash does not
have an autonomous LLM or guardian retry loop; use `resume` after inspection.

```sh
# Graceful stop: current bounded jobs finish their receipts; no new tasks start.
python3 research/elliptic-curves/cas/run_broad_rank_search.py stop --folder "$RUN"
# Inspect before clearing the deliberate stop marker.
python3 research/elliptic-curves/cas/run_broad_rank_search.py status --folder "$RUN"
rm "$RUN/STOP"
python3 research/elliptic-curves/cas/run_broad_rank_search.py resume --folder "$RUN"
# Check retained input/output integrity; not a fresh arithmetic replay.
python3 research/elliptic-curves/cas/run_broad_rank_search.py verify --folder "$RUN"
```

`REPORT.json` includes per-parent/per-arm baseline and current lower-bound
histograms, calls, unknowns, aliases and CPU. Detailed gain timelines live in
`runtime/research/broad-cases/*/batch-*/broad-state.json`; rank>=27 events are exported
under `events/`. `COMPLETE.json` means all finite scheduled stages are terminal,
not that every case is a resolved mathematical negative.

## Relationship to the other broad wrapper

The concurrently committed `run_broad_mw17_search.py` and its snapshots are left
untouched. That wrapper defaults to the two X1092 inputs and delegates to the
24-call class1 policy. This successor has a separate directory/schema, explicitly
includes the productive native R17 seed path, gives fuller first-stage coverage,
and retains its own resume/integrity semantics. Do not combine their arms or costs
as if they were one frozen experiment. Do not run both unintentionally.
