# Maintained warm-transfer runner

This is the operational entry point for the already-authorized `11952` warm
roster (41, 72, 186). The original generic-start control remains the independently
replayed 17-to-17 finite no-gain experiment. No mathematical status is promoted
by this repair, and finding a larger rank is not guaranteed.

## Run

From the repository root:

```sh
git pull --ff-only
cd research
python3 -m unittest discover -s elliptic-curves/tests -p 'test_v3_warm_pipeline.py' -q
python3 elliptic-curves/cas/run_v3_warm_start_overnight.py resume --hours 10
```

`resume` returns after launching one detached controller. The wall budget is a
limit, not an estimate. The original per-search/per-replay limits and numerical
V3 selector stay unchanged. A new session records the exact implementation and
input hashes. A case-local wall/RSS stop is retained and permits the remaining
cases to be considered; a proof, replay, input-integrity, or software error
stops the roster. There is no automatic retry loop or budget escalation.

```sh
python3 elliptic-curves/cas/run_v3_warm_start_overnight.py status
python3 elliptic-curves/cas/run_v3_warm_start_overnight.py diagnose
```

Status includes the actual phase (`PREFLIGHT`, `SEARCHING`, or `REPLAYING`),
PID/start-token identity, heartbeat, supervisor record, and exact current
`worker_log` path. `diagnose` prints a bounded tail of that log. The controller
log remains `artifacts/local/elliptic-curves/v3-warm-start-overnight-v1/autorun.log`.

To stop an owned new controller safely:

```sh
python3 elliptic-curves/cas/run_v3_warm_start_overnight.py stop
```

The stop command verifies PID identity and uses a Linux pidfd; it never sends
signals to an arbitrary stale PID. Resume later with the same command above.
A legacy live controller must be inspected explicitly; it is never auto-killed.

## Recovery of case 41

Leave every original file in place. Do not run the legacy `repair-resume`,
remove directories, regenerate the parent bank, rename chart files, or edit
certificates. The maintained runner reuses the existing seed/orbit inputs and
checks their hashes.

Its sequence is: real Sage seed/map preflight, structurally validate the old
sealed terminal, independently replay it in numeric chart order, then start
72 and 186. The 1,508 search charts are not executed again. Old failed logs and
supervisor records remain evidence; new attempts use unique session directories.

New independent replay work is checkpointed under each case's
`warm-verification/`. A resumed replay reuses a check only after its input and
checker hashes match. Completion creates `warm-replay.json` and
`warm-verified.json`; original `full-replay.json`, `verified.json`, the native
control and the 302 calibration records are not overwritten.

## Reproduced faults and fixes

1. `:03d` is minimum padding, not a fixed-width number. Sorting 1,508 chart
   filenames lexically makes index 101 read `chart-1000.json`. The checker then
   correctly rejects its centre against schedule entry 101. Numeric ordering,
   contiguous-index checks and explicit filename/payload/schedule checks fix
   traversal without dropping any assertion or changing any search evidence.
2. The old seed comparison equated serialized string coordinates with
   `Fraction` coordinate tuples. That predicate is false even when their exact
   rational values agree. Normalize both sides, retaining order and sign.
   The earlier cache/dispatch explanations did not establish this root cause.
3. Search transcripts bind an entire `MWState` content key, not just a list of
   mathematically equal points. Replay restores the historical state record,
   verifies its exact finite-reduction witnesses and checks its marked basis.
4. The controller is now standard-library-only. Numerical imports occur only
   inside explicit `sage -python` workers. No source-string replacement,
   `exec`-patched `run_case`, or transfer-wrapper monkeypatch chain remains in
   this active path. Frozen numerical modules themselves are reused.
5. A lifetime file lock, PID start tokens, zombie checks, live-worker checks,
   atomic status updates and explicit replay phases replace stale status and
   launch races. SIGTERM cleanup uses the existing bounded process supervisor.
6. A seeded epoch commits its exact state before publishing any chart. Resume
   retains completed charts, rebuilds only missing audit work, preserves partial
   artifacts, and accounts for completed chart slots. A sealed search always
   advances to replay rather than starting over.
7. Independent replay checks mod-3 and mod-5 on no-gain epochs as well as gains,
   binds them to the actual final cloud, rejects incomplete audits, and checks
   the frozen chart height, timeout and backend hash.

## Tests and limits of validation

`test_v3_warm_pipeline.py` contains 32 deterministic Python tests, including the
exact 1,508-file ordering failure, 82/1,000/1,001/1,508/4,096 boundaries, altered
centres, missing and duplicate indices, seed representations, immutable writes,
Sage subprocess arguments, inherited locks, reused/zombie PIDs, stale status,
crashes between chart/cloud/stage publication, and sealed-search recovery.
The arithmetic in controller fault-injection tests is explicitly mocked; these
are not mathematical certificates. The tests and compilation pass in the
coding sandbox. Full native Sage/PARI and the user's retained 1,508-chart
replay were not executed there; the real seed-to-chart preflight and full
independent replay run on the research machine before any result is accepted.

## Search meaning

The numerical method remains the frozen generic-shell-8/10, full-extension,
diversified V3 selector. New independent points cause immediate subgroup
rebuilding; returned points require exact certificates. The current campaign
extends existing curves, not newly chosen parameters. A gain is exported as
`warm-discovery.json` with the curve, certified points, and proof binding, but
is not automatically called a catalogue record or submitted anywhere.

No taller boxes, extra anchors, fresh parameter sweep or retrospective target
points have been added to disguise an engineering repair as a research gain.
Any later deepening/population experiment needs its own explicit fixed plan;
the null results here remain bounded exposures, not rank upper bounds or a
proof that the parent-specific method cannot transfer.
