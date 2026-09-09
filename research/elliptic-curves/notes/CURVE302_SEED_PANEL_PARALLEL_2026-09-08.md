# Parallel scheduling of the frozen twelve-seed panel

**Completed 2026-09-09:** all twelve cases independently replayed, eleven
reached31, and recovered-local-01 stopped at29 after a complete no-gain
epoch. Total13206 charts, with zero censored charts. No worker remains active.
See the [canonical results and evidence](CURVE302_SEED_UNIVERSALITY_RESULTS_2026-09-09.md).
The takeover and process references below document the completed run.

The original controller spent each case's search and independent replay on one
process before starting the next case. The first completed seed's search phase
took 3937.978 seconds, of which 271.019 seconds were backend point searches;
independent final replay was additional. These measurements are recorded in
[the performance audit](V3_FUTURE_SEARCH_PERFORMANCE_2026-09-08.md).

The scheduling-only continuation ran up to **four independent Sage case
processes**. Available capacity at takeover was 32 logical CPUs and about
29 GiB available memory. Dispatch follows the original bottleneck ranking;
completion order may differ. Each case calls the original panel context,
`det1092_v3_worker.run_search` and the original `verify_case`. The seeds,
V3 sources, numerical policy, per-chart limits, full replay and case evidence
paths retain their original bindings. A separate 68-file manifest binds the
new scheduler and its inputs. It does not amend the old policy hashes.

Parallelism can affect wall-time censoring through resource contention. Four
workers is a conservative initial concurrency cap; every original timeout
and censor record is retained. A fourfold end-to-end speedup is not measured
or promised. The separately benchmarked integer-CVP optimization is not
installed in these frozen runs.

## Takeover and evidence

The serial controller PID3972817 was stopped during independent replay of
`recovered-strict-01`, after its search sealed rank31 in1513 charts.
Its partial replay restarted with the original checker; no search charts were
repeated for that sealed case. The verified `residual-strict-03` result
(rank31,385 charts) was reused with protocol/terminal hash checks.

Only the exact serial worker was signalled, using a Linux PID handle after
checking its command, process start token and sealed-search state. Its old
controller state and log are retained. The new detached controller started as
PID3191781. Subsequent PIDs are operational state, not durable identities.

All new scheduling evidence is under
`artifacts/local/elliptic-curves/curve302-seed-universality-panel-v1/parallel-v1/`:
`manifest.json`, `serial-state-before-takeover.json`, `state.json`,
`controller.log`, and one log per active seed. Original case certificates and
checkpoints remain in their original seed directories. Final summary results
are ordered by the original roster, regardless of finish order.

Exclusive controller and case locks prevent duplicate parallel workers.
The original status file advertises the live parallel controller, so the old
launcher refuses a concurrent resume. A failed child stops dispatch and the
remaining owned children; all finite search/checkpoint evidence is retained.

## Commands

The run is complete; no resume is needed. From `research/`, read its final status:

```sh
python3 elliptic-curves/cas/run_curve302_seed_panel_parallel.py status
```

Historical recovery reference: after reviewing an interrupted run, resume
with the same manifest and concurrency:

```sh
python3 elliptic-curves/cas/run_curve302_seed_panel_parallel.py resume --workers 4
```

Do not use the old serial resume command while the parallel continuation is
active. Its status tail still points at the historical serial log.

## Validation

Fourteen focused tests pass: seven original panel regressions and seven
scheduler checks covering dispatch order, bounded concurrency, out-of-order
completion, exclusive locks, PID reuse, receipt tampering, sealed-search reuse
and refusal to take over an unsealed search. A foreground Sage context check
validated the first four roster contexts at rank18 with zero charts and all
68 bound files. Live status then confirmed four separate processes and their
first landscape/replay output. The frozen mathematical experiment remains
retrospective known-seed amplification on302.
