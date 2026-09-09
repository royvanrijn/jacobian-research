# Unattended norm12 search

The campaign completed normally and is stopped. All four configured parent
banks exhausted their finite search coverage, with no rank gains. Work is
paused at the user’s request; no restart or new campaign is scheduled.

It finished at 09:28:09 UTC on 2026-09-09, before the current Linux boot
at approximately 09:40:19 UTC. The controller is no longer running.

| Curve | Final certified lower bound | Completed calls |
| --- | ---: | ---: |
| 188 / ICARM619 | 28 | 1,568 |
| 48 | 27 | 1,558 |
| 71 | 27 | 1,552 |
| 40 | 27 | 1,518 |

All 53 automatic batches completed and replayed in 101.3 minutes. Final
cumulative coverage is 6,196 calls, including 5,096 added by the controller.
The post-reboot integrity check verified all 53 exported-result hashes and
424 bound-receipt hashes; all 6,196 expected final chart files are present.
No loss was detected within these checks. This was a file-integrity check,
not a fresh arithmetic replay or a repository-wide data-loss audit.

[Completion and integrity receipt](../../artifacts/generated-results/elliptic-curves/unattended_norm12_v1/completion-and-integrity.json).
Exhaustion applies only to these parent banks and finite bounds; it does
not prove exact ranks or exclude unseen points.

## Controller and retained starting inputs

The pure Python controller made no model/API calls and required no AI turn
between batches. It ran as a detached local process, not a recurring Codex
automation, and does not automatically restart after a reboot.

| Job | Curve | Starting certified subgroup | Initial completed calls |
| --- | --- | ---: | ---: |
| job1 | 188 / ICARM619 | M28 | 600 |
| job2 | 48 | M27 | 200 |
| job3 | 71 | M27 | 200 |
| job4 | 40 | M27 | 100 |

The controller rotates jobs sequentially, one bounded worker at a time.
Each batch adds at most100 point calls, independently replays the result,
exports sealed evidence, and advances only that job's verified parent.
Previously completed coverage is inherited from exact bound receipts.
The retained checkpoint records the completed scheduling ledger. The old
suffixes have been consumed; do not relaunch them.

A new certified gain stops the campaign for mathematical review. If a
complete-cloud audit requests reconciliation, the controller replays the
saved witnesses and standalone certificate twice before stopping. It also
stops on any failed worker, changed source binding, unknown terminal
condition, exhausted coverage, a requested stop, or a resource limit.
No bounded miss proves an upper bound. Mathematical status is not promoted
automatically; new gains require review and registration in MATH_STATUS.

Configured limits:24 hours before starting another batch,100 batches total,
10GiB minimum free disk before each batch,3GiB process-tree RSS per worker,
1800 seconds per search/replay phase and120 seconds per reconciliation phase.
A current batch finishes its replay even if the campaign time limit or STOP
request arrives mid-batch. Existing map and point limits are unchanged.

## Status and stop

From the repository root:

```sh
python3 research/elliptic-curves/cas/run_unattended_norm12.py status --folder research/artifacts/local/elliptic-curves/unattended-norm12-v1
python3 research/elliptic-curves/cas/run_unattended_norm12.py stop --folder research/artifacts/local/elliptic-curves/unattended-norm12-v1
```

Status checks the PID's process-start token, not just a status file. Stop is
graceful, after the current batch and replay. The state file records the
current phase, accepted batches, result paths and latest verified parents.
If the process dies unexpectedly, status reports interruption; partial
outputs remain intact and are not automatically retried or overwritten.

- [Controller](../cas/run_unattended_norm12.py).
- [Final checkpoint](../../artifacts/local/elliptic-curves/unattended-norm12-v1/state.json).
- [Frozen configuration](../../artifacts/local/elliptic-curves/unattended-norm12-v1/config.json).
- [Controller log](../../artifacts/local/elliptic-curves/unattended-norm12-v1/controller.log).
- [Parent evidence](PRODUCTIVE_PARENT_SPAN_REASSESSMENT_2026-09-09.md).

Four focused tests pass: unchanged-budget versus exhausted scheduling,
gain/reconciliation precedence, unknown/decreasing-rank rejection and process
identity. The [first real batch](../../artifacts/generated-results/elliptic-curves/unattended_norm12_v1/integration-check.json)
passes search and independent replay at700 cumulative curve188 calls, M28
unchanged. The same live controller automatically advanced to curve48,
confirming the unattended handoff. Controller
sources are preserved in `controller-sources.zip`; every batch archives its
bound arithmetic sources and saves supervised search/replay receipts.

The existing Codex Goal loop is separate from this process. Use `/goal pause`
to prevent further AI continuation turns; pausing that goal does not stop
this detached Python worker. This slash command is documented in the
[official goal guide](https://learn.chatgpt.com/use-cases/follow-goals).
