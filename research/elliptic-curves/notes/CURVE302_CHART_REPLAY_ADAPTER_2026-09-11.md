# Curve302 historical chart replay adapter

This adapter completed a local historical replay and three retrospective
chart-exposure analyses. It has no current claim in
[MATH_STATUS.json](../../MATH_STATUS.json), and the ignored local raw transcript
and outputs are not portable evidence. It is not a current work queue or a
request to rerun the adapter.

The recorded replay normalizes the historical chart ledger without inferring
exposure from quotient height. It binds seed and epoch to the complete path,
retains only the maximal cumulative chart snapshot at each stage, deduplicates
recorded points, and promotes an exposure only after an exact group-equality
check. Unbound point evidence remains `UNKNOWN`.

The full historical interface and its commands are preserved byte-for-byte in
the [elliptic archive](../../archive/elliptic-curves/notes/CURVE302_CHART_REPLAY_ADAPTER_2026-09-11.md.txt).
The retained regression code is
[`run_curve302_chart_replay_adapter.py`](../cas/run_curve302_chart_replay_adapter.py).
Any new replay requires a separately scoped reason, immutable inputs, and a
portable evidence plan; it must not automatically continue into downstream
experiments.
