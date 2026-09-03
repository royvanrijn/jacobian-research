# AGENTS.md — elliptic curves over `Q`

This programme inherits the repository rules. `../MATH_STATUS.json` remains the sole mathematical-status authority.

## Programme state

**PAUSED since 2026-09-03.** Compute is reserved for prime-gap calculations.

Do not launch rank-32 searches, Selmer/descent campaigns, family sweeps, high-height point searches, or K3 specialization scans unless this programme is explicitly resumed. Maintenance, literature work, and small exact regressions are fine.

## Frozen targets

- Rank-at-least-31 is certified for ICARM curve 302; no unconditional exact-rank upper bound is known.
- Rank-at-least-30 is certified for curve 273.
- Rank-at-least-21 / low-conductor branches of the original operational target are replayed.
- The next record target remains rank `>=32` or a stronger unconditional exact-rank/conductor result.
- A point list is not a rank lower bound until independence is certified. Selmer dimension, analytic estimates, Nagao scores, and bounded-search survival are not rank lower bounds.

## If work resumes

1. start from [`README.md`](README.md), `../MATH_STATUS.json`, and the canonical curve notes;
2. require a genuine residual descent/Selmer gate before expensive point search;
3. for alternate Q80, build native calibration fibres first — published R17 rank-25--28 controls do not occur at rational alternate parameters;
4. minimize every specialization and re-prove point independence;
5. pin all external CAS inputs, outputs, software versions, limits, and failure semantics.

## Evidence discipline

- Exact rank requires matching unconditional lower and upper bounds.
- Conductor means the exact conductor of a global minimal model, not a discriminant radical.
- Preserve raw search/checkpoint data under ignored local-artifact paths and compact certificates under `../artifacts/generated-results/elliptic-curves/`.
- Bounded misses remain bounded experiments.
- Do not promote Selmer classes to Mordell–Weil directions without the global argument.
- Preserve scripts and tests even when their search campaign is archived; they may remain useful regressions.

For K3 construction context use [`../elkies-k3/README.md`](../elkies-k3/README.md) and [`../elkies-k3/AGENTS.md`](../elkies-k3/AGENTS.md). The historical raw q12/Q80 routes are not active starting points.
