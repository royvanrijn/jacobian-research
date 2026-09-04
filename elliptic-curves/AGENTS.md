# AGENTS.md — elliptic curves over `Q`

This programme inherits the repository rules. `../MATH_STATUS.json` remains the sole mathematical-status authority.

## Programme state

**ACTIVE.** The programme is open for theorem-directed breakthrough work.

Rank-32 searches, Selmer/descent campaigns, family sweeps, high-height point
searches, and K3 specialization scans must have an explicit mathematical gate,
declared limits, checkpoints, and a reproducible certificate plan.

## Current targets

- A target-free A1/MW16 parameter experiment is active alongside the other
  rank-jump experiments.  Known record equations, parameters, points, ranks,
  target `j`-invariants, and jump labels stay out of its selection and
  execution.
- Rank-at-least-31 is certified for ICARM curve 302; no unconditional exact-rank upper bound is known.
- Rank-at-least-30 is certified for curve 273.
- Rank-at-least-30 is certified for curve 398.  Its two recovered norm-eight
  A1/MW16 survivor labels on `X948` are exact base-equivalent presentations of
  one fibration; their specialized integral MW16 bases generate the same
  subgroup and do not provide two generic structures.
- Rank-at-least-21 / low-conductor branches of the original operational target are replayed.
- The next record target remains rank `>=32` or a stronger unconditional exact-rank/conductor result.
- R17/MW17 specialization, the A1/MW16 parameter experiment, and curve-302
  parent reconstruction are simultaneous peer paths.  Do not select a family
  by generic rank alone.
- A point list is not a rank lower bound until independence is certified. Selmer dimension, analytic estimates, Nagao scores, and bounded-search survival are not rank lower bounds.

## Breakthrough workflow

1. start from [`README.md`](README.md), `../MATH_STATUS.json`, and the canonical curve notes;
2. keep proof gates separate from search-budget gates: an unconditional
   certified rank/Selmer upper bound below the target excludes a fibre, while
   incomplete or conditional descent data only schedule finite, checkpointed
   search; a certified independent point subgroup proves its lower bound
   without waiting for descent;
3. use curve 398's duplicate survivor pair only as a mandatory
   `PGL2`/Weierstrass deduplication regression; run the A1/MW16 parameter
   experiment on one representative and do not put its known target or any
   other record in the experiment loop;
4. use the exact native alternate-Q80 calibration fibres from the complete ICARM sweep — especially rank-at-least-29 curve 12 on `11952` — while keeping the published-R17 rank-25--28 transfer misses as negative controls;
5. minimize every specialization and re-prove point independence;
6. pin all external CAS inputs, outputs, software versions, limits, and failure semantics.

## Evidence discipline

- Exact rank requires matching unconditional lower and upper bounds.
- Rank lower bounds from exactly verified independent rational points do not
  require a Selmer upper bound or a completed descent.
- Conductor means the exact conductor of a global minimal model, not a discriminant radical.
- Preserve raw search/checkpoint data under ignored local-artifact paths and compact certificates under `../artifacts/generated-results/elliptic-curves/`.
- Bounded misses remain bounded experiments.
- Do not promote Selmer classes to Mordell–Weil directions without the global argument.
- Preserve scripts and tests even when their search campaign is archived; they may remain useful regressions.

For K3 construction context use [`../elkies-k3/README.md`](../elkies-k3/README.md) and [`../elkies-k3/AGENTS.md`](../elkies-k3/AGENTS.md). The historical raw q12/Q80 routes are not active starting points.
