# Curve302 historical chart replay adapter

## Purpose

The real seeded-V3 transcript schema does not store quotient words in each chart.
Instead, cumulative `charts` snapshots contain chart definitions with
`centre/index/mapping/search`, while separate historical files include MW-state
and recorded rational-point evidence. The original chart-exposure controller
therefore stopped correctly at `UNKNOWN_UNSUPPORTED_OR_AMBIGUOUS_CHART_SCHEMA`.

This adapter reconstructs the normalized `curve302-chart-exposure-ledger.v1`
without inferring exposure from quotient height.

## Exact reconstruction

1. Seed and epoch are bound from the complete relative path, including
   `.../<seed>/replay-M17/epoch-NN/...`.
2. The raw tree is parsed once. For every historical stage only the maximal
   cumulative `charts` snapshot is retained.
3. `elliptic-curves.mw-state.v1` snapshots are used to derive a common set of
   seventeen rational points and exactly one additional seeded point for each of
   the fourteen exceptional direction ids.
4. The Weierstrass coefficients are inferred exactly from those rational points
   and every basis point is checked against the resulting equation.
5. `elliptic-curves.recorded-point-mod2-rank.v3` point records are joined to
   charts by stage plus chart/index evidence. Points embedded in or referenced
   by a chart's `search` payload are also accepted as positive evidence.
6. Each distinct recorded rational point is recognized in the displayed rank-31
   basis by a high-precision Neron--Tate coordinate solve. The proposed integer
   coordinates are then verified by exact Sage elliptic-curve group equality.
   Only after that exact check is its primitive `D/M17` word promoted to an
   exposure.
7. The generated ledger is normalized and validated by the existing consumer.
   All 180 actual acquisitions must occur in their own historical gain stages or
   the adapter stops.

Unbound recorded-point files do not become zero exposures. Their stages have
completion downgraded to UNKNOWN, so absence-based multiplicity and ordering
claims remain gated. Positive exact hits are retained.

## Runtime shape

The raw transcript tree is read once during `run`. The expensive arithmetic is
then proportional to the number of *distinct recorded rational points*, not the
1.27 million cumulative chart records. Neron--Tate coordinates are attempted at
128, 256, 512, then 1024 bits and exact group equality is required.

Progress is emitted as:

```
CURVE302_REPLAY_SCAN|...
CURVE302_REPLAY_RECOGNITION|unique=25|state=IN_D
...
```

## Commands

First use the cheap structural audit if desired:

```bash
sage -python research/elliptic-curves/cas/run_curve302_chart_replay_adapter.py audit
```

The real pipeline is:

```bash
sage -python research/elliptic-curves/cas/run_curve302_chart_replay_adapter.py run
```

A successful replay automatically launches the existing exposure census,
norm-matched multiplicity test, and stage-local ordering counterfactuals. Outputs
are placed under:

```
research/artifacts/local/elliptic-curves/curve302-chart-replay-v1/
```

The exact normalized ledger is `chart-exposure-ledger.json`; downstream results
are under `experiments/`.

To stop after producing the exact ledger:

```bash
sage -python research/elliptic-curves/cas/run_curve302_chart_replay_adapter.py run --no-experiments
```

After a complete run:

```bash
sage -python research/elliptic-curves/cas/run_curve302_chart_replay_adapter.py check
```

The downstream check recomputes the three analyses byte-for-byte from the frozen
ledger. The adapter itself never starts a V3 point search or modifies an existing
campaign.
