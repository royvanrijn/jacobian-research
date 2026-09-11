# Curve302 historical chart-exposure controls

Date: 2026-09-11

## Question

The seeded Curve302 V3 trajectories acquire many different quotient vectors, yet their
integral prefix lattices recover the same growing core. Static quotient height is a
real bias, but complete short-vector controls do not reproduce that convergence.
This experiment asks whether the missing bias is already visible in the historical
pointed-quartic chart transcripts.

This is a **read-only retrospective analysis**. It launches no point search and never
propagates a counterfactual gain into a later historical V3 stage.

## Evidence boundary

A chart is counted as exposing a quotient direction only when its historical audit or
replay contains an **explicit 14-dimensional quotient word** attached to a result
record (or when the same information is supplied in the normalized ledger schema
below). A nominal chart direction, a quotient-height guess, or a point with no exact
quotient coordinate is not promoted to exposure evidence.

This boundary is intentional. If the existing raw transcript is point-only, `prepare`
stops with `UNKNOWN`-style failure evidence and writes `SCHEMA_PROBE.json`; it does not
turn the absence of an adapter into zero exposure. In that case the next action is a
narrow, independently checked replay/export adapter that binds those points to fixed-D
quotient coordinates.

A second gate concerns **absence**. A chart may contribute positive exposure evidence
whenever an exact quotient word is recorded, but a missing candidate is only usable as
a negative comparison when that chart explicitly carries `complete=true` (or an
equivalent complete/exhaustive flag recognized by the adapter). Norm-matched
multiplicity and stage-local ordering controls are promoted only when every chart in
that historical stage has complete coverage. Otherwise the stage remains
`UNKNOWN_INCOMPLETE_CHART_COVERAGE`; recorded positive hits are still retained in the
census.

## Candidate population

For every historical stage the frozen candidate population is

1. the first 1,000 primitive directions in the completed exact static-height
   enumeration; plus
2. every actually acquired primitive quotient direction outside that prefix.

Directions already rationally contained in the current historical prefix are removed.
The actual acquisition is required to remain in the population and to occur in the
historical exposure evidence for that stage.

Static-height matching uses the same frozen decade rank bands as the previous control,
with complete norm shells retained by the sealed vocabulary implementation.

## Experiment A — exact historical exposure census

For every historical stage and candidate direction `v`, record:

- number of historical charts exposing `v`;
- first chart order exposing `v`;
- minimum explicit parameter/slope height, when present;
- minimum quartic coefficient bit size among exposing charts;
- exposing chart IDs.

`L2`, `L1`, and `L4` are additionally recorded even when they lie outside the actual
acquisition's static-rank control band.

The output is `exposure-census.json`.

## Experiment B — norm-matched representation multiplicity

At each stage compare the actual acquisition with all still-unknown candidates in its
frozen static-rank band. Report tied descriptive percentiles for:

- exposure count (higher is better);
- minimum parameter height (lower is better);
- minimum quartic coefficient bit size (lower is better).

These are retrospective descriptive ranks. They are not p-values and are not
prospective rank-success probabilities.

The output is `multiplicity.json`.

## Experiment C — stage-local scheduling counterfactuals

Reorder only the already-generated charts of each historical stage under:

- original order;
- reverse order;
- ascending quartic coefficient bit size;
- 256 frozen random permutations;
- 256 frozen random permutations within original score bands, when every chart has a
  score-band label.

The first chart in that ordering that exposes any still-unknown frozen candidate is the
stage-local counterfactual gain. For that gain record, exactly over `Z`:

- whether it equals the historical acquisition;
- whether adjoining it is a primitive lattice extension;
- whether it contains the observed next common integral core;
- whether it contains `L2`, `L1`, or `L4`.

The calculation is deliberately stage-local. A counterfactual gain does not define a
valid later historical state, so it is never propagated.

The output is `counterfactuals.json`.

## Normalized ledger escape hatch

If the historical raw JSON layout is not recognized but an exact replay can export
quotient words, pass `--ledger FILE` with schema:

```json
{
  "schema": "curve302-chart-exposure-ledger.v1",
  "direction_ids": ["recovered-local-01", "... 13 more ..."],
  "runs": [
    {
      "seed": "recovered-local-01",
      "stages": [
        {
          "epoch": 0,
          "charts": [
            {
              "chart_id": "...",
              "order": 0,
              "score_band": "...",
              "quartic_coefficients": ["..."],
              "search_bound": 100000,
              "exposures": [
                {
                  "quotient_word": [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  "parameter_height": 37,
                  "within_bound": true
                }
              ]
            }
          ]
        }
      ]
    }
  ]
}
```

Every run and epoch must match the sealed trajectory roster exactly. The total chart
count must match the historical trajectory audit. An actual acquisition missing from
its stage exposure evidence is a hard failure.

## Run

From the repository root:

```bash
python3 research/elliptic-curves/cas/run_curve302_chart_exposure.py probe
```

Use `probe` first if the raw transcript schema is uncertain. Then run:

```bash
python3 research/elliptic-curves/cas/run_curve302_chart_exposure.py run
python3 research/elliptic-curves/cas/run_curve302_chart_exposure.py check
```

Optional explicit sources:

```bash
python3 research/elliptic-curves/cas/run_curve302_chart_exposure.py run \
  --source /path/to/curve302-short-vector-core-results \
  --structure /path/to/curve302-closure-structure-results \
  --raw-root /path/to/historical-seeded-v3-evidence
```

Or use `--ledger /path/to/chart-exposure-ledger.json`.

The default output folder is
`research/artifacts/local/elliptic-curves/curve302-chart-exposure-v1/`.

## Fail-closed behavior

- Inputs, selected raw evidence and code are hash-frozen before stages run.
- Historical roster, epoch count and chart count must agree exactly.
- Point-only/unsupported evidence is diagnostic, never silently treated as no exposure.
- Interrupted or unsealed stages require a fresh output folder.
- Each stage has CPU and address-space limits.
- `check` recomputes all three analyses from the frozen inputs and requires byte-for-byte
  output agreement.

A completed run is still only retrospective model checking on one known displayed
`D/M17`. It is not a propagation theorem and not yet a prospective rank selector.
