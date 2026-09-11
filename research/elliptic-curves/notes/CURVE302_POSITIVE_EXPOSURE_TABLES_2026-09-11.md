# Curve302 positive-evidence tables

This post-processing analysis consumes the sealed exact chart replay and extracts four tables without using chart-completeness assumptions.

The analysis never interprets an absent chart hit as non-exposure. All comparisons are among explicitly recorded positive exposures.

## Targets

The named axis targets are:

- `L2 = recovered-local-02`;
- `L1 = recovered-local-01`;
- `L4 = recovered-local-04`.

It also identifies the first observed common integral cores of ranks five and six from the sealed intrinsic/observed filtration bundle (`C5`, `C6`) rather than hard-coding their generators.

For a current historical prefix `S` and positively exposed primitive direction `v`, hypothetical saturation uses

`Sat_Z(S + Z v) = span_Q(S,v) intersect Z^14`.

Thus an integral target is contained after saturation exactly when its generators lie in the rational span of `S` and `v`.

## Tables

1. `lead-times.json`: first direct positive exposure, first saturation-enabling positive exposure, first positive core-progress exposure, and first actual integral containment for each target in each run.
2. `precontainment-multiplicity.json`: chart/epoch lower bounds for positive target exposure before actual containment.
3. `coexposure-choices.json`: every explicitly exposed rationally-new direction at each historical gain stage, with the actual gain batch marked inside that positive set.
4. `saturation-impact.json`: exact saturation impact of every positively exposed candidate on the next observed common core and on `L2/L1/L4/C5/C6`.

The output `SUMMARY.md` prints compact aggregate versions of all four tables; CSVs contain the detailed row data.

These are retrospective descriptive tables on one known displayed quotient. They are not significance tests, complete choice-set comparisons, prospective rank selectors, or propagation theorems.
