# Curve302: completed twelve-seed universality panel

The frozen V3 policy reached certified rank31 from **11 of the 12 additional
known rank18 seed subgroups**. All cases completed independent replay.
The exception, `recovered-local-01`, reached rank29 and then exhausted its
final schedule without gain. Including the two previously verified winners
gives **13/14 tested seeds reaching31, and all14 reaching at least29**.

The panel is complete; no worker remains active. The
[exported result](../../artifacts/generated-results/elliptic-curves/curve302_seed_universality_panel_v1.json)
retains the twelve original replay summaries and all stage records.

## Frozen order and results

Each case starts from generic M17 plus one supplied known exceptional point,
independently certified as rank18. The two earlier winners,
`recovered-strict-02` and `recovered-strict-03`, are excluded from this panel.
The twelve seeds follow the already-frozen retrospective bottleneck ranking;
search outcomes do not select or reorder them. Each case retains the original
V3 numerical policy, finite stopping conditions and independent replay.

| Frozen position | Supplied seed | Final certified lower bound | Charts |
| ---: | --- | ---: | ---: |
| 1 | residual-strict-03 | 31 | 385 |
| 2 | recovered-strict-01 | 31 | 1513 |
| 3 | recovered-local-04 | 31 | 975 |
| 4 | residual-strict-05 | 31 | 422 |
| 5 | residual-strict-07 | 31 | 1034 |
| 6 | recovered-local-01 | 29 | 2430 |
| 7 | recovered-local-02 | 31 | 705 |
| 8 | recovered-local-03 | 31 | 568 |
| 9 | residual-strict-01 | 31 | 975 |
| 10 | residual-strict-02 | 31 | 1618 |
| 11 | residual-strict-04 | 31 | 525 |
| 12 | residual-strict-06 | 31 | 2056 |
| | Panel total | | 13206 |

All eleven rank31 cases stop with `TARGET_LOWER_BOUND_REACHED`.
All twelve new runs have zero censored charts.
The [earlier two-seed runs](CURVE302_SEEDED_V3_RESULTS_2026-09-08.md)
used855 and637 charts; the combined fourteen-seed total is **14698 charts**.

## The exception is a completed bounded stall

For `recovered-local-01`, the search certifies 18→19→…→29, then stops
with `COMPLETE_FINITE_NO_GAIN`. Its final rank29 epoch executes all
**1434 scheduled charts**, after scoring **131072 retained-anchor extension
cosets** and performing **1672 exact-CVP refinements**. Its final independent
mod-3/mod-5 certificates both give lower bound29; the full run uses2430 charts.

The rank29 endpoint describes the recovered subgroup. The ambient curve
is still the known rank-at-least31 curve302. This frozen seed/policy combination
does not close to31 within its allowed schedule. It does not exclude other
representatives or further directions.

## Interpretation

Across all fourteen tested labels, recovered-local seeds reach31 in3/4 cases,
recovered-strict seeds in3/3, and residual-strict seeds in7/7.
These counts support broad amplification on302 once an exceptional direction
is supplied. They do not prove that every exceptional point, representative,
or starting subgroup behaves this way. One local exception is insufficient
to establish a general local/strict distinction. Predictive correlation
between the bottleneck metric and realized chart cost is not established here.

The seeds are retrospectively supplied known points. This does not solve
prospective first-seed acquisition, establish transfer to other curves,
produce a new curve or rank record, compute a conductor, or prove exact rank.

## Scheduling history and evidence

The run began sequentially. After one case completed and the second search
sealed, the [four-worker continuation](CURVE302_SEED_PANEL_PARALLEL_2026-09-08.md)
parallelized independent cases in the original dispatch order. The second
case's partial independent replay restarted; its sealed search was reused.
Each case used the original search and full replay functions. Completion
order could differ; this table remains in the frozen roster order.

The [provenance export](../../artifacts/generated-results/elliptic-curves/curve302_seed_universality_panel_provenance_v1.json)
pins63 retained evidence files and copies the complete68-file scheduling
manifest. On2026-09-09, every manifest binding matched; each exported result
equalled its independent replay receipt, protocol and terminal hashes matched,
stage chart totals agreed, and all stored mod-3/mod-5 stage ranks agreed
with the certified stage bounds. These receipt checks do not repeat the
completed arithmetic replays.

Raw evidence remains under
`research/artifacts/local/elliptic-curves/curve302-seed-universality-panel-v1/`:
`summary.json`, one seed directory per case, and `parallel-v1/` scheduling
records. Each seed retains its input/proof/protocol, chart transcripts,
terminal record and `seeded-verified.json`.

- [Original panel](../cas/run_curve302_seed_universality_panel.py).
- [Original independent case replay](../cas/run_curve302_seeded_v3_amplifier.py).
- [Parallel scheduler](../cas/run_curve302_seed_panel_parallel.py).

From `research/`, read the completed status with:

```sh
python3 elliptic-curves/cas/run_curve302_seed_panel_parallel.py status
```
