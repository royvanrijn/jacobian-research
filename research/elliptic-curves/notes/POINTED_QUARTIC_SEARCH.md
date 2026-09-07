# Shared pointed-quartic search

`PointedQuarticSearch` is the active half-lattice search API for every rank.
It uses the existing GMP modular worker and exact denominator/lattice
arithmetic. The [control regression](../../artifacts/generated-results/elliptic-curves/universal_pointed_control_regression_v1.json)
re-sieves all **1,034** frozen calibrated MW16 boxes, reproduces every square
hit and mapped point, and rechecks the integral relations and finite-reduction
certificates for **14/12/11/10/8 = 55** quotient directions.

This consolidates the implementation. It does not calibrate sensitivity on
MW17 or MW18 and does not add a prospective gain. The previously completed
104-fibre, 856-box prospective null result remains the result recorded in
[MATH_STATUS.json](../../MATH_STATUS.json) and the
[sensitivity proof note](MW16_SENSITIVITY_RECOVERY_2026-09-05.md).

## API and exact coordinates

For Python imports, add `elliptic-curves/cas` to `PYTHONPATH`.

```python
from pointed_quartic_search import PointedQuarticSearch

search = PointedQuarticSearch(
    curve=[0, 0, 0, -1, 1],
    subgroup=[{"x": "0", "y": "1"}],
    centre={"coefficients": [1]},
    coordinate_policy={"kind": "metric", "weight": "16",
                       "matrix": [1, 0, 0, 1]},
)
result = search.search(height=100, seconds=2, checkpoint_dir="/tmp/pointed-demo")
```

The implementation is [pointed_quartic_search.py](../cas/pointed_quartic_search.py).
The input curve is `[A,B]` or the five rational Weierstrass coefficients.
Rational completion of the square and cube transports general models to
short form, with an exact inverse for returned points. No minimal-model
computation, quartic reduction, or factorization is required.

Subgroups have arbitrary finite length, including zero with an explicit
centre. Their generators need not be independent or saturated. A centre is
an integral coefficient vector, an explicit finite point, or both; supplying
both requires exact equality. Centre selection remains in the caller:
generic CVP, adaptive recovered-lattice CVP, binary residues, and rescue
classes all enter through the same argument.

`metric:16` is the default inherited from the MW16 calibration. `gauss`
selects weight one; `raw` selects the original chord slope. A rational
invertible 2-by-2 matrix acts after the selected chart. For example,
`gauss:1,2,1` selects `z=(s+1)/2`. Every coefficient and ordinate scaling is
checked exactly. Different coordinate policies define different finite
boxes. PARI's `hyperellred` is available only in the explicitly selected
historical MW16 regression comparator.

The finite box is `|n| <= H`, `1 <= d <= H`, `gcd(n,d)=1`, with the transformed
point at infinity checked separately and both square-root signs transported.
The known original pointed endpoints `O,Q` are recorded separately from
discoveries. `denominator_start` and `denominator_end` permit disjoint shards.
The worker cap is in seconds per shard; centre construction and exact replay
are outside that cap. Heights are bounded by `10^6`.

Checkpoints bind the curve, every subgroup point, centre, coordinate policy,
height, denominator interval, worker budget, and source hashes. A checksum
and exact chart/hit/map replay validate cached results. Only a completed
interval is reused; interrupted intervals are retried. A timeout remains
incomplete and cannot imply absence of points or an upper bound on rank.
Completed shards survive interruption of a larger campaign. Completeness of
a saved transcript relies on the pinned worker; map replay is not a second
enumeration. Quotient independence is certified by the campaign, outside
this search service.

## Active callers

| Caller | Shared route |
|---|---|
| MW16 sensitivity and Nagao finalists | `checkpoint` / `CheckpointedBackend` adapters |
| MW17 jump-v2 and R17 refresh/CRT | `run_quartic_search` adapter; the old raw implementation is removed |
| MW17 zero-gain rescue | same jump-v2 adapter for generic rescue and adaptive slots |
| Curve 385, curve 398, curve 400, parent ladder and A1 parameter searches | common lattice helper's search slot; sparse restart semantics retained |
| Anchored MW18 specializations | MW18 input adapter in `run_pointed_quartic_search.py` |
| Future families or explicit points | the same API or generic JSON job manifest |

The [integration canary](../../artifacts/generated-results/elliptic-curves/universal_pointed_integration_v1.json)
executes the curve-specific, MW17, zero-gain and CRT adapters on a real MW17
chart with PARI search calls forbidden. Unit checks additionally run a real
MW18 specialization through the generic job runner. These are small
integration checks, not sensitivity calibrations.

`run_pointed_quartic_search.py --input jobs.json --output results.json`
executes a declared manifest with schema `elliptic-curves.pointed-quartic-jobs.v1`.
Each entry in `jobs` supplies `id`, `curve`, `subgroup`, `centre`,
`coordinate_policy`, `height`, and `seconds`; denominator interval fields
are optional. `--max-new` limits new chart calls. The output contains full
replay records, with no inferred quotient score.

For MW18, use `--mw18-candidate ID --centres centres.json --height H
--seconds S --export-jobs jobs.json`; `centres.json` is a list of centre
specifications. This reads the existing exact eighteen-section specialization
certificate. Adding `--output results.json` executes the declared jobs.
Neither mode automatically starts the 178-fibre campaign.

Frozen MW17 and sparse protocols still pin their historical source bytes.
The [runtime amendment](../cas/pointed_quartic_migration.py) verifies those
pins at the preserved revision, then records the current implementation and
coordinates separately. Candidate populations, centre-selection rules and
rank certificates retain their checks. New checkpoint directories and
runtime validation prevent historical/current search results from mixing.

## Regression controls and replay

The older implementations remain solely for historical controls:

- `half_lattice_fake_descent_replay.sage`: historical PARI minimize/reduce/search;
- `half_lattice_direct_reduction.py`: historical direct-reduction comparator;
- `half_lattice_pointed_sieve.py`: immutable v1 regression entry points **and
  the shared exact arithmetic/kernel primitives** used by the new API;
- `mw16_sensitivity_backend.py`: historical coordinate/search comparator,
  selected only by `run_mw16_sensitivity.sage --regression-controls`;
- old calibration, height-compression, ablation and Nagao/Fermigier point
  scripts referenced by historical certificates: regression material.

Unpointed Selmer-cover solubility and higher-genus parameter searches solve
different input problems; without a known rational centre they cannot use
this API. They are not alternative half-lattice backends.

The complete pre-migration source tree is pinned at
`d30a742133f0658185c3bd4c99f0b0f815f2f74b`. Historical full campaign commands
use that revision, including `replay_mw16_sensitivity.sh`. Current runners
use the shared backend; their finite boxes must not be described by the old
PARI protocols. MW16's evidence bundles also retain the executed source bytes.
`replay_pointed_quartic_snapshot.py` extracts them into an isolated directory
and runs the original checker without weakening source equality checks.

```sh
python3 -m unittest discover -s elliptic-curves/tests -p test_pointed_quartic_search.py
sage -python elliptic-curves/cas/check_pointed_quartic_integration.sage
python3 elliptic-curves/cas/replay_pointed_quartic_controls.py \
  --output artifacts/local/elliptic-curves/pointed-quartic-search/controls.json
python3 elliptic-curves/cas/replay_pointed_quartic_snapshot.py \
  --bundle artifacts/generated-results/elliptic-curves/mw16_sensitivity_controls_v1.json.gz \
  --summary artifacts/generated-results/elliptic-curves/mw16_sensitivity_controls_summary_v1.json
```

The control command checks all 1,034 boxes and recovered groups, using
content-addressed checkpoints. It does not rerun the 104 prospective fibres.
