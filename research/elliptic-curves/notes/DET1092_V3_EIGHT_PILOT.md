# Same-parent V3 pilot: eight frozen determinant-1092 fibres

## Decision and scope

Pause further adaptation to 11952. Its clean generic-start 17-to-17 control and
three independently replayed warm 27-to-27 exposures remain unchanged. They do
not prove saturation, exact rank, or impossibility of transfer.

This experiment tests the intermediate level: V3 execution on other fibres of
the same determinant-1092 parent as the successful 302 calibration. It reuses
the original 48 score-selected record-scale fibres. They are **new to this V3
policy**, not untouched by all previous searches. No fresh parameter population,
old point-search results, public exceptional points or catalogue ranks enter
selection or execution.

## Run

From the repository root:

```sh
git pull --ff-only &&
cd research &&
python3 -m unittest discover -s elliptic-curves/tests -p 'test_det1092_v3_trial.py' -q &&
python3 elliptic-curves/cas/run_det1092_v3_trial.py launch --hours 10
```

The script launches one detached controller and returns. Preparation and all
numerical work run under `sage -python`, not the controller's system Python.
`--hours` is a session resource budget, not an estimated running time. No agent
needs to poll continuously or decide the next chart.

```sh
python3 elliptic-curves/cas/run_det1092_v3_trial.py status
python3 elliptic-curves/cas/run_det1092_v3_trial.py diagnose
python3 elliptic-curves/cas/run_det1092_v3_trial.py stop
python3 elliptic-curves/cas/run_det1092_v3_trial.py resume --hours 10
```

`status` reports PREPARING / PREFLIGHT / SEARCHING / REPLAYING, PID start-token
identity, supervisor result, heartbeat and the exact worker-log path. `stop`
uses a verified Linux pidfd. `resume` preserves previous sessions and reuses
checked charts/replay checkpoints. It never erases or recomputes a completed
search just because replay was interrupted. A resource stop is not a negative
rank observation; a proof/replay/integrity failure stops the controller.

New data live under:

- `artifacts/local/elliptic-curves/det1092-v3-eight-pilot/`
- `artifacts/local/elliptic-curves/det1092-v3-eight-controller/`

All 11952 and 302 source files, certificates, logs and numerical modules are
left unchanged. The working warm-transfer checkpoint contracts, certified-state
constructor, finite audits and independent landscape checker are reused through
explicit function calls, without a wrapper monkeypatch chain.

## Exactly which eight?

The original `selection-result.json` is pinned to SHA-256
`1dbce6d7bc287176afeae198a60fe9cae0cc8c1db9c42ffe97f71cdeaeb6b869`, taken from the
retained intake manifest at `d440079b`. A byte-identical decompressed copy of the
tracked intake package is accepted when the local raw file is absent. A changed
local file is rejected, not silently replaced by a different selection.

The rule is fixed in code before any V3 outcome:

1. Four strongest across the 48 by the existing integer `score_units`, with
   denominator, numerator and ID as the original deterministic tie-breakers.
2. One central moderate-stratum candidate from each j-height band 10 and 11
   (index 2 of the four score-ordered moderate entries in that band).
3. One SHA-selected remaining candidate from each band, excluding the six
   already selected. SHA domain: `det1092-v3-eight-pilot-v1`.

The forty remaining rows receive a frozen SHA order too. The eight rows, their
roles and reserve order are published in `plan.json` **before any seed-rank
check**. There is no rank-dependent replacement. All eight seeds are certified
before `roster.json` releases any point searches.

These SHA controls are conditioned on the original score-selected population,
not uniform random elliptic curves. The experiment cannot estimate the overall
frequency of rank jumps from this sample.

## Exact inputs and unchanged numerical choices

The reduced parent and full generic orbit table are pinned to their retained
Git object identities. Preparation checks generic section polynomial identities,
the determinant-1092 Gram matrix, the orbit-table certificate binding, every
specialized curve equation, all point memberships, and finite independence of
its 17 ordered generic sections. Each worker repeats exact specialization and
certificate checks before installing its own-case artifact-read guard.

The numerical V3 module, low-shell anchor rule, Babai scoring, exact-CVP
shortlist, diversified finalists, tie handling, factor-free maps, chart height,
per-chart timeout, and chart/epoch budgets are not tuned. The experiment's target
is rank at least 32 rather than the calibration's rank-31 goal. Every certified
gain triggers the same immediate subgroup rebuild. The first complete no-gain
epoch ends that fibre. A timeout stays censored.

At M17 there is only **one extension per retained anchor**. The diversified
extension mechanism grows only after a new direction is found. Consequently the
first epoch measures bootstrap sensitivity; one should not interpret a generic
no-gain run as a direct test of all later cascade behavior. Coverage is complete
over extensions of the retained anchors, not over the whole M_r/2M_r chart atlas.

The full 1092 orbit table contains the zero class, unlike the filtered 11952
bank. Preflight deliberately picks a nonzero norm-8/10 entry and checks its
exact quartic map; it does not mistake the first table row for a usable centre.

## Results and decisions

Completed on the research host: **all eight cases independently verified
17→17**, with 82 charts each (656 total), no resource stops and no eligible
reserve expansion. The summary is retained at
`artifacts/local/elliptic-curves/det1092-v3-eight-pilot/summary.json`.
See the [completed comparison](CURVE302_SEEDED_V3_RESULTS_2026-09-08.md).

Every case requires independent landscape, numeric chart-order, exact historical
state, map/witness and mod-2/3/5 replay, including no-gain epochs. Reports are
`trial-replay.json` and `trial-verified.json`. A gain additionally produces
`trial-discovery.json` containing the equation, rational points and proof binding;
it is not automatically claimed as a new catalogue entry or a record.

After eight verified results, `summary.json` reports gains, censorship and whether
an unchanged forty-fibre follow-up is eligible. **This launcher never releases
those forty automatically.** The intended follow-up requires a positive verified
gain and completion of this initial panel; its ordering is already frozen.

The eight completed 17-to-17 results mean no recovery in these eight finite
exposures. They do not distinguish low jump incidence from inadequate
bootstrap visibility. Extra Selmer classes, if investigated later, would not
by themselves certify extra rational points either. Keep incidence, visibility
and global solubility separate, as in `RANK_JUMP_REASSESSMENT_2026-09-05.md`.

## Validation boundary

46 Python regression tests and compilation pass in the coding environment.
Tests cover selection-order invariance, exclusion of outcome annotations,
original population/hash enforcement, exact homogeneous evaluation, unchanged
numerical policy fields, zero-orbit preflight, 1,508-chart numeric ordering,
immutable checkpoints, duplicate/hole and centre-mismatch rejection, Sage argv,
PID reuse/zombies, lifetime locks, resource/failure classification, recovery of
sealed searches, and a mocked eight-case dispatch that never runs the reserves.
Mocked worker tests are **not** mathematical certificates. Native Sage/PARI and
the actual eight-fibre experiment were not executed in that environment; those
run with real preparation, preflight and independent replay on the research host.

No mathematical status entry is promoted by adding this experiment.
