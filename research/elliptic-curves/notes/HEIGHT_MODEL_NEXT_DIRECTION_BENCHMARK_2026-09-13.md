# Three-arm test of certified coordinate coverage

The nine-arm run has finished. The intended search comparison remains
incomplete: all six bounded-model arms hit the cold-support preparation limit
before making a point call. Three baseline gains and twelve individual target
locations have independent exact replays. The
[completion packet](../../artifacts/generated-results/elliptic-curves/height_model_next_direction_v1/completion.json)
retains the costs and separate failure stages.

This is the bounded runtime gate following the
[first-chart height certificate](POINTED_CHART_HEIGHT_BOUNDS_2026-09-13.md).
It tests whether selecting equivalent models by a certified bound improves
CPU to the first independently certified direction. It changes neither the
bound algorithm nor the local-neighbour vocabulary.

The frozen inputs are the retained Curve302 M27, M29 and M30 subgroups and
their original first 256 centres, in order. These are retrospective recovery
controls for 27→28, 29→30 and 30→31. They are not three independent curves.
No withheld point or winning-chart index enters construction or selection.

## Results and residual

| Control | Baseline CPU through independent gain | Best-single arm | Pair arm |
|---|---:|---|---|
| Curve302 M27→M28 | 49.965 s | Cold-support timeout | Cold-support timeout |
| Curve302 M29→M30 | 30.588 s | Cold-support timeout | Cold-support timeout |
| Curve302 M30→M31 | 59.033 s | Cold-support timeout | Cold-support timeout |

Every timeout hit the declared 30-second support limit, inside the five-minute
arm allowance. Including startup and starting-subgroup replay, the six failed
arms cost 32.634–32.829 CPU seconds each. All starting and successful enlarged
subgroups passed both finite-group implementations. The baseline made 143
point calls; the bounded arms made zero. The nine arms used **336.107 CPU
seconds** and **355.894 elapsed seconds**. Diagnostic attempts, including their
failures, used a further **9.545 CPU seconds** and **9.997 elapsed seconds**.
Input/source sealing used 0.441 CPU seconds separately. These measurements do
not include all development and test execution and are not a whole-project cost.

Thus this cold implementation has not paid its preparation cost. This does
not establish that the single model or pair searches slowly: their search
stage was never reached. Full discriminant factorisation is this frozen
implementation's preparation policy, not a mathematical prerequisite for
every possible valid height bound. The safe unfactored Bezout fallback from
the original theorem remains valid. No larger factorisation allowance, warm
search replacement, revised selector or follow-up point search was launched.

The post-seal individual-model calculation did complete. At the **baseline's
winning anchor in each case**, the literal winning point has these exact
parameter heights:

| Control | Factor-free | Full minimal | First neighbour | Second neighbour |
|---|---:|---:|---:|---:|
| M27→M28, centre index 50 | 40,774 | 33,301 | 33,301 | 47,231 |
| M29→M30, centre index 28 | 15,550 | 4,692 | 12,467 | **1,781** |
| M30→M31, centre index 62 | 39,979 | **35,191** | 36,388 | 41,176 |

The neighbour primes are respectively `(7,7)`, `(5,7)` and `(5,11)`.
The best uniform single bounds select the full minimal, full minimal and
second-neighbour models respectively. They select a smallest actual target
height only in the first row, with a tie there. In particular, the smallest
M29 target coordinate is on a model that the global bound does not select.
Three selected points on one curve do not estimate a general predictor's
accuracy; they do show that the uniform bound's ranking need not coincide
with the ranking for a particular new direction.

All twelve individual boxes at height 125,000 contain their target parameter,
but **none of the twelve certified sufficient elliptic-height ranges contains
the target**. The exact slack decomposition is retained. If `g` is actual
cancellation and `S=max(|N|,|D|)/H(parameter)^4`, then

\[
B_2-\bigl(h(t)-\tfrac14h_x(2P-Q)\bigr)
=\tfrac14\log(C/g)+\tfrac14\log(S/L).
\]

Here the finite contribution is approximately 9.85–13.99, versus 0.38–2.51
from the real contribution. The worst-case cancellation cap is the dominant
source of slack at these targets. This explains the limited pointwise
predictiveness of the guaranteed range in this panel; it does not change the
validity of the coverage theorem or determine CPU to an unknown direction.

The first diagnostic stopped when the two first neighbours were both at 7:
the original product-of-two-prime partition formula does not cover that case.
That joint result remains **UNKNOWN**. A diagnostic-only successor retains
four independently verified individual models instead of suppressing the
row. Its initial path-binding error and corrected replay are both retained
and charged; no experimental arm was altered. On the other two winning
charts the unchanged joint formula replays and allocates pair heights 58,751
and 46,367. Those pairs have **not** been point-searched in this experiment.

The [portable data](../../artifacts/generated-results/elliptic-curves/height_model_next_direction_v1/target-accessibility.json)
include the rational points, coordinates, model identities and local/real
slack factors. The portable checker replays all twelve bounds, exact
`2P-Q` transports, parameter heights and finite-rank certificates without
search, reduction or factorisation:

```sh
sage -python research/elliptic-curves/cas/verify_height_model_accessibility.py
```

## Arms and fixed allocation

| Arm | Preparation and search |
|---|---|
| Current cheaper V3 baseline | Factor-free chart reduction; height 125,000. |
| Best single bounded model | Construct the full minimal model and the first two neighbours from the original fixed rule; independently replay all bounds; select the smallest uniform `C/L`; height 125,000. |
| Certified complementary pair | Same model bank; select the pair with smallest joint `max_cell min_model D`; search each at the least integer height matching the selected single model's guaranteed elliptic-height coverage. |

The baseline is the cheaper representation from the
[completed September 12 comparison](NEXT_DIRECTION_RANK32_BENCHMARK_2026-09-12.md).
The older dual-map arm is not silently described as this baseline. Ties use
model indices. Within a pair, search the individually stronger bound first.
The pair height is the exact integer ceiling of
`125000 * (D_pair/D_single)^(1/4)`. The original first-chart regression gives
models 0 and 2 at height **84,441**. This matches a guaranteed height range;
it does not assert that the two boxes contain every point in the old box.

No new minimisation choices, larger prime range, additional centres or adaptive
height increases are admitted. The original joint formula requires neighbours
at distinct primes and its stated local partition hypotheses. An unsupported
chart is recorded as incomplete and the unchanged bank order continues.

## Complete cost at the supplied interface

Every arm starts a fresh process from its equation, subgroup and supplied
centre bank. Metering includes imports, starting-rank certification, new model
construction and reduction, bound production and independent replay, point
search, exact maps, duplicate rejection, finite admission and a separate
process replay. The final replay uses both the existing finite-quotient
implementation and independently enumerated Sage finite elliptic groups.

The bound prototype used an existing discriminant-prime certificate. The
bounded search arms must instead obtain their support from the equation by
fresh integer factorisation and unconditional prime proofs, with no supplied
factor hints. This gets at most **30 seconds** inside the arm allowance. A
miss is a preparation timeout, not free cached preparation or a mathematical
failure of the covering. The retained support is permitted only for the
separately metered post-seal diagnostics.

Historical selection of the input centre bank is outside the supplied
`(E,M,Q)` interface. Its full cold landscape cost remains UNKNOWN; this
experiment is not a comparison from an unprepared elliptic curve. Copying and
sealing inputs and sources is metered separately. No map or point cache is
shared across arms. Source snapshots preserve the actual implementation even
if other agents change the checkout.

One worker uses one fixed CPU with numerical thread counts set to one. Linux
child-subreaper accounting collects descendant CPU separately from elapsed
time. Arm order rotates over the three controls; each cell has one repetition.
Order effects and system contention therefore cannot be eliminated statistically.

## Stops and diagnostic boundary

Each of the nine arms has **300 seconds total**, including its independent
replay. The search worker is externally capped at 240 seconds and stops
starting new work after its 210-second internal threshold. Individual map,
model-bank and search limits are 5, 20 and 10 seconds. RSS is capped at 3 GiB
for an arm and 1 GiB for map and preparation workers. All limit hits remain
visible. Infrastructure or endpoint-replay failure stops the controller.

An arm stops at its first independently certified quotient gain. There is no
production follow-up, retry or automatic enlargement. The nine-arm allowance
is at most 45 minutes; final diagnostics get at most five more minutes. The
controller has an aggregate 51-minute ceiling including scheduling allowance.

Only after all nine endpoints are sealed does a separate evaluator receive
their winning points. For each literal point at its winning anchor it records:

- the exact rational parameter and height in all four model presentations;
- `H(parameter)^4 / (D * H_x(2P-Q))` and the logarithmic height distortion;
- membership in the requested search boxes and the certified sufficient range;
- the selected single and pair, their exact allocation and joint bound.

The evaluator independently checks `2P-Q` using Sage and exact rational group
arithmetic. Different arms can win on different points and anchors. Such rows
must not be presented as timing three representations of one identical target.
Diagnostic work is separately charged and cannot guide construction, allocation
or a second search. A failed search does not prove that every possible new
direction was inaccessible.

## Follow the run

The complete local packet is
`research/artifacts/local/elliptic-curves/height-model-next-direction-v1/`.
Its immutable `protocol.json` and `source-lock.json` were sealed before launch.
The controller updates `state.json`; every arm keeps its inputs, maps, bounds,
logs, supervisors, point witnesses, independent replay and `seal.json`.
`all-arms-sealed.json` gates the target evaluator. Process liveness is checked
using the PID and Linux start token; an old status file alone is not liveness.

From the repository root:

```sh
python3 research/elliptic-curves/cas/run_height_model_benchmark.py status \
  --folder research/artifacts/local/elliptic-curves/height-model-next-direction-v1

tail -f research/artifacts/local/elliptic-curves/height-model-next-direction-v1/controller.log

python3 research/elliptic-curves/cas/run_height_model_benchmark.py stop \
  --folder research/artifacts/local/elliptic-curves/height-model-next-direction-v1
```

The completed controller ran detached and checkpointed. These status and stop commands
do not restart it. The launch command refuses an existing launch receipt.

The allocation regression passes three tests, including exact integer rounding
and rejection of an expanding pair bound. An externally bounded preparation
replay also reproduced the original pair before launch, without point search.
The frozen experimental endpoint, rather than this launch description, decides
whether there was any measured benefit.
