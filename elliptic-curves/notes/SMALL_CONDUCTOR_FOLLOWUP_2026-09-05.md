# Bounded follow-up of the new small-conductor rank22 curve

The [curve and exact conductor proof](NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md)
remain at **rank at least 22**. The follow-up search retained 127 of its
301 planned charts before the declared 1,800-second / 1.5-GiB limit.
No extra independent direction was admitted. Its entire retained point cloud
contains 7,753 points up to sign, including the starting independent basis.
An independent finite-quotient audit through prime 997 also proves lower
bound 22, with 131 available quotient rows. This does not prove rank 22 is
an upper bound or that another rational point does not exist.

The [complete-cloud certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_followup_recorded_mod2_v1.json)
has a separate exact replay. The
[diagnostic bundle](../../artifacts/generated-results/elliptic-curves/small_conductor_followup_evidence_v1.json)
retains the chart prefix, numerical centre proposals, fixed protocol,
supervision outcome, point proof and replay logs. No chart exhausted its
whole height-100,000 box; the mean retained denominator fraction is about
0.7552. Neither the remaining 174 classes nor the unsearched parts of the
retained boxes have an absence certificate.

## Scheduling and the catalogue correction

The candidate was selected only after completing the twelve-address batch,
its independent rank proofs and the exact conductor calculation. Its
conductor lies below the smallest recorded rank-at-least-23 conductor.
The worker uses only the new curve's 22 certified points, crossing 301
generic parity classes with cyclic nonzero six-bit exceptional quotient
words and ordering their specialized representatives by numerical height.
The Gram computation explicitly requests 384-bit precision. An extra
independent point would stop the worker for verification.

A prepared predecessor used the older catalogue placement of second. The
refreshed 586-row gate found the newly recorded conductor of existing entry
575 and rejected that launch. The corrected worker's frozen gate records
the smaller-conductor IDs `[376,575]`. The predecessor source, protocol
and rejection remain historical evidence. No predecessor worker ran.

The [refreshed inventory replay](../../artifacts/generated-results/elliptic-curves/refreshed_new_curve_inventory_replay_v1.json)
independently compares all 36 new equations with all 586 pinned catalogue
equations. There are no Q-isomorphism matches. The earlier twelve-address
evidence ZIP separately passed extracted-copy replays of the conductor,
twelve lower bounds and 36-curve inventory. The first isolated inventory
attempt reached a 120-second cap; the retained 300-second retry passed.

## Exact-rank diagnostic remains incomplete

All ten bad-prime local root numbers were computed from the already proven
factorization; their product, including the infinite place, gives global
sign **+1**. This proves even analytic order of vanishing, not algebraic
rank parity or exact rank. See the
[PARI root-number documentation](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html#ellrootno).

A separately bounded PARI 2-descent probe used the integral equation and
all 22 known points, with a 180-second / 1.5-GiB limit and a 512-MB PARI
stack. It reached the time limit during `ellrankinit`, before a BNF
checkpoint or unconditional `bnfcertify` stage. It supplies **no upper
bound**, including no completed conditional bound. The source, protocol,
initial-stage checkpoint and timeout log are retained. No repeated descent
campaign or L-series summation was launched.

## The remaining admission cost was measured and reduced

The frozen `MWState` validates every old point observation after each new
admission, including already-known basis points. At a pinned 108-chart
prefix it had 11,018 observations, compared with 3,122 initially. Repeated
construction entails 55,828,668 old/new observation membership checks over
that prefix. Profiling one known-point admission attributes most time to
this repeated validation, followed by serialization for the state hash.

[`CachedObservationMWState`](../cas/research_runtime/cached_observation_state.py)
is an optional subclass whose validator is bound to the hash of the frozen
original source. It caches exact curve-membership booleans by immutable
curve and point values, with at most 65,536 entries. All other validation,
admission, observations, state keys and serialized records are preserved.
It does not cache rank conclusions or suppress observations. The running
point worker was not changed.

On three known-point admissions into the retained long history, the
original took 0.1707–0.1808 seconds per admission and the cached validator
0.0380–0.0464 seconds, with identical complete state records each time.
These small measurements are not a general speed guarantee. Separate tests
cover known points, ambiguity, a new independent direction, portable replay,
off-curve observations and changing the curve behind a repeated coordinate.
The [validation artifact](../../artifacts/generated-results/elliptic-curves/cached_observation_state_validation_v1.json)
binds the retained input and both implementations.

The new follow-up replayer uses this exact membership cache and the already
validated quotient-only reduction cache. The original runtime and replay
sources remain unchanged. State serialization still grows with the full
observation history; this cache does not solve that separate storage and
hashing cost.
