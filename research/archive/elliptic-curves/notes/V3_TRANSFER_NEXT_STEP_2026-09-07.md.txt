# After 302: one controlled V3 transfer, then three warm rank-27 jobs

**Execution handoff.** Do not retune V3 or continue optimizing against curve302.
First close its ongoing independent replay. Then run the bounded experiment
below. This handoff adds no mathematical-status promotion and no rank claim.

Input code snapshot: `305d4e788b370bd9fdce9321dde7441fc68ea7c3`.
The reported V3 17→31 search success is a recovery of a known curve, not a new
record. Its final packaging/replay must finish before this experiment opens.

## Why this is the next step

Test whether the calibrated algorithm transfers to another **parent surface**
and improves previously certified high-rank inputs. All four jobs below use
native/compact `11952` on the determinant948 surface, not302's determinant1092
surface. The complete generic norm8/10 bank must therefore be rebuilt in the
correct ordered generic basis. Do not reuse302's orbit numbers or its Gram.

| Fixed order | Input | Starting subgroup | Role |
|---|---|---:|---|
|1|Redacted native11952/curve12 fixture|17|Known-positive recovery control; exceptional points withheld|
|2|`11952`, `-2448/11`, inventory41|27|Warm inventory extension|
|3|`11952`, `2012/211`, inventory72|27|Warm inventory extension|
|4|`11952`, `4286/1881`, inventory186|27|Warm inventory extension|

The order is fixed, not selected from new outcomes. All seeds and the parent
bank are frozen before any of these point searches. No new parameter sweep,
MW16 adapter, other parent, or automatic website submission is authorized.
The public-data rank28 reproduction ID188/#619 is not a warm discovery seed.

## Commands for the executing agent

Work from `research/`. Finish the already running V3 replay rather than
starting a duplicate heavy checker. Its final JSON must be the output of
`check_visibility_cascade_v3.sage --start 17 --output ...`, **not** a progress
checkpoint. Finish the independent metric replay too:

```sh
sage -python elliptic-curves/cas/check_visibility_metric_v3.sage --start 17
```

Then prepare the new campaign, replacing only the replay filename:

```sh
sage -python elliptic-curves/cas/v3_transfer_campaign.sage prepare \
  --v3-replay artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3/FINAL_FULL_REPLAY.json
```

Preparation verifies the completed31 gate, reconstructs the compact11952
sections through their existing exact transport, enumerates its norm≤10 generic
lattice ellipsoid with **exact rational LDL**, validates every exported
parity/word/norm, and checks the rational specialization and independent seeds.
The existing Sage point exports are parsed as literal data, never executed.
Exact Q-isomorphisms transport them to the corresponding family models.
The first17 basis points must be the actual specialized generic sections.

Now execute one case and its full independent replay:

```sh
sage -python elliptic-curves/cas/v3_transfer_campaign.sage next
sage -python elliptic-curves/cas/v3_transfer_campaign.sage status
```

Repeat `next` for the remaining three cases **only after a completed,
independently replayed gain on the positive control**. The program enforces
that gate and stops after the four-case roster. Each invocation owns at most
one bounded worker at a time; it does not launch detached work itself.
An interrupted/censored or failed attempt is retained and requires review,
not an automatic retry with new limits. Do not run internal `*-worker` actions.

The orbit preparation has a20-million-node and one-hour wall cap. Each case
has a4096-chart cap, at most16 epochs, a rank≥32 stop, and a two-hour wall cap;
its replay has a separate two-hour cap. These are experimental limits, not
completion estimates. The shared supervisor enforces one worker and3GiB group
RSS, including descendant cleanup. All per-chart bounds and V3 shortlist
parameters remain those of the frozen calibration. No full minimization or
conductor computation is inserted into point search.

## What is reused and what changes

`adaptive_visibility_cascade_v3.sage::landscape`, `visibility_selection_v3.py`,
exact parity CVP, the factor-free mapper, and the existing finite-rank checkers
are imported unchanged. The new driver changes **inputs/parent bank** and the
rank32/budget stop only; it does not tune scores, quantiles, or preferred orbits.
After every certified gain it cancels the remaining epoch, preserves the
basis prefix, and rebuilds the complete extension landscape, just as V3 does.

The independent V3 checker is reused with an explicit seed/bank adapter.
Its historical `replay-M17/` subdirectory name is retained for compatibility;
**it does not imply that a warm case starts at17**. Every job/report stores the
actual starting rank. The additional replay verifies cumulative-cloud bytes
against precisely the seed and executed chart witnesses, recomputes height
metrics, and checks terminal no-gain clouds modulo3 and5 as well as2.

Generic shell enumeration is deliberately bounded and exact. A node/wall
failure publishes no complete bank and launches no point search. There is no
unproved completeness assertion from floating `qfminim`. For comparison, PARI's
own [qfminim documentation](https://pari.math.u-bordeaux.fr/dochtml/html/Vectors__matrices__linear_algebra_and_sets.html#qfminim)
explicitly qualifies its floating enumeration. LLL here is only a checked
unimodular change of basis before exact enumeration.

## Decision rules / report

- **Control has no certified gain:** stop. Report the finite exposure and
  determine whether parent-bank coverage, shortlisting or coordinate exposure
  failed. Do not interpret this as low true rank or retune on warm outcomes.
- **Control gains, warm cases do not:** retain the transfer calibration and
  bounded negatives. No larger scan is automatically justified.
- **Any warm gain:** replay, retain the new point/certificate and report
  `old lower bound → new lower bound`. Check the current catalogue afterward
  for an existing higher-rank submission or a size record; distinguish a rank
  improvement from a newly discovered curve.
- **Rank≥32:** preserve all sources and raw evidence; obtain a second
  standalone lower-bound proof before any record announcement. Do not stop
  at an incremental rank counter.

Report points/initial rank, final lower bound, chart count, exact-CVP/landscape
costs, resource failures, stale cancellations and source/checker hashes. Keep
all bounded misses and all preparation costs. The campaign's frozen roster
and `verified.json` files are the operational authority; existing
`MATH_STATUS.json`, V1/V2/V3 certificates and running jobs remain untouched.

A successful control is **not** population-level validation. A warm inventory
is selected experimental data, not an unbiased sample. Full coset enumeration
and semantic fingerprints are not a proof of full-policy basis invariance:
LLL/Babai, anchors and finite shortlisting still matter. A recovery cascade
explains point exposure, not why a fibre has its arithmetic rank.

## Cheap checks and remaining validation boundary

```sh
python3 -m unittest discover -s elliptic-curves/tests -p test_v3_transfer_contract.py
python3 -m py_compile elliptic-curves/cas/v3_transfer_contract.py \
  elliptic-curves/cas/v3_transfer_orbits.py elliptic-curves/cas/v3_transfer_campaign.sage
```

The added Python contracts and exact low-dimensional shell enumeration were
unit-tested. Sage/PARI is required for the actual17-dimensional parent adapter,
seed certification and solver/replay. That full transfer has **not** been run
by the author of this handoff; preparation and the positive control are
mandatory integration gates. Do not promote preparation itself as a result.
