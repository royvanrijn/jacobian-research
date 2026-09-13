# Extreme-anchored MW18 exact certification — historical continuation boundary (2026-09-04)

<!-- status-consumer: EC-K3-R17-EXTREME-ANCHORED-MW18-SPECIALIZATIONS 48b1e4d97ce68fc4 -->

The 178-fibre independence certificate below remains canonical. The
operational stop state and its follow-ups are historical: the
[MW18 deep-centre calibration](../elliptic-curves/notes/MW18_DEEP_CENTRE_CALIBRATION_2026-09-05.md)
was a completed bounded control, and the full cover-census replay did not
finish. The [elliptic-curve programme map](../elliptic-curves/README.md)
names current work. This source does not authorize a census, Selmer, or
point-search restart.

## Historical stop state (2026-09-04)

The campaign is intentionally stopped.  No census, Selmer, specialization,
or point-search process remains active.

The useful new result is complete: all 178 `H<=1000` Nagao finalists were
specialized exactly, including the five projective `r=infinity` rows.  On
every fibre the seventeen direct R17 sections and the cover section satisfy
the raw short Weierstrass equation exactly.  For every row their images have
full column rank 18 in a product of finite quotients
`E(F_p)/2E(F_p)`, and a separate good-reduction group order proves
`E(Q)[2]=0`.  Infinite descent therefore proves that the eighteen points are
integrally independent.

The exact artifact is
[`../artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json`](../artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json),
SHA-256
`912d44bcbf9570d22a77300602875aaf4d1791e9e782c6f78f084f4a8f2fe562`.
It records 178 successful specializations, 178 certified rank-at-least-18
subgroups, no structural failures, and no unresolved independence checks.

This is a lower-rank result only.  It proves neither exact specialized rank
18 nor exact generic rank 18.  It does not supply a 2-Selmer upper bound.

## Arithmetic-size warning

The raw model coefficient sizes across the 178 rows are 290 bits minimum,
8,479.5 bits median, and 8,913 bits maximum.  Point-coordinate sizes are 158
bits minimum, 4,273 bits median, and 4,493 bits maximum.  The small end is the
historical rank-28 anchor; the refreshed direct-chart finalists are generally
at the large end. Any separately scoped future use should minimize or otherwise
reduce a fibre before committing a long descent or point-search budget.

## Partial residual-Selmer pilot

[`../elliptic-curves/cas/run_r17_extreme_anchored_mw18_residual_selmer.py`](../elliptic-curves/cas/run_r17_extreme_anchored_mw18_residual_selmer.py)
defines a deterministic thirteen-model cohort: the top Nagao row on each of
the nine covers, plus the four distinct extreme-anchor raw fibres.  Exact
duplicate anchor presentations are collapsed.

Three genuine PARI `ellrank` workers were run with all eighteen points, a
60-second wall limit, 3 GB RSS limit, and 1 GB PARI stack:

| candidate | outcome | wall seconds |
|---|---|---:|
| `curve-536-08234-19188-p405d124` | strict timeout, no Selmer result | 60.120 |
| `curve-531-08234-12f61-m901d737` | strict timeout, no Selmer result | 60.119 |
| `curve-531-08234-0a9bf-p808d403` | strict timeout, no Selmer result | 60.120 |

Candidate four was interrupted before it produced a checkpoint. Earlier local
checkpoint rows with input key beginning `c1200e60` are wrapper-syntax failures
and have no arithmetic content. The three real timeouts have input key
beginning `6adc7ca0`; a later exact-artifact regeneration changed the current
input hash. The historical supervisor's proposed rerun is not a current queue.

No complete residual 2-Selmer dimension was obtained, hence there are no
authorized survivors. The serious point-search stage was not run. A timeout
is not a Selmer bound or evidence against high rank. The invalid
all-wrapper-failure aggregate was removed; no valid residual-Selmer aggregate
exists. Any replacement study needs its own objective, inputs, limits, and
checkpoint plan.

## Fresh no-cache replay

A fresh invocation of the cover compiler with `--check` was started without
reading the local cache.  It completely replayed the `07ca9` sieve, its two
exact survivors, and their curve-545 quotient span.  On `08234` it completed
the primes 1009 and 1013, leaving 51,267 unresolved bisection--fibre pairs and
one not-yet-validated bisection, then was stopped during the next prime.

This partial run did not reach the final byte comparison and is not an
independent full replay. The canonical exact cover artifact is unchanged.
The `08234` local sieve checkpoint records only the two completed fresh
primes. A fresh no-cache audit would begin from raw inputs and needs separate
scope; it is not a cleanup action.

## Historical replay interfaces

The retained exact specializer,
`elkies-k3/scripts/specialize_r17_extreme_anchored_mw18_finalists.sage`,
checks the completed 178-row certificate. The historical residual pilot is
`elliptic-curves/cas/run_r17_extreme_anchored_mw18_residual_selmer.py`; the
raw cover compiler is
`elkies-k3/scripts/certify_r17_extreme_anchored_mw18_covers.sage`.
Their pinned inputs, outputs, and original invocation metadata are preserved
with the generated artifacts. Running any of them requires a separately
scoped purpose; no result in this note is contingent on a new replay.

## Retained implementation files

The implementation files are:

- `elkies-k3/scripts/specialize_r17_extreme_anchored_mw18_finalists.sage`;
- `elliptic-curves/cas/run_r17_extreme_anchored_mw18_residual_selmer.py`.

The shared `mod_l_reduction_independence.py` remains at its repository version;
the specializer uses its existing exact finite-quotient implementation.

The worktree contains substantial unrelated concurrent changes.  Do not
normalize, revert, or fold them into this lane.
