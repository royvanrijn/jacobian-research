# Full11952 retained-score trial

The completed short-score scan covers **20,888,422,894 primitive nonzero
parameters** n/d with |n|,d at most131072 in the compact11952 family. Both signs
and all1024 denominator residue classes are included. Each of2048 slices keeps
its512 strongest short scores, giving1,048,576 retained candidates. Every
retained short sum and good-prime count replays exactly against the projective
cache. The full scan took838.583 seconds; its eight-slice timing gate passed
before the remaining blocks ran.

All1,048,576 candidates then receive the unchanged extended S1 score through
32749. Cached integer scoring takes11.437 seconds; the separate short-score
replay takes2.278 seconds. The completed eleven-stage controller checks the
bank, both scores, fixed selection and fresh scalar traces. These are recorded
elapsed times, not guarantees for other workloads.

The [frozen64 selection](../../artifacts/generated-results/elliptic-curves/full11952_h131072_selection_v1.json)
orders by combined S1, good-prime count, denominator and signed numerator. It
excludes789 previously measured equations and earlier selected rational
isomorphism classes. All64 pass fresh scalar GP agreement for the selection
band and a disjoint validation band;256 direct character sums also pass.
Validation does not influence ordering. Extended scoring covers the retained
million rows, **not all20.89 billion parameters**; no global extended-score
optimality or completeness of high-rank discovery is asserted.

The mathematical motivation is the separately replayed
[known29 retention diagnostic](EXTENDED_CACHE_RETENTION_GATE_2026-09-06.md):
the same control moves from position435 to2 against the same967 saved outer
candidates when the extended score becomes available. This motivated retaining
every slice through the stronger score. The control is not an input to the
prospective ordering or point maps, and its conditional position is not a
population quantile or a measured discovery probability.

After the64 selections and all point maps were frozen, the
[full-bank control audit](../../artifacts/generated-results/elliptic-curves/full11952_known29_retention_v1.json)
confirmed that89074/31895 survived: position291 in its512-row slice,
478285 under short S1 among all1048576 survivors,30 under extended S1,
and29 in the final64 after equation exclusions. Its cached values agree with
the earlier independent scalar control call. Thus a global short prefix below
478285 within this retained population would discard this known strong curve.
The actual policy preserves it. This is one retrospective positive selection
control, not a measured success rate. The789 exclusions are previous cohort
equations; they are not the full public catalogue or every separate control.
The selected known curve cannot count as a novel discovery.

## Fixed point exposure

The point protocol is frozen after all score checks pass. It gives each of the
64 equations all49 generic maximum-parity charts, at height125000 and ten
seconds per chart: at most3136 initial boxes. All64 map files must precede any
point search. There are two workers,600 seconds per curve, and a20000-second
batch cap. A provisional lower bound28 stops that curve pending exact proof.
There is no adaptive wave, refill or automatic budget extension.

Streaming verification replays every search history and exact point-cloud
independence. Any provisional bound at least28 also receives immediate rational
geometry and modulo3/5 checks. After the fixed cohort finishes, the finalizer
checks all geometry and clouds modulo2,3,5 and compares the equations against
the pinned593-curve catalogue. All3136 point boxes,64 exact histories and
full-cloud checks pass. The cohort adds a
[new27-point curve and a new23-point curve](FULL11952_NEW_RANK27_2026-09-06.md),
raising the inventory to187. The known control reaches the28-point stop and
matches ICARM12; it is excluded from the additions.

Sources are `scan_full11952_h131072.py`, `score_full11952_retained.py`,
`finish_full11952_retained.py`, `full11952_64_r17_pari_batch.py` and its named
preparer/verifier/finalizer under `../cas/`. Frozen protocols, raw outputs and
supervisor records are retained under `artifacts/local/elliptic-curves/` in
`full11952-h131072-short-v1`, `full11952-h131072-retained-v1`,
`full11952-retained-controller-v1` and `full11952-64-r17-pari-v1`.
All194 isolated point-history, geometry, cloud and equation-comparison stages
pass. The point-only bundle binds the selection as context; it does not rerun
the full scan, million-row scores or trace-cache construction. Scores, bounded
misses and finite point clouds prove no exact rank, upper bound, saturation,
point absence or universal novelty.
