# Rank32 next-direction benchmark

The completed September 12 protocol, its inputs, worker history, commands, and
failure records are preserved byte-for-byte in
[the archive](../../archive/elliptic-curves/notes/NEXT_DIRECTION_RANK32_BENCHMARK_2026-09-12.md.txt)
(`sha256: d4632d18e48ac0d801e2b0a4d99ca29c693bd2cce5c3a2321d0dc156f782c76c`).
It does not authorize a new search, a larger height box, or a controller restart.

## Completed result

Within the retained centre banks, factor-free reduction at height 125,000 was
the cheaper known-control representation: four development recoveries used
151.019 CPU seconds versus 239.016 for dual-map V3, and the independent
Curve302 M30 validation used 83.570 versus 168.420. These are fixed-bank arm
costs, including preparation, map construction, search, and replay; they omit
historical cold landscape construction and do not prove a general speed law.

The frozen follow-up made 2,164 completed point-search calls on three rank-27
inputs (709, 729, and 726) with no new certified direction. Curve90 had all
376 candidate maps in the sealed prior-exposure exclusion list, so it produced
zero new coverage and is not a fourth failed exposure. Every completed call's
map, point, and rank replay passed. This bounded no-gain result neither raises
a rank lower bound nor supplies an upper bound or a dependence proof for every
returned point.

The [completion receipt](../../artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json),
[preflight](../../artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/preflight.json),
and [launch receipt](../../artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/launch.json)
bind the terminal state, input seals, replays, and box accounting. The local
packet named by the archive remains the full reproducibility record.

## Reuse boundary

Use factor-free reduction only as a cheaper representation in a separately
authorized, compatible experiment. Do not re-search these completed boxes or
interpret the no-gain endpoint as an exclusion. A rank-32 attempt needs a
declared coverage or admission gap beyond the completed banks; the current
[programme map](../README.md) chooses any future work. The
[centre-policy lesson](../../knowledge/ALGORITHMS.md#method-ec-centre-policy-keep-shallow-enumeration-deep-centres-and-adaptive-coverage-distinct)
contains the reusable rule.
