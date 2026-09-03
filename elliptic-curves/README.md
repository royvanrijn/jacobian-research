# Elliptic curves over `Q` — ACTIVE

This programme is open for theorem-directed breakthroughs in exceptional
Mordell--Weil rank, low conductor, and the elliptic-K3 constructions behind
them. `../MATH_STATUS.json` is the sole status authority; this page is the
active navigation map.

## Current milestone

- ICARM curve 302: certified `rank E(Q) >= 31`, trivial torsion, global minimality, exact conductor/local data, and two independent point-independence implementations. No unconditional rank upper bound.
- ICARM curve 273: independently replayed `rank E(Q) >= 30`.
- ICARM curve 356: certified `rank E(Q) >= 29` with exact conductor/local data.
- ICARM curves 285/286 and curve 394: certified rank-at-least-21 results; curve 394 is the compact Elkies `t=3/8` specialization with exact conductor replay.
- The pinned K3 now has two explicit rootless arithmetic MW17 charts over `QQ`: published R17 and the direct degree-two alternate-Q80 chart from `norm12-orbit-11952`.
- The complete 43-chart norm-twelve atlas excludes curves 273 and 302, but proves that curves 351, 356, 376, 377, and 385 are untwisted fibres of one eight-chart published-R17 class. Their displayed-subgroup quotients by generic MW17 are `Z^8,Z^12,Z^5,Z^6,Z^12`.
- Rank `>=32`, unconditional exact rank for curve 302, and sharper conductor records remain open.

See [`../elkies-k3/README.md`](../elkies-k3/README.md) for the current K3 milestone.

## Compute gates

- Gate broad rank-32 Nagao, point, two-cover, or Selmer searches with a declared
  residual arithmetic target.
- Gate residual descent campaigns on the low-conductor near misses with exact
  local and quotient data.
- Give new family sweeps and expensive K3 specialization scans fixed limits and
  checkpointed outputs.
- Build native calibration fibres before alternate-Q80 specialization work.

Existing scripts, tests, local checkpoints, and generated certificates are retained for reproducibility.

## Canonical entry points

- [`notes/ICARM_CURVE302_RANK31.md`](notes/ICARM_CURVE302_RANK31.md) — rank-at-least-31 certificate.
- [`notes/ICARM_CURVE273_RANK30.md`](notes/ICARM_CURVE273_RANK30.md) — rank-at-least-30 certificate.
- [`notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md`](notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md) — rank-at-least-29 record/fingerprint.
- [`../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md`](../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md) — exact 43-chart record sweep and common five-fibre R17 construction.
- [`notes/ICARM_7FFF_ZIP_SEQUENCE.md`](notes/ICARM_7FFF_ZIP_SEQUENCE.md) — rank-at-least-21 curves 285/286.
- [`notes/ICARM_CURVE394_RANK21.md`](notes/ICARM_CURVE394_RANK21.md) — compact R17 rank-at-least-21 specialization.
- [`notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md`](notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md) — preserved low-conductor descent inputs.
- [`notes/ELKIES_RANK_JUMP_FINGERPRINTS.md`](notes/ELKIES_RANK_JUMP_FINGERPRINTS.md) — published-R17 specialization controls and quotient fingerprints.
- [`../archive/elliptic-curves/`](../archive/elliptic-curves/) — bounded-search history and superseded command surfaces.

## Active fronts

The useful gates are still:

1. use the exact pinned-ICARM fibre inventory—especially the `074d9` controls and native alternate-Q80 curve 12—to target a rank-32 neighbourhood, while separately pursuing an unconditional upper bound for curve 302;
2. complete genuinely residual Selmer/descent calculations before authorizing expensive point search;
3. use curve 12 and the other exact native alternate-Q80 fibres from the complete pinned-ICARM sweep to compare its specialization tail to published R17;
4. pursue low-conductor survivors only after exact quotient/descent gates justify them.

A heuristic score, point list without independence, incomplete Selmer calculation, or bounded miss is not a rank theorem.

## Reproduction

Use [`REPRODUCE.md`](REPRODUCE.md) and the exact checker paths recorded in `../MATH_STATUS.json`. The normal regression suite remains `make verify-elliptic-curves`; long CAS/search jobs remain separate targeted replays.

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 291a539d07b842b9 -->
<!-- status-consumer: EC-CF-NEARMISS-DESCENT-INPUTS 25c9f212e5162216 -->
<!-- status-consumer: OP-EC-NEXT e135b23ef9910845 -->
