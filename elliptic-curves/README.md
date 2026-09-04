# Elliptic curves over `Q` — ACTIVE

This programme is open for theorem-directed breakthroughs in exceptional
Mordell--Weil rank, low conductor, and the elliptic-K3 constructions behind
them. `../MATH_STATUS.json` is the sole status authority; this page is the
active navigation map.

## Current milestone

- ICARM curve 302: certified `rank E(Q) >= 31`, trivial torsion, global minimality, exact conductor/local data, and two independent point-independence implementations. No unconditional rank upper bound.
- Curve 302 point-cloud reconstruction: exact mod-2 and mod-3 finite Kummer
  codes have rank 31 with no visible first-17 boundary; elementary
  squareclass, degree-six held-out interpolation, and fixed-`X` deformation
  probes are negative in their declared bounded models.
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
- [`notes/ICARM_CURVE302_POINT_CLOUD_RECONSTRUCTION.md`](notes/ICARM_CURVE302_POINT_CLOUD_RECONSTRUCTION.md) — direct 31-point reconstruction probes and calibrated claim boundary.
- [`notes/ICARM_CURVE273_RANK30.md`](notes/ICARM_CURVE273_RANK30.md) — rank-at-least-30 certificate.
- [`notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md`](notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md) — rank-at-least-29 record/fingerprint.
- [`../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md`](../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md) — exact 43-chart record sweep and common five-fibre R17 construction.
- [`notes/ICARM_7FFF_ZIP_SEQUENCE.md`](notes/ICARM_7FFF_ZIP_SEQUENCE.md) — rank-at-least-21 curves 285/286.
- [`notes/ICARM_CURVE394_RANK21.md`](notes/ICARM_CURVE394_RANK21.md) — compact R17 rank-at-least-21 specialization.
- [`notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md`](notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md) — preserved low-conductor descent inputs.
- [`notes/ELKIES_RANK_JUMP_FINGERPRINTS.md`](notes/ELKIES_RANK_JUMP_FINGERPRINTS.md) — published-R17 specialization controls and quotient fingerprints.
- [`notes/HALF_LATTICE_FAKE_DESCENT_REPLAY_2026-09-04.md`](notes/HALF_LATTICE_FAKE_DESCENT_REPLAY_2026-09-04.md)
  — blind rank-28 half-lattice replay, exact productive-class ledger,
  equal-budget deep/random/shallow ablation with sealed +12 holdouts, and the
  failure of the pointed quartics to supply a prospective local-solubility
  predictor.
- [`notes/R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md`](notes/R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md) — corrected residual (2/4/8)-Selmer image filtration for curves 356 and 385.
- [`../archive/elliptic-curves/`](../archive/elliptic-curves/) — bounded-search history and superseded command surfaces.

## Active fronts

The useful gates are still:

1. reverse-engineer curve 302 from its complete 31-point configuration,
   calibrated on the actual transported generic subgroup of the known
   Fermigier--Mestre curve-245 control and without dimensioning the search at
   17;
2. use the exact pinned-ICARM fibre inventory—especially the `074d9` controls and native alternate-Q80 curve 12—to target a rank-32 neighbourhood, while separately pursuing an unconditional upper bound for curve 302;
3. accumulate proved residual-Selmer constraints monotonically, rejecting below
   the required residual dimension 15; permit only explicitly bounded point
   search while the full descent is open, and still require a complete
   unconditional descent for every Selmer or exact-rank claim;
4. pursue low-conductor survivors only after exact quotient/descent gates justify them.

A heuristic score, point list without independence, incomplete Selmer calculation, or bounded miss is not a rank theorem.

<!-- status-consumer: EC-ICARM-CURVE302-POINT-CLOUD 1e1eb37dd6d4350f -->

## Reproduction

Use [`REPRODUCE.md`](REPRODUCE.md) and the exact checker paths recorded in `../MATH_STATUS.json`. The normal regression suite remains `make verify-elliptic-curves`; long CAS/search jobs remain separate targeted replays.

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 8a4c932153e2bb2d -->
<!-- status-consumer: EC-CF-NEARMISS-DESCENT-INPUTS 25c9f212e5162216 -->
<!-- status-consumer: OP-EC-NEXT e135b23ef9910845 -->
