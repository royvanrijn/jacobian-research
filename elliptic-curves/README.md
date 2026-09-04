# Elliptic curves over `Q` — ACTIVE

This programme is open for theorem-directed breakthroughs in exceptional
Mordell--Weil rank, low conductor, and the elliptic-K3 constructions behind
them. `../MATH_STATUS.json` is the sole status authority; this page is the
active navigation map.

## Current milestone

The target-free A1/MW16 parameter experiment now runs alongside the other
active fronts.  It samples new rational parameters on one deduplicated exact
pencil and measures exact quotient rank beyond specialized MW16 without
loading known-record targets, parameters, public points, ranks, or target
`j`-invariants.

- ICARM curve 302: certified `rank E(Q) >= 31`, trivial torsion, global minimality, exact conductor/local data, and two independent point-independence implementations. No unconditional rank upper bound.
- Curve 302 point-cloud reconstruction: exact mod-2 and mod-3 finite Kummer
  codes have rank 31 with no visible first-17 boundary; elementary
  squareclass, degree-six held-out interpolation, and fixed-`X` deformation
  probes are negative in their declared bounded models.
- ICARM curve 273: independently replayed `rank E(Q) >= 30`.
- ICARM curve 398: independently replayed `rank E(Q) >= 30`, trivial torsion,
  a singleton rational isogeny class, and exact semistable conductor/local
  data.  Its hidden `A1`/MW16 fibration has now been recovered on the
  `norm12-orbit-11952` chart, including the exact parameter and a saturated
  sixteen-section specialization.  The apparent second norm-eight survivor is
  exactly the same fibration after an affine base change and Weierstrass
  scaling.  From a redacted MW16 input, the generic half-lattice plus adaptive
  quotient search blindly recovers all fourteen held-out directions and the
  full displayed rank-30 subgroup.
- ICARM curve 356: certified `rank E(Q) >= 29` with exact conductor/local data.
- ICARM curves 285/286 and curve 394: certified rank-at-least-21 results; curve 394 is the compact Elkies `t=3/8` specialization with exact conductor replay.
- The pinned K3 now has two explicit rootless arithmetic MW17 charts over `QQ`: published R17 and the direct degree-two alternate-Q80 chart from `norm12-orbit-11952`.
- The refreshed complete 43-chart norm-twelve atlas decides every equation in
  the hash-pinned 573-curve ICARM response: 86 hits and 3,352 class misses,
  with all 479 native comparisons untwisted.  The new priority cohort includes
  curve 543 with displayed quotient `Z^12`, six rank-at-least-28 fibres with
  quotient `Z^11`, and four further exact priority quotients.  Five lower-rank
  hits also have exact quotients; on curve 499, adjoining the non-contained
  generic MW17 subgroup enlarges the displayed subgroup by `Z/3Z`.  Curves 542
  and 548 miss all six classes despite independently replayed rank lower bounds.
  The preserved high-rank misses include curves 273, 302, and 398; these are
  six-class atlas exclusions, not exclusions from other K3 fibrations.
- A redacted, fixed-policy blind ladder on the sixteen quotient-eligible new
  hits recovers exact quotient ranks
  `6,6,11,3,11,10,11,10,5,6,8,8,12,0,11,8` in curve-id order.  Exact
  tied-margin Kendall association with the displayed jumps is `tau_b=.7503`
  (`p=60852/2421619200`), while score at least ten selects 7/8 of the
  `+10/+11/+12` tail and 0/8 below it (`p=1/1430`).  This passes the frozen
  extreme-jump detector gate but does not replace residual Selmer.  The blind
  curve-478 basis also improves its unconditional lower bound to rank at
  least 23.
- Rank `>=32`, unconditional exact rank for curve 302, and sharper conductor records remain open.

The rank-32 roadmap is parallel.  The target-free A1/MW16 experiment searches
directly for rank jumps while the R17/MW17 and other experiments continue.
Curve 398 remains an exact historical calibration, not an input target: its
generic MW16 and fourteen displayed quotient directions are already recovered
blindly.  A rank-32 fibre in the same family needs sixteen quotient directions.

See [`../elkies-k3/README.md`](../elkies-k3/README.md) for the current K3 milestone.

## Compute gates

- Gate broad rank-32 Nagao, point, two-cover, or Selmer searches with a declared
  residual arithmetic target.
- Make the residual target family-relative: a certified generic rank `r`
  requires at least `32-r` independent quotient directions, hence 15 for MW17
  and 16 for A1/MW16.
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
- [`notes/ICARM_CURVE398_RANK30_AND_CONSTRUCTION.md`](notes/ICARM_CURVE398_RANK30_AND_CONSTRUCTION.md)
  — rank-at-least-30 certificate, exact base equivalence of the two recovered
  A1/MW16 survivor presentations, equality of their specialized integral MW16
  groups, and blind rank-14 quotient rediscovery from the first presentation.
<!-- status-consumer: EC-K3-CURVE398-A1-MW16-RECOVERY 75978a18cc26690f -->
<!-- status-consumer: EC-K3-CURVE398-TWO-PARENT-COLLISION 626a440519ff77f3 -->
- [`notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md`](notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md) — rank-at-least-29 record/fingerprint.
- [`notes/ICARM_573_CURVE_REFRESH_OVERVIEW_2026-09-04.md`](notes/ICARM_573_CURVE_REFRESH_OVERVIEW_2026-09-04.md)
  — exact 573-curve atlas refresh, complete appended-row intake, sixteen new
  specialization quotients, and the curve-499 commensurability obstruction.
- [`notes/R17_REFRESH_BLIND_JUMP_LADDER_2026-09-04.md`](notes/R17_REFRESH_BLIND_JUMP_LADDER_2026-09-04.md)
  — redacted 16-fibre adaptive half-lattice ladder over observed jumps `+3`
  to `+12`, exact pre-complement rank responses, passing frozen ordinal and
  `q>=10` upper-tail tests, post-freeze fibration/`j`-class sensitivity, and
  the new curve-478 rank-at-least-23 lower bound.  High score may schedule
  within the calibrated norm-twelve R17 setting; low score may not veto, and
  alternate-Q80 extreme-tail transfer remains unvalidated.
<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER b7518cc41268489a -->
- [`notes/MW17_JUMP_V2_2026-09-04.md`](notes/MW17_JUMP_V2_2026-09-04.md)
  — the 2,239-fibre immutable-population evaluation now ranks the old
  bounded-box, balanced CRT, six atlas-family, and alternate-Q80 candidates
  only by exactly certified quotient rank recovered beyond MW17.  It runs
  `07ca9` and `08234` first, never uses initial 43-chart gain as a hard filter,
  checkpoints every fibre, and stops globally on a certified `+15`.
- [`../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md`](../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md) — exact 43-chart record sweep and common five-fibre R17 construction.
- [`notes/ICARM_7FFF_ZIP_SEQUENCE.md`](notes/ICARM_7FFF_ZIP_SEQUENCE.md) — rank-at-least-21 curves 285/286.
- [`notes/ICARM_CURVE394_RANK21.md`](notes/ICARM_CURVE394_RANK21.md) — compact R17 rank-at-least-21 specialization.
- [`notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md`](notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md) — preserved low-conductor descent inputs.
- [`notes/ELKIES_RANK_JUMP_FINGERPRINTS.md`](notes/ELKIES_RANK_JUMP_FINGERPRINTS.md) — published-R17 specialization controls and quotient fingerprints.
- [`notes/R17_SMALL_FIELD_CLASS_QUOTIENT_LAB.md`](notes/R17_SMALL_FIELD_CLASS_QUOTIENT_LAB.md)
  — the frozen 100-fibre rank-blind R17 cohort and fail-closed two-phase test
  of whether the unconditional class quotient predicts later certified
  half-lattice escape; Phase 0 is frozen and the BNF feature campaign is open.
- [`notes/HALF_LATTICE_FAKE_DESCENT_REPLAY_2026-09-04.md`](notes/HALF_LATTICE_FAKE_DESCENT_REPLAY_2026-09-04.md)
  — blind rank-28 half-lattice replay, exact productive-class ledger,
  equal-budget deep/random/shallow ablation with sealed +12 holdouts, and the
  failure of the pointed quartics to supply a prospective local-solubility
  predictor; it also records the frozen two-stage replacement detector now
  running on the pre-existing 2,560-fibre CRT cohort and the fail-closed rule
  that its binary escape endpoint cannot promote a rank-32 candidate.  Its
  chart-order policy gives legacy depth/old-deep/Hamming fields search-order
  meaning only, invalidates them on every lattice or basis change, and forbids
  absence or Selmer inference from a miss.
<!-- status-consumer: EC-K3-R17-074D9-HALF-LATTICE-PROMOTION-GATE 9a1f080523c9ecae -->
- [`notes/HALF_LATTICE_HEIGHT_COMPRESSION_MECHANISM_2026-09-04.md`](notes/HALF_LATTICE_HEIGHT_COMPRESSION_MECHANISM_2026-09-04.md)
  — exact midpoint and height identities for the pointed quartic charts, a
  3,865-chart coordinate-map audit plus 394 compact earlier-control records,
  and the separation between deep-hole old-point exclusion and target-relative
  height compression.  It replaces mask count by current-lattice covering
  radius plus reduced-coordinate distortion as the proposed rank-32 scheduling
  geometry; the existence of a new point inside an empty ball remains a
  separate arithmetic question.  Its first curve-385 builder pilot samples the
  full current `M29/2M29` parity space rather than another quotient-weight
  shell, commits 16 fresh deep charts, and obtains 0 finite points and no group
  growth at the historical bound; this is a clean old-point-exclusion result,
  not a saturation or rank conclusion.
<!-- status-consumer: EC-HALF-LATTICE-HEIGHT-COMPRESSION 3baeaf370aec751c -->
<!-- status-consumer: EC-CURVE385-HEIGHT-COMPRESSION-BUILDER-PILOT c0d6f2d67018def4 -->
- [`notes/CURVE385_ITERATED_HALF_LATTICE_RECOVERY_2026-09-04.md`](notes/CURVE385_ITERATED_HALF_LATTICE_RECOVERY_2026-09-04.md)
  — the quotient-bit iteration from blind `M20` to blind `M29`; post-freeze
  mutual integral coordinates prove equality with the displayed public
  rank-29 subgroup. Its exact weight profile shows that weight one spans seven
  of the nine new directions and weight at most two spans all nine; a frozen,
  checkpointed sparse-mask rank-32 protocol replaces monolithic enumeration.
  Its v2 operational amendment replaces the insufficient combined four-state
  allowance with independent limits for three rank-changing and four
  saturation-only group changes, while preserving the completed v1 campaign.
  Stability and exact rank remain open.
<!-- status-consumer: EC-K3-R17-CURVE385-INDEPENDENT-RESTART-BUDGETS 39cfce110e3e494f -->
- [`notes/R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md`](notes/R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md) — corrected residual (2/4/8)-Selmer image filtration for curves 356 and 385.
- [`../archive/elliptic-curves/`](../archive/elliptic-curves/) — bounded-search history and superseded command surfaces.

## Active fronts

The useful gates are now:

1. execute the A1/MW16 parameter experiment without any known-record target in
   candidate selection or execution, and certify every quotient gain exactly;
2. continue the peer R17/MW17 path using the exact refreshed ICARM inventory—
   especially curve 543, the six new rank-at-least-28 fibres, the `074d9`
   controls, and native alternate-Q80 curve 12—with the passing blind
   extreme-jump detector as a scheduling signal, never as a substitute for
   the residual-Selmer promotion gate;
3. accumulate proved residual-Selmer constraints monotonically, rejecting
   below the family-relative requirement `32-r`; permit only explicitly
   bounded point search while the full descent is open, and still require a
   complete unconditional descent for every Selmer or exact-rank claim;
4. continue curve-302 parent reconstruction independently of the parameter
   experiments;
5. pursue an unconditional upper bound for curve 302 and low-conductor
   survivors only after exact quotient/descent gates justify them.

A heuristic score, point list without independence, incomplete Selmer calculation, or bounded miss is not a rank theorem.

<!-- status-consumer: EC-ICARM-CURVE302-POINT-CLOUD 1e1eb37dd6d4350f -->

## Reproduction

Use [`REPRODUCE.md`](REPRODUCE.md) and the exact checker paths recorded in `../MATH_STATUS.json`. The normal regression suite remains `make verify-elliptic-curves`; long CAS/search jobs remain separate targeted replays.

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 8a4c932153e2bb2d -->
<!-- status-consumer: EC-K3-R17-NORM12-ICARM-573-REFRESH a93ce35de34fde21 -->
<!-- status-consumer: EC-CF-NEARMISS-DESCENT-INPUTS 25c9f212e5162216 -->
<!-- status-consumer: OP-EC-NEXT 086b745cb3ca078e -->
