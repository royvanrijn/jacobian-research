# Elliptic curves over `Q` — ACTIVE

The [new small-conductor curve](notes/NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md)
has 22 independently certified points and an exact 76-digit conductor. It
would be third among rank-at-least-22 entries with recorded conductors in
the pinned 586-curve ICARM catalogue; four such entries lack conductors.
The completed searches now give thirty-six distinct curves with certified
lower bounds 22–25, including three at least 25, and no match in that snapshot.
The [consolidated JSON](../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v2.json)
contains all equations, points and rank certificates; an
[equation CSV](../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v2.csv) is also available.
The [machinery audit](notes/ELLIPTIC_BREAKTHROUGH_AUDIT_2026-09-05.md)
records the selection and checkpoint fixes; new rank-at-least-28/32 targets
and exact ranks remain open. Catalogue absence does not prove universal novelty.
The [retained-point audit](notes/RECORDED_POINT_ADMISSION_AUDIT_2026-09-05.md)
checks 202 transcripts and recovers a missed rank-26 certificate on a known
control at prime 257; the new curves’ lower bounds remain unchanged.
The [cross-family incidence proof](notes/COMPACT_CROSS_FAMILY_INCIDENCE_2026-09-05.md)
checks all 384 pairs of the first 32 curves and twelve family presentations; it
finds no additional generic directions beyond the original families.

The [compact six-family atlas](notes/COMPACT_SIX_R17_ATLAS_2026-09-05.md)
reduces the compiled R17 equations to 141–169 coefficient bits and exactly
transports all 102 generic sections. Its balanced pilot and all retained
adaptive follow-up transcripts have been replayed.

The [five compact MW16 families](notes/COMPACT_FIVE_MW16_ATLAS_2026-09-05.md)
add 80 exactly transported generic sections and 141–181-bit equations,
broadening the compact input base to eleven family models.

The [completed rank-jump diagnostics](notes/RANK_JUMP_DIAGNOSTICS_2026-09-05.md)
explain retained MW18 visibility, recover 31/31 masked ordinary-fibre
directions, compare disjoint-prime selectors, and measure the height cost of
one explicit MW18-to-MW19 construction. No selector is promoted.

The [fibre-height population experiment](notes/FIBRE_HEIGHT_POPULATION_2026-09-05.md)
compares Nagao selection with actual arithmetic height and measured chart cost
on fresh bounded MW16 and MW18 populations.

New arithmetic and search work uses the [shared runtime](notes/SHARED_RESEARCH_RUNTIME.md):
cached labelled fields, subspace descent, lazy MWState search and retained-witness replay.

The [MW18 deep-centre calibration](notes/MW18_DEEP_CENTRE_CALIBRATION_2026-09-05.md)
uses the exact generic height geometry: deepest, diverse deep, and nearest-first
recover 22, 21, and 0 directions across five anchor presentations. The frozen
prospective gate fails.

This programme is open for theorem-directed breakthroughs in exceptional
Mordell--Weil rank, low conductor, and the elliptic-K3 constructions behind
them. `../MATH_STATUS.json` is the sole status authority; this page is the
active navigation map.

All active half-lattice searches now use the rank-agnostic
[`PointedQuarticSearch`](notes/POINTED_QUARTIC_SEARCH.md): MW16, MW17, MW18,
curve-specific and zero-gain routes share one GMP sieve. The migration
replays all 1,034 calibrated control boxes and all 55 quotient directions.

The [2026-09-04 external audit](notes/EXTERNAL_AUDIT_2026-09-04.md) records
mathematical corrections, bounded certificate replays, and remaining proof
and reproduction gaps.

The [rational-solubility theorem package](../elkies-k3/RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md)
connects residual Selmer classes to the point-search charts. An exact
eleven-fibre replay supplies 110 independent soluble 2-cover classes and
their forced-zero Cassels--Tate rows; unknown complements remain unclassified.
The [soluble-versus-Sha comparison](notes/EXCEPTIONAL_SOLUBLE_VS_SHA_PANEL_2026-09-05.md)
adds 62 marked quartic presentations on 356/385/398/400/543 and a certified
rank-16 Sha obstruction control. Pairing rank detects Sha; the surviving
radical and coefficient size do not certify solubility.

## Current milestone

The target-free A1/MW16 parameter experiment has completed alongside the other
active fronts.  It sampled 104 new rational-parameter fibres across five exact
fibrations in nine anonymous coordinate charts, without loading known-record
targets, parameters, public points, ranks, or target `j`-invariants.

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
- The complete A1/MW16 atlas now supplies a five-curve blind ladder rather
  than nine independent hits: exact base-change audits collapse the repeated
  labels on curves 398, 400, and 548 to one fibration each.  The initial
  maximum-depth wave recovers 38 of 55 demonstrated quotient directions; the
  historical curve-398 and curve-400 adaptive waves raise that blind
  total to 54 of 55.  Curve 400 is an exact `M16 -> M21 -> M28` recovery.
  The completed target-free height-300 parameter experiment retains 104 exact,
  pairwise nonisomorphic fibres across all nine anonymous coordinate charts.
  Direct exact searches exhaust all 856 maximum-depth MW16 quartic charts
  through height 100,000, with no timeout or structural failure, but return no
  affine point and no quotient direction beyond MW16.  This is a bounded null
  result, not a rank upper bound; zero candidates advance to Selmer or
  unrestricted point search.
- The [specialized pointed-quartic sieve](notes/POINTED_QUARTIC_SIEVE.md)
  completes all 856 frozen prospective charts at height 10,000 without generic
  minimization or reduction. Its exact denominator/lattice transforms yield
  1,537--1,789-bit quartics; all modular searches together take 53.1 seconds.
  That first run finds no prospective point and recovers twenty initial
  control directions. The new [sensitivity calibration](notes/MW16_SENSITIVITY_RECOVERY_2026-09-05.md)
  varies height, horizontal coordinates, rational slopes, and centres. Its
  frozen metric policy recovers **55/55** control directions, including the
  missing curve-401 direction. The subsequent 104-fibre rerun completes all
  856 height-100,000 boxes with no finite point, quotient gain, or timeout;
  this is a bounded null result, not a rank upper bound.
  The [model-size audit](notes/ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md#exact-arithmetic-model-audit)
  now supplies global minimal models, all sixteen transported sections and
  arithmetic-selected quartic coordinates for every finalist.
- ICARM curve 356: certified `rank E(Q) >= 29` with exact conductor/local data.
- ICARM curves 285/286 and curve 394: certified rank-at-least-21 results; curve 394 is the compact Elkies `t=3/8` specialization with exact conductor replay.
- The [fixed-field comparison](notes/FIXED_FIELD_COMPARISON_2026-09-05.md)
  is complete. At `u=-2,1,2`, the profiles
  `(local dimension, CT rank, radical dimension, certified realized dimension)`
  are all `(13,12,1,0)`. The frozen extension gives `(17,16,1,0)` at `u=-3`
  and `(15,14,1,0)` at `u=3`; the `u=0` control is `(20,0,20,20)`.
  All five remaining classes received bounded point searches, with no hit.
  Their point-or-Sha status stays unknown; the success criterion was not met.
  The [u=-1 baseline](notes/FIXED_CUBIC_U_MINUS1_CASSELS_TATE_2026-09-05.md)
  retains its rank-16 pairing, two-dimensional radical, and three unresolved
  inherited classes. Its known rank-one point lies outside the inherited span.
  Earlier minimal-model, empty-box and tangent-conic evidence remains indexed
  in that proof note.
  The bounded relative-norm and local-reconstruction conic-solver comparison
  also leaves the first auxiliary point unconstructed.
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-SOLVER-COMPARISON 6a178bc3a4ada43b -->
  The longer norm/reconstruction follow-up checked 12.48 million further
  candidate vectors without a point; no higher cover is available.
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-LONG-SEARCH 825fb4cd6ed84cb1 -->
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-MINIMAL-MODELS 90216b8c456edd20 -->
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-SEARCH-GEOMETRY 678f7beb805a4530 -->
<!-- status-consumer: EC-FIXED-CUBIC-TANGENT-CONIC-GATE 26a49e30ff3128d3 -->
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-CASSELS-TATE df45391a84f0e3c9 -->
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-RANK1 7e488a894d136732 -->
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

The rank-32 roadmap is parallel.  The completed target-free A1/MW16 wave
searched directly for rank jumps while the R17/MW17 and other experiments
continued.
Curve 398 remains an exact historical calibration, not an input target: its
generic MW16 and fourteen displayed quotient directions are already recovered
blindly.  A rank-32 fibre in the same family needs sixteen quotient directions.
The target-informed five-parent ladder is a separate prospective lane: retain
all nine bounded-height coordinate charts, count five fibration-level
observations, and recover every quotient jump exactly.  Residual Selmer may
exclude or prioritize a fibre, but incomplete descent does not veto a finite,
checkpointed point search.

See [`../elkies-k3/README.md`](../elkies-k3/README.md) for the current K3 milestone.

## Proof and compute gates

- Separate proof exclusions from resource authorization.  Only an
  unconditional certified rank/Selmer upper bound below 32 excludes a fibre;
  incomplete or conditional arithmetic changes scheduling only.
- Give every rank-32 Nagao, point, two-cover, or Selmer search declared finite
  limits, checkpoints, and an exact candidate/model binding.
- Make the residual target family-relative: a certified generic rank `r`
  requires at least `32-r` independent quotient directions, hence 15 for MW17
  and 16 for A1/MW16.
- Gate residual descent campaigns on the low-conductor near misses with exact
  local and quotient data.
- Accept 32 exactly verified independent rational points as an unconditional
  rank-at-least-32 success without requiring descent to finish; exact rank
  still needs a matching upper bound.
- Treat a large residual 2-Selmer quotient only as the combined envelope of
  new Mordell--Weil directions and `Sha[2]`.  The constructive gate is explicit
  locally soluble 2-coverings with rational points, followed by exact
  independence modulo the growing known subgroup; no scalar class-group
  statistic substitutes for that step.
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
- [`notes/ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md`](notes/ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md)
  — exact reduction of nine atlas labels to five fibrations, complement-blind
  `38/55` initial and `54/55` best ladder recovery, the complete blind
  curve-400 `+12` recovery, and the 104-fibre prospective
  Nagao-to-half-lattice attempt.  All 856 prospective quartic charts time out,
  so exact model/section size reduction is the next engineering gate and no
  Selmer or expensive-search promotion occurs.
<!-- status-consumer: EC-K3-ICARM-MW16-BLIND-LADDER acfa3bdcebb18137 -->
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
<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER a2d7034fb8977c18 -->
- [`notes/MW17_JUMP_V2_2026-09-04.md`](notes/MW17_JUMP_V2_2026-09-04.md)
  — the 2,239-fibre immutable-population evaluation ranks the old
  bounded-box, balanced CRT, six atlas-family, and alternate-Q80 candidates
  only by exactly certified quotient rank recovered beyond MW17.  It runs
  `07ca9` and `08234` first, never uses initial 43-chart gain as a hard filter,
  checkpoints every fibre, and stops globally on a certified `+15`.  The run
  was stopped on a negative signal after all 66 completed leading-`07ca9`
  fibres returned measured `+0`: all 2,838 charts completed with no timeout or
  backend failure.  This is detector-budget evidence, not a rank upper bound.
  A separately frozen one-in-eight zero-gain rescue arm retains the same 2,239
  addresses and exact-gain ranking.  Assigned clean zeros receive the next 301
  generic classes in seven batches, switching unused slots to the established
  adaptive policy after the first escape; the total remains 344 charts.  The
  arm is frozen but unrun.
<!-- status-consumer: EC-K3-MW17-JUMP-V2-ZERO-GAIN-RESCUE 39ac93b60152bf88 -->
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
- [`notes/QUOTIENT_GEOMETRY_TABLE_2026-09-04.md`](notes/QUOTIENT_GEOMETRY_TABLE_2026-09-04.md)
  — the complete 30-presentation quotient-Gram, regulator, successive-minimum,
  and 230-direction chart-geometry table.  All thirteen nonempty strictly
  partial initial recoveries reject a scalar intrinsic-quotient-height cutoff;
  the pointwise projection CVPs separate optimal half-lattice position from
  actual chart phase and coordinate distortion.
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
- [`notes/FIXED_CUBIC_FIELD_VARYING_CURVE_EXPERIMENT_2026-09-04.md`](notes/FIXED_CUBIC_FIELD_VARYING_CURVE_EXPERIMENT_2026-09-04.md)
  — explicit fixed-2-division-field family, certified rank-20 Kummer anchor,
  exact whole-span local kernels for five curves, and the covering equations
  for the still-open point-realization phase.
<!-- status-consumer: EC-FIXED-CUBIC-VARYING-CURVE-LOCAL-KUMMER 46ca45db3e702eb6 -->
- [`../archive/elliptic-curves/`](../archive/elliptic-curves/) — bounded-search history and superseded command surfaces.

## Active fronts

The useful gates are now:

1. use the completed target-free A1/MW16 parameter experiment as the baseline:
   its 104-fibre, 856-chart wave found no bounded quotient gain, so any later
   wave starts from this null result rather than from historical target fibres;
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
5. use the completed [fixed-field comparison](notes/FIXED_FIELD_COMPARISON_2026-09-05.md)
   to select any next bounded experiment before enlarging a point budget;
   all five tested deformations leave one unresolved inherited class, versus
   three at `u=-1`; restrict any further point solving to certified radicals;
6. pursue an unconditional upper bound for curve 302 and low-conductor
   survivors only after exact quotient/descent gates justify them.

A heuristic score, point list without independence, incomplete Selmer calculation, or bounded miss is not a rank theorem.

<!-- status-consumer: EC-ICARM-CURVE302-POINT-CLOUD 1e1eb37dd6d4350f -->

## Reproduction

Use [`REPRODUCE.md`](REPRODUCE.md) and the exact checker paths recorded in `../MATH_STATUS.json`. The normal regression suite remains `make verify-elliptic-curves`; long CAS/search jobs remain separate targeted replays.

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 8a4c932153e2bb2d -->
<!-- status-consumer: EC-K3-R17-NORM12-ICARM-573-REFRESH a93ce35de34fde21 -->
<!-- status-consumer: EC-CF-NEARMISS-DESCENT-INPUTS 25c9f212e5162216 -->
<!-- status-consumer: OP-EC-NEXT 80385bab71bd299c -->

<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-SELMER-PANEL 539bd8ec36b36c44 -->

<!-- status-consumer: EC-K3-ICARM-MW16-POINTED-SIEVE cb83c1afae1d0141 -->

<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-VS-SHA-COMPARISON f37417a9fda3ee3f -->

<!-- status-consumer: EC-K3-ICARM-MW16-SENSITIVITY f88886c066d6cb45 -->

<!-- status-consumer: EC-FIXED-FIELD-COMPARISON 02c49a8120aeb7bd -->
