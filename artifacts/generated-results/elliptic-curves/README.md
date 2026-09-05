# Active elliptic-curve artifacts

The shared runtime also retains the [MW18 census witness bundle](runtime_mw18_census_witnesses_v1.zip)
and [replay summary](runtime_mw18_census_witnesses_v1.json), with
[positive MW states and finite facts](runtime_mw18_complete_replay_v1.json). The default
`certify_r17_extreme_anchored_mw18_covers.sage --check` verifies its polynomial
group-law/chord witnesses and the eight positive covers. See the
[engineering ledger](../../../elliptic-curves/notes/SHARED_RESEARCH_RUNTIME.md).

- [Runtime chart-policy sweep](runtime_chart_policy_sweep_v1.json) and
  [portable exact witnesses](runtime_chart_policy_sweep_v1.zip): 21 controls,
  four policies and 1,008 completed charts; [method and replay](../../../elliptic-curves/notes/SHARED_RESEARCH_RUNTIME.md).

This directory is the complete active artifact surface for
[`elliptic-curves/`](../../../elliptic-curves/README.md). Exploratory outputs and superseded results live in the
[`archive`](../../../archive/elliptic-curves/README.md); large or resumable
runs belong under the ignored `artifacts/local/elliptic-curves/` directory.

[`CATALOG.tsv`](CATALOG.tsv) is the index. Every row gives the exact SHA-256,
an evidence label, the governing `MATH_STATUS.json` identifier when one
exists, and a one-line scope statement. The evidence labels are deliberately
strict:

- `theorem-certificate`, `exact-rank-certificate`, and
  `exact-lower-bound-certificate` are proof-bearing;
- `conditional-bound` states its hypothesis;
- `exact-computation` proves only the calculation described;
- `bounded-experiment` never promotes a negative search to a theorem;
- `partial-reproduction` records precisely what remains dependent on a public
  source or missing local replay;
- `source-transcription`, `reproducibility-fixture`, and `search-plan` carry no
  mathematical conclusion by themselves.

Important distinctions made explicit by the catalogue:

- `universal_pointed_control_regression_v1.json` replays all 1,034 calibrated
  control boxes through the shared API and rechecks all 55 quotient directions.
  `universal_pointed_integration_v1.json` records bounded active-caller canaries.
  See the [shared API and claim limits](../../../elliptic-curves/notes/POINTED_QUARTIC_SEARCH.md).

- `mw16_sensitivity_*v1.json.gz` and their summaries retain the exact
  coordinate/centre calibration and the subsequent gated prospective replay.
  The selected policy recovers 55/55 control directions; all 856 prospective
  height-100,000 boxes complete with no gain or timeout. See the
  [canonical sensitivity note](../../../elliptic-curves/notes/MW16_SENSITIVITY_RECOVERY_2026-09-05.md).

- `exceptional_soluble_vs_sha_comparison_v1.json` extends the five primary
  target rows with 62 marked quartics and a same-curve soluble/Sha control.
  The inherited residual block has dimension 18 and CT rank 16; full
  target complements remain unknown. `exceptional_selmer_feasibility_v1.json`
  preserves the bounded setup failures. See the
  [comparison and proof bounds](../../../elliptic-curves/notes/EXCEPTIONAL_SOLUBLE_VS_SHA_PANEL_2026-09-05.md).
- `exceptional_soluble_selmer_panel_v1.json` certifies 110 independent
  residual point classes across eleven fixed exceptional fibres, with
  exact cubic square characters, explicit quadrics, and rational witnesses.
  Their Cassels--Tate rows vanish against every Selmer class. No complete
  Selmer complement or insoluble Selmer control is certified; see the
  [canonical theory](../../../elkies-k3/RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md).

- `fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json` fixes the
  Fermigier rank-20 cubic field and certified 20-dimensional global Kummer
  span, then computes the complete local intersection on the whole span for
  `u=-2,-1,0,1,2`.  The resulting dimensions are `13,18,20,13,13`; all newly
  bad primes are checked and `u=0` is the full-dimension positive control.
  The artifact uses no class group and records explicit cover inputs, but
  certifies no rational point, new rank lower bound, or full Selmer group.
<!-- status-consumer: EC-FIXED-CUBIC-VARYING-CURVE-LOCAL-KUMMER 46ca45db3e702eb6 -->

- `fixed_field_u_minus1_*v1.json` retain bounded point searches and exact
  conic/quartic model maps. The point `(A+1,A-B+1)` certifies rank at least
  one, with a valuation-parity witness separating its Kummer class from
  the entire inherited space. These search artifacts alone do not classify
  the inherited directions. `fixed_field_point_realization_positive_controls_v1.json`
  replays ten mapped points at `u=0`, including parameter-infinity cases.
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-RANK1 7e488a894d136732 -->

- `fixed_cubic_u_minus1_cassels_tate_v1.json` and its compressed arithmetic
  evidence certify pairing rank 16 on the eighteen inherited classes.
  Exactly three nonzero radical combinations remain point-solving
  candidates; all other classes have a proved Sha obstruction. The summary
  includes the full matrix, symplectic pairs, radical masks and quartics.
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-CASSELS-TATE df45391a84f0e3c9 -->

- `fixed_field_radical_models_v1.json` and its compressed evidence retain
  six globally minimal degree-four models, exact maps, completed p-adic
  lattice searches without hits, and six incomplete four-descent attempts.
  All three radical classes remain **UNKNOWN**.
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-MINIMAL-MODELS 90216b8c456edd20 -->

- `fixed_field_radical_search_geometry_v1.json` certifies that elementary
  real inequalities already excluded all six nominal height boxes. This
  corrects the interpretation of the searches without deciding a class.
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-SEARCH-GEOMETRY 678f7beb805a4530 -->

- [`fixed_field_tangent_conics_v1.json`](fixed_field_tangent_conics_v1.json)
  and its [evidence](fixed_field_tangent_conics_evidence_v1.json.gz) retain
  three exact cubic tangent conics, eight invertible reductions and incomplete
  construction attempts. No genuine lift or new target obstruction is claimed.
<!-- status-consumer: EC-FIXED-CUBIC-TANGENT-CONIC-GATE 26a49e30ff3128d3 -->

- `quotient_geometry_table_v1.json` is the complete 30-presentation join of
  five usable R17 controls, sixteen refreshed R17 ladder fibres, and nine
  A1/MW16 parent presentations.  It stores full displayed-quotient Grams,
  regulators, successive minima, and 230 recovered-direction
  projection/phase/distortion decompositions.  All 230 pointwise projection
  CVPs agree at rounding scales `10^5` and `10^6`.  All thirteen nonempty,
  strictly partial initial recoveries fail the necessary scalar
  successive-minimum-prefix test; seven remain failures at the final stage and
  three containment comparisons remain null.  This is an exact
  rational-subspace comparison built on numerical height data, not a rank upper
  bound or an interval height certificate.

- `icarm_mw16_pointed_sieve_h10000_summary_v1.json` pins the full compressed
  856-chart specialized-sieve replay and the separate initial control ledger.
  Every prospective height-10,000 box completes; 28,134 exact square tests
  yield no finite point. The new coordinates recover twenty initial control
  directions, without a claim to reproduce the historical adaptive 54/55.
  `verify_icarm_mw16_pointed_sieve.py --check --replay-charts` checks source
  hashes, every exact chart map, and the control-group certificates.

- `icarm_mw16_parent_presentation_audit_v1.json` proves that the nine
  complete-A1 hit labels are exactly five fibration classes; repeated labels
  remain nine bounded-height coordinate charts but only five statistical
  observations.  `icarm_mw16_blind_ladder_calibration_v1.json` records exact
  complement-blind recovery of 38/55 demonstrated directions initially and
  54/55 after the completed curve-398 and curve-400 adaptive waves.  The
  prospective height-300 chain retains 104 exact fibres, while
  `icarm_mw16_nagao_finalist_half_lattice_h300_v1.json` records 856/856 chart
  timeouts at its first declared budget.  This is a censored bounded
  experiment, not a negative rank result; no Selmer or expensive-search gate
  opens.

  The distinct target-free replay is
  `a1_mw16_target_free_parameter_candidates_h300_v1.json` followed by
  `a1_mw16_target_free_parameter_search_h300_v1.json`.  Its sampler reads the
  anonymous family template rather than target-bearing parent rows.  It
  produces 104 pairwise nonisomorphic fibres, and direct exact searches finish
  all 856 maximum-depth charts through height 100,000 with zero affine points,
  timeouts, failures, or quotient gains.  This is a completed bounded null
  experiment, not a rank upper bound.
  `mw16_short_models_h300_v1.json.gz` stores all global minimal models,
  explicit section transports, renewed independence certificates and 856
  quartic maps; `mw16_short_models_h300_summary_v1.json` summarizes sizes.
  `mw16_short_models_chart_benchmark_v1.json` records the tiny fixed-budget
  benchmark. The canonical note explains why the intrinsic `j` sizes rule
  out an orders-of-magnitude bit-length reduction on these same fibres.
<!-- status-consumer: EC-K3-ICARM-MW16-BLIND-LADDER acfa3bdcebb18137 -->

- `r17_refresh_jump_ladder_protocol_v1.json` and
  `r17_refresh_jump_ladder_blind_v1.json` preserve the stopped cross-class
  assertion and the already sealed curve-478 response.  The v2 protocol
  transparently corrects the initial class set to the top 43 exact generic
  depths before the other fifteen outcomes.  Its blind artifact contains only
  exact pre-complement rank responses; the verifier opens public complements
  afterward, and `r17_refresh_jump_ladder_analysis_v2.json` records passing
  exact ordinal and `+10/+11/+12` tail endpoints plus post-freeze
  fibration/`j`-class and `q>=11` sensitivity analyses.  The latter limit
  scheduling use to the calibrated norm-twelve R17 setting and do not validate
  alternate-Q80 extreme-tail transfer.  These are fixed-panel detector results,
  not a full-rank, saturation, Selmer, or population theorem.
<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER a2d7034fb8977c18 -->

- `half_lattice_height_compression_analysis_v1.json.gz` reconstructs 3,865
  detailed half-lattice chart maps and their presearch lattice, coefficient,
  invariant, reduction, and known-basis distortion data.  Its explicitly
  posthoc rank-28 audit exactly replays all 2,560 chart/target source decisions
  and every prefix quotient gain from reduced-coordinate height, while the
  target-free scalar comparison remains negative across eleven positive chart
  orders after adding 394 compact `+8` through `+12` control records.  This
  certifies the midpoint/old-point-exclusion mechanism, not a prospective
  success probability or any conclusion from a bounded miss.
<!-- status-consumer: EC-HALF-LATTICE-HEIGHT-COMPRESSION 3baeaf370aec751c -->

- `curve385_height_compression_pilot_protocol_v1.json` freezes 16 charts from
  32 fresh scale-stable local maxima found by deterministic ascent in the full
  current `M29/2M29` parity space.  The bound-100,000 result
  `curve385_height_compression_pilot_blind_v1.json` completes all sixteen with
  zero finite points and exact rank `29 -> 29`.  This supports old-point
  exclusion only; it is not a point-absence or saturation result.
<!-- status-consumer: EC-CURVE385-HEIGHT-COMPRESSION-BUILDER-PILOT c0d6f2d67018def4 -->

- `icarm_curve398_two_parent_collision_v1.json` compiles the second exact
  A1/MW16 survivor for curve 398 and proves that it is the first fibration after
  an affine `PGL2(Q)` base change and constant `Q`-Weierstrass scaling.  Its
  specialized generic group equals the first presentation's group integrally.
  The intersection and sum both have rank 16; the quotient in public `M30` is
  torsion-free `Z^14`, so the index is infinite.  This is a deduplication and
  basis-consistency certificate, not a transversality or rank-upper-bound
  result.

<!-- status-consumer: EC-K3-R17-TRAINING-EXACT-ARITHMETIC-GROUP-GATE 427bf822e774c81e -->

- `r17_training_arithmetic_group_inputs_v1.json.gz` is an outcome-free compact
  freeze of the R17 development parameters, split labels, selection/holdout
  membership, and controls.  Its exact audit
  `r17_training_arithmetic_group_audit_v1.json` finds 100,000 distinct rational
  `j`-classes, no labelled/holdout twist overlap, and no control/development
  twist overlap.  The learned-score laboratory entry now fails closed unless
  this audit explicitly authorizes the unchanged v1 score artifact.
- `curve385_sparse_restart_budget_v2.json` is an operational amendment for
  future curve-385 sparse rank-32 searches.  It preserves the frozen v1
  protocol and primary no-growth evidence while replacing their combined
  four-state allowance with independent limits of three rank-changing and
  four saturation-only group changes.  Its adverse-path regression admits two
  saturation-only changes before the three unit gains from rank 29 to 32.  A
  budget stop still proves no rank upper bound or saturation statement.
<!-- status-consumer: EC-K3-R17-CURVE385-INDEPENDENT-RESTART-BUDGETS 39cfce110e3e494f -->
- `r17_mw17_only_selmer_control_inputs_v1.json` freezes the strict prospective
  replay for record fibres 356 and 385.  Its source-hash-pinned Magma programs
  contain only each minimal curve and exactly seventeen MW17 points; mechanical
  audits exclude the twelve held-out coordinate rows, labels, half-ideals,
  external reads, and cover/point searches.  The associated
  `r17_mw17_only_selmer_control_run_v1.json` ledger records zero completed
  replays and leaves the operational Selmer candidate gate false.
- `elkies_2026_record_pair_relative_2selmer_inputs_v1.json` is the older exact
  fixture-sequenced Magma input for record fibres 356 and 385, including the
  specialized MW17 controls and twelve held-out points.  It is a pinned input
  for post-discovery diagnostics, not the prospective control and not a
  completed descent or Selmer certificate.  Its whole-file SHA-256 is
  `a0492c02910c035c9702a10224132ddbccd47236089d1ef7c647108b132b9e92`.
<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-RANK-ESCAPE-DETECTOR-V2 eda7a0053b31b7c9 -->

- `latent_lattice_calibration_truth_v1.json` contains exact withheld control
  embeddings; it is not selector input. `latent_lattice_calibration_v2.json`
  is the active corrected-semantics gate artifact.  It still rejects the
  bounded selector before any target search.  The superseded `v1` selector
  bytes are retained only for provenance and are not a nonexistence theorem
  for a common wgxli lattice.
  `latent_lattice_finite_calibration_v1.json` adds disjoint development and
  held-out finite-code ensembles.  It passes R17 proposal recall but fails the
  Fermigier control and finite-profile selection, so it also forbids target
  use.
  `latent_lattice_shape_calibration_v1.json` adds exact multiplicity-preserving
  cross-bound enclosure intersections.  It puts the primitive Fermigier
  rank-12 truth at rank 65 in a bounded top-128 ledger and exactly recovers the
  held-out R17 rank-25 control, but symmetric R17 selection still fails.  Its
  status is `PASS_PROPOSAL_CALIBRATION_SELECTOR_FAIL`, so it likewise forbids
  target use.
  `latent_lattice_relation_consensus_v1.json` is a supervised exact-signal
  benchmark: every R17 leave-one-out split retains a rank-17 coefficient core
  on the held-out fibre, while exact rational-ray normalization records the
  fibre-specific Fermigier saturation denominators.  The published embeddings
  supply the alignments, so its `PASS_CONTROL_EXACT_RELATION_SIGNALS` status
  does not authorize blind target use.
  `latent_lattice_hypergraph_matcher_v1.json` certifies exact equal-core
  rebasing and primitive rectangular `17 x r` lifts for all four supervised
  held-out controls.  Full-cloud replay sees 238, 266, 291, and 304 rays.
  `latent_lattice_metric_relation_search_v1.json` is the complementary blind
  bounded experiment: 39,714 states and 500 lift attempts reach at most 49
  replayed rays against a 100-ray gate.  Its
  `FAIL_BLIND_R17_RECOVERY_GATE_CLOSED` status keeps all targets gated.
  `latent_lattice_partial_replay_v1.json` adds exact proper-subspace
  saturation and replay.  A supervised 103-ray rank-16 path lifts to a
  primitive target image and replays 194 rays and 318 relations; its six-block
  finite signature is exact.  The independently bounded oracle-center beam
  makes 400 partial audits but sees at most 29 replayed rays and no full
  embedding.  Its `PASS_EXACT_PARTIAL_REPLAY_SELECTOR_FAIL` status therefore
  keeps the same gate closed.
  `latent_lattice_star_component_v1.json` tests a whole center star with exact
  mod-2/mod-3 rank pruning.  The supervised visible truth star has rank 11 and
  a primitive 32-ray replay, while the 512-state, 500-audit bounded ledger
  reaches withheld overlap only 9/11.  Its
  `FAIL_STAR_COMPONENT_RECALL_GATE_CLOSED` status is another proposal-generator
  failure and does not authorize targets.

- `icarm_273_282_302_family_discovery_v1.json` screens 2,334 generated
  one-parameter families exactly.  It rediscovers curve 282 in both the
  canonical Fermigier coordinate `u=11671/42` and the generated six-root
  Mestre coordinate `T=11671/21`, verifies `Q`-isomorphism in both models, and
  finds no match for curves 273 or 302 in this bounded construction space.
  Its SHA-256 is
  `5d19571dc74f9e8c270bcaaa943f19e3876bb2a13ffc33c37d479a3206fd8770`;
  the producing script's pinned SHA-256 is
  `dbb4ae2580c45d4a25800225e11301c51158c97253b3c3631054a15c60769101`.
  The related fixed-root fingerprint ledgers are
  `icarm_construction_fingerprints_v1.json` for curves
  273/281/282/285/286 and `icarm_construction_fingerprints_v2.json` with curve
  302 added. Their historical repository-model diagnostic covers only
  uncompressed JSON beneath this generated-results tree, not gzip, archive,
  or all repository models. Run
  `python3 elliptic-curves/cas/audit_icarm_construction_recognition_artifacts.py`
  to check these committed bytes and scope flags without regenerating a
  family, sieving, or factoring.
- `icarm_curve282_conductor_parameter_recovery_v1.json` verifies global
  minimality, selected local Tate data, two-chart discriminant-root profiles,
  bounded CRT/Gauss recovery of `u=11671/42`, and the exact Fermigier/target
  `j`-identity.  Its exact leaf residues are replay inputs; the valuations
  alone determine only coarser p-adic balls.
- `record_first17_subgroups_v1.json` exactly computes the first-seventeen
  coordinate, quotient, finite-Kummer, and bad-component codes for curves 273
  and 302, while its canonical-height/theta layer is numerical at 100 digits.
  Its saturation index one is only inside each displayed subgroup, not in the
  full Mordell--Weil group.
- `icarm_curve273_rank30_v1.json` is the pinned `ECR30` lower-bound
  certificate. It records 30 exact point checks and a stored `31 x 30`
  finite-quotient certificate, while explicitly excluding numerical heights
  from the proof and making no exact-rank claim. Its SHA-256 is
  `e2a7a322fbd4703af4239f497749a69a68f9d5149aa8a1f696b39ab3941a3284`.
- `icarm_curve302_rank31_v1.json.gz` is the deterministic compressed `ECR31`
  certificate. It proves rank at least 31, not exact rank 31; its public
  BSD/GRH statement and numerical regulator are not used. Its compressed
  SHA-256 is
  `fc50b4b9ec5fe1dd1fe31aa299f13d8bc3476d43f3ed98e2ade5a4fc8972aa04`
  and its decompressed JSON SHA-256 is
  `3be0d6fe82c58e0f9284df5d9340332944a1d906508ea986d4abe00357036991`.
  `python3 elliptic-curves/cas/audit_icarm_rank_lower_bound_artifacts.py`
  checks both pinned files and their source provenance without performing
  curve arithmetic, finite-group enumeration, or matrix-rank computation.
- The ICARM 285/286 analysis exactly proves independence of 21 displayed
  points and now independently replays global minimality and every local
  conductor exponent.
- `icarm_curve394_rank21_v1.json` specializes the compact Elkies R17 family at
  `t=3/8`, proves a generic-17 plus public-4 basis independent, and replays the
  exact conductor locally.  It proves rank at least 21 at
  `log(N)=166.252098...`, not exact rank 21.
- `conductor_first_near_miss_descent_targets_v1.json` pins exact known-subgroup
  Kummer inputs for four fixed fibres; it does not claim a complete Selmer
  group or rank upper bound.
- `conductor_first_family_anchor_pilot_v1.json` is the closed 27-fibre
  Fermigier/Mestre anchor-neighborhood ledger. Nine sieve survivors have
  exact global/local Tate data and full-dimensional known mod-2 Kummer
  subgroups. The `u=481` fibre has an exact rank-at-least-14 point certificate;
  no complete Selmer or residual-cover classification is claimed.
- The Elkies compact-`t` artifacts certify the rank-25--28 positive controls
  and the complete height-10000 three-block calibration. The PARI and eclib
  rank-28 residual-descent artifacts are strict timeouts with no Selmer bound;
  both explicitly forbid point search.
- `r17_frozen_nagao_shell_h10001_30000_v1.json` exhaustively scores all
  972,697,152 primitive parameters in the disjoint `10000 < H <= 30000` shell
  with that unchanged three-block rule and pins pooled-Nagao and deterministic
  random controls. It is a complete heuristic ranking, not rank evidence.
- `r17_frozen_nagao_shell_search_v1.json` records the exact preexisting-
  bisection-atlas outcomes on the three matched 128-row lanes. Seven distinct
  fibres have one certified quotient direction beyond the generic 17; one
  split row is censored. No unrestricted point search or Selmer result is
  claimed.
- `elkies_2026_relative_2selmer_suite_inputs_v1.json` pins raw-basis
  unconditional Magma jobs for the rank-21 and rank-25--28 controls and the
  first ten frozen high-Nagao candidates. The paired
  `elkies_2026_relative_2selmer_suite_run_v1.json` records that Magma is not
  available on this host. These are exact inputs and a backend audit, not a
  completed Selmer computation.
- `elkies_2026_relative_2selmer_open_rank21_300s_v1.json` is the first run of
  the open-source Sage/PARI `ell2cover` replacement on the rank-21 control.
  With twelve proved factor hints, a 2 GB PARI stack, and a 4 GB RSS envelope,
  it reaches the strict 300-second limit inside `ellrankinit` at 440,283,136
  bytes peak observed RSS. It returns no BNF certificate, cover basis, or
  Selmer dimension and is explicitly incomplete.
- `elkies_2026_relative_2selmer_open_nagao0001_120s_v1.json` applies the same
  frozen open method to the top high-Nagao candidate `t=-5643/6760`. It reaches
  the strict 120-second limit inside `ellrankinit` at 230,608,896 bytes peak
  observed RSS, so it contains no Selmer dimension or candidate-promotion
  evidence.
- `elkies_2026_relative_2selmer_checkpointed_rank21_c1p1_600s_v1.json`
  calibrates the separated BNF stage on the rank-21 control with PARI
  `tech=[0.1,4,20]`. Relation collection reaches the strict 600-second limit
  at 256,798,720 bytes peak observed RSS. Its status is `INCOMPLETE_BNF`; the
  later Selmer, embedding, and quotient-cover stages were not run and no rank
  or Selmer bound is claimed.
- `elkies_2026_relative_2selmer_global_engine_benchmarks_v1.json` records two
  further 600-second rank-21 global-field failures: a reduced-field,
  locally tuned PARI 2.17.4 build and exact Hecke 0.40.2 method 2. Both stop
  before the class group is available, so neither run reaches the Selmer,
  embedding, or covering stages.
- `elkies_2026_relative_2selmer_open_bottleneck_benchmarks_v2.json` extends
  that fail-closed benchmark with four PARI technical-parameter trials, two
  Hecke bound-240 trials, an aggressive PARI restart build, and three exact
  bounded relation ledgers. It also pins three archimedean Hecke variants, an
  exact factor base augmented by all 25 S-prime ideals, and a direct
  S-multiplier relation collector. It now also records the official PARI 2.19
  development branch's six-parameter threaded engine: the best 300-second
  run reaches a 1,996-ideal factor base but retains a 1,635-relation request.
  The third ledger uses exact multi-large-prime sparse-hypergraph elimination
  at factor-base bound 1,000; its 666 closed dependencies give zero rank gain
  modulo the canonical S-span.
  None completes the global class/unit group; relation deficits and residual
  dimension 34 are explicitly not Selmer bounds.
- `elkies_2026_known_kummer_quotients_controls_v1.json` and
  `elkies_2026_known_kummer_quotients_suite_v1.json` are class-group-free
  exact lower-bound audits. Residue squareclasses certify known mod-2 ranks
  `21,25,26,27,28`, hence exceptional dimensions `4,8,9,10,11`, on the five
  controls and generic rank 17 on all ten high-Nagao candidates. They do not
  compute the full Selmer quotient, unknown classes, or blind recovery.
- `elkies_2026_known_exceptional_quotient_covers_v1.json` pins all 42 basis
  covers and hashes the full local enumeration of 3,851 nonzero classes in
  the known exceptional control subgroups. Every explicit intersection of
  quadrics has a verified rational witness. These are realized positive
  controls; the artifact does not assert that the known subgroup exhausts the
  relative Selmer quotient.
- `elkies_2026_bisection_specialization_controls_v1.json` is the complete
  195,600-test evaluation of the 39,120 equation-level bisections at the four
  rank-25--28 controls and ICARM curve 394. It finds split counts
  `6,3,2,1,25` and known-complement class-span dimensions `5,3,2,1,4`, with
  no finite-quotient escape. Full-rank relation blocks verified by exact group
  addition prove generated-subgroup ranks `25,26,27,28,21`; these do not give
  upper bounds for the full curves.
- `elkies_2026_bisection_visibility_record_curves_v1.json` row-reduces those
  classes into deterministic visible and complementary quotient bases.  It
  records the ten-dimensional rank-28 target packet and the exact mechanism
  boundary that translated trace shells repeat existing bisection classes.
  It also proves, through degree-24 irreducibility witnesses, that the 2024
  rank-29 curve and ICARM 273, 302, and 398--400 are not rational fibres of the
  published `R17` fibration, including after quadratic twisting.  The rank-28
  control recovers `5471*t+9529`; other fibrations, families, and isogeny
  constructions remain open.
- `elkies_2026_deep_cover_exceptional_quotients_v1.json` gives exact integral
  `L_t/M_t` coordinates for every split degree-two point at the rank-25--28
  controls and combines them with 69 sampled norm-20 trisections, all 160
  norm-26 deep trisections, and 53 sampled norm-34 quadrisections.  None of
  the 282 new covers has a rational component at a control, so the captured
  ranks remain `5,3,2,1` through degree four.  The norm-20 and degree-four
  layers are deterministic samples, not complete non-splitting theorems.
- `elkies_2026_rank28_bad_place_kummer_ledger_v1.json` proves the complete
  factorization of the rank-28 2-division cubic discriminant and contains all
  thirteen finite/2-adic/real local blocks for the generic seventeen points.
  Their combined coordinate rank is 15; this is a known-image calculation,
  not an ambient `K(S,2)`, local-solubility, or Selmer certificate.
- `elkies_2026_rank28_residual_2selmer_pari_factored_8g_v1.json` is the first
  PARI descent supplied with that factor certificate. Its strict 600-second
  run reached 5,698,514,944 bytes peak observed RSS but returned no Selmer
  dimension, so it remains search-forbidden. The longer supervised artifact
  has suffix `_8g_30min_v1`; its strict 1,800-second run reached
  6,040,723,456 bytes peak observed RSS and likewise returned no dimension.
- `elkies_2026_rank28_s_class_pari_v1.json` isolates the next PARI stage. The
  exact factor-supplied maximal-order setup completes, but the strict
  120-second run stops inside class-group relation generation before
  `bnfcertify`; it is not an `S`-class or Selmer bound.
- `elkies_2026_rank28_s_class_pari_polredabs_v1.json` repeats that envelope on
  PARI's exact reduced cubic with an explicit original-generator map. The
  polynomial order index drops by 27, but the run reaches the same 153-request
  random-relation plateau and returns no class-group bound. It remains
  search-forbidden.
- `elkies_2026_rank28_bnf_free_s_class_pilot_v1.json` audits 172 exact
  canonical principal rows after a paired-special-ideal pilot. Its displayed
  factor-base quotient dimension 141 is explicitly uncertified because bound
  1,000 is below the 1,202,640 Bach/ERH generation threshold; the pilot found
  no noncanonical relation and forbids point search.
- `elkies_2026_rank28_generic17_local_signature_v1.json` recomputes all 53
  bad-place coordinates for the generic 17 and public complement 11. The
  generic and full rank-28 local-signature ranks are both 15, with zero
  incremental rank from every exceptional direction; this is exact evidence
  that known Kummer signatures are not a Mordell--Weil quotient or Selmer
  bound.
- `elkies_2026_rank28_generic17_local_coverage_v1.json` certifies full known-
  point coverage at four of eleven odd bad primes and at the real place. The
  other seven odd places and the two-adic place remain unresolved.
- `elkies_2026_rank28_norm_one_local_pilot12_v1.json` is a bounded selected-
  place audit of 12 of 49 norm-one cover candidates. Of 84 cover/place tasks,
  60 have certified local points and 24 remain inconclusive; none has a
  certified local obstruction. This is not class-group completion, everywhere
  local solubility, Selmer membership, or search authorization.
- `elkies_2026_rank28_public11_selmer_candidates_v1.json` records the exact
  cubic Kummer classes `X(Q)-theta` of the eleven certified public complement
  points. Their norms are displayed squares and their quotient independence
  proves a residual 2-Selmer lower bound of 11.
- `elkies_2026_rank28_public11_two_cover_controls_v1.json` materializes those
  eleven intersections of quadrics and verifies `[1:0:0:1]` on each one.
  `elkies_2026_rank28_public11_global_cover_witness_audit_v1.json` rechecks the
  witnesses and records local solubility at every place. These are positive
  controls, not an ambient Selmer enumeration or upper bound.
- `elkies_2026_rank28_construction_fingerprints_v1.json` reverse-engineers the
  eleven public complement directions relative to specialized R17.  It records
  stable rounded-Gram closest representatives, quotient-height profiles, exact
  Kummer and bad-place codes, irreducible `[2]`/`[3]` preimage fields and their
  Galois/Frobenius data, the trivial rational-isogeny graph, symbolic generic
  R17 multiplication covers, and an explicitly heuristic consensus clustering.
  The cubic Kummer map exposes all eleven known directions but is a
  representation of known points, not yet a construction of unknown ones.
- `newfamily_rank14_t83_6_v1.json` proves only rank at least 14; the separate
  `newfamily_rank14_t83_6_pari_exact_rank_v1.json` supplies the PARI interval
  `[14,14]` used for the exact-rank statement.
- Bounded Fermigier searches are indexed as experiments even when every
  calculation inside the declared box is exact.

Run the catalogue/archive integrity check with:

```sh
python3 elliptic-curves/scripts/audit_artifact_catalog.py
```

The pre-cleanup bytes and every provenance-only refresh are recorded under
[`archive/elliptic-curves/`](../../../archive/elliptic-curves/README.md).

<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-SELMER-PANEL 539bd8ec36b36c44 -->

<!-- status-consumer: EC-K3-ICARM-MW16-POINTED-SIEVE cb83c1afae1d0141 -->

<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-VS-SHA-COMPARISON f37417a9fda3ee3f -->

<!-- status-consumer: EC-K3-ICARM-MW16-SENSITIVITY f88886c066d6cb45 -->
