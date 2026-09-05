# Active CAS modules

This directory is intentionally narrower than the computational archive. It
contains exact status checkers and their shared arithmetic, plus code for the
current rank-32, low-conductor, residual-Selmer, and K3-construction gates.
Stable user-facing commands are listed in [`../scripts/`](../scripts/) and
[`../REPRODUCE.md`](../REPRODUCE.md).

## Start here

- [Shared arithmetic and search runtime](../notes/SHARED_RESEARCH_RUNTIME.md):
  cached contexts, subspace-first descent, lazy MWState search, regulator gates,
  common worker supervision and portable proof replay. New searches start with
  `run_mw_search.py`; full BNF is an explicit upper-bound option.

- [`pointed_quartic_search.py`](pointed_quartic_search.py) is the single active
  half-lattice backend. MW16/MW17, curve-specific and zero-gain callers use
  its adapters; `run_pointed_quartic_search.py` accepts future-family manifests
  and exact MW18 specializations. See [API, migration and regression controls](../notes/POINTED_QUARTIC_SEARCH.md).

- `compare_exceptional_soluble_vs_sha.sage`: replay the five large-jump
  soluble subspaces, 63 marked quartic transports, and the rank-16 Sha
  control; `--check` includes exact CT arithmetic and probe-log hashes.
  `run_exceptional_selmer_feasibility.py` is the separate three-curve,
  45-second-per-curve cached arithmetic probe with known factor support;
  complete Selmer is explicit and point search is disabled.
- `certify_exceptional_soluble_selmer_panel.sage`: exact known soluble
  residual subspaces and 110 witnessed 2-covers on eleven exceptional
  fibres, without a point search or full descent. Use `--check` to replay.
- `check_icarm_curve302_rank31_pinned.py`: deterministic rank-at-least-31
  replay and compressed-artifact hash check.
- `verify_icarm_curve356_rank29.py`: independent rank-at-least-29, minimal
  model, conductor, and local-reduction replay for the new rank-29 size record.
- `analyze_icarm_curve356_lineage.py`: hash-pinned curve-351/356 ordered
  denominator and numerical height-Gram comparison; this does not identify a
  family.
- `verify_icarm_curve273_rank30.py`: independent rank-at-least-30 replay.
- `verify_icarm_curve398_rank30.py`: independent rank-at-least-30, trivial
  torsion, singleton rational-isogeny-class and complete semistable local-data
  replay, bad-node fingerprint, and exact exclusion from the one
  equation-explicit A1/MW16 control family.
- `verify_icarm_curve398_two_parent_collision.sage`: independent Sage replay
  of the exact base/Weierstrass equivalence, curve-398 parameter transport,
  group laws, integral basis transition, ranks, and Smith quotient for the two
  A1/MW16 survivor presentations.
- `prepare_icarm_mw16_parent_ladder_inputs.sage` and
  `run_icarm_mw16_parent_ladder_blind.sage`: reconstruct all nine complete-A1
  hit presentations and run the complement-blind exact maximum-depth MW16
  calibration.  Responses are nested within five target curves.
- `audit_icarm_mw16_parent_presentations.sage`: exact base-change,
  Weierstrass-scaling, target-transport, and integral-subgroup audit proving
  that the nine labels are five fibrations.
- `run_icarm_mw16_curve400_adaptive_calibration.sage` and
  `verify_icarm_mw16_blind_ladder_calibration.py`: complete the curve-400
  124-chart adaptive wave and compare the five blind responses with atlas
  jumps only after the search.  Best recovery is 54 of 55 demonstrated
  directions; this is purposive detector calibration, not population
  inference.
- `sieve_icarm_mw16_parent_presentations_nagao.py`,
  `specialize_icarm_mw16_nagao_finalists.sage`,
  `run_icarm_mw16_nagao_finalist_half_lattice.sage`, and
  `merge_icarm_mw16_nagao_finalist_half_lattice_shards.py`: enforce the
  prospective local-ordering to exact-half-lattice gate on 104 fibres.  The
  first 856 chart attempts are wholly timeout-censored on very large
  coefficients, so none advances to Selmer or unrestricted point search.
- The finalist runner defaults to `PointedQuarticSearch` with calibrated
  weight-16 coordinates and the shared GMP worker. The new API reproduces
  all 1,034 frozen control boxes and 55 quotient directions. Historical
  `replay_mw16_sensitivity.sh` commands use the pinned pre-migration revision;
  `replay_pointed_quartic_snapshot.py` checks its self-contained bundles.
  See [sensitivity recovery](../notes/MW16_SENSITIVITY_RECOVERY_2026-09-05.md).
- `extract_a1_mw16_family_template.py`,
  `build_a1_mw16_target_free_parameter_candidates.sage`,
  `run_a1_mw16_target_free_parameter_search.sage`, and
  `merge_a1_mw16_target_free_parameter_search_shards.py`: run the actual
  target-free experiment on anonymous family data.  The frozen height-300
  sampler produces 104 pairwise nonisomorphic fibres; direct exact searches
  complete all 856 deepest MW16 charts through height 100,000 with no affine
  point, timeout, failure, or recovered quotient direction.  This is bounded
  search evidence only.
- `prepare_mw16_short_models.sage`, `mw16_model_size.py`, and
  `verify_mw16_short_models.sage`: global minimal models, all sixteen section
  maps and fresh independence certificates for the 104 finalists; arithmetic
  coordinate selection on all 856 quartics; a nine-chart benchmark at two
  seconds per call. See the
  [model-size audit](../notes/ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md#exact-arithmetic-model-audit).
<!-- status-consumer: EC-K3-ICARM-MW16-BLIND-LADDER acfa3bdcebb18137 -->
- `analyze_record_first17_subgroups.py`: exact first-seventeen coordinate,
  quotient, finite-Kummer, and bad-component comparison for curves 273 and
  302, plus a 100-digit canonical-height/theta profile.
- `analyze_icarm_7fff_zip_sequence.py`: exact independence replay for the
  public curves 281, 282, 285, and 286, plus repository-local global
  minimality and complete local conductor reconstruction for 285 and 286.
- `verify_icarm_curve245_rank20.py`: fully local low-conductor rank-at-least-20
  certificate.
- `verify_icarm_curve394_rank21.py`: compact-`t=3/8` specialization identity,
  exact generic-plus-public rank-at-least-21 certificate, and complete local
  conductor replay for the current ICARM rank-21 conductor anchor.
- `certify_mestre_dsquare_rank19_frontiers.py`: exact rank-at-least-19
  low-conductor frontiers and conditional fixed-fibre diagnostics.
- `run_conductor_first_pari_diagnostic.py`: bounded PARI `ellrank` supervisor
  for the four conductor-first near misses; provisional upper endpoints stay
  GRH-conditional, while returned points receive independent exact mod-2
  certification.
- `fixed_cubic_field_curve_family.py` and
  `run_fixed_cubic_field_curve_family.sage`: keep the Fermigier rank-20 cubic
  field and its certified 20-dimensional Kummer span fixed while varying the
  curve through `alpha_u=theta+u*theta^2`.  The bounded `|u|<=2` integer run
  computes exact whole-span local kernels `13,18,20,13,13`, checks every new
  bad prime, and emits explicit covering inputs without a class-group
  computation.  It proves no point realization or rank on a new curve.
<!-- status-consumer: EC-FIXED-CUBIC-VARYING-CURVE-LOCAL-KUMMER 46ca45db3e702eb6 -->
- `run_fixed_field_point_realization.py`: minimizes equivalent conic/quartic
  presentations of the surviving `u=-1` covers, searches basis combinations
  and translations by the exact point `(A+1,A-B+1)`, and replays every map
  and actual Kummer identity. Valuation parity above 19 proves this point
  independent of the inherited span; bounded misses alone do not classify
  inherited classes. The following pairing certificate now selects the
  remaining candidates.
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-RANK1 7e488a894d136732 -->
- `build_conductor_first_s_class_envelopes.py`: exact four-target comparison
  of cubic-field discriminants and materialized Bach/ERH factor-base sizes;
  this orders the BNF-free relation collectors but is not a class-group or
  Selmer computation.
- `run_fermigier_rank20_fixedfb_quadratic_specialq.py`: reusable
  quadratic-in-theta special-q collector for a declared cubic/S-prime set;
  its optional sparse-hypergraph engine cancels multi-large-prime columns
  while retaining exact principal-generator witnesses, and its bounded
  adaptive mode can reuse residual degree-one ideals as new special ideals.
- `run_fermigier_rank20_minkowski_specialq.py`: full-ideal Minkowski
  special-q collector.  Its hybrid path trial-divides beyond the factor base,
  proves retained probable-prime cofactors, can batch-GCD unresolved composite
  cofactors using a product/remainder tree with a narrow exact pairwise
  fallback, and its sparse-hypergraph path retains every exact partial edge as
  well as closed dependencies. Primitive projective normalization rejects
  rational-multiple fake cycles. It supplies proved declared discriminant
  factors to PARI before maximal-order construction. This is relation
  collection, not class-group completion or a Selmer calculation.
- `merge_bnf_free_minkowski_relation_ledgers.py`: replays the retained sparse
  edges from compatible independent runs in one hypergraph, deduplicates
  projective generators across runs, and stores exact generator witnesses for
  any genuinely cross-run cycles. A forest remains a negative checkpoint, not
  evidence of class-group completion.
- `select_bnf_free_minkowski_feedback_specials.sage`: ranks repeated residual
  prime-ideal vertices in a mergeable ledger, reconstructs their degree-one
  residues exactly, and emits a reproducible adaptive special-q seed list.
  This is graph scheduling only; it proves no relation or arithmetic bound.
- `refine_r17_unresolved_ideal_vertices.sage`: reuses the record-pair
  Minkowski caches without completing their norm factorizations. It divides
  only proved factor-base/shared prime ideals, retains the exact reduced
  residual ideal as a sparse vertex, and accepts only replayable incidence
  dependencies. Exact matches with certified MW29 half-ideal classes are
  killed, but the bounded factor base supplies no global upper bound.
- `close_r17_residual_ideal_vertices.sage`: targets those residual vertices
  directly. Directional or projectively deduplicated Minkowski samples
  `beta in I` factor the smaller quotient `(beta)/I`; certified MW29 source
  ideals are killed before elimination. The same run adds principal `(p)`
  rows for every used outside rational prime.
- `augment_r17_targeted_closure_canonical.sage`: cache-only replay for target
  ledgers made before the outside `(p)` rows were added. It rechecks every
  stored target factorization and adds those exact canonical rows without
  rerunning a lattice search or integer factorization.
- `certify_nagao_rank20_t5081.py`: exact Nagao rank-at-least-20 certificate.
- `newfamily/certify_rank_t83_6.py`: exact-rank-14 Sage/PARI replay.
- `elkies_residual_selmer_gate.py`: fail-closed rank-32 residual-dimension
  policy. Complete descent controls theorem claims; an incomplete monotone
  sieve may reject on a proved upper bound below 15 or authorize only a search
  with explicit finite limits. Missing BNF data mean no finite bound yet.
- `run_elkies_2026_rank28_residual_selmer.py`: resource-bounded genuine PARI
  or Selmer-only eclib 2-descent on the public rank-28 positive control. Its
  `pari-factored` backend consumes the separately proved complete
  2-division-discriminant factorization and records a strict memory/wall cap.
- `build_elkies_2026_rank28_bad_place_ledger.py`: proves the complete
  2-division-discriminant factorization and computes the generic-seventeen
  Kummer images at every bad finite place, at 2, and at infinity. The complete
  known-image ledger is exact input to descent, not a Selmer upper bound.
- `run_elkies_2026_rank28_s_class_pari.py`: isolates factor-supplied maximal-
  order certification, class-group relation generation, and `bnfcertify` in
  an owned resource-bounded worker. Its `--field-model polredabs` path records
  an exact reduced polynomial and generator map before the same certified
  maximal-order calculation. A completed class quotient would still precede
  the separate norm and local-solubility gates.
- `run_fermigier_rank20_minkowski_specialq.py --elkies-rank28`: feeds the same
  proved factor support into the exact BNF-free principal-relation collector.
  Its factor-base-1000 paired-cover pilot is explicitly uncertified and does
  not authorize search.
- `build_elkies_2026_rank28_local_coverage.py`: recomputes the bad-place
  signatures of the generic 17 and public complement 11. The eleven globally
  independent exceptional directions add zero local-signature rank, an exact
  positive control showing that signature rank is not a Mordell--Weil or
  Selmer quotient dimension.
- `audit_bnf_free_local_kummer_coverage.py`: compares known-point projections
  with exact odd-prime and real local-dimension bounds. Equality certifies a
  full local image at that place; every strict inequality and the two-adic
  place remain unresolved.
- `run_bnf_free_two_cover_local_supervisor.py`: replaces monolithic bounded
  cover audits with owned, resumable `(cover,p)` workers. Each cache block is
  input- and limit-bound; smooth lifts and obstructions are retained, while
  timeouts and state caps stay mathematically inconclusive.
- `build_elkies_2026_rank28_public_selmer_controls.py`: converts the certified
  public complement eleven into exact classes `X(Q)-theta`, square norms, and
  rational two-cover witnesses. Together with the finite-reduction quotient
  certificate this proves a residual 2-Selmer lower bound of 11 and validates
  the cover layer on genuine classes; it supplies no upper bound.
- `replay_half_lattice_search_ablation.sage`: fixture-blind equal-budget
  comparison of generic/specialized deep, union, five SHA-random, median, and
  shallow half-lattice chart sets. It checkpoints development and sealed
  holdout phases separately, uses identical per-cover reduction/search limits,
  and records parent-plus-child CPU cost.
- `verify_half_lattice_search_ablation.sage` and
  `summarize_half_lattice_search_ablation.py`: load public points only after
  hashing the blind artifact, compute exact quotient ranks over `Q` and
  mod 2 for every arm, and emit the compact cross-phase comparison. Bounded
  misses remain bounded-search misses.
- `half_lattice_chart_policy.py`: executable state-binding and interpretation
  contract for all half-lattice chart orders. It records that the quartics are
  birational point-search charts, gives legacy depth/old-deep/Hamming fields
  ordering meaning only, rejects cached orders after a basis, lattice,
  height-form, quotient-coordinate, or chart-universe change, and forbids all
  absence, rank-upper-bound, and Selmer inference from misses.
- `production_search_gates.py`: keeps theorem exclusions and resource
  authorization independent.  A certified upper bound below the target blocks
  production search; incomplete descent is scheduling information, while a
  separately bounded search may proceed and certified independent points give
  an unconditional lower bound.
- `run_mw17_jump_v2_zero_gain_rescue.sage`: executes the outcome-blind
  one-in-eight rescue assignment over the unchanged MW17-jump-v2 population.
  Assigned clean zeros search generic class ranks 44--344 in seven batches and
  switch unused slots to the existing adaptive policy after first escape,
  retaining a 344-chart total cap.
<!-- status-consumer: EC-K3-MW17-JUMP-V2-ZERO-GAIN-RESCUE 39ac93b60152bf88 -->
- `analyze_half_lattice_height_compression.sage`: reconstructs the exact
  reduced horizontal maps and presearch height data for 3,865 completed
  detailed charts.  It proves the chart midpoint/empty-ball identities,
  verifies constant quartic invariants, audits target-free predictors across
  eleven positive chart orders (including 394 compact earlier-control
  records), and on the sealed rank-28 ledger exactly replays all source
  visibility and prefix quotient gains from reduced-coordinate height.
  Target-relative results are explicitly posthoc and bounded misses retain no
  arithmetic conclusion.
<!-- status-consumer: EC-HALF-LATTICE-HEIGHT-COMPRESSION 3baeaf370aec751c -->
- `build_quotient_geometry_table.sage`: joins the five usable R17 controls,
  sixteen refreshed R17 ladder fibres, and nine A1/MW16 presentations to their
  displayed quotient lattices and exact blind-recovery traces.  It records
  every quotient Gram, regulator, successive minimum, intrinsic quotient
  energy, projection coefficient, optimal and actual half-lattice phase,
  reduced-coordinate distortion, and first recovery stage, then tests exact
  rational recovered subspaces against the deterministic successive-minimum
  flags.  Projection CVPs are checked at two Gram-rounding scales; numerical
  heights are not interval certificates.
- `run_curve385_height_compression_pilot.sage`: two-phase current-lattice
  builder and bounded pilot.  Its build phase samples all 29 parity bits by
  deterministic local ascent, reduces and calibrates 32 fresh charts, then
  freezes a diverse 16-chart order without calling `hyperellratpoints`.  The
  separate search phase completes all sixteen charts at height 100,000 with no
  finite points and no M29 group growth.  The default invocation verifies the
  protocol/result binding; the miss has no saturation or rank meaning.
<!-- status-consumer: EC-CURVE385-HEIGHT-COMPRESSION-BUILDER-PILOT c0d6f2d67018def4 -->
- `run_curve385_iterated_half_lattice_search.sage`: starts from the three
  directions blindly recovered by curve 385's generic-deepest 43 charts,
  saturates inside the discovered group, and searches the 301 height-prioritized
  lifts with a nonzero new quotient word. It reaches blind rank 29 and stops at
  the predeclared next-round limit rather than claiming stability.
- `verify_curve385_iterated_half_lattice_search.sage`: hashes the frozen blind
  ledger before loading the public fixture, then proves by mutually inverse
  determinant-one integral coordinate matrices that the blind and displayed
  public rank-29 subgroups are equal.
- `analyze_curve385_quotient_weight_profile.py`: uses only the frozen blind
  ledger's exact integral coordinates to certify quotient-rank gains `7,9,9`
  at weights `1,<=2,<=3`, and audits all 28 bases of the three-bit quotient.
- `build_curve385_sparse_quotient_rank32_protocol.py` and
  `curve385_sparse_quotient_policy.py`: freeze the staged 12-bit search,
  deterministic alternate quotient bases, restart-on-growth rule, exact
  acceptance gates, and fail-closed limits before any new search outcome.
- `run_curve385_sparse_quotient_rank32_search.sage`: checkpointed rank-32
  v1 runner starting from the blind `M29`. It is retained byte-for-byte for
  the completed primary campaign and its frozen combined four-state budget.
- `curve385_sparse_restart_policy.py` and
  `build_curve385_sparse_restart_budget.py`: classify exact group changes as
  rank-changing or saturation-only and freeze independent limits of three and
  four.  Their regression proves that two saturation-only changes do not
  prevent the three unit rank gains from 29 to 32.
- `run_curve385_sparse_quotient_rank32_search_v2.sage`: source-pinned future
  runner with the independent counters, strict v2 checkpoint validation, and
  target-rank precedence. It searches complete sparse stages, deduplicates
  exact base-point charts, and recomputes from weight one after every exact
  group enlargement.
<!-- status-consumer: EC-K3-R17-CURVE385-INDEPENDENT-RESTART-BUDGETS 39cfce110e3e494f -->
- `build_elkies_2026_rank28_relative_descent_magma.py`: replays the certified
  generic Kummer image and emits an unconditional basis-level Selmer job whose
  rejection gate precedes all residual-cover construction.
- `parse_elkies_2026_rank28_relative_descent.py`: accepts only a complete log
  from the exact generated Magma source and emits the common fail-closed gate.
- `run_elkies_2026_relative_2selmer_open.py`: open-source Sage/PARI replacement
  for the generalized R17 suite. Its blind worker uses
  `ellrankinit`/`bnfcertify`/`ell2cover`, stores the complete binary-quartic
  basis and maps bounded-search points back to the elliptic curve. The parent
  process introduces generic and held-out control points only afterward and
  labels recovered classes with exact finite-reduction coordinates. Resource
  stops remain incomplete artifacts.
- `run_elkies_2026_relative_2selmer_checkpointed.py`: resumable open descent
  which persists a fully certified cubic BNF, transports it through an exact
  `polredbest` field isomorphism, applies Simon's norm/sign/local conditions,
  embeds MW17 by exact squareclass tests, and builds/searches explicit
  intersections of quadrics for quotient classes before loading held-out
  exceptional points. Its v3 Selmer record also retains every local allowed
  subspace and the exact local-condition matrix rank after deleting each
  place. Enumeration limits and resource stops are recorded fail-closed.
  Its supervisor loads the authoritative control object, so this staged worker
  isolation is not the strict prospective record replay.
- `quotient_rank_escape_detector_v2.py`: exact backend-independent `F_2`
  measurement layer for a completed all-place descent.  It canonicalizes the
  global condition row space, quotients by the actual MW17 image, computes
  summed/independent local codimensions and every leave-one-place-out residual
  dimension, then loads the twelve held-out record directions.  It refuses
  incomplete place sets and makes no pairing claim on a coordinate complement.
  The present record certificate is Outcome D because its global squareclass
  domains have not completed.
<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-RANK-ESCAPE-DETECTOR-V2 eda7a0053b31b7c9 -->
- `build_r17_mw17_only_selmer_replay.py`: builds the strict prospective
  controls for curves 356 and 385.  Each generated Magma source contains only
  the minimal curve and exactly seventeen specialized generic points, requests
  the full unconditional 2-Selmer group, quotients only by MW17, and exits at
  `blind_freeze`.  A source audit forbids held-out coordinate rows, labels,
  MW29 tokens, half-ideal shortcuts, external reads, and point/cover searches.
- `run_r17_mw17_only_selmer_replay.py`: resource-bounded supervisor and strict
  transcript parser for those two sources.  It rejects a transcript whose
  final protocol record is not `blind_freeze`; both controls must certify
  MW17 rank 17 and residual dimension at least 12 before the Selmer candidate
  gate is called operational.  The committed ledger has zero completed runs.
- `build_mw29_relative_2selmer_matrix.py`: backend-independent proof gate for
  post-discovery closure on the record fibres. It row-reduces the certified
  29-dimensional Kummer image
  first, forms every local equation only on its complement, emits the actual
  residual kernel basis, greedy rank-gain order, pairwise local intersections,
  every leave-one-place-out rank, and an exact minimum annihilating place cut
  within a declared search budget. A zero matrix kernel is promoted only when
  the supplied global envelope and local maps are certified. Certified Selmer
  parity is applied to upper bounds, so an even residual bound of at most one
  closes at zero.
- `audit_mw29_relative_selmer_witness_bound.py`: non-enumerative upper-bound
  gate. From `dim V <= D` and exact global witnesses whose combined norm/local
  obstruction syndromes have rank `r`, it certifies
  `dim(Sel_2/im(MW29)) <= D-r-29`. This lets an F2-only class/ray-class bound
  and partial relation collection prune monotonically without pretending that
  the explicit witnesses span the anonymous class-group remainder. It likewise
  supports certified residual parity and rejects auxiliary fingerprints as
  condition blocks.
- `run_mw29_relative_2selmer_from_bnf.sage`: quotient-native certified-BNF
  post-discovery backend for curves 356 and 385. It constructs only the global norm envelope,
  embeds and exactly verifies all 29 point squareclasses, quotients them before
  requesting a local image, checkpoints after every place, and stops as soon
  as a certified local subset annihilates the residual envelope. It constructs
  neither the full Selmer basis nor any cover before this gate. The same JSON
  ambient manifest can instead be supplied by an F2-only class-relation or
  ray-class upper-bound backend.  Because all 29 point classes are quotient
  inputs, this backend cannot calibrate the prospective MW17-only gate.
- `run_elkies_2026_pari219_bnf_benchmark.py`: owns and benchmarks the six-
  parameter threaded BNF collector introduced on PARI's 2.19 development
  branch. It retains a binary checkpoint only after `bnfcertify`; timeouts
  retain factor-base and relation-deficit telemetry but no class-group or
  Selmer claim.
- `run_elkies_2026_record_pari219_bnf.py`: record-pair specialization with
  exact discriminant-factor hints. Its `class-quotient-upper` mode starts with
  `bnfinit(...,0)` and certifies only that the true class group is a quotient
  of the computed group via `bnfcertify(...,1)`. Thus a completed mod-two
  dimension is an unconditional upper bound without fundamental-unit
  certification; its checkpoint metadata explicitly says that it is not a
  full-BNF Selmer input.
  Resource stops preserve compact per-strategy relation-search telemetry and
  make no class-group, Selmer, or rank claim.
- `run_elkies_2026_pari219_selmer_from_bnf.py`: reloads such a certified BNF
  in the same GP build, applies the shared Simon norm/sign/local-condition
  implementation, and records algebraic Selmer representatives plus the
  exact local-condition matrix rank after deleting each place. This avoids
  importing a PARI 2.19 binary checkpoint into Sage's older libpari.
- `audit_elkies_2026_known_kummer_quotients.py`: fast class-group-free exact
  audit of the *supplied* Kummer subgroup. It evaluates `4*x(P)-zeta` in
  squarefree residue factors of an integral two-division polynomial, certifies
  MW17 plus the known exceptional directions on all five controls, and checks
  MW17 on the ten high-Nagao inputs. It neither computes a Selmer upper bound
  nor searches for unknown directions.
- `build_elkies_2026_known_quotient_covers.py`: enumerates all 3,851 nonzero
  classes in the five certified known exceptional quotient subgroups and
  constructs an exact intersection of quadrics with a verified rational point
  for every class. The compact manifest hashes the full local ledger; this is
  a realized-class control corpus, not a full Selmer quotient.
- `run_elkies_2026_s_class_hecke_monitor.jl`: checkpoints Hecke's exact
  principal relation rows modulo S columns while its open-source class-group
  collector runs. It can add every S-prime ideal to the factor base and method
  4 directly multiplies target ideals by S-ideals before short-element
  enumeration. Closing its bounded factor-base quotient would still require a
  factor-base generation proof, units, and local Selmer conditions.
- `r17_kummer_quotient_search.py` and
  `run_r17_kummer_quotient_sclass_collector.sage`: dependency-free policy
  helpers plus the exact Sage arithmetic driver.  Candidate lattices cycle
  through single, paired, and sparse products of all certified Kummer
  half-ideals, preferentially sampling the exceptional block.  The driver
  rotates through nonzero products of unresolved columns of configurable
  width and records every exact row's rank gain both modulo generic MW17 and
  modulo the full known subgroup. Before attempting any norm factorization it
  also caches reduced-ideal HNFs: an exact collision, together with the two
  retained principal multipliers, is immediately a quotient relation. Its
  optional `idealredmodpower2` engine reduces modulo ideal squares and serves
  as a diagnostic control; current pilots show that neither it nor ordinary
  multi-target reduction replaces a batch special-q sieve.
- `run_r17_kummer_quotient_sclass_suite.py`: gives curves 351, 356, 376, 377,
  and 385 identical bounded budgets in both quotient objectives, then records
  relation structure and descriptive displayed-MW-gain correlations in a
  local JSON summary.  Neither script's materialized quotient dimension is a
  global upper bound without a separate factor-base generation proof.

## Shared arithmetic

The principal reusable modules are `fermigier_mestre.py`, `mestre_root_tuples.py`,
`nagao_1994.py`, `nagao_linear_sections.py`, `multiple_root_lifting.py`,
`crt_lattice.py`, `finite_quotient_escape.py`,
`mod2_reduction_independence.py`, `mod_l_reduction_independence.py`, and
`pari_bridge.py`. Prefer extending one of these over creating another copy of
the same arithmetic.

## Current research code

- `prepare_r17_refresh_jump_ladder_inputs.sage`,
  `run_r17_refresh_jump_ladder_blind_v2.sage`, and
  `analyze_r17_refresh_jump_ladder_v2.py` implement the redacted sixteen-fibre
  `+3` through `+12` detector validation.  The runner cannot read the public
  complements or jump labels; the analyzer replays the frozen exact Kendall
  and upper-tail tests after the blind hash is sealed.  The v1 cross-class
  deepest-count failure is preserved rather than rewritten.
<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER a2d7034fb8977c18 -->
- Files containing `bnf_free`, `residual_selmer`, or `curve273` implement the
  unfinished residual 2-Selmer chain. Intermediate success is not a rank
  theorem.
- For Elkies compact-`t` candidates, residual dimension below 15 rejects rank
  32. Only a completed unconditional global/local descent of dimension at
  least 15 may unlock a same-curve two-cover or expensive point search.
- The retained `fermigier_rank20`, `mixed_small_prime`, and
  `six_root_low_conductor` drivers support improved-conductor and residual
  descent work now that the original rank-21/cutoff branch is closed.
- `newfamily/` has its own [workflow index](newfamily/README.md).
- The exact exceptional-transport and Mestre two-section checkers remain
  active because they have entries in `MATH_STATUS.json`; their old surrounding
  search campaigns are archived.

## Archive boundary

Superseded versions, completed bounded searches, negative scans, and their
tests/artifacts are indexed under
[`archive/elliptic-curves/`](../../archive/elliptic-curves/). Do not move an
archived bounded result back into the active tree merely because it is
interesting; promote only a compact reproducible result with the correct
evidence label and a canonical note.

<!-- status-consumer: EC-K3-ICARM-MW16-POINTED-SIEVE cb83c1afae1d0141 -->

- `run_fixed_field_comparison.py` runs and replays the
  [frozen six-curve comparison](../notes/FIXED_FIELD_COMPARISON_2026-09-05.md):
  665 CT entries, the zero control, and five radical-only bounded point searches.
  All five new deformations leave one unresolved nonzero inherited class.

- `run_fixed_cubic_cassels_tate.sage` computes the restricted `u=-1`
  Cassels--Tate matrix and searches its three nonzero radical classes.
  `verify_fixed_cubic_cassels_tate.sage` independently checks exact cover
  maps, cubic identities, local Hilbert symbols and support: pairing rank
  16, radical dimension 2, with point realization still unknown there.
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-CASSELS-TATE df45391a84f0e3c9 -->

- `run_fixed_field_radical_covers.py` replays six globally minimal
  quadric-intersection models for the three radical masks and their
  Q-translations. `--prepare` writes standalone lattice-search and
  four-descent jobs. Exact maps are certified; no class is decided.
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-MINIMAL-MODELS 90216b8c456edd20 -->

- `audit_fixed_field_radical_search_geometry.py` certifies that all six
  nominal height boxes were empty before square testing, using rational
  Sturm sequences and integer denominator bounds. No global class is decided.
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-SEARCH-GEOMETRY 678f7beb805a4530 -->

- `run_fixed_field_tangent_conics.py` certifies three cubic tangent conics
  and eight exact reductions. `--prepare --lattice` replays the bounded local
  solver with Sage class-group calls disabled. No genuine lift is constructed.
<!-- status-consumer: EC-FIXED-CUBIC-TANGENT-CONIC-GATE 26a49e30ff3128d3 -->
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-SOLVER-COMPARISON 6a178bc3a4ada43b -->
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-LONG-SEARCH 825fb4cd6ed84cb1 -->

- `solve_fixed_field_conic.py --verify` replays 320 local-reconstruction cells
  and 25,920 exact candidate misses, plus separate positive controls. Its
  relative-norm run timed out; no target conic point or higher cover resulted.

- `solve_fixed_field_conic.py --long-verify` audits the retained 45-minute
  norm timeout, 131,072-cell deep reconstruction transcript, and partial
  exhaustive `W=1` prime-37 transcript. No target point or higher cover resulted.

<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-VS-SHA-COMPARISON f37417a9fda3ee3f -->

<!-- status-consumer: EC-K3-ICARM-MW16-SENSITIVITY f88886c066d6cb45 -->

<!-- status-consumer: EC-FIXED-FIELD-COMPARISON 02c49a8120aeb7bd -->
