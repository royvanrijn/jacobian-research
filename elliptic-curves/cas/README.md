# Active CAS modules

This directory is intentionally narrower than the computational archive. It
contains exact status checkers and their shared arithmetic, plus code for the
current rank-32, low-conductor, residual-Selmer, and K3-construction gates.
Stable user-facing commands are listed in [`../scripts/`](../scripts/) and
[`../REPRODUCE.md`](../REPRODUCE.md).

## Start here

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
<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-RANK-ESCAPE-DETECTOR-V2 f07ee569c95bf3a1 -->
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
