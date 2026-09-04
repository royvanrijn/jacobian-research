# Retrospective discovery audit worklist — 2026-09-04

This is a maintenance and execution ledger, not a mathematical-status source.
`MATH_STATUS.json` remains the sole authority. A checked item means that the
named audit action was performed; it does not strengthen the scope of any
theorem or bounded experiment.

## Goal

Re-audit every recorded discovery and computational route against current
knowledge; recover missed failures and discarded valid candidates; repair
unsound pruning, stale assumptions, scripts, certificates, and navigation;
and leave every remaining unknown, bounded result, and expensive continuation
with an explicit fail-closed gate and reproducible certificate plan.

## Snapshot and notation

The latest strict working-registry pass has 1,163 entries: 1,004 proved, 102
partial, 17 open, 18 parked, 16 archived, and 6 falsified. The first
recent-history pass compared the 85 entries committed after `dada420` through
the then-current September 4 head and found 77 distinct checker programs.
Concurrent September 4 work has since added further status entries; those are
included in the live per-entry trace and must still be re-snapshotted at audit
close.

- `[x]` means completed and evidenced in this pass.
- `[ ]` means actionable work remains.
- **Blocked-current** means the item touches an actively changing file or
  certificate and must not be “fixed” by blessing an intermediate hash.
- **Long replay** means mathematically in scope but intentionally not hidden
  inside the cheap default check.

## P0 — authority, workspace, and fail-closed release gates

- [x] Read the root, K3, and elliptic-curve repository instructions and keep
  `MATH_STATUS.json` authoritative.
- [x] Inventory the dirty worktree before editing and preserve all unrelated
  September 4 work.
- [x] Compare the August 31 registry, the committed September 4 registry, and
  the live working registry rather than assuming one static head.
- [x] Confirm that canonical-source and checker paths resolve for the initial
  1,145-entry snapshot.
- [x] Validate the process ledger schema and generated-document markers.
- [x] Run repository hygiene, Markdown-link, JSON-schema, and whitespace
  checks on the first repaired snapshot.
- [x] Re-run the full registry path/dependency/consumer validation on the live
  1,163-entry working registry; all current paths, edges, and source hashes
  pass strictly.
- [x] Wait for the concurrent inherited-cover/complete-character-closure edit
  to settle, then re-audit rather than blessing its intermediate hash.  The
  final live checker hash now matches its registry row and no blocker remains.
- [x] Regenerate `STATUS.md` only through `scripts/render_status.py` after the
  checker hashes stabilized, and byte-check the rendered result.
- [x] Re-audit
  `EC-K3-R17-NORM12-11952-SINGLETON-PO0-TOP150` after its untracked checker
  settles. Its final checker SHA-256
  `084611e538ed995a6d9ce2b60c593bc891ac5977170dcb98a8ed5463e767fa8e`
  matches the registered hash, so the intermediate fail-closed blocker is
  cleared.
- [ ] Re-run `make check`, `git diff --check`, Markdown-link validation, process
  ledger validation, and repository hygiene as the final release gate.
- [ ] Recount all states and re-diff every status entry against the audit-start
  commit so concurrently committed discoveries cannot fall between snapshots.

## P0 — concrete repairs already found by the retrospective

- [x] Strengthen the singular-pencil shard merger so interval counts cannot
  conceal a duplicated parity class and an omitted good class.
- [x] Require exact equality between the merged parity-class set and the
  priority-table domain.
- [x] Validate shard direct-model hashes instead of checking only that a model
  path key is present.
- [x] Validate legacy-local or current-global sequential trace indices and
  emit global indices in all future direct singular-search shards.
- [x] Select the correct default priority table separately for `11952` and
  `103b2`, and reject unknown source labels unless a table is explicit.
- [x] Add byte-identical `--check` mode to the singular merger.
- [x] Replay all 63,917 `11952` and 63,925 `103b2` classes byte-identically
  under the stronger merger.
- [x] Mutate one copied shard adversarially and verify that the duplicate-class
  guard fails closed.
- [x] Pin the strengthened merger source hash in the singular genus-one status
  entry without changing the theorem scope.
- [x] Detect that the determinant-720 physical-corridor certificate pinned an
  obsolete hash of its route generator.
- [x] Replay the six-edge corridor, confirm that only provenance changed, and
  refresh the stored generator hash.
- [x] Follow that hash change into the integral rank-transfer glue census,
  refresh the dependent hash, and replay the census exactly.
- [x] Detect that the projective first-jet positive control silently discarded
  its intentionally degenerate mod-17 witness because a fixed control fibre is
  singular modulo 17.
- [x] Restore the degenerate algebraic control while retaining nonsingular,
  53-equation positive controls at 53 and 67.
- [x] Replay all 210 active mod-17 target charts with no solution and no
  timeout, and pin the repaired first-jet checker source hash.
- [x] Replay the active 2,550 mod-53 and 4,160 mod-67 artifacts under the
  repaired checker with exact-input-hash-validated unit caches: both retained
  zero solutions and zero timeouts. The mod-53 comparison used a disposable
  normalization of cache-policy metadata and did not rewrite the certificate.
- [x] Detect that the committed prescribed-root shard summarizer rejected
  duplicate ambient labels but did not require the documented complete set of
  13 ambients, so an omitted whole shard could still print `PASS`.
- [x] Pin the ordered NS0001--NS0048 universe, exact 13 ambient labels, 16 D5
  anchors, and current source hashes in that merger; verify the existing
  four-shard artifact byte-identically and adversarially confirm that omitting
  the `D16_E8`/`D24` shard now fails closed.
- [x] Detect and document that the proved corrected H3 q8 checker is not
  portable from a clean checkout: its q6 child-Jacobian, transported-zero, and
  E7-infinity JSON prerequisites are ignored and absent.  Correct the canonical
  note, historical frontier, registry scope, and failure ledger so the old
  one-line command is no longer presented as standalone.
- [x] Replace the finite-etale checker's tautological target-normalization
  assertion (`-0/2 == 0`) with an exact quotient-algebra reconstruction of all
  three determinant-one normalized output coordinates.  The complete existing
  symbolic regression passes, and `FEKC1` now pins the strengthened source.
- [x] Correct two active SIC narratives that still treated the full mixed
  bidegree-`(3,3)` moment--nullcone equality as open.  `SIC2B33RDS` already
  falsifies that equality with a full-rank semistable Rodrigues survivor; its
  separate all-order integration-by-parts cutoff makes that orbit SIC-safe.
  Preserve the actual open problem: whether every bidegree-`(3,3)` all-order
  pure-zero orbit is SIC-safe.
- [x] Repair the bidegree-`(3,3)` sparse-census validation so stored totals
  cannot hide a duplicated support, omitted support, or mislabeled outcome.
  Exact cheap audits now reconstruct all 560 size-three, 1,401 mixed size-four,
  3,864 mixed size-five, 7,588 mixed size-six, 11,200 mixed size-seven, and
  12,780 mixed size-eight keys and require each exactly once; they also pin the
  two size-six survivor indices/supports and the separately resolved size-eight
  parity timeout.
- [x] Replace the final size-nine checker's count-only predecessor complement
  with the explicit union of rectangles, fringes, crosses, both three-line
  classes, the full-line class, and the `(3,2,2,2)^2` class.  A cheap audit
  proves this union is exactly the other 4,370 mixed supports and that the
  final 7,050 supports complete all 11,420 supports in 2,924 symmetry orbits.
- [x] Harden the reusable sparse-support shard merger itself: every record in
  every half-open shard interval must now equal the reconstructed ordered
  mixed-support key at that global index, rather than merely matching an
  aggregate count and contiguous interval.
- [x] Correct the SIC2C4 and three-pair Image-Mathieu checkers' console
  wording: checking 99 instances is a bounded replay of formulas proved for
  all orders in their canonical notes, not itself an “all-order” computation.
  Preserve the independent written all-order proofs as the theorem sources.
- [x] Expand the underspecified `O1` registry scope.  It now records the
  Keller and local-Artin hypotheses, uniqueness of the source trivializer,
  the determinant-one clause, and the dimension-three finite-jet statement,
  while explicitly denying bounded-degree orbit triviality, reduced-base
  algebraization, or disappearance of global stable moduli.  Retain the
  rank-four translation as the exact nonalgebraization control.

## P0 — discovery and process ledgers

- [x] Extend `RESEARCH_TIMELINE.md` through the direct alternate-Q80 equation,
  noncyclic closure, complete public atlas, singular/genus-two boundary,
  quotient-class exclusions, and arithmetic marking obstructions.
- [x] Update `KNOWLEDGE_BASE.md` with target-free compilation, arithmetic
  marking, finite-quotient slicing, sharded coverage, formal-local boundaries,
  and mutable-source projection lessons.
- [x] Expand the K3 process ledger from 55 to 75 events and from 17 to 24
  mechanisms, then regenerate all process-atlas tables.
- [x] Correct the route/failure ledger’s obsolete “alternate equation open” and
  historical q24-frontier language.
- [x] Add the false-assumption records for target-fitted prediction, geometric
  versus arithmetic rank, formal versus rational sources, ranked prefixes,
  quotient-class overreach, and incomplete shard coverage.
- [x] Record the seven exact native quotient/visibility audits, 51 fitted
  norm-eight directions, and the 556-curve live-source projection while
  retaining 57 quotient rows as literal `UNKNOWN`.
- [x] Record the NS0031 split-Clifford/`X_0(37)` rational-marking obstruction
  and replace the obsolete NS0031 algebraization handoff by arithmetic-first
  reranking of the remaining lattices.
- [x] Record the 66-row rank-19 arithmetic-marking classifier in the timeline,
  knowledge base, and K3 process atlas: 1 existing control possible, 3
  excluded, 62 literal `UNKNOWN`, and no equation-work handoff. Distinguish the
  coarse even-Clifford/norm-one curve from the full stable
  discriminant-kernel marking curve.
- [x] Promote the exact incremental-CRT beam false negative into the registry,
  timeline, knowledge base, and failure ledger; document both active beam
  helpers as heuristic proposal mechanisms rather than completeness sieves.
- [x] Add the already-committed 2,134-row prescribed-root MW1 census to the
  mathematical registry with its exact 13-ambient/16-anchor/48-NS boundary and
  its non-arithmetic, non-isometry-classification caveats.
- [x] Add a short repository-root navigation link to this worklist and keep it
  out of the mathematical-status hierarchy.
- [ ] At close, verify that every new September 4 status entry is represented
  in the research timeline or is deliberately omitted as an implementation
  detail.

## P1 — per-entry registry audit protocol

Structural coverage completed in this pass:

- [x] Add `scripts/audit_discovery_registry.py` and execute it across all 1,163
  live rows, tracing 2,570 internal dependency edges, 690 external dependency
  edges, 484 generated-artifact locks, 32 local-artifact locks, inverse
  consumers, update edges, canonical notes, checkers, and current hashes.
- [x] Extend the per-entry trace with the literal title and scope, proof type,
  assurance flags, canonical-note digest, forbidden attack classes,
  supersession notes, and status IDs mentioned in the scope, including whether
  each mention has a declared registry edge.
- [x] Write and byte-check the full per-entry trace at
  `artifacts/local/discovery-retrospective-registry-audit-v1.json`; keep this
  dynamic ledger local so it cannot become a second status authority.
- [x] Add duplicate guards for dependencies, software locks, and replacements;
  remove the sole duplicated lock found in the live registry.
- [x] Link the two obsolete IV-star marking assumptions reciprocally to the
  proved corrected D13 equation result, rather than incorrectly using the
  still-partial marking row as their formal invalidator.
- [x] Promote already-stated failure boundaries to machine-readable forbidden
  attack classes for `OP-EC-NEXT`, `OP-SIC2-B33`, `OP-GVC3-MIN`,
  `OP-HC4-D5`, the parked NS0024 route, and the active different-NS programme.
- [x] Review and encode evidence-backed forbidden routes for eleven further
  open/parked rows: the bypassed binary-GVC promotion, closed GVC3 tagged
  search, binary-GVC polarization, global sparse census, both legacy LR faces,
  restricted minima, minimal-boundary suspension, resolved GMC(2), closed
  rank-ten matroid gate, and closed elliptic rank/conductor disjunction.
- [x] Add the missing `UCUT3` dependency to `OP-SUSP`; its scope already used
  that exact counterexample to forbid universal cubic gradedness.
- [x] Review the remaining 16 open/parked rows. Encode exact forbidden routes
  for `OP-AHQ`, `OP-CCDM`, `OP-CR`, `OP-GS`, `OP-KCOMP`, `OP-KDESC`,
  `OP-KMOD`, `OP-KMON`, `OP-LR-NE`, `OP-PLANE`, `OP-QO`, and `OP-RITT`.
  For `OP-ARITH`, `OP-FC`, `OP-PF`, and `OP-QP-PROV`, record explicit reviewed
  omission reasons: their current dependencies do not prove a route-specific
  no-go beyond their short open or absorbed scopes, so adding prohibitions
  would be boilerplate rather than evidence.
- [x] Teach the renderer and registry audit to distinguish an unreviewed open
  problem from a reviewed row with no evidence-backed forbidden attack. The
  live report now has zero unreviewed open/parked attack-boundary rows.
- [x] Add checker source-audit markers and shared-checker users to the local
  per-entry trace. The live registry uses 922 distinct checker paths, 74 of
  them shared; markers locate early loop exits, pruning vocabulary, bounds,
  timeouts, saturation, and denominator handling without treating their mere
  presence as a defect.
- [x] Inspect all 18 distinct partial-result checkers containing literal loop
  `break` exits (25 status rows). Their exits are monotone completion/unit
  detection, exact linear rejection, Weyl or hull termination, or explicitly
  recorded failure/timeout paths; none silently upgrades its partial scope.

Apply this checklist to every non-archived registry row, recording exceptions
by ID rather than silently normalizing them:

- [ ] Confirm the title and scope still match the strongest current canonical
  proof and do not omit a later counterexample or narrowing theorem.
- [ ] Confirm every dependency is logically necessary, correctly directed,
  and not merely a historically adjacent calculation.
- [ ] Confirm the canonical source is unique, active, and contains the theorem
  boundary rather than only an experiment narrative.
- [ ] Confirm every named checker is the program that proves the recorded
  payload, not a generator, downstream summarizer, or unrelated smoke test.
- [ ] Recompute each checker-source hash after all live edits settle.
- [ ] Confirm each software lock exists and distinguishes source programs,
  raw replay inputs, compact certificates, and optional local checkpoints.
- [ ] Validate document consumers and their status markers after every scope
  change.
- [ ] Check reciprocal `supersedes`, `replaced_by`, `closes_problems`,
  `narrows_problems`, and `invalidates_assumptions` edges.
- [ ] Ensure a proved bounded reproduction is described as a proved execution
  of a bounded protocol, never as the unbounded mathematical conclusion.
- [ ] Ensure `independent_replay`, `formal_verification`, and `external_review`
  are evidence flags, not inferred from proof type or a passing local checker.
- [ ] Preserve falsified and archived evidence needed to explain why a route
  or assumption must not be reused.
- [ ] For entries without checkers, identify whether the proof is genuinely
  prose-only or is missing a compact deterministic witness that should exist.
- [ ] For shared checkers, replay every parameter/source mode named by status;
  a default-mode pass is not coverage of alternate modes.
- [ ] For all partial entries, name the first exact missing gate and remove any
  prose that implies the partial calculation is “almost” a theorem.
- [ ] For open and parked problems, verify forbidden attack classes and
  witnesses reflect every later negative result.

## P1 — generated-certificate provenance graph

- [x] Scan 981 tracked generated JSON certificates recursively for recognizable
  repository-local SHA-256 input edges: 1,164 edges were checked.
- [x] Classify 24 owners with 60 historical/current-code hash differences
  rather than treating every immutable old-code hash as a stale theorem.
- [x] Repair the active determinant-720 corridor and glue-census chain found by
  exact checker replay.
- [x] Confirm that 22 mismatch owners are inactive historical/exploratory
  outputs whose old program hashes should remain evidence of the executed run.
- [x] Replay the generic-linear and specialized Brown-PRS third-`q12`
  certificates whose embedded worker hashes predated the current worker:
  exact mathematical payloads were unchanged apart from explicit unused-mode
  nulls, and their provenance now names the current worker.
- [ ] Refresh the downstream generic-quartic certificate after its local exact
  first-jet prerequisite. The prerequisite was rebuilt successfully, but the
  full multi-million-bit discriminant/division replay was stopped when the
  user clarified that this pass is cleanup-only. Do not resume it without
  explicit compute scope. **Long replay.**
- [x] Add `scripts/audit_generated_provenance.py`, a reusable read-only scan
  that reports matching, changed, and missing inputs separately for active
  co-locked dependencies and historical/unclassified pins. Its first live
  pass found 1,192 recognizable edges in 983 tracked JSON artifacts and
  isolated exactly the two already-known active third-`q12` worker drifts.
- [ ] Never bulk-refresh embedded hashes: regenerate through the owning checker
  and diff the mathematical payload first.
- [ ] For every mutable external source, preserve a minimal exact projection or
  raw snapshot sufficient for the claim; a URL plus old hash is not replayable
  once upstream bytes change.
- [x] Repair the historical `wgxli` lineage replay so it can validate its
  September 1 projection offline.  The default checker now joins two later
  committed, hash-pinned sufficient projections (the original ids 1--474
  equation inventory and all thirteen `wgxli` point records), cross-checks
  every claim-relevant field, and reproduces the original 27,432-byte artifact
  byte-for-byte at SHA-256 `f875f391...c1c03`; the mutable live-URL audit is an
  explicit optional mode.
- [ ] Add output hashes to canonical notes when a certificate is primary and
  currently lacks a documented digest.
- [ ] Rebuild the full `EC-K3-H3-Q6` intermediate chain and preserve a compact
  sufficient projection of the three inputs needed by
  `derive_h92_q6_child_q8_corrected2cover_qq.sage`; this is an explicit replay-
  packaging calculation and remains outside the cleanup-only execution scope.

## P1 — recent K3 and elliptic-curve discoveries

- [x] Replay the clean Python checkers for ICARM bounded rebasing, rootless-J1
  bound, deep-cover quotients, visibility complexity, rank-28 character glue,
  NS0024 relative-`U`, both inherited-product modes, and the V4 rank screen.
- [x] Replay the 43-chart/474-curve exact lineage atlas and confirm six
  `PGL2` classes, 69 hits, and curve 12 quotient `Z^12`.
- [ ] Replay every checker behind the seven-fibre native quotient/visibility
  audit after the current `08f72` extensions settle. **Blocked-current.**
- [ ] Compile saturated chart transports for the 54 recognized fibres whose
  quotient dimensions remain unknown; preserve nulls until each exact replay.
  The highest-rank tranche is complete: curves 11, 391, and 423 each have an
  exact displayed `Z^11` quotient transport.
- [x] Group predictor training/test splits by the six family classes.  The new
  1,536-row ordinary cohort locks two complete rational-`PGL2` families until
  a predictor trained on the other four complete families is hash-frozen;
  no fibre-level random split is permitted.
- [ ] Record actual search exposure denominators before comparing public hit
  counts across families; do not infer effort from successes.  The historical
  denominator remains unknown; the prospective protocol instead commits a
  new, fully enumerated 1,536-row scheduled denominator.
- [x] Replay the rootless J2 frame classifier modes that complete within the
  cheap window and retain the full Niemeier-first enumeration as a long replay.
- [x] Replay the direct noncyclic `4A1/MW13` equation and reverse target-free
  hop.
- [x] Replay arithmetic rank transfer, exact E6/E6+A1 controls, NS0024 route
  checks, NS0031 physical corridor, and the core finite transfer certificates
  that complete inside the declared limit.
- [ ] Complete long replays for core generation, prime-local bridge mutation,
  103b2 MW glue, determinant-78 frame classification, rootless low-degree
  census, and full J2 Niemeier enumeration.
- [ ] Replay the NS0031 rational-marking checker and independently inspect the
  two external theorem inputs: Vélu’s `X_0(37)(Q)` classification and the
  marked-K3 period/Clifford implication.
- [x] Run the same rational-marking/moduli screen on determinant 720 before any
  coefficient search. The stable curve is `X_0(60)` and has no rational
  noncuspidal points; the known rational `3A5` point has determinant 20 and is
  only a negative source control.
- [x] Reverse the global foundry to rank all 827 transcendental rows by their
  full marked arithmetic gate before NS complement, rootlessness, equation
  incidence, or planner cost. The current new NS/equation handoff is empty.
- [ ] Keep the NS0024 and NS0031 local/geometric artifacts as controls while
  preventing either from re-entering the rational-source queue.
- [x] Replay the rank-28 pilot, mixed-trace and simultaneous-splitting bounded
  searches, 103b2 hard-fibre products, alternate laboratory, control-j
  preimages, rationalized-D6 chart, and V4 shortlist.
- [ ] Complete the long 103b2 MW-lattice, 63,917-class product-inversion, and
  high-throughput splitting replays with checkpoints.
- [x] Audit every `PASS_EXHAUSTIVE_BOUNDED` status string. The only two live
  instances are the rank-28 simultaneous- and mixed-trace splitting artifacts;
  both record the primitive height-10,000 box, exact visit/test counts, the
  selected cover or trace domain, and a proof boundary denying any global
  rational-point or rank conclusion.
- [ ] Revisit the two V4 base-rank timeouts separately; the stored `62 complete,
  2 timeout` screen must never be summarized as 64 completed ranks.
- [ ] Compute nonzero product quotient classes only after the exact full
  involution/2-Selmer lattice is available; zero-class exclusion is not an
  existence or rank theorem.
- [ ] Keep the genus-two normalization result bounded until infinity charts,
  simultaneous bad reduction, and parameters beyond the CRT box are covered
  or a theorem removes them.
- [ ] Add independent replay implementations for the highest-risk complete
  character and singular-normalization exclusions.

## P1 — historical K3 failure modes to re-test everywhere

- [x] Search all active K3 scripts for unqualified ADE/MW equality used as an
  object identity; require marked `F`, `O`, components, chamber, and integral
  transports.  The 404 live-checker hits use ADE/MW tuples as classification
  fields or as one gate inside a marked construction, not as the sole object
  key.  Route composition checks primitive `U`, nef/horizontal-wall data,
  mutually inverse determinant-one NS transports, literal frame Grams, and a
  terminal integral isometry.  The NS0024 completion selects within a fixed
  core/bridge/glue construction and exports each full child Gram; it explicitly
  denies a marked elliptic-neighbour identity.  The noncyclic R17 control also
  records that three `4A1/MW13` frames with different norm-four counts are
  distinct.  No active script equates two fibrations from ADE/MW labels alone.
- [x] Search for pseudo-zero or chamber-zero reuse and require equation-side
  effectivity before any route score or continuation.  The 111 live-checker
  hits and 14 directly affected status rows retain the distinction.  The
  q4/orbit230 audit compares the chamber pseudo-zero with the exact effective
  section and withdraws the `4199` cost; the q6/orbit1307 route withdraws its
  component-10 continuation; active q4/q8/q10 routes consume explicit
  equation-effective-zero certificates or remain partial where effectivity is
  missing.  The noncyclic R17 lattice control proves its physical zero by the
  primitive nef fibre, `(-2)` square, fibre intersection one, and component
  intersections, while still denying an equation.  No active route score or
  equation continuation reuses a rejected chamber pseudo-zero.
- [x] Search binary-quartic and cover maps for missing degree-two accounting;
  regression-test heights and collision degrees.  The 195 live-checker hits
  preserve the distinction between old-fibre degree, cover degree, geometric
  branch degree, and specialization degree.  The shared extension checker
  cancels numerator/denominator factors before counting odd branch support,
  accounts for infinity, rejects split covers and branch degree four, requires
  exact positive-definite anti-invariant height data, and adds rank only from
  the rational deck eigenspaces.  The rank-28 glue checker explicitly applies
  `height_after_pullback = 2 * height_before_pullback`; the norm-twelve glue
  certificate independently records `2 * 12 = 24`.  Its adversarial self-test
  passes under the repository virtual environment.  No missing factor of two
  or collision-degree promotion required repair.
- [x] Search rational normalization and CRT code for incomplete denominator
  clearing, especially omitted `Dx`-type factors.  The 266 focused hits across
  live K3 checkers retain common denominators before coefficient rows, require
  polynomial quotients by exact remainder/denominator tests, and saturate
  denominator-created affine components only when the chart denominator is
  declared.  The corrected q8 compiler literally enforces
  `R*h*Dy == Ny*Dx (mod Nx)`, while the genus-two reconstruction rejects
  nonpolynomial residual products/lifts before characteristic-zero
  substitution.  The generic neighbour compiler's spurious-component
  adversarial regression passes.  No active formula still uses the withdrawn
  `Ny`-without-`Dx` congruence.  The separate missing-q6-intermediate issue is
  now recorded as replay packaging, not mistaken for a denominator failure.
- [x] Search resolved-RR code for row-per-blow-up-centre condition counts and
  replace them with saturated connected-component modules.  The 102 focused
  live-checker hits preserve the repaired model: q24 distinguishes 12 centres
  from 13 geometric components, transports the selected connected effective
  cluster, recomputes the exact local-surface quotient after each chart
  transform, and checks every dimension loss against its modular replay.  The
  generic compiler represents local constraints as named quotient blocks with
  declared trivializations; its regression rejects a non-kernel pencil and
  requires explicit saturation for a denominator-created component.  A lone
  row in the fixed reverse `I2` calculation is a genuine single connected
  nonidentity component with a unit chord denominator, not one row per centre.
  No active RR compiler retains the withdrawn centre-count heuristic.
- [x] Search specializations for inherited generic ranks, roots, pole orders,
  or MW coordinates; force a new typed specialization state and recomputation.
  Across 171 live K3 checker paths this produced 274 broad textual hits, 40
  rank/root/pole/coordinate-focused hits, and 59 affected non-archived K3
  status rows.  Exact specialization claims reconstruct the specialized
  curve and points, verify the specialized group law, and certify independence
  by finite-quotient or exact rank checks; inherited generic sections are
  recorded as a specialized subgroup, not as the full specialized group.
  The point-factory control keeps specialized `MW13` coordinates explicitly
  `NOT_YET_COMPUTED`, fibre-root calculations re-factor the specialized
  discriminant, and all bounded or numerical specialization probes deny rank
  upper bounds.  No untyped generic-to-specialized promotion required repair.
- [x] Search active checker paths and consumers for first-hit or beam-pruned
  output used as a completeness claim. The retrospective bridge predictor is
  already non-prospective, retains only 4 of 5 rootless controls, and has zero
  ranking power on a held-out 277-candidate shell. Separately, the elliptic CRT
  regression proves that a width-one height beam misses a height-53 completion
  and keeps height 1409; the analogous q12/orbit5867 projective constructor is
  explicitly bounded. Both active beam helpers now forbid mathematical
  exclusion from beam survival, and the open problems carry explicit no-go
  metadata.
- [x] Search ranked-prefix consumers for mathematical rejection of candidates
  outside the prefix. The alternate-Q80 1,024-cover laboratory, genus-one
  trace sample, frozen Nagao cohort, Fermigier control tables, A11 cost display,
  genus-two/three selected traces, and Golay size-five prefix all retain their
  literal population or cohort boundaries. No active consumer upgrades an
  outside-prefix omission to a mathematical rejection; the open-problem and
  failure ledgers explicitly prohibit that upgrade.
- [x] Audit the committed shard aggregators for the singular norm-eight search,
  prescribed-root foundry census, 24A1 completion/canonicalization, NS0007
  fixed-case census, degree-three spectrum, and bounded SIC scouts. The first
  two needed and now have adversarial omission/duplication guards; the others
  already pin contiguous ranges or literal expected shard/key manifests and
  keep their bounded domains explicit.
- [ ] Re-audit the alternate-bisection chunk merger after its concurrent
  `1183a`/`098fc` expansion settles; do not bless its intermediate hashes.
  **Blocked-current.**
- [x] Review the 56 non-archived registry rows whose title or scope invokes
  modular, finite-field, Hensel, p-adic, or formal-lift evidence. Local-only
  rows explicitly deny algebraization and rationality; rows promoted to `QQ`
  include exact rational reconstruction plus literal characteristic-zero
  substitution. The D5 low-slice eliminants correctly reclassify the two
  tempting p-adic lifts as nonrational algebraic points, and the NS0031 formal
  branch remains blocked by the rational-marking obstruction.
- [x] Review all 64 non-archived registry rows mentioning geometric Picard or
  Mordell--Weil data. Arithmetic-rank rows either display exact `QQ` sections
  and an independence/saturation certificate or compute the Galois action;
  geometric-only neighbour and foundry rows explicitly deny arithmetic
  descent. NS0024 and NS0031 remain preserved as geometric controls but are
  excluded as fully rational rank-19 markings, and the new 2,134-row census is
  explicitly non-arithmetic.
- [x] Search quotient-sliced conclusions for accidental promotion of one class
  to all quotient classes.  The 11 direct `quotient class`/`Tate class`/`zero
  class` registry matches and the broader 49-row live class/coset/squareclass
  review preserve their declared slices.  In particular,
  `EC-K3-R17-NORM12-11952-PRODUCT-ZERO-TATE-CLASS-EXCLUSION` excludes only the
  zero Tate class under its height-eight direct-polynomial and local-component
  gates; its canonical note and the deep-trace follow-up both leave the full
  quotient, every nonzero class, section existence, and any Mordell--Weil rank
  consequence `UNKNOWN`.  No ledger repair was needed.

## P2 — elliptic curves over `Q`

- [x] Re-audit every claimed rank lower bound by literal point substitution,
  minimal-model transport, torsion handling, and independent finite-quotient
  rank.  The 16 active rational-specialization controls use exact point and
  model transport plus finite-quotient independence and torsion exclusion;
  generic-family ranks use separate Shioda--Tate/Picard arguments, never
  numerical heights.
  **Completed record-curve provenance subset:** the `ECR30` and `ECR31`
  status rows now register their pinned JSON/gzip certificates and all direct
  proof helpers, including the shared exact mod-2 quotient implementation and
  the separate curve-273 Sage replay. Static review confirms that a bounded
  failure to reach full quotient rank would be reported only as an
  inconclusive nonsaturation possibility, never dependence. A stdlib-only
  auditor checks bytes, source provenance, stored matrix shapes, and claim
  flags without curve arithmetic or rank recomputation. The numerical
  regulators, public BSD/GRH statement, exact-rank conclusion, and unproved K3
  provenance remain excluded. No PARI, Sage, finite-group, or matrix-rank
  calculation was run.
- [x] Reconcile the fixed-root and generated-family ICARM construction ledgers.
  The v2 fingerprint row now includes curve 302 and narrows its formerly broad
  “repository-wide” diagnostic to what the historical script actually scanned:
  uncompressed generated-results JSON only. The separately committed
  2,334-family discovery is now registered with its 2,329-family census,
  bounded Fermigier-generator additions, deduplication, exact-survivor counts,
  curve-282 isomorphisms, and bounded negative scope. A stdlib-only auditor
  checks hashes and stored scope flags without generating families, running a
  modular sieve, or factoring. Static source review found the modular exits
  one-sided and fail-closed: degree drops are retained and only a good-prime
  no-root certificate rejects. No recognition calculation was run.
- [x] Re-audit every exact-rank phrase for a matching unconditional upper bound.
  Curve 302 remains rank at least 31, not exact rank 31; `EC-MD235-GENERIC-
  RELATION` concerns only the displayed subgroup; and the rank-19 Mestre upper
  bounds remain explicitly conditional on GRH and BSD.  `EC-NF-R14` is the
  sole active `E/Q` exact-rank promotion and requires both its exact eclib lower
  bound and PARI's unconditional `ellrank` interval `[14,14]`.
- [x] Re-audit conductor claims through global minimization and local Tate data,
  never a discriminant radical.  Every active exact-conductor claim reaches a
  PARI `ellglobalred` computation on the globally minimal model or explicitly
  reconstructs the conductor from `elllocalred` exponents; discriminant-degree
  comparisons in the Mestre family remain labelled conductor geometry only.
- [x] Check all Selmer/descent outputs for the distinction between Selmer
  dimension, Mordell--Weil quotient, and possible Tate--Shafarevich classes.
  Corrected the PARI adapter's misleading `pari_sha_two_dimension` name: its
  third `ellrank` field is the even Cassels-pairing quotient rank
  `dim(Sha[2]/2Sha[4])`, and a tested helper now reconstructs the total Selmer
  dimension from `C=r2+T+s` before applying the residual gate.
- [x] Keep timeout, `norm-one`, and incomplete local-class outputs fail-closed
  in the rank-32 gate.  Forty exact policy tests cover incomplete backends,
  provisional S-class data, norm-one local inconclusives, model/parameter
  binding, residual thresholds, and all five expensive-search entry points.
- [x] Require the family-relative residual dimension `32-r` on the same
  minimal curve before any expensive rank-32 cover or point campaign, where
  `r` is the certified generic rank.  The retained MW17 entry points enforce
  threshold 15; a reconstructed A1/MW16 campaign must enforce threshold 16.
  Every search gate must be bound to the exact parameter and global minimal
  model.
- [ ] Treat curve-398 A1/MW16 reconstruction and curve-302 parent-family
  reconstruction as a first-class rank-32 path parallel to R17/MW17.  For 398,
  recover the exact fibration, parameter, sixteen-section map, and rank-14
  displayed quotient before neighbourhood search.  For 302, work from the
  complete 31-point configuration without assuming a generic rank or
  `17+14` split.  Keep curve 273 and low-conductor near misses behind completed
  residual descent or a sharper unconditional upper-bound route.
- [x] Preserve exact native K3 fibres as calibration controls without assuming
  rational parameter transport between different fibrations.  The compact
  published-R17 rank-25--28 fibres remain bound to their literal parameters,
  models, and specialized sections; the native alternate-Q80 curve-12 control
  is admitted only through an exact rational `j`-preimage, twist test, and
  model isomorphism.  The six distinct norm-twelve `j`-map classes retain
  separate base coordinates, with `PGL2(Q)` transport used only inside a
  proved class.
- [x] Audit all public-database checkers for mutable-source preservation,
  deterministic point ordering, and snapshot-specific claims.  The active
  ICARM consumers either embed/hash-pin their model and ordered point list or
  now default to two sufficient committed projections.  Offline regressions
  replayed the 69-fibre/1,545-point snapshot, all 2,844 original curve/class
  decisions, the five exact section/quotient certificates, the 65,536-sign
  rebasing control, and the earlier curve-351/356 fingerprint.  Network access
  is now an explicit drift/refresh mode rather than a theorem dependency.

## P2 — foundational Keller, cancellation, and stable boundary programmes

- [x] Re-read every core theorem from `F1` through the finite-etale and
  multiplicity layers against the current falsification graph.  The scoped
  chain `F1`, `W1`, `C1`, `RQG1`, `URKC1`/`URK1`, `S1`, `B1`, `P1`, `T1`,
  `M1`, `RQG2`--`RQG5`, `CAF1`/`LPCF1`, `FEKF1`/`FEKC1`, and
  `UQFM1`--`UKFM2` keeps explicit construction families, presentation
  descent, canonical-boundary invariance, complete fibres, and stable-class
  conclusions separate.  The falsified `OP-UG3` does not invalidate any
  multiplicity row: `UKFM2` limits the cubic collapse to three controlled
  mechanisms, while `UCFM1` proves multiplicity by fiber-invisible phantom
  boundary components.  Repair the one vacuous `FEKC1` normalization test and
  the stale `OP-SUSP` formal-tail narrative found during this review.
- [x] Recheck marked-root reconstruction for omitted boundary points and
  accidental promotion of an affine inverse to a global one.  The foundational
  model uses both projective-root charts and sends the simple infinity root on
  `c=0` to the finite `x=0` source divisor; its inclusion and exceptional-fibre
  checkers explicitly split that case.  Weighted, cancellation, and
  quadratic-gauge full-fibre claims require squarefreeness plus their literal
  reconstruction-open condition.  The finite-etale quotient reconstruction
  checks both compositions and now also checks the determinant-one target
  scaling scheme-theoretically.
- [x] Recheck valuation/cancellation scripts for polynomiality after all
  denominators and for full source--core--target determinant balance.  `BCI1`
  keeps polynomiality as a separate required datum and checks the complete
  source--core--target ledger for all three established families; the family
  notes supply weighted admissibility, the finite cancellation jet congruence,
  or the quadratic-gauge coefficient-weight identity before using the
  cancelled Jacobian.  No rational zero--pole identity is promoted by itself.
- [x] Recheck stable-multiplicity claims against the universal-cubic gradedness
  counterexample and keep only the boundary-minimal theorem where applicable.
  The multiplicity rows are scoped to explicit weighted, quadratic-gauge,
  power-shift, or fiber-invisible cubic families; none infers universal cubic
  gradedness.  `UCFM1` instead uses the growing phantom-boundary count as its
  stable separator.  Remove the stale active `OP-UG3` roadmap language and
  retain `UCUT3` as the exact falsification.
- [x] Re-audit the global low-degree sparse Keller census without launching a
  new census or Gröbner calculation.  The singleton-bucket closure is an exact
  branch-complete rule on the declared coefficient torus; valuation, sign,
  and modular stages are organizers rather than exclusion sieves; and all 913
  orbit representatives reach the stored characteristic-zero stage.  Add a
  solver-free audit mode that checks the eight manifest-pinned artifacts,
  reconstructs each residual two-element orbit from its representative,
  validates stored determinant triples, and requires identical ordered IDs in
  every intermediate and exact ledger.  The committed data pass, while
  support-seven attainment and the cardinality-unbounded census remain open.
- [x] Separate the plane sparse-support theorem from its two large bounded
  regressions in the maintenance interface.  Static review confirms that the
  support-at-most-six result is sourced by the written chain argument, exact
  unbounded Z3 exponent formulas, explicit singleton Rabinowitsch identities,
  and the four exceptional saturated shear charts; the 14,653,584-row `2+2`
  and 5,290,000-row `3+3` enumerations are checks only.  Add a subsecond audit
  mode that hash-pins the committed JSON and Singular program and validates
  all claim boundaries, masks, digests, survivor labels, and shear formulas
  without starting either enumeration or any solver.
- [x] Recheck the adjacent affine-support/Newton bridge specifically against
  accidental reuse of the support-six theorem on the five-term F2 terminal
  block.  Its two component supports occupy characters `{1,4}` and `{0,1,3}`
  modulo five, all six character-pair brackets are present in the committed
  ledger, and the non-target sectors cancel rather than descend.  Add a
  solver-free artifact audit that keeps the generic triangular proposition
  separate from its degree-four-through-twelve witnesses and requires the
  explicit “does not exclude `(75,125)`” boundary.
- [x] Recheck finite normalization and reconstruction scripts for saturation,
  nilpotent thickness, and component removal.  Foundational and cancellation
  boundary completeness use the canonical finite normalization, both affine
  and boundary primes, and exact local degree sums; reduced incidence and the
  `mr(m+1)` nilpotent contact thickness stay distinct.  The cubic frontend
  separately computes support-bidual and cotangent saturation, radical support,
  special-fibre lengths, and six-generator different modules.  The audit found
  stale prose reopening the already-closed squarefree formal-tail saturation
  queue, not a missing saturation in the checkers; that route is now corrected
  in `OP-SUSP` and its canonical notes.
- [x] Re-audit cubic normalization, minimal-boundary, and controlled-suspension
  notes for parallel stale status narratives.  The active cubic-closure note
  incorrectly left the six singular 24-parameter saturation families open
  despite the later `KDSQ6`/`SSADPALL` all-orders closure; the normalization
  frontend and universal-cubic testbed retained the same obsolete boundary in
  earlier sections, and the quartic-kernel note omitted the completed
  squarefree continuation from its status and reproduction surface.  Correct
  all affected notes plus the `OP-SUSP`
  dependency/scope/forbidden-route ledger: the live gates are
  boundary-geometric `S2`/local-CM or Cartier input, the three non-squarefree
  leading symbols, global Keller compatibility, and coefficient base-change
  rigidity—not another finite quartic-axis sweep.
- [ ] Identify high-risk theorems with no independent replay and add compact
  low-dependency witnesses where feasible.
- [ ] Keep presentation-free stable descent, global minima, and full stable
  moduli open unless an exact new certificate closes them.

## P2 — GVC, SIC, Image/Mathieu, Ritt, and extended geometry

- [ ] Recheck every finite-degree or finite-support GVC/SIC census for exact
  domain enumeration, symmetry quotienting, and no missing support orbit.
  **Completed subset:** the complete bidegree-`(2,2)` 501-support regression
  and the bidegree-`(3,3)` size-three through size-nine chain now have
  exact key/coverage checks.  The bidegree-`(4,4)` two-row layer now also has
  solver-free exact-artifact audits for the dense degree-604/mu8 chart, all
  135 separated-row boundary orbits, all six dense off-diagonal
  representatives covering ten row pairs,
  all 1,174 off-diagonal boundary orbits, all 942 stored rank-two unit opens,
  both delayed-fibre records, all 60 single-shear charts, and the exact
  150-chart/78-orbit double-shear reversal cover.  Other extended-geometry
  censuses remain queued.
- [ ] Recheck moment and Fitting ideals for the intended saturation and marked
  scheme thickness.
  **Completed completed-moment subset:** `SIC2MDF35`, `SIC2MSP35`,
  `SIC2MIA36`, and `SIC2MCR` now pin their committed output ledgers and direct
  helper sources. The canonical note's automatic-ledger command was missing
  `--ladder-beta-check 32` and cited a superseded output hash; both are
  corrected. The older bounded relation ledger is now identified as the
  `c91498bbca85` historical schema instead of being presented as a
  byte-identical output of the extended current producer. Static review keeps
  all finite-field relation exclusions support-bounded, all Hilbert tests
  necessary-only, and every fixed-field conclusion restricted to its diagonal
  or single-phase slice. A stdlib-only auditor checks bytes and stored scope
  flags; no SymPy invariant scan or Singular process was run. Fitting and
  remaining moment programmes stay queued.
  **Completed degree-four field subset:** `SIC2C4I`, `SIC2C4D`, `SIC2C4SP`,
  `SIC2C4P1M`, and `SIC2MA` now pin all six committed ledgers and the helper
  sources that actually prove their composite scopes. In particular,
  `SIC2C4D` no longer relies on its diagonal checker alone for the separate
  22-even-parameter and weight-16 nonrelation assertions. Static review found
  no sieve promotion: full modular column rank excludes only the declared
  relation support, rank defects remain candidates, diagonal and fixed-locus
  fibers remain slice results, and the four-point phase-one completeness
  remains only over `F_101`. Its rational branches and orbit identity are
  exact, but additional characteristic-zero components and generic degree are
  still open. The pinned phase-one path asserts solver success before writing;
  timeouts produce no certificate. No SymPy, Singular, or `msolve` run was
  performed.
  **Completed quartic `q_2` normal-jet subset:** `SIC2Q2NJ4` and
  `SIC2Q2C3D4` now pin their two committed artifacts and the shared producer.
  The maintenance audit preserves all three logical boundaries: dimensions
  six and four are finite-field quadratic/cubic jet statements at one
  normalized point, the quartic dominant-sheet radicals do not prove formal
  isolation, and the `F_2=0` boundary plus characteristic-zero/global
  integrality remain open. The 300-second quartic row is a checked timeout;
  the separate 240-second off-axis timeout is a retained historical field
  serialized by the cubic checker rather than rerun by it, with inference
  explicitly `none`. No SymPy or Singular calculation was run.
- [ ] Recheck polarization steps for division by a characteristic-dependent
  scalar or an implicit genericity assumption.
- [ ] Recheck the two-pair all-order and bidegree results against their explicit
  counterexamples and status replacements.
- [ ] Recheck Ritt/Hurwitz deformation computations for formal-to-algebraic
  promotion and for omitted obstruction stages.
- [ ] Recheck tangent-space conclusions: a vanishing ordinary quotient does
  not erase filtered or boundary-decorated moduli.
  **Status-layer repair complete:** `O1` now distinguishes Artin/formal source
  triviality, order-dependent reduced representatives, and failed global
  algebraization; checker-by-checker Ritt cotangent review remains queued.
- [ ] Re-run cheap exact coefficient and independent-Python/Singular/Lean
  controls before any large Gröbner refresh.
- [ ] Keep minimum GVC(3) degree, efficient polarization, and complete
  bidegree-(3,3) SIC classification open.
- [x] Re-audit `RSMCENSUS`: all eleven declared `(d,r)` rows are present;
  modular full rank is used only for characteristic-zero algebraic
  independence, negative Hilbert coefficients only for finite-prefix
  semistable existence, and every nonnegative or palindromic observation is
  explicitly bounded through degree 85 rather than promoted to an hsop or
  all-order theorem.
- [x] Replace count-only coverage in the retained cubic Gaussian null-cone
  checker with exact reflection-orbit comparisons for its two-, six-, and
  seven-weight domains.  Its cheap audit now also requires the exact three
  unresolved keys handed off by the five-weight frontier and matches all 31
  stored input hashes to literal characteristic-zero unit outcomes.
- [x] Make the predecessor chain for the cubic four-weight closure explicit:
  reconstruct all 33 mixed-sign supports, require the earlier
  four-support/24-chart complement, remove the separately certified symmetric
  four-chart support, and check that the resulting three supports and 20
  charts are exactly the seven stored rational unit systems.
- [x] Reconcile the retained cubic theorem with the later all-degree `G2T`
  result: `G2T` subsumes its GMC consequence, while `G2N` correctly remains a
  separate finite radical-containment and moment-ten strengthening rather
  than being marked wholly replaced.

## P2 — plane Jacobian programme

- [x] Re-audit the closed `(72,108)` coefficient row from its Proposition-4.3
  tail census, quotient first block, exact Case-1 systems, and archive
  certificates rather than a navigation summary.  The compact determinantal
  checker previously generated a Singular program but did not itself compose
  the exact unit decision, two special-fibre certificates, or sign-branch
  transport named by `PJ72108D1`.  Pin its three direct archive inputs and add
  a true top-level verifier.  The cleanup mode validates 34 manifest entries,
  three additional transport inputs, all four reconstructed residual rows,
  and the special-fibre Bezout identity without running Singular or the 89 MB
  multiplication.  The complete replay remains available as an explicit
  separate mode.  The alternative no-vertical Belyi closure now also has a
  solver-free audit that pins the certificate and quotient graph, requires
  five distinct 21-sheet permutation representatives, checks all three
  rank/nullity ledgers, and preserves the saturated `B_8 != 0` terminal unit
  record.  Both rows still depend on the published general reduction.
- [x] Re-audit `PWB1--PWB6` and the two committed plane wild-boundary ledgers
  without replaying normalization, point counting, packet enumeration, or
  symbolic identities.  The status rows now pin the ledgers and the directly
  imported boundary-lattice helper.  A stdlib-only auditor checks their exact
  hashes, row partitions, fail-closed source paths, and scope boundaries.  In
  particular, the empty post-support balanced reconstruction queue does not
  discard the 20 original prescribed-cover and 6 comparison rows still marked
  `needs_reconstruction`; `packet_gate_only` rows remain abstract necessary
  templates; and the six-row `F_3` degree-seven scan does not exclude other
  plane-cover architectures.  No SymPy or Singular calculation was run.
- [ ] Recheck every local conductor, cusp, node, and connector calculation for
  saturation and compatible global gluing assumptions.
- [ ] Recheck pole and valuation ledgers for omitted infinity components or
  incorrect common-power cancellation.
- [ ] Recheck finite-cover and monodromy claims for actual braid/meridian data,
  not abstract endpoint matchings.
- [ ] Verify the F2 carrier and outgoing-tail handoffs do not present local
  compatibility as a global compactification.
- [ ] Keep the programme parked where actual global carrier/compactification
  data are missing; no bounded chart census proves `JC(2)`.

## P2 — HC4, DC2, and Hessian programme

- [ ] Recheck scalar, rank-one, rank-two, repeated-factor, and Schur packets for
  a complete pivot stratification.
  **Completed clean split-linear subset:** static review of `HC4NHM4--12`
  confirms the exact denominator partition chain `3+1`, `2+2`, `2+1+1`,
  and `1+1+1+1`.  `HC4NHM4` remains a nonempty leading
  Hessian--Schur packet rather than an exclusion; `HC4NHM5` excludes only its
  genuinely two-component prolongation.  The squarefree domain is partitioned
  exactly into 16 all-concurrent, 8 triple-concurrent/transverse-fourth, 16
  no-three-concurrent, and 8 triple-concurrent/tangent-fourth flag patterns.
  The notes and status scopes retain the positive-defect, lower-Smith,
  nonlinear-denominator, and rank-at-most-two frontiers.  `HC4NHM10--11` now
  carry their missing `narrows_problems` edges to `OP-HC4-D5`.  This cleanup
  did not rerun the `HC4NHM12` Singular saturations.
  **Completed scalar relative-pencil subset:** the staged degree-six frontier
  in `HC4RSD24` was still being repeated in the introduction to
  `HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md` even though the same note later
  closes degrees six and seven and `HC4MR1` closes the complete auxiliary
  relative-nilpotent pencil branch in all degrees.  The headline now names
  that current endpoint while preserving the theorem-local historical
  frontier and the essential warning that neither `HC4` nor `JC2` follows.
  The primary `HC4MR1` checker proves only the final local moving-frame step;
  its previously unregistered committed prolongation artifact is now a
  software lock with a no-rewrite audit that requires the artifact's explicit
  written-proof boundary.  No prolongation calculation was rerun.
  The shared `HC4RSD17--23,25--28` artifact is now registered on all eleven
  status rows as well.  Its apparent degree-seven `open_frontier` is retained
  as stage-local historical provenance and explicitly marked superseded by
  `HC4RSD40`/`HC4MR1`; its audit-only mode checks that distinction without
  regenerating the symbolic ledger.
  The same stage/current distinction is now enforced for the shared
  `HC4RSD11--16` scalar-dichotomy ledger: its higher-degree pencil handoff is
  preserved as historical `HC4RSD16` data, registered on all six consumers,
  and explicitly superseded by `HC4MR1`.  No quadratic or higher-degree
  symbolic identity was rerun.
  The three quadratic-pivot predecessor ledgers now have the same fail-closed
  treatment.  Their `HC4RSD8`, `HC4RSD9`, and `HC4RSD10` frontier fields are
  preserved byte-for-byte as historical stage data, but their audit modes and
  canonical notes point to `HC4RSD9--10`, `HC4RSD11--16`, and `HC4MR1` as the
  corresponding closures.  All three artifacts are newly registered on their
  theorem rows; no block determinant or normal-form calculation was rerun.
  **Completed affine-pivot scope split:** the `HC4RSD6` constant-span-deficient
  frontier is still meaningful for classifying affine Schur representations,
  but `HC4RSD7` already bypasses it for collisions inherited at a common pivot
  value.  The latter does not cover nonlinear pivots or collisions joining
  different affine fibers.  Both committed ledgers are now registered and
  have audit-only modes that enforce this distinction without rebuilding the
  determinant identities.
  **Completed scalar-kernel predecessor subset:** the five `HC4RSD1--5`
  ledgers and their shared equation helper are now registered and
  cleanup-auditable without importing SymPy. `HC4RSD1--2` also register the
  projective-polar atlas they actually consume. The committed artifact bytes
  remain unchanged, but their stage-local frontier text is no longer treated
  as current: `HC4RSD2--5` close all affine-in-`x` primitive kernel lines and
  fixed primitive two-component constant-support kernels, while `HC4MR1`
  closes the nonzero-corner auxiliary constant-Hessian-pencil branch. Fixed
  three/four-component and parameter-moving nonlinear singular kernels,
  nonlinear zero-corner exact remainders, and moving matrix pivots remain
  explicit. No symbolic identity or search was rerun.
  **Completed pure-sextic/septic artifact subset:** all twelve committed
  `HC4RSD29--40` stage ledgers were present but absent from their status-row
  software locks. They are now individually registered and covered by one
  stdlib-only batch auditor that checks exact bytes and status/format mapping.
  Their generators were not changed or run: `HC4RSD32` and `HC4RSD40` remain
  the degree-six/seven endpoints, and `HC4MR1` remains the later all-degree
  endpoint only inside the auxiliary relative-nilpotent pencil branch.
  **Completed nonreduced predecessor navigation subset:** static review of
  `HC4NHM1--3` found no hidden artifact or permissive CAS fallback. The
  `HC4NHM2` checker fails closed when Singular is unavailable, while each
  checker explicitly separates its executable identities from the written
  normalization/DVR proof. The stale `HC4NHM2--3` “next” text is now joined
  to current knowledge: `HC4NHM4--12` close every split-linear clean
  denominator partition. Clean nonlinear denominator components,
  positive-defect packets, and lower-Smith strata remain open. No checker was
  executed and no mathematical claim was strengthened.
  **Completed projective-polar provenance subset:** the `HC4PPG1--9` rows
  previously named only the Python environment even though their checkers
  share a Segre helper and form an atlas -> Rees/strata -> vertex and
  conditional-sieve artifact chain. Each row now registers the exact
  artifacts it produces or consumes, and one stdlib-only auditor validates
  all five committed ledgers and the helper by hash and semantic scope. It
  explicitly preserves the 624 surviving rows as necessary numerical
  configurations, not realized candidates, and retains higher torsion and
  exceptional codimension-three packets. No symbolic or Macaulay2 replay was
  run.
  **Completed all-dimensional Segre provenance subset:** `PGS1--PGS3` now pin
  their three committed ledgers, shared implementation helper, and four
  Macaulay2 calibration sources.  Their Python verifiers previously rewrote
  generated artifacts on every ordinary run; default behavior now compares
  exact committed bytes and requires explicit `--write` for replacement.  A
  stdlib-only auditor checks all 84 transform regressions, the exact 270- and
  2,160-row parameter-key domains, and the distinctions among complete Segre
  vectors, top-degree-only controls, uncomputed families, smooth generic
  first-Segre laws, and singular DVR-profile laws. No SymPy or Macaulay2
  calculation was run.
  **Completed mixed-canonical artifact subset:** all ten committed
  `HC4MCP1--10` ledgers are now registered on their status rows together with
  the directly imported source files needed to interpret them. A stdlib-only
  auditor checks their exact bytes and semantic scope without importing
  SymPy, invoking Singular, launching a search, or rewriting an artifact.
  Static inspection confirms that the early modular breaks discard a chart
  only after unequal values give a valid characteristic-zero nonconstancy
  witness; agreements remain survivors for exact handling. `HC4MCP1--4` and
  `HC4MCP7--9` retain their declared finite boxes. `HC4MCP10` removes the
  coefficient and constant-direction bounds only for the 54 `HC4MCP6`
  resonance families; moving pivots, other supports, and longer words remain
  open.
  **Completed nonlinear clean-denominator static subset:** `HC4NHM13`,
  `HC4NHM15`, `HC4NHM18`, and `HC4NHM19` have no generated ledgers to
  register. Their checkers are self-contained; the two Singular-backed entry
  points require the executable and fail closed if it is missing. The notes
  already preserve the exact chain from conic-divisible exclusion through
  at-most-four-root double-conic closure and the unsaturated many-root warning,
  while the smooth-cubic degree packets `(2,7)`, `(3,6)`, and `(4,5)` remain
  open. Only the stale “next action” heading was relabeled historical. No CAS
  command was run.
  **Completed relative rank-two navigation subset:** the seven local proof-map
  notes `HC4RSD56--63` still exposed their discovery-time “what remains” and
  “next target” sections without saying that the targets were later closed.
  Each now identifies its historical stage and exact successor: `HC4RSD60`
  closes `[2,2]`, `HC4RSD63` closes `[3,1]`, and `HC4MR1` subsumes both in the
  all-degree auxiliary relative-nilpotent branch.  This is navigation cleanup,
  not a stronger claim; unrestricted `HC4` and `JC2` remain untouched.
  **Completed repeated-linear direct subset:** the `HC4DIR27` headline stopped
  at a three-packet lower-rank sextuple reduction even though later sections
  in the same canonical note synchronize those packets to one degree-five,
  order-one pure-cube scalar-parent resonance and `HC4DIR28` excludes that
  family as a collision source.  The current headline now matches the proved
  endpoint while retaining nonlinear factors, multiple repeated factors,
  multiplicity at least seven, and top-Hessian rank at most two as separate
  frontiers.  Its previously documented but unregistered generated artifact
  is now a software lock with a no-rewrite audit of the exact endpoint and the
  written UFD/DVR proof boundary.  No symbolic identities were rerun.
- [ ] Recheck denominator extraction and localization before interpreting a
  Schur complement or moving kernel.
  **Completed historical-failure guard:** `HC4NHM21` already records the
  explicit lower-Smith counter-witness to unsaturated discriminant membership
  and correctly leaves the clean invariant elimination open.  Its checker now
  hash-pins the normal-layer helper it imports and has a cleanup-only mode that
  validates that provenance without replaying symbolic algebra.  This prevents
  helper drift from silently reviving the false unsaturated route; the wider
  denominator/localization packet review remains queued.
  **Completed committed Fitting-denominator subset:** the previously
  unregistered cube-torsion and fourth-power ledgers are now the explicit
  partial record `HC4QSE5`. Static review confirms that the fourth-power
  scanner's early exit is monotone: it fires only once all fifteen targets
  already lie in the current column span. The status and a stdlib-only auditor
  enforce the real boundaries instead: primitive rows are interpreted only on
  `D(nu)`, the scans do not cover proper finite-field extensions or prove a
  characteristic-zero support theorem, and all timed-out symbolic routes have
  mathematical conclusion `none`. The additional rational point is only a
  cube-to-fourth-power nilpotence jump, not a reduced Schur component. No
  Singular calculation or finite-field scan was rerun.
- [ ] Recheck all “all degree” frontends for the exact hypothesis under which
  induction or recurrence closes.
  **Completed squarefree direct-filtration subset:** `HC4DIR2` is genuinely
  degree-free only under characteristic zero, top-Hessian generic rank three,
  and squarefree nonzero ternary Hessian determinant.  Its universal step is
  the written factorwise divisibility/negative-degree proof, not the checker's
  two sample determinant faces.  The stale status claim that an off-diagonal
  term was removed by a triangular change is corrected: the proof forces the
  assumed first nonzero term to vanish.  The documented committed identity
  artifact is now registered and has a no-rewrite, no-symbolic-replay audit
  path.  Rank at most two and the nonsquarefree locus remain open.
  **Completed diagonal/Meng--Yang subset:** `HC4FSD1--3` and `HC4MYGJ2`
  correctly separate degree-free written arguments from the committed bounded
  degree-4--8/order-1--12 discovery tables.  The latter are now source- and
  artifact-locked on the three results they actually regress, and their scope
  grid has a cleanup-only audit.  The artifact's old recommendation to
  control arbitrary lower layers is retained as historical provenance but is
  explicitly superseded by `HC4FSD3`; formal Meng--Yang solvability still
  does not imply polynomial termination or collision preservation.
- [ ] Recheck modular searches for good-prime coverage and characteristic-zero
  reconstruction.
  <!-- status-consumer: HC4FF1 c659c03cbc63d15e -->
  **Completed collision-first finite-field subset:** the formerly
  narrative-only `GF(11)`/`GF(13)` campaign is now registered as `HC4FF1`.
  A cleanup-only audit pins all five artifacts and three generating scripts,
  reconstructs the exact 45,181,194-choice nonempty one/two-direction domain,
  checks the 96+32+128+144 selected dense-family keys, and requires all 800
  support-prime outcomes to be unit with no timeout/error.  The open-problem
  ledger now explicitly forbids promoting this fixed-normalization,
  fixed-prime, sampled-support evidence to characteristic zero or unrestricted
  `HC4`; no search was rerun.
- [x] Recheck smooth quartic polar and exceptional-slice decompositions for
  omitted tangent or coisotropic strata.
  **Completed static scope/provenance audit:** `HC4NHM14`, `HC4NHM16--17`,
  `HC4NHM20`, and `HC4NHM22--24` retain the correct theorem boundaries.  The
  first visible divisor is the polar of the reducible binary resultant, not a
  tangent/bitangent/flex discriminant of the prospective quartic, so no such
  contact stratum was silently removed.  The generic smooth polar conic and
  every generic two-line component type are excluded, but parametrization,
  linear-pivot, witness-determinant, secondary, and hidden-denominator strata,
  the complementary residual-line chart, and all other reciprocal boundary
  types remain open.  Four imported equation-builder dependencies are now
  authoritative software locks; their checkers have provenance-only modes
  that do not construct equations, exact fields, or Singular jobs.  The
  Fermat-symmetry note's corrupted inline mathematics was also repaired.  No
  symbolic, exact-field, or Singular certificate was rerun.
- [ ] Keep mixed/coisotropic Schur descent and direct degree-five residual
  packets open; do not summarize the packet library as unrestricted `HC(4)`.

## P2 — formalization, papers, and external review

- [ ] Map each Lean theorem to the exact registry scope and list all arithmetic
  or geometric hypotheses that remain external.
- [ ] Re-run the dependency-minimal formal targets before describing a result
  as formally verified.
- [ ] Confirm external certificate URLs, prover identities, and refereed flags
  remain accurate.
- [ ] Recheck paper theorem statements against current status and generated
  artifacts before rebuilding PDFs.
- [ ] Keep a publication’s theorem, the repository’s independent replay, and a
  new derived corollary as three distinct evidence records.
- [ ] Refresh literature claims only from primary sources and record the exact
  version/date when the distinction affects scope.

## P2 — archive, supersession, and navigation

- [ ] Verify every active README points to current canonical notes rather than
  a superseded handoff.
- [ ] Replace obsolete operational instructions by short tombstones linking to
  the new route; preserve scripts and evidence.
- [ ] Verify archived narratives are not canonical sources or active checkers.
- [ ] Recheck archive manifests and hashes after any move.
- [ ] Keep historical script-input hashes unchanged when they document the code
  actually run; do not confuse them with current-source replay claims.
- [ ] Add explicit historical labels to notes that still say “current,” “next,”
  or “open” after a later theorem closed the route.
- [ ] Ensure `RESEARCH_TIMELINE.md`, `KNOWLEDGE_BASE.md`, programme READMEs, and
  process/failure ledgers agree with `MATH_STATUS.json` without becoming
  competing status authorities.

## Completion criteria

- [ ] Every registry entry has a recorded semantic review outcome and every
  active checker has either a passing exact replay or a named bounded/blocked
  reason.
- [ ] Every complete negative search has exact domain, coverage, exception,
  and merge certificates plus at least one adversarial omission/duplication
  control.
- [ ] Every partial or bounded result exposes its first unproved gate in both
  canonical prose and registry scope.
- [ ] Every mutable external input needed for replay is preserved by a
  sufficient exact projection or immutable source.
- [ ] All repaired checker and certificate hashes propagate through dependent
  artifacts and documented output digests.
- [ ] The generated status, discovery timeline, knowledge base, process atlas,
  failure ledger, and programme navigation agree with the final registry.
- [ ] The final cheap suite passes; long calculations are either replayed with
  checkpoints or left as explicit actions without theorem promotion.
- [ ] The closeout report distinguishes confirmed discoveries, corrected
  scripts/claims, historical immutable drift, current-work blockers, bounded
  misses, and genuinely new theorem-directed opportunities.
