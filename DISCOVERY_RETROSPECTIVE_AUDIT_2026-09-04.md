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

The working registry currently has 1,150 entries: 995 proved, 98 partial, 17
open, 18 parked, 16 archived, and 6 falsified. The first recent-history pass
compared the 85 entries committed after `dada420` through the then-current
September 4 head and found 77 distinct checker programs. Concurrent September
4 work has since added further status entries; those belong to the live-current
queue below and must be re-snapshotted at audit close.

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
- [ ] Re-run the full registry path/dependency/consumer validation after the
  live September 4 files stop changing.
- [ ] Resolve the current checker-source hash drift for the actively edited
  inherited-cover and complete-character-closure scripts only after their new
  `08f72`/`08f72`-adjacent artifacts are final. **Blocked-current.**
- [ ] Regenerate `STATUS.md` through `scripts/render_status.py` only after the
  preceding checker hashes are stable; never hand-edit it. **Blocked-current.**
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

## P0 — discovery and process ledgers

- [x] Extend `RESEARCH_TIMELINE.md` through the direct alternate-Q80 equation,
  noncyclic closure, complete public atlas, singular/genus-two boundary,
  quotient-class exclusions, and arithmetic marking obstructions.
- [x] Update `KNOWLEDGE_BASE.md` with target-free compilation, arithmetic
  marking, finite-quotient slicing, sharded coverage, formal-local boundaries,
  and mutable-source projection lessons.
- [x] Expand the K3 process ledger from 55 to 73 events and from 17 to 23
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
- [x] Add a short repository-root navigation link to this worklist and keep it
  out of the mathematical-status hierarchy.
- [ ] At close, verify that every new September 4 status entry is represented
  in the research timeline or is deliberately omitted as an implementation
  detail.

## P1 — per-entry registry audit protocol

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
- [ ] Repair the historical `wgxli` lineage replay so it can validate its
  September 1 projection offline; the live database has grown from 474 to 556
  records although the five target values used by the fingerprint remain
  unchanged.
- [ ] Add output hashes to canonical notes when a certificate is primary and
  currently lacks a documented digest.

## P1 — recent K3 and elliptic-curve discoveries

- [x] Replay the clean Python checkers for ICARM bounded rebasing, rootless-J1
  bound, deep-cover quotients, visibility complexity, rank-28 character glue,
  NS0024 relative-`U`, both inherited-product modes, and the V4 rank screen.
- [x] Replay the 43-chart/474-curve exact lineage atlas and confirm six
  `PGL2` classes, 69 hits, and curve 12 quotient `Z^12`.
- [ ] Replay every checker behind the seven-fibre native quotient/visibility
  audit after the current `08f72` extensions settle. **Blocked-current.**
- [ ] Compile saturated chart transports for the 57 recognized fibres whose
  quotient dimensions remain unknown; preserve nulls until each exact replay.
- [ ] Group any predictor training/test split by the six family classes so
  fibre-level leakage cannot inflate performance.
- [ ] Record actual search exposure denominators before comparing public hit
  counts across families; do not infer effort from successes.
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
- [ ] Run the same rational-marking/moduli screen on determinant 720 before any
  coefficient search; its known rational `3A5` point has determinant 20 and is
  only a negative source control.
- [ ] Rerank all remaining different-NS candidates using arithmetic
  admissibility before equation incidence or planner cost.
- [ ] Keep the NS0024 and NS0031 local/geometric artifacts as controls while
  preventing either from re-entering the rational-source queue.
- [x] Replay the rank-28 pilot, mixed-trace and simultaneous-splitting bounded
  searches, 103b2 hard-fibre products, alternate laboratory, control-j
  preimages, rationalized-D6 chart, and V4 shortlist.
- [ ] Complete the long 103b2 MW-lattice, 63,917-class product-inversion, and
  high-throughput splitting replays with checkpoints.
- [ ] Audit every `PASS_EXHAUSTIVE_BOUNDED` status string so “exhaustive” is
  always visibly qualified by its finite height/trace/source box.
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

- [ ] Search all active K3 scripts for unqualified ADE/MW equality used as an
  object identity; require marked `F`, `O`, components, chamber, and integral
  transports.
- [ ] Search for pseudo-zero or chamber-zero reuse and require equation-side
  effectivity before any route score or continuation.
- [ ] Search binary-quartic and cover maps for missing degree-two accounting;
  regression-test heights and collision degrees.
- [ ] Search rational normalization and CRT code for incomplete denominator
  clearing, especially omitted `Dx`-type factors.
- [ ] Search resolved-RR code for row-per-blow-up-centre condition counts and
  replace them with saturated connected-component modules.
- [ ] Search specializations for inherited generic ranks, roots, pole orders,
  or MW coordinates; force a new typed specialization state and recomputation.
- [ ] Search neighbour enumerators for first-hit or beam-pruned output used as
  a completeness claim.
- [ ] Search ranked-prefix consumers for mathematical rejection of candidates
  outside the prefix.
- [ ] Search sharded campaigns for coverage by counts alone; require exact key
  sets or independently certified disjoint domains.
- [ ] Search modular/Hensel outputs for language implying algebraization or a
  rational characteristic-zero point.
- [ ] Search geometric Picard/MW claims for promotion to arithmetic rank without
  Galois and full-marking evidence.
- [ ] Search quotient-sliced conclusions for accidental promotion of one class
  to all quotient classes.

## P2 — elliptic curves over `Q`

- [ ] Re-audit every claimed rank lower bound by literal point substitution,
  minimal-model transport, torsion handling, and independent finite-quotient
  rank.
- [ ] Re-audit every exact-rank phrase for a matching unconditional upper bound;
  curve 302 remains rank at least 31, not exact rank 31.
- [ ] Re-audit conductor claims through global minimization and local Tate data,
  never a discriminant radical.
- [ ] Check all Selmer/descent outputs for the distinction between Selmer
  dimension, Mordell--Weil quotient, and possible Tate--Shafarevich classes.
- [ ] Keep timeout, `norm-one`, and incomplete local-class outputs fail-closed
  in the rank-32 gate.
- [ ] Require residual dimension at least 15 on the same minimal curve before
  any expensive rank-32 cover or point campaign.
- [ ] Revisit curve 273, curve 302, and low-conductor near misses only through
  completed residual descent or a sharper unconditional upper-bound route.
- [ ] Preserve exact native K3 fibres as calibration controls without assuming
  rational parameter transport between different fibrations.
- [ ] Audit all public-database checkers for mutable-source preservation,
  deterministic point ordering, and snapshot-specific claims.

## P2 — foundational Keller, cancellation, and stable boundary programmes

- [ ] Re-read every core theorem from `F1` through the finite-etale and
  multiplicity layers against the current falsification graph.
- [ ] Recheck marked-root reconstruction for omitted boundary points and
  accidental promotion of an affine inverse to a global one.
- [ ] Recheck valuation/cancellation scripts for polynomiality after all
  denominators and for full source--core--target determinant balance.
- [ ] Recheck stable-multiplicity claims against the universal-cubic gradedness
  counterexample and keep only the boundary-minimal theorem where applicable.
- [ ] Recheck finite normalization and reconstruction scripts for saturation,
  nilpotent thickness, and component removal.
- [ ] Re-audit cubic normalization, minimal-boundary, and controlled-suspension
  notes for parallel stale status narratives.
- [ ] Identify high-risk theorems with no independent replay and add compact
  low-dependency witnesses where feasible.
- [ ] Keep presentation-free stable descent, global minima, and full stable
  moduli open unless an exact new certificate closes them.

## P2 — GVC, SIC, Image/Mathieu, Ritt, and extended geometry

- [ ] Recheck every finite-degree or finite-support GVC/SIC census for exact
  domain enumeration, symmetry quotienting, and no missing support orbit.
- [ ] Recheck moment and Fitting ideals for the intended saturation and marked
  scheme thickness.
- [ ] Recheck polarization steps for division by a characteristic-dependent
  scalar or an implicit genericity assumption.
- [ ] Recheck the two-pair all-order and bidegree results against their explicit
  counterexamples and status replacements.
- [ ] Recheck Ritt/Hurwitz deformation computations for formal-to-algebraic
  promotion and for omitted obstruction stages.
- [ ] Recheck tangent-space conclusions: a vanishing ordinary quotient does
  not erase filtered or boundary-decorated moduli.
- [ ] Re-run cheap exact coefficient and independent-Python/Singular/Lean
  controls before any large Gröbner refresh.
- [ ] Keep minimum GVC(3) degree, efficient polarization, and complete
  bidegree-(3,3) SIC classification open.

## P2 — plane Jacobian programme

- [ ] Re-audit the closed `(72,108)` row from canonical carrier equations and
  exact boundary certificates, not navigation summaries.
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
- [ ] Recheck denominator extraction and localization before interpreting a
  Schur complement or moving kernel.
- [ ] Recheck all “all degree” frontends for the exact hypothesis under which
  induction or recurrence closes.
- [ ] Recheck modular searches for good-prime coverage and characteristic-zero
  reconstruction.
- [ ] Recheck smooth quartic polar and exceptional-slice decompositions for
  omitted tangent or coisotropic strata.
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
