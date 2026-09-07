# Research timeline

This is a curated chronology of the repository's mathematical development.
Dates describe when work entered or was corrected in this repository, not
priority or publication dates.  Commit messages and artifact mtimes establish
order; mathematical status comes only from [`MATH_STATUS.json`](MATH_STATUS.json).
The detailed, time-resolved Elkies--K3 chronology lives in
[`elkies-k3/ELKIES_K3_PROCESS_ATLAS.md`](elkies-k3/ELKIES_K3_PROCESS_ATLAS.md).

## 2026-07-20 to 2026-07-21: foundational map and first geometry

- The repository began from an explicit three-dimensional constant-Jacobian
  map with a rational three-point collision.
- The inverse-pencil, incidence, normalization, image, omitted-value,
  repeated-root, discriminant, and finite-field frameworks were separated
  into exact statements and checkers.
- The weighted marked-root construction emerged as the first all-degree
  extension.  Early overclaims were downgraded or archived as proof
  boundaries became clearer.

Durable lesson: start from the finite marked-root cover and distinguish its
normalization, affine reconstruction open, collision fibre, and omitted
boundary.

## 2026-07-22 to 2026-07-23: stable moduli, symplectic lifts, and cancellation

- Degree-five stable-moduli and marked-point dimension results were added,
  followed by symplectic/Weyl lifts and exact factorization bridges.
- Cancellation and plane-boundary programmes were reorganized around
  controlled suspensions, contact resultants, normalized factorization, and
  intrinsic boundary data.
- Hurwitz/Ritt calculations grew from individual intersections into
  deformation and synchronization problems.

Correction retained: coarse or presentation-dependent data do not by
themselves classify stable polynomial left--right equivalence.

## 2026-07-24 to 2026-07-26: arithmetic fibres and formalization

- The finite-etale Keller-fibre theorem was made constructive, including the
  rank-two exclusion, quotient reconstruction, explicit quintics, and
  local-to-global arithmetic compilation.
- Lean projects were consolidated and expanded for finite-etale Keller maps,
  GMC(2), and support saturation.  The formalized layers were explicitly
  separated from arithmetic, geometric, or literature inputs outside the
  prover.
- The paper registry and publication boundaries were established.

Durable lesson: generated examples become reliable research objects when a
compact certificate, independent arithmetic replay, written implication, and
formal layer have distinct interfaces.

## 2026-07-27 to 2026-07-31: GVC, SIC, Image, and high-dimensional searches

- Exact low-dimensional GVC and Gaussian-moment work produced the two-real
  theorem, homogeneous three-variable failures, radial propagation, and
  factorial/prime-power structure.
- Two-pair SIC and Image calculations found explicit counterexamples and
  sharp phase distinctions, while increasingly large modular and Gröbner
  searches mapped the unresolved bidegree-(3,3) boundary.
- Several broad searches were retained as bounded experiments after they
  failed to yield proof-sized systems.

Correction retained: moment vanishing, nullcone membership, and the Mathieu
property are distinct; unsaturated moment ideals can carry the wrong scheme
thickness.

## 2026-08-01 to 2026-08-05: plane-JC and GVC consolidation

- The audited plane `(72,108)` route was closed and the `F2` boundary work
  advanced through Kummer, residue-cover, Laurent, conductor, and logarithmic
  packages.
- The homogeneous GVC(3) counterexample and its spillovers were certified,
  while factorial and Laplacian frontiers were incorporated into papers and
  formal projects.
- Repository structure was consolidated around canonical notes, one status
  authority, and targeted reproduction commands.

Correction retained: local conductor length and affine singularities are not
automatically logarithmic point defects on a resolved Keller boundary.

## 2026-08-06 to 2026-08-10: HC4 and Hessian obstruction programme

- Scalar reverse-Schur and minimal-excess branches were closed across all
  degrees, followed by Wronskian, developable-image, moving-kernel, Krylov,
  Frobenius, and direct homogeneous-filtration reductions.
- Squarefree and repeated top-Hessian packets, affine-plane bridges, and
  projective polar strata substantially reduced the direct degree-five
  frontier.
- The surviving questions were explicitly split between direct `HC4`
  packets, mixed/coisotropic Schur descent, and lower-rank `DC2` or HN
  problems.

Correction retained: a method obstruction is not a collision exclusion, and
closure of many normal-form packets is not an unrestricted Hessian
conjecture theorem.

## 2026-08-11 to 2026-08-18: consolidation and elliptic-curve search

- The repository underwent a broad canonical-source and verification audit.
- A separate elliptic-curve programme introduced exact family arithmetic,
  Hensel--CRT--lattice searches, finite-quotient independence, conductor
  replay, and candidate-promotion gates.
- Exact rank-at-least-17, 19, and 20 frontiers and several family components
  were certified; negative scans were moved out of the active command
  surface.

Durable lesson: point lists, analytic scores, Selmer dimensions, and raw
discriminant radicals are different evidence types.  Exact rank and exact
conductor require separate upper-bound and local-minimality arguments.

## 2026-08-19 to 2026-08-21: rank-17 lattice recovery and H3 source

- The determinant-948 rootless `MW17` lattice and its 1,311 height-four pairs
  became the endpoint fingerprint.
- Direct large section solves proved poorly sized.  Reverse elliptic-neighbour
  searches traded MW rank for reducible-fibre roots and exposed low-rank
  equation charts.
- Shimura/Humbert/Kumar data identified the actual level-474 H3 source
  `E7+E8/MW2`; the first `q6` neighbour was executed exactly over
  characteristic zero.

Correction retained: an elliptic neighbour changes the embedded `U` on one
fixed Neron--Severi lattice; it does not create rank.  A CM specialization is
a different operation and may raise Picard rank.

## 2026-08-22: q8 repair and parallel Q80 route

- The generic Q80 low-q lattice corridor reached rootless `MW17`, and its
  CM24 specialization began an equation-level compiler route.
- The H3 q8 calculation was repaired after finding a duplicated 2-cover
  multiplier and a missing denominator factor.  The corrected exact
  `13 -> 2` Riemann--Roch calculation produced `D13/MW4`.
- Sixteen degree-46 diagnostics were archived and two obsolete IV-star
  markings were explicitly falsified.

Durable lesson: a full-rank modular obstruction only applies to the exact
marked, saturated module supplied to it.

## 2026-08-23: characteristic-zero Q80 closeout and H3 q24 route

- The difficult terminal Q80 section was recovered as a group-law difference
  of easier exact sections, closing the CM24 characteristic-zero route.
- The H3 `D13 --q24/orbit85--> D12/MW5` horizontal and resolved geometry were
  recovered, while a complete marked lattice chain was identified with the
  pinned R17 endpoint.
- Physical component resolution replaced blow-up-centre counting and exposed
  the correct local condition modules.

Correction retained: matching a sequence of ADE/MW labels does not establish
that two routes use the same marked fibrations.

## 2026-08-24: rank 31 and the A11/2A5 equation frontier

- ICARM curve 302 was independently certified to have rank at least 31,
  trivial torsion, a global minimal model, and exact conductor/local data.
  The retrospective now pins both this compressed certificate and the curve
  273 rank-at-least-30 certificate on their `ECR31`/`ECR30` rows, together
  with the shared finite-quotient helper and direct proof inputs. A
  dependency-free auditor enforces that these remain unconditional lower
  bounds only: heights, the public BSD/GRH statement, exact rank, and K3-family
  recognition are not inferred. No curve arithmetic was replayed.
- The H3 q24 section and `D12/MW5` equation were closed, then orbit42 gave an
  exact `A11/MW6` child.
- Construction-fingerprint and target-coset audits found the correct A11 q8
  divisor.  The split-I12 recurrence method then produced the exact
  component-9-zero `2A5/MW7` equation.
- Cheap changed-zero routes were re-audited against equation-effective
  sections.  Scores 4,199, 10,334, and 13,518 were withdrawn even though
  portions of their abstract lattice paths remain exact.

Durable lesson: route cost is meaningful only after the proposed zero and
horizontal curves are actual effective curves on the equation model.

## 2026-08-25: physical suffix and direct endpoint strategy

- The physical `q4/orbit208` equation to `3A3/MW8` was completed.
- A q323-free marked lattice route was selected:

  ```text
  3A3/MW8
    --q4/orbit1584--> D4+A3+3A1/MW7
    --q4/orbit164--> 2A3+2A1/MW9
    --q8/orbit376--> 4A1/MW13
    --q12/orbit5867--> rootless/MW17.
  ```

- By the end of the day, the first two edges of that suffix had exact characteristic-zero
  Riemann--Roch planes, quartics, minimal Jacobians, effective zeros, and
  equation markings. The `q8/orbit376` horizontal was also exact over
  `QQ(t)`; its resolved pencil and the preferred final `q12/orbit5867`
  equation were the remaining equation gates. `q12/orbit4484` was retained as
  a certified fallback, and both marked lattice transports and the endpoint
  isometry were already exact.
- Historical sources were re-audited.  They confirm the existence and broad
  modular/Hensel construction strategy for the rank-17 surface but do not
  provide the missing final equation, section basis, or neighbour maps.
  Therefore intermediate work is now minimized to what the next edge needs;
  the final rootless equation will be tested by a direct 17-section endpoint
  certificate.
- The status index was reconciled with eight committed checker-source hashes.
  This was a provenance-only repair: no Sage result was regenerated and no
  mathematical claim was promoted. See
  [`STATUS_HASH_REFRESH_2026-08-25.md`](archive/provenance-audits/STATUS_HASH_REFRESH_2026-08-25.md).

## 2026-08-26: rootless endpoint theorem and reverse-route audit

- Exact chord saturation closed the `q8/orbit376` equation and pointed its
  `4A1/MW13` child at P1229.
- A complete polynomial section-shell audit replaced the nominal q12 word by
  an equation-effective word. The resulting `q12/orbit5867` Jacobian has
  `24I1`, and seventeen exact sections reproduce the determinant-948 R17
  height lattice.
- A rational point on the stored quartic proved source identity. Counts at
  two good primes proved geometric Picard rank 19, and the unique possible
  index-two even enlargement was excluded. Thus the full geometric
  Mordell--Weil group is saturated R17 of exact rank 17.
- Reverse fixed-corridor work made the final historical arrows
  equation-explicit, then parked the route after preserving a corrected
  physical `4A1` lift. Arithmetic use of the endpoint became the priority.

Durable lesson: route, equation, section basis, Picard upper bound, saturation,
and source identity are separate proof layers; a direct endpoint package can
close them without completing every historical intermediate presentation.

## 2026-08-27 to 2026-08-31: published chart, calibrated arithmetic, and bisections

- The compact model and seventeen sections from Elkies's 2026 paper were
  replayed exactly and identified with the reconstructed q12 endpoint by a
  rational base change and Weierstrass scaling.
- The four published rank-25--28 fibres were locked as exact positive controls,
  with quotient gains `8,9,10,11` over the generic rank 17. A complete
  height-10000 scan showed that weakest-block scoring across three disjoint
  prime ensembles ranks all four inside the declared top one percent.
- A fail-closed residual 2-Selmer policy was implemented: rank 32 requires
  residual dimension at least 15 on the same minimal curve before covers or
  expensive point searches. Resource-bounded PARI and eclib attempts on the
  rank-28 control timed out without a Selmer result; an unconditional Magma
  path is generated but not completed on the current host.
- The first conic cover acquired an explicit eighteenth section and rational
  parameter. The paired cover acquired two explicit sections, exact maps from
  a rank-at-least-four elliptic curve, and a bounded heuristic sieve.
- Curves 285 and 286 received repository-local global-minimality and complete
  Tate-algorithm conductor replays, closing the low-conductor branch of the
  original operational target as well as the already-closed rank branch.
- Exact lattice enumeration found 39,120 rational-bisection translation
  orbits and 8,895,801 disjoint-priority pairs. All 39,120 classes were then
  constructed at equation level and gave distinct quadratic squareclasses.
  This proves injectivity on the complete survivor set, produces 39,120
  generic-rank-at-least-18 covers, and excludes this mechanism as a source of
  a rank-two collision or generic-rank-19 family.

Durable lesson: calibrate heuristics on exact controls, make descent gates
fail closed, and quotient bisections by section translation before spending
equation-level work on extension collisions.

## 2026-09-01: rootless completeness, public lineages, and source foundry

- Niemeier-first enumeration proved that the determinant-948 Neron--Severi
  lattice has exactly two rootless J2 frame classes.  This turned the alternate
  class from an observed candidate into a complete finite geometric target.
- Exact coefficient and first-jet fingerprints isolated the five-member wgxli
  public-curve lineage.  Signed-permutation and one-shear rebasing failures were
  retained as bounded model rejections rather than a nonexistence theorem.
- The q12/orbit5867 point factory received a complete two-primary boundary
  closeout, while small-prime cyclic-isogeny exclusions removed several simple
  explanations for exceptional published-R17 fibres.
- Equation-first lattice foundry work produced primitive determinant-720 and
  determinant-1184 candidate lattices, marked source data, and physical
  corridors.  At this stage these were geometric and finite-field objects, not
  arithmetic K3 sources over `QQ`.
- The complete declared prescribed-root source cover across thirteen rooted
  Niemeier ambients and sixteen D5 anchors produced 2,134 reduced-Gram MW1
  records covering all 48 foundry NS classes and no MW0 record.  The audit now
  pins this previously prose-only result in the registry and makes its merger
  reject a missing ambient shard, a changed NS-id universe, a changed input
  hash, or an anchor-count mismatch.

Correction retained: a complete geometric frame classification, a public-fibre
fingerprint, a bounded coordinate rejection, and a rational source equation are
four different evidence types.

## 2026-09-02: rational-source controls and predictor auditing

- Rational quadratic base changes supplied exact Picard-rank-19 E6 and E6+A1
  controls, complete first neighbor shells, explicit orbit-103 and orbit-96
  equations, and examples where geometric Mordell--Weil directions split under
  Galois.
- Target-fitted genus-one pencils explained all eleven exceptional directions
  at the published rank-28 fibre, while declared height-bounded simultaneous
  and mixed-trace splitting searches found no further split.  Degree-three and
  sampled degree-four visibility experiments remained partial.
- The determinant-720 Golay source and determinant-1184 NS0031 source acquired
  exact physical same-NS corridors to MW17 frames.  Their source-equation and
  arithmetic-marking gates remained separate.
- Early experiments and stale command surfaces were archived with manifests;
  exact negative results and regression inputs were preserved.

Durable lesson: target-fitted visibility is explanatory, not predictive.  A
geometric route to MW17 cannot be promoted to arithmetic MW17 until the source
and all required divisor classes are defined over the base field.

## 2026-09-03: integral rank transfer and the direct alternate-Q80 equation

- Integral involution eigensublattices, discriminant glue, theta convolution,
  bridge mutation, reverse masks, and root-system witnesses were consolidated
  into an exact rank-transfer and inverse-ADE layer.  Several cheap linear masks
  were shown to be nonselective and demoted from elimination gates to ranking
  data.
- The arithmetic Shioda--Tate theorem made the rational-source condition
  explicit: geometric rank transfer alone does not preserve arithmetic rank.
- A universal marked degree-two chord compiler and bounded relative-`U`
  completeness theorem connected target-free lattice planning to exact
  equation compilation without identifying planner reachability with arithmetic
  realization.
- Minimum-incidence search on the compact published R17 equation found the
  norm12/orbit11952 divisor.  Its direct degree-two hop produced a polynomial
  `24 I1` alternate-Q80 equation over `QQ`, seventeen saturated rational
  sections, and arithmetic generic rank 17.  The historical degree-11511 Q80
  transport and million-bit third-q12 reconstruction became provenance rather
  than operational routes.
- Exact preimage polynomials proved that the four published rank-25--28 controls
  are not rational fibres of alternate Q80, forcing native calibration.
- The alternate arithmetic laboratory certified 121 inherited covers, 7,260
  products, a cost-ranked 1,024-cover native prefix, 64 rational genus-one V4
  bases, and seventeen rank-one base Jacobians.  Complete smooth-bisection
  character maps on alternate Q80 and hidden `103b2` were injective; the native
  prefix remained explicitly bounded.

Correction retained: a completeness theorem for one component of a mixed
pipeline does not transfer to a ranked prefix, and an exact preimage miss on
one chart does not exclude a high-rank fibre on another chart.

## 2026-09-04: noncyclic closure, complete public atlas, and different-NS pivot

- The direct noncyclic chain `published R17 -> 4A1/MW13 -> published R17`
  became equation-explicit with thirteen saturated rational sections, maximal
  `Z/4+Z/8` bridge data, and target-free reverse selection.  This closed the
  determinant-948 construction control without creating a new surface.
- All 43 norm-twelve shared-zero equations and all 474 curves in the pinned
  ICARM snapshot were compared exactly.  The six rational `PGL2` j-classes give
  69 hits and 2,775 misses; all 376 native chart/fibre comparisons are
  untwisted.  Curve 12 supplies the first native alternate-Q80
  rank-at-least-29 control and has displayed quotient `Z^12` over the generic
  rank-17 subgroup.
- The point-forced cubic class-group calculation was normalized by the
  specialized generic `MW17` and extended to all five published-R17 controls
  plus native alternate-Q80 curve 12.  Every exceptional block adds zero
  bad-valuation rank; the certified residual 2-class-image lower-bound strata
  are `+5 -> 3`, `+6 -> 5`, `+8 -> 6`, and `+12 -> 10..11`.  This explains
  the full-BNF pressure across both frames by localizing the already-known
  exceptional Kummer information in global 2-class directions.  It does not
  explain why the rational points appear: the separation is largely forced by
  the zero valuation-rank identity, and prospectively a residual Selmer space
  still mixes Mordell--Weil directions with `Sha[2]`.
<!-- status-consumer: EC-K3-R17-KUMMER-CLASSGROUP-PRESSURE-COMPARISON 74b1dae24470b531 -->
- A later projection from the live 556-curve database preserved the 69
  recognized fibres and 1,545 displayed points.  Seven priority fibres now
  have exact displayed quotients and complete fixed-cover visibility spans;
  five alternate-Q80 fibres have exact fitted norm-eight incidence for all 51
  preferred quotient directions.  The other 57 quotient rows remain literal
  `UNKNOWN` pending chart-specific saturated section transports.
- Complete singular arithmetic-genus-one searches on alternate Q80 and hidden
  `103b2` found no nonsplit rational quadratic normalization.  Multi-prime and
  CRT genus-two searches also missed, but remain bounded; global genus-two
  injectivity is still unknown.
- Exact involution cohomology and minimum-norm enumeration reduced a possible
  height-eight zero Tate class to norm-eight and norm-twelve trace carriers.
  Complete inversion and 833 residual trace/target tests reject those carriers.
  This excludes the zero class only; existence, nonzero quotient classes, and
  product-twist ranks remain open.
- The determinant-950 NS0024 arithmetic route was closed negatively.  A full
  rational NS0024 marking would force a rational non-CM point on the relevant
  degree-475 Fricke quotient, contradicting the modular-curve obstruction and
  the geometric MW-rank-one frame.  Geometric NS0024 work and larger fields are
  not excluded.
- Determinant-1184 NS0031 briefly became the preferred different-NS
  replacement.  Its model-157 point has an exact one-parameter formally smooth
  `ZZ_7` branch and a five-edge physical corridor, but the split-Clifford
  modular curve maps to `X_0(37)` and neither noncuspidal rational point lifts
  through the required level-four non-split Cartan condition.  Thus a full
  rational NS0031 marking is impossible over `QQ`; its local and geometric
  certificates remain valid controls, and the remaining frames must be
  reranked through the arithmetic-marking gate before equation work.
- The arithmetic-marking classifier then applied that gate to all 66 exact
  rootless-MW17 candidate lattices in the 827-surface catalogue.  It classified
  one candidate as possible (the existing determinant-948 control), excluded
  determinant 720, NS0024, and NS0031, and left 62 literal `UNKNOWN`;
  consequently it emitted no equation-agent handoff. For determinant 720 the
  primitive coarse curve is `X_0(15)`, its stable marked curve is `X_0(60)`,
  and the latter has no rational noncuspidal point. The known rational `3A5`
  point instead saturates to determinant 20.
- Lane B was then reversed globally. The arithmetic-first planner orders all
  827 rank-three transcendental rows without rootless-frame data and admits
  `NS=T^perp`, rootlessness, or equation work only after an exact rational
  non-CM point on the full stable curve. Its 823-row research queue currently
  has 24 coarse genus-at-most-two diagnostics and no new positive handoff. The
  construction milestone was strengthened to require MW17 plus a certified
  positive-rank low-genus carrier and independent section, with an integral
  `V4` lattice of character ranks `17+1+1+1` as the stretch target.
- The retrospective registry audit promoted the already-tested incremental-CRT
  beam failure into a first-class status boundary.  On four prime groups, the
  width-one survivor has height 1409 while a discarded branch completes to
  `48/53`, of height 53.  Both active beam helpers now state that finite beams
  rank proposals and cannot exclude omitted CRT classes.
- The same audit removed the live ICARM endpoint from the theorem replay path.
  Committed sufficient projections now reproduce the historical `wgxli`
  artifact byte-for-byte, all 2,844 pinned curve/class decisions, the
  69-fibre/1,545-point snapshot, and the five exact section/quotient
  certificates; live retrieval is an explicit drift audit only.
- The residual rank-32 adapter retained its correct Selmer formula but renamed
  PARI `ellrank`'s third field from an apparent full `Sha[2]` dimension to the
  actual Cassels-pairing quotient rank.  A tested helper now enforces
  `dim Sel_2(E)=r2+dim E(Q)[2]+s`, and all 40 fail-closed policy regressions
  pass without running a descent or point search.
- The retrospective also reconciled the cubic cancellation narratives with
  their later all-orders formal results.  `KDSQ6` and `SSADPALL` already close
  the squarefree formal-tail saturation queue; the active `OP-SUSP` gates are
  now the boundary-geometric `S2`/local-CM or Cartier implication, the three
  non-squarefree leading symbols, global Keller compatibility, and coefficient
  base-change rigidity.  The falsified universal `OP-UG3` route remains only
  as a counterexample boundary.
- The finite-etale reconstruction checker no longer treats the identity
  `-0/2=0` as evidence that determinant-one target normalization preserves the
  distinguished fibre.  It now reconstructs in the exact quotient algebra and
  checks all three normalized target coordinates; the existing symbolic suite
  passes with the stronger assertion.
- Two SIC navigation passages were corrected to match the existing
  falsification graph.  The semistable Rodrigues point already disproves the
  bidegree-`(3,3)` moment--nullcone equality, but is SIC-safe; the remaining
  open target is the complete bidegree-`(3,3)` SIC classification, not revival
  of the stronger nullcone equality.
- The sparse SIC census chain gained fail-closed coverage checks.  The stored
  size-three through size-eight records must now match every exact relevant
  support key once, and the final size-nine checker reconstructs the union of
  all seven predecessor families instead of inferring “earlier closed” from a
  count-preserving set complement.  Cheap audit modes verify the complete
  560/1,401/3,864/7,588/11,200/12,780 and 11,420-support domains without
  rerunning msolve.
- The quartic two-row SIC atlas gained the same solver-free integrity layer.
  Exact-key replays now verify the dense degree-604/mu8 record, all 135
  separated-row boundaries, all six dense off-diagonal representatives and
  their ten-pair reversal cover, and all 1,174 off-diagonal boundaries exactly
  once.  They also rederive the stored strata, minor-open covers, 942 unit
  outcomes, and both delayed-fibre records without replacing the original
  characteristic-zero eliminations.
  The same pass pins all 60 single-shear labels and makes the double-shear
  150-chart/78-orbit reversal coverage explicit instead of count-only.
- The retained cubic Gaussian null-cone theorem now has an exact-key audit
  across its two-, five-, six-, and seven-weight closure systems.  It matches
  the three five-weight leftovers to their predecessor artifact, checks the
  complete reflection cover rather than aggregate chart totals, and binds all
  31 committed inputs to their QQ unit records without rerunning `msolve`.
  Its finite radical/cutoff content remains distinct from the later all-degree
  `G2T` proof even though that proof subsumes the GMC consequence.
- Its four-weight predecessor chain is fail-closed as well: the 33-support
  universe, the earlier four-support/24-chart complement, the separately
  removed symmetric four-chart support, and the final three-support/20-chart
  seven-system certificate now agree by exact keys rather than totals.
- The global low-degree Keller census received a solver-free integrity path.
  It validates the eight pinned artifacts and reconstructs the residual
  support orbits and determinant triples, then requires the same ordered 913
  support IDs through valuation, sign, all three modular, and exact rational
  ledgers.  This confirms the committed support-at-most-six routing without a
  new algebra run and leaves support-seven attainment and the unbounded census
  explicitly open.
- The plane sparse-support theorem now has a separate maintenance path as
  well.  It pins the existing JSON and Singular bytes and checks the exact
  arbitrary-degree claim records, bounded-regression labels, mask digests,
  exceptional shear formulas, and unique quadratic survivor without running
  the 14.6-million/5.29-million enumerations or a solver.
- Its affine-support/Newton bridge now has the same guardrail.  The committed
  audit requires all six Kummer character-pair brackets of the five-term F2
  terminal block and its nontrivial character profile, preventing the
  constant-Jacobian support-six theorem from being applied without the missing
  lower bands; it explicitly preserves the non-exclusion of `(75,125)`.
- The compact Case-1 `(72,108)` determinantal route now names a real verifier
  rather than its Singular-input generator.  Three directly consumed replay
  files are hash-pinned, and the top-level checker composes the adjacent-minor
  decision, both special-fibre unit certificates, and exact sign-branch
  transport.  Its cleanup mode verifies 34 manifest entries plus the compact
  reconstruction without invoking Singular or multiplying the 89 MB identity.
- The alternative no-vertical `(72,108)` Belyi closure also gained a
  committed-only guard: the certificate and quotient graph are pinned, all
  five 21-sheet permutation labels are distinct and well formed, the three
  deformation rank/nullity ledgers balance, and the terminal `B_8`-saturated
  unit record is required without reconstruction or a Singular run.
- The HC4 double-conic invariant gate was rechecked at the historical
  unsaturated-ideal failure.  Its theorem boundary was already correct: the
  lower-Smith witness invalidates unsaturated discriminant membership, clean
  saturation by `Phi2` is mandatory, and a discriminant still closes only the
  squarefree open.  The witness checker now hash-pins its harmonic-layer helper
  and offers a provenance-only cleanup mode, so no symbolic calculation is
  needed to detect replay-source drift.
- Static review of the clean split-linear HC4 quartic-denominator chain found
  the partition and flag coverage intact: `3+1`, `2+2`, `2+1+1`, followed by
  squarefree rows split as `16+8+16+8`.  The nonempty `HC4NHM4` leading packet
  is not confused with its failed prolongation, and all nonlinear,
  positive-defect, lower-Smith, and rank-at-most-two boundaries remain open.
  Missing open-problem narrowing edges on `HC4NHM10--11` were restored without
  rerunning their symbolic or Singular calculations.
- The `HC4DIR2` all-degree squarefree filtration was reconciled with its
  current proof.  The status ledger no longer invokes an obsolete triangular
  coordinate change: squarefree factorwise divisibility makes the putative
  first off-diagonal quotient polynomial of negative degree, so it vanishes.
  The checker is now explicitly only a two-face identity replay, while the
  written proof carries the universal quantifier; its committed artifact is
  registered and can be hash-audited without symbolic replay or rewriting.
- The diagonal/Meng--Yang all-degree frontend now treats its committed
  degree-4--8/order-1--12 tables explicitly as bounded regressions.  Their
  generating source and artifact are registered and auditable without CAS;
  the artifact's discovery-time lower-layer recommendation is preserved but
  marked superseded by `HC4FSD3`.  Formal all-order Meng--Yang recursion still
  carries no polynomial-termination or marked-collision conclusion.
- The historical HC4 collision-first finite-field campaign is now a first-
  class bounded record `HC4FF1`.  A no-search audit pins five artifacts and
  three sources, verifies all 45,181,194 nonempty one/two-direction choices
  and the exact `96+32+128+144` selected dense families, and requires all 800
  stored support-prime ideals to be unit.  The direct degree-five problem now
  explicitly forbids treating these two-prime, fixed-normalization records as
  a characteristic-zero or unrestricted exclusion.
- Static review of the smooth-quartic reciprocal chain confirmed that the
  polar decomposition is not being used as a tangent/contact classification.
  The generic polar conic and all generic two-line orbit types are closed, but
  every lower denominator, secondary-pivot, complementary-chart, and other
  reciprocal boundary remains explicit.  The four checkers that consume the
  shared 81-equation builders now register those source locks and provide
  provenance-only maintenance modes; no Singular or exact-field calculation
  was rerun.  Corrupted inline mathematics in the Fermat-symmetry note was
  repaired at the same time.
- The scalar relative-pencil narrative was brought forward from its obsolete
  degree-six handoff: its own later theorems close scalar directions through
  degree seven, and `HC4MR1` closes the complete auxiliary relative-nilpotent
  pencil branch in all degrees while leaving unrestricted `HC4` and `JC2`
  untouched.  The final affine-plane checker remains only a local certificate,
  not an end-to-end replay; its committed prolongation artifact is now
  registered and can be checked, including that proof boundary, without
  recomputing or rewriting it.
  The earlier shared `HC4RSD17--28` ledger is also registered and auditable;
  its degree-seven `open_frontier` is preserved as historical stage data and
  explicitly superseded by `HC4RSD40`/`HC4MR1`, rather than silently rewritten.
  Likewise, the `HC4RSD11--16` scalar-dichotomy ledger now identifies its
  higher-degree pencil handoff as historical and registers the committed
  artifact on every consumer without replaying it.
  The same protection now covers the three preceding quadratic-pivot ledgers:
  their rank-one/rank-two and cancellation handoffs remain immutable stage
  evidence, while maintenance output names the later closing theorems instead
  of presenting those fields as current work.
  The affine-pivot pair is now guarded too: `HC4RSD6`'s residual locus remains
  a representation-classification problem, but `HC4RSD7` bypasses it for
  inherited collisions and still leaves nonlinear or different-fiber
  collisions open.
  The five earlier scalar-kernel ledgers are now guarded in the same
  cleanup-only way. `HC4RSD1--5` pin their committed artifacts and shared
  equation helper, while `HC4RSD1--2` also pin the projective-polar atlas they
  consume. Their historical frontier text is interpreted against current
  knowledge: the affine and fixed two-component subcases are closed by later
  `HC4RSD` steps, the nonzero-corner auxiliary pencil is closed by `HC4MR1`,
  and the larger nonlinear singular kernels, nonlinear zero-corner exact
  remainders, and moving matrix pivots remain. No symbolic calculation was
  replayed.
  The twelve `HC4RSD29--40` pure-sextic/septic stage artifacts are now
  individually registered too. A single stdlib-only auditor verifies their
  exact hashes and status mappings without importing SymPy, invoking Singular,
  or regenerating any ledger. This preserves the degree-six/seven closure and
  the later `HC4MR1` auxiliary-pencil endpoint without suggesting a proof of
  unrestricted `HC4`.
  The earlier `HC4NHM1--3` nonreduced-Hessian route map was also reconciled
  statically. Its checkers have no generated ledgers to register, and the
  sole external Singular step already fails closed rather than skipping.
  Historical `3+1` and later split-linear handoffs now point to their
  `HC4NHM4--12` closures, while nonlinear clean denominators, positive defect,
  and lower-Smith strata remain current. No CAS checker was run.
  The `HC4PPG1--9` projective-polar chain now registers its previously implicit
  shared Segre helper and all produced/consumed atlas, Rees, vertex, and
  conditional-sieve ledgers. A stdlib-only auditor verifies those committed
  bytes and requires the intended boundary: 624 necessary numerical rows are
  not existence results, and the higher-torsion and exceptional
  codimension-three strata remain. No symbolic or Macaulay2 calculation was
  replayed.
  The ten `HC4MCP1--10` mixed-canonical ledgers are now first-class status
  inputs as well. Their exact files and directly imported sources are pinned,
  and a dependency-free maintenance auditor preserves the crucial scope
  split: unequal modular values are valid nonconstancy witnesses, equal values
  are retained rather than sieved out, and only `HC4MCP10` is coefficient-
  uniform—and then only on the 54 short-word `HC4MCP6` resonance families.
  Moving pivots, other supports, and longer words remain open. No search,
  SymPy calculation, or Singular calculation was run.
  Static review of the nonlinear clean-denominator predecessors
  `HC4NHM13/15/18/19` found no omitted artifacts or permissive CAS fallback:
  the Singular-backed scripts fail closed and the other two are self-contained.
  Their notes already distinguish the at-most-four-root double-conic closure,
  the unsaturated many-root frontier, and the three surviving smooth-cubic
  degree packets; a stale “next action” heading was relabeled historical. No
  checker was run.
  The committed HC4 cube-torsion/Fitting-denominator and fourth-power support
  experiments had no mathematical-status row. They are now the partial record
  `HC4QSE5`, with a dependency-free auditor for their exact files and source.
  It preserves the chart restriction `nu != 0`, the finite-field and extension-
  field limits, and the fact that every timed-out symbolic route has no
  conclusion. The reconstructed point `(-5/3,-1/6)` is only a nilpotence-order
  jump, not a reduced exceptional Schur component. No Singular or scan replay
  was run.
- The ICARM construction-recognition records were reconciled without replaying
  their searches. The fixed-root v2 ledger includes curve 302, and its old
  repository-model diagnostic was narrowed to uncompressed generated-results
  JSON rather than all repository models. The already-committed 2,334-family
  generated-space ledger is now registered separately: it exactly recognizes
  curve 282 in two coordinates and excludes curves 273 and 302 only inside
  that declared bounded space. A stdlib-only auditor pins the sources and
  artifacts; no family generation, modular sieve, or factorization was run.
- The completed-moment artifacts were reconciled without rerunning invariant
  or Gröbner calculations. Four status rows now pin the automatic, diagonal,
  and single-phase ledgers and their direct helpers. The automatic-ledger note
  had omitted its 32-step beta-regression option and retained an obsolete
  hash; the older bounded relation JSON also came from an earlier committed
  producer rather than the extended current schema. Both provenance errors
  are now explicit. Slice fixed-field theorems, necessary-only Hilbert tests,
  and support-bounded modular nonrelations remain separate claim types.
- The adjacent degree-four moment-field chain now pins the ledgers behind its
  composite scopes: the quartic witness, 22 even parameters, bounded
  weight-16 nonrelations, diagonal and single-phase fields, and phase-one
  four-point chart. The review preserves the essential separations between
  full-rank modular support exclusion and field equality, raw slice degree and
  invariant-quotient degree, and `F_101` fiber completeness versus exact
  rational branch identities. Solver timeouts remain non-results. No SymPy,
  Singular, or `msolve` calculation was rerun.
- The two quartic `q_2` normal-jet rows now pin their stored finite-field
  ledgers. The review keeps the positive-dimensional quadratic/cubic jets,
  dominant-sheet quartic radicals, and unresolved off-axis branch distinct;
  none proves formal isolation, treats `F_2=0`, or lifts the global result to
  characteristic zero. It also records that one 240-second timeout is
  historical metadata serialized by the checker, not a current calculation.
  No SymPy or Singular calculation was rerun.
- The `PWB1--PWB6` wild-boundary chain now pins its two committed result
  ledgers and direct lattice helper. A stdlib-only audit preserves the
  distinction between the empty post-support balanced queue and the 26
  unresolved reconstruction rows in other stored architectures, and keeps
  bounded `packet_gate_only` survivors separate from actual covers. The
  first odd degree-seven scan remains complete only inside its six-row `F_3`
  retained-polynomial family. No symbolic identity, packet enumeration,
  normalization, or point-count calculation was rerun.
- The three all-dimensional projective-gradient Segre ledgers are now
  first-class inputs to `PGS1--PGS3`, together with their shared helper and
  independent Macaulay2 calibration sources. Their Python producers no
  longer overwrite committed evidence on an ordinary verification run;
  replacement is explicit through `--write`. A dependency-free audit checks
  exact parameter-key coverage and preserves the separation between complete
  vectors, aggregate top-degree information, and unresolved singular-profile
  data. No SymPy or Macaulay2 calculation was rerun.
- The seven precursor notes for the relative rank-two Jordan packets now mark
  their “remaining” and “next” sections as historical.  Their exact closure
  chain is explicit: `HC4RSD60` closes square-zero `[2,2]`, `HC4RSD63` closes
  length-three `[3,1]`, and `HC4MR1` subsumes both.  This prevents old local
  handoffs from being mistaken for current open packets.
- The repeated-linear direct-HC4 note no longer advertises the obsolete
  three-packet exact-sextuple handoff.  Its later internal reductions leave one
  degree-five, order-one pure-cube scalar-parent resonance, which `HC4DIR28`
  excludes as a collision source.  The generated identity artifact is now
  registered and can be checked for that endpoint and its written UFD/DVR
  proof boundary without reconstructing or rewriting it.
- The formal-orbit ledger was narrowed to its actual theorem.  Local-Artin
  Keller deformations are uniquely source-trivial, and individual
  determinant-one jets in dimension three have reduced representatives, but
  neither statement algebraizes the compatible formal family or erases
  filtered/global stable moduli; the existing rank-four translation remains
  an exact counter-control.

Durable lesson: exclusion of one quotient class is not a rank theorem, formal
smoothness is not algebraization, incremental CRT height is not a monotone
pruning invariant, and a geometric Neron--Severi lattice must pass a separate
Galois/rational-marking gate before arithmetic equation work.

## Current handoff

The repository-wide continuation queue is
[`STATUS.md`](STATUS.md#active-open-problems). The principal K3 foundry gate is
the 823-row transcendental-first arithmetic queue, beginning with the exact
coarse genus-zero rows at determinants 378, 256, and 512. Each still needs its
literal stable kernel and rational-point decision before `NS=T^perp` is
constructed. The determinant-948 route is the single possible classifier
control; determinant 720, NS0024, and NS0031 are closed over `QQ`. A positive
row must ultimately supply arithmetic MW17 plus the certified carrier section,
or the stronger saturated `V4` character lattice.

The principal rank-32 gate remains a completed residual 2-Selmer quotient on
the same minimal curve, followed only on a passing fibre by cover or point
search.  Parallel exact questions are an unconditional upper bound for curve
302, residual descent for curve 273 and the low-conductor near misses, nonzero
alternate-Q80 product quotient classes, and the unresolved rational visibility
directions.  Smooth and singular-genus-one character collisions are closed on
the two direct norm-twelve charts; genus two remains open beyond the recorded
bounded screen.

## 2026-09-05: GVC finite certificates, HC4 repair, and a JC2 polynomial gap

- [GVC2SC](extended-geometry/BINARY_GVC_FINITE_CERTIFICATE.md) extracts a
  unique Hall direction over the original field, a finite rational-input
  decision procedure, and the mixed cutoff `m > (deg P + deg Lambda) deg Q`.
  The active GVC manuscript includes the written theorem and exact replay.
- [HC4MRA1](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md) corrects the earlier claim
  that the adapted motion determinant is constant: the frozen normalization
  controls `pq/a^2`. A new prolongation excludes the positive sign, while
  the negative sign retains a compatible finite jet. The subsequent
  [HC4MRA2 proof](HC4_NEGATIVE_MOTION_POLYNOMIAL_OBSTRUCTION.md) excludes its
  polynomial realization by restricting N to an affine leaf and comparing
  degrees in `2nn''-3(n')^2=0`. This restores HC4MR1/2 through a global
  polynomial argument. Earlier closure descriptions are superseded by the
  corrected proof; the original notes and conditional certificate survive.
- [PF2D6O1](plane-jc/F2_DEGREE_6_10_POLYNOMIAL_GAP.md) distinguishes a
  numerical infinity semigroup from a polynomial parametrization in the
  prescribed degrees. A literal `b^12` identity proves the sharp odd gap
  bound 21 for degree-`(6,10)` normalization pairs and excludes the F2
  normal terminal row `r=9`. A birational polynomial target attains
  `r=7`; normal rows `5,7` and nonnormal conductor slices remain open.
