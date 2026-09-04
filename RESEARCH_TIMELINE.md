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
  [`STATUS_HASH_REFRESH_2026-08-25.md`](STATUS_HASH_REFRESH_2026-08-25.md).

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

Durable lesson: exclusion of one quotient class is not a rank theorem, formal
smoothness is not algebraization, and a geometric Neron--Severi lattice must
pass a separate Galois/rational-marking gate before arithmetic equation work.

## Current handoff

The repository-wide continuation queue is
[`STATUS.md`](STATUS.md#active-open-problems). The principal K3 foundry gate is
an arithmetic-first reranking of the remaining different-NS frames, followed
only for a lattice whose full rational rank-19 marking survives by target-free
marked-`U` selection and exact endpoint compilation. The determinant-948 route
is a complete control; NS0024 and NS0031 arithmetic MW17 are both closed over
`QQ`. Determinant 720 is the strongest remaining lattice/corridor control, but
its known rational `3A5` point saturates to determinant 20 and is not a valid
source.

The principal rank-32 gate remains a completed residual 2-Selmer quotient on
the same minimal curve, followed only on a passing fibre by cover or point
search.  Parallel exact questions are an unconditional upper bound for curve
302, residual descent for curve 273 and the low-conductor near misses, nonzero
alternate-Q80 product quotient classes, and the unresolved rational visibility
directions.  Smooth and singular-genus-one character collisions are closed on
the two direct norm-twelve charts; genus two remains open beyond the recorded
bounded screen.
