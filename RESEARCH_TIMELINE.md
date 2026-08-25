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

- The first two edges of that suffix now have exact characteristic-zero
  Riemann--Roch planes, quartics, minimal Jacobians, effective zeros, and
  equation markings.  The `q8/orbit376` horizontal is also exact over
  `QQ(t)`; its resolved pencil and child equation remain open.  The preferred
  final `q12/orbit5867` equation remains open, with `q12/orbit4484` retained
  as a certified fallback.  Both marked lattice transports and the endpoint
  isometry are exact.
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

## Current handoff

The repository-wide continuation queue is
[`STATUS.md`](STATUS.md#active-open-problems).  The nearest concrete K3 task is
the resolved `q8/orbit376` equation lift, followed by preferred
`q12/orbit5867` (or fallback `q12/orbit4484`) and direct R17 endpoint
certification.  The nearest elliptic arithmetic tasks are an unconditional
upper bound for curve 302 or a rank-at-least-32 example, exact repository-local
conductors for curves 285/286, and the residual 2-Selmer closure for curve 273.
