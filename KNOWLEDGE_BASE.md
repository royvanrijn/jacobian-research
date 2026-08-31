# Research knowledge base

This page is the repository-wide guide to reusable methods, failed methods,
and proof discipline.  It is a synthesis, not a second mathematical-status
registry.  For the current state of every named result, use
[`MATH_STATUS.json`](MATH_STATUS.json) or its generated view
[`STATUS.md`](STATUS.md).  For the order in which the programmes developed,
use [`RESEARCH_TIMELINE.md`](RESEARCH_TIMELINE.md).

Snapshot date: **2026-08-31**.

## Evidence ladder

The repository repeatedly succeeded by separating five questions that are
easy to conflate:

1. **Was a finite search completed?**  This gives a bounded experiment only.
2. **Does an exact candidate satisfy its identities?**  This gives an exact
   symbolic or modular computation in the stated model.
3. **Does the computation imply the mathematical claim?**  This needs the
   written hypotheses, reduction, and proof boundary in a canonical note.
4. **Can another implementation replay it?**  Independent replay is recorded
   separately from proof type.
5. **Was it formalized or externally reviewed?**  These are additional
   assurance signals, never inferred from a passing checker.

The most common historical error was to skip from step 1 or 2 to step 3.
`MATH_STATUS.json` prevents that by typing every entry as proved, partial,
open, parked, archived, or falsified and by recording its proof type,
dependencies, checker, and software lock.

## Methods that worked

| Method | Where it paid off | Reusable lesson |
|---|---|---|
| Marked-root reconstruction and exact inverse identities | Foundational, weighted, quadratic-gauge, finite-etale-fibre, and arithmetic constructions | Keep the finite normalization, its affine reconstruction open, and the omitted boundary distinct.  A primitive root presentation is computational data, not the invariant itself. |
| Valuation ledgers and controlled-boundary cancellation | Weighted and cancellation Keller maps; nonproperness and stable multiplicity | Prove polynomiality and the full source--core--target determinant balance together.  A formal cancellation alone does not produce a polynomial map. |
| Saturation before elimination | Collision ideals, Fitting modules, resolved Riemann--Roch modules, conductor and Cox calculations | Remove boundary or wrong-thickness components before interpreting rank, length, or membership.  Unsaturated equations routinely created false obstructions. |
| Exact modular discovery followed by characteristic-zero certification | Keller searches, GVC/SIC calculations, elliptic curves, and Elkies--K3 neighbours | Use finite fields to discover sparse structure, then lift or reconstruct and verify literal identities over the intended field.  A good-prime result is evidence until that final step. |
| Independent low-dependency replay | Foundational map, proof-carrying arithmetic compiler, elliptic-curve independence, and several certificates | Freeze a compact certificate and rederive its decisive identities with a genuinely separate arithmetic path where risk warrants it. |
| Finite-quotient independence | High-rank elliptic-curve specializations | Verify points exactly, prove trivial relevant torsion, and obtain full column rank in products of good-reduction quotients.  This proves a lower bound, not an upper bound. |
| Positive-control calibration plus fail-closed descent gates | Elkies rank-25--28 fibres and the rank-32 search | Require a heuristic to rank every exact positive control across disjoint prime blocks, then require the actual residual Selmer quotient on the same minimal curve before any expensive point search.  A timeout, signature, or incomplete local-class list authorizes nothing. |
| Neron--Severi lattice navigation plus marked equation lifts | Elkies--K3 programme | Search primitive embedded hyperbolic planes cheaply in the lattice, but carry the full marking and then realize a selected divisor by an exact Riemann--Roch/quartic/Jacobian calculation.  ADE/MW labels alone do not identify a state. |
| Exact group-law decomposition of difficult sections | Q80 terminal closeout and H3 A11-to-2A5 lift | Recover a hard target from a small exact section shell and modular marking rather than solving one large nonlinear system from scratch. |
| Recurrence and branch-value interpolation | H3 A11-to-2A5 and q4/orbit164 equation lifts | Look for forced boundary coefficients, rational branch values, and denominator-square cancellation before nonlinear elimination; these reduced large-looking RR problems to exact linear calculations. |
| Local-to-global compilers | Arithmetic fibres, boundary packages, Ritt strata, and plane-JC work | Make local obligations explicit and emit proof-carrying global data.  This exposes precisely which local checks are theorem inputs and which remain experiments. |
| Deformation, cotangent, and Postnikov compression | Ritt/Hurwitz and formal-boundary calculations | Replace rapidly growing global Gröbner systems by finite obstruction modules or cellular stages when the geometry supplies compatible filtrations. |
| Formalization of stable algebraic layers | GMC(2), support saturation, and finite-etale Keller fibres | Formalize reusable algebraic implications and explicit certificates.  Do not describe an entire headline theorem as formal when arithmetic or geometric inputs remain outside Lean. |

## Methods that failed or were too coarse

These are route exclusions or cautions.  Unless `MATH_STATUS.json` says
otherwise, they are not general impossibility theorems.

| Attempt | What went wrong | Replacement |
|---|---|---|
| Large direct Gröbner solves before structural reduction | K3 section systems with dozens or hundreds of variables timed out; several GVC/SIC and deformation systems grew past a useful certificate boundary. | Compress first by lattice transport, group law, symmetry, saturation, recurrences, or a finite obstruction module.  Record timeout as a bounded experiment. |
| Treating a bounded miss as a lower bound, upper bound, or classification | Sparse-support, section-height, modular, Hensel, and candidate scans only cover their declared box. | State the box, seed, primes, and stopping rule.  Promote only an exact surviving candidate or a separately proved exhaustive reduction. |
| Using ADE/MW type as a fibration identifier | Distinct K3 markings with the same type led to invalid section transports and route comparisons. | Store `F`, `O`, physical components, horizontal classes, chamber, and forward/inverse integral NS transports. |
| Counting blow-up centres as independent components or conditions | The H3 q24 resolved calculation overcounted local Riemann--Roch conditions. | Build the physical component graph and saturated local module, then impose quotient conditions on connected exceptional geometry. |
| Reusing a chamber pseudo-zero as an equation-effective section | K3 q4/q6 route scores 4,199 and 10,334, and the q104 comparator 13,518, were invalidated by physical intersection tests. | Identify the actual rational zero on the equation before promoting a changed-zero cost or continuation. |
| Losing a cover degree or denominator | A duplicated binary-quartic 2-cover and an omitted `Dx` factor inflated the H3 q8 collision degree from 10 to 46. | Track every cover/isogeny degree and clear the complete rational expression before local or CRT normalization. |
| Unsaturated moment/Fitting ideals | They retained the wrong marked-collision thickness and suggested false scheme structure. | Saturate by the geometrically named boundary ideal and verify both retained nilpotents and removed components. |
| Ordinary tangent-space quotients for stable moduli | The unrestricted infinitesimal left--right quotient vanishes even where stable parameters are known. | Use filtered, boundary-decorated, or formal-to-algebraic invariants; tangent dimension alone does not see the moduli. |
| Coarse conductor or class-group accounting as a plane-JC obstruction | Local cusp/connector normalizations and conductor classes remained compatible, and abstract endpoint matchings were not finite. | Recover actual global carrier equations, braid/meridian data, marked valuations, and finite-cover bounds. |
| Analytic scores, discriminant radicals, or Selmer dimensions as elliptic rank/conductor claims | Nagao scores are heuristics, raw radicals are not exact conductors, and residual Selmer classes may lie in Tate--Shafarevich groups. | Re-minimize, run exact local conductor analysis, verify points, and prove independence; label upper bounds and hypotheses separately. |
| Universal cubic gradedness and similar overstrong extrapolations | Explicit counterexamples show that useful minimal-boundary classifications do not extend to every cubic representative. | Retain the valid boundary-minimal theorem and record the universal claim as falsified rather than silently narrowing its statement. |

The complete machine-readable replacement and falsification graph is in
[`STATUS.md`](STATUS.md#falsified-claims).  In particular, the unmodified
Hurwitz--LL quotient, the all-degree two-pair moment--nullcone equality, a
uniform inverse-degree-nine bound, universal cubic gradedness, and two old K3
IV-star markings are explicitly falsified.  Sixteen further K3 degree-46
diagnostics are archived because they used the wrong q8 normalization.

## Programme map and present boundary

| Programme | Durable result | Current boundary |
|---|---|---|
| Verified Keller core | Explicit characteristic-zero threefold collision; marked-root, weighted, finite-etale-fibre, monodromy, atomicity, arithmetic, and stable-multiplicity theorems | Presentation-free stable descent, global target/degree minima, and the full stable-moduli object remain open. |
| Gaussian/GVC/SIC/Image | GMC(2) proved; GMC fails from dimension 3; GVC has its exact dimensional phase; two-pair Image/Mathieu failures and extensive SIC frontiers are explicit | Minimum GVC(3) degree, efficient polarization, and complete bidegree-(3,3) SIC classification remain open. |
| Cancellation and stable boundary | Controlled suspensions, cubic normalization packages, and many intrinsic boundary invariants are exact | Extract a canonical minimal-boundary package or prove the remaining point-flatness/curvilinearity statements. |
| Plane Jacobian programme | The audited `(72,108)` row is closed and large parts of the logarithmic/conductor compiler are exact | The active geometric frontier is parked pending actual global branch-carrier and compactification data; no proof of `JC(2)` is claimed. |
| HC4/DC2/Hessian programme | Large scalar, rank-two, repeated-factor, Schur, Wronskian, and polar packets are closed | Mixed/coisotropic Schur descent and the direct degree-five residual packets remain open; none of this is an unrestricted `HC(4)` proof. |
| Elliptic curves over `Q` | Curve 302 has unconditional rank at least 31 with exact local arithmetic; curves 285/286 have fully replayed rank-at-least-21 and sub-cutoff conductor certificates | Exact rank 31 or rank at least 32; residual descent for the four rank-19/20 conductor near misses and curve 273; sharper exact conductor/rank records. |
| Elkies--K3 reconstruction | The complete physical H3 equation route reaches the rootless `24I1/MW17` endpoint, which has exact source identity, geometric Picard rank 19, and full saturated determinant-948 lattice R17; the compact published chart and all seventeen sections are identified exactly | Compute genuine residual 2-Selmer quotients before rank-32 point search. The complete 39,120-class rootless-bisection map is injective on quadratic squareclasses and gives that many generic-rank-at-least-18 covers, so this mechanism cannot produce a rank-two collision or generic rank 19. |

## Working rules distilled from the repository

- Read the canonical note and status entry before reusing a result.
- Preserve the distinction between a search frontier and a mathematical
  frontier.
- Keep generated outputs reproducible: command, parameters, versions, seed,
  and whole-file hash.
- Prefer the smallest exact intermediate record that is sufficient for the
  next construction step; require a stronger certificate at a final endpoint.
- When a route fails, preserve the reason and proof boundary.  Negative work
  is valuable only when later work can tell exactly what it excluded.
- When a stronger result appears, update the canonical source and
  `MATH_STATUS.json`, regenerate `STATUS.md`, and then repair summaries and
  consumer markers.  Do not maintain a second status queue in a narrative
  document.

## Where to continue

The active and parked queues are generated in [`STATUS.md`](STATUS.md#active-open-problems).
More focused entry points are:

- [`UNIFYING_THESIS.md`](UNIFYING_THESIS.md) for the conceptual Keller-map
  spine;
- [`cancellation/RESEARCH_ROADMAP.md`](cancellation/RESEARCH_ROADMAP.md) for
  minimal-boundary continuation;
- [`extended-geometry/README.md`](extended-geometry/README.md) for GVC, SIC,
  Ritt, moduli, and boundary programmes;
- [`plane-jc/README.md`](plane-jc/README.md) for the plane-JC compiler;
- [`elliptic-curves/README.md`](elliptic-curves/README.md) for arithmetic
  rank and conductor;
- [`elkies-k3/ELKIES_K3_PROCESS_ATLAS.md`](elkies-k3/ELKIES_K3_PROCESS_ATLAS.md)
  for the detailed K3 chronology and method ledger.
