# Research knowledge base

This page is the repository-wide guide to reusable methods, failed methods,
and proof discipline.  It is a synthesis, not a second mathematical-status
registry.  For the current state of every named result, use
[`MATH_STATUS.json`](MATH_STATUS.json) or its generated view
[`STATUS.md`](STATUS.md).  For the order in which the programmes developed,
use [`RESEARCH_TIMELINE.md`](RESEARCH_TIMELINE.md).

Snapshot date: **2026-09-04**.

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
| Generic-normalized Kummer class pressure | Published-R17 and alternate-Q80 rank-jump controls | Compare the ideal-square-root image forced by known point classes only after quotienting the specialized generic MW17 contribution.  On the six certified controls the lower-bound strata `3,5,6,10..11` strictly follow the known jumps `5,6,8,12`, while the unnormalized full `Cl[2]` bounds do not.  This explains descent pressure but cannot predict a jump prospectively because the exceptional points are inputs. |
<!-- status-consumer: EC-K3-R17-KUMMER-CLASSGROUP-PRESSURE-COMPARISON 74b1dae24470b531 -->
| Neron--Severi lattice navigation plus marked equation lifts | Elkies--K3 programme | Search primitive embedded hyperbolic planes cheaply in the lattice, but carry the full marking and then realize a selected divisor by an exact Riemann--Roch/quartic/Jacobian calculation.  ADE/MW labels alone do not identify a state. |
| Galois-equivariant arithmetic marking gates | Different-NS K3 foundry and rational-surface controls | A geometric Picard-rank-19 lattice or MW17 frame is not an arithmetic rank-17 fibration.  Carry the Galois action on the trivial lattice, MW lattice, and discriminant glue, and require the full rational marking before equation work. |
| Full-marking classification against coarse arithmetic curves | Rank-19 different-NS reranking | Distinguish the coarse even-Clifford/norm-one curve from the stable discriminant-kernel marking curve, and pin every classification decision to exact lattice, order, and lift data.  Of 66 exact candidates, one existing control remains possible, three are excluded, and 62 stay `UNKNOWN`; only a proved full-marking pass may enter equation compilation. |
| Transcendental-first arithmetic foundry | Lane B K3 construction | Start from the full rank-three `T` ledger, compute the literal stable marked curve and a rational non-CM point, then construct and saturate `NS=T^perp`; inspect rootless frames and equations only after that pass. Coarse genus is a priority diagnostic, not a marking certificate. |
| Prescribed character blocks and graph glue | Low-genus K3 carriers and `V4` rank transfer | A stronger construction target fixes invariant rank 17, one rank-one block for each nontrivial `V4` character, and every 2-primary half-sum before equation search. Two characters and a positive-rank base do not imply the product-character section or saturated rank 20. |
| Target-free planning plus universal degree-two compilation | Direct alternate-Q80 and noncyclic R17 bridges | Select a primitive marked `U` from an intrinsic predicate, minimize its equation incidence, and compile it with one exact chord/quartic/Jacobian interface.  Keep planner reachability, compiler correctness, and arithmetic realization as separate gates. |
| Exact finite-quotient slicing | Alternate-Q80 product twists and integral rank transfer | Derive the invariant/anti-invariant or discriminant quotient first, enumerate its exact classes and minimum representatives, and solve class-sliced systems.  Excluding one class never excludes the others. |
| Coverage-certified sharded enumeration | Rootless bisections, prescribed-root foundry censuses, and singular rational-normalization searches | A distributed result needs exact expected key sets or half-open domain ranges, source hashes, gap/overlap checks, exact exception accounting, and a fail-closed merged certificate.  Aggregate counts and individual no-hit chunks do not carry theorem scope. |
| Formal local certification with an explicit algebraization boundary | NS0031 marked source | A unit Jacobian minor can prove a compatible formally smooth `ZZ_p` branch.  Algebraization, a rational characteristic-zero point, Picard control, and a rational NS marking remain separate obligations. |
| Minimal projections of mutable public data | ICARM lineage and native rank-jump calibration | A live-URL hash detects database drift but does not preserve historical bytes.  Pin the exact claim-relevant records, source count, selection rule, null fields, and projection hash; treat later fetches as new snapshots. |
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
| Analytic scores, discriminant radicals, or Selmer dimensions as elliptic rank/conductor claims | Nagao scores are heuristics, raw radicals are not exact conductors, and residual Selmer classes may lie in Tate--Shafarevich groups.  In PARI's four-field `ellrank` output, the third field is the Cassels-pairing quotient rank `dim(Sha[2]/2Sha[4])`, not all of `Sha[2]`. | Re-minimize, run exact local conductor analysis, verify points, and prove independence; derive `dim Sel_2` from the documented backend fields and label upper bounds and hypotheses separately. |
| Treating target-fitted visibility as a predictor | Low-genus pencils could be made to contain known exceptional directions, but only after those targets were supplied. | Label post-hoc visibility separately, freeze any predictor before inspection, and require held-out exact controls before it may rank new candidates. |
| Pruning incremental CRT states by their current rational height | Exact four-prime regression gives height 1409 for the width-one survivor while a discarded branch later completes to `48/53`, of height 53. Shortest-representative height is not monotone as congruences are added. | Treat finite beams as proposal/ranking devices only. A completeness claim needs the full Cartesian product or a separately proved lossless reduction; never interpret an omitted state as an excluded solution. |
| Using enriched lattice scores as hard sieves | “Maximize bridge minimum” retained only four of five rootless fixed-core controls, then had zero ranking power on a held-out 277-candidate shell while costing more than direct root classification. | A heuristic may order work but may not discard candidates.  Use only mathematically monotone lower bounds for rejection, and distinguish fixed-core exact theta convolution from prospective core generation. |
| Promoting geometric MW rank to arithmetic MW rank | Exact NS0024 frames and same-NS corridors did not imply that all nineteen divisor classes were rational over `QQ`. | Compute the Galois action and require a full rational marking.  The NS0024 marking obstruction parks the arithmetic route while preserving its geometric evidence. |
| Promoting a modular or formal lift to a rational source | A smooth GF(7) point and compatible Hensel data for NS0031 did not give an algebraic `QQ` source; the later split-Clifford/`X_0(37)` obstruction proves that the required full rational marking cannot exist. | Certify the formal branch at its exact scope, but run the arithmetic moduli and rational-marking gate before algebraization or coefficient search. Preserve local/geometric branches after a rational obstruction instead of treating them as failed computations. |
| Treating a coarse norm-one curve as the full marking curve | Low genus or a rational point on the even-Clifford/norm-one quotient can still fail the stable discriminant-kernel lift required to make all nineteen divisor classes rational. | Compute the full marking group and exact lift obstruction.  Keep unresolved candidates `UNKNOWN`, and emit no equation-work handoff until a full rational marking is proved possible. |
| Letting a ranked prefix inherit exhaustive scope | The cheapest 1,024 alternate-Q80 bisections are a useful laboratory, but cost optimization does not cover the remaining native classes. | Keep complete inherited/smooth-character enumerations and cost-ranked native prefixes separately typed in data, prose, and status. |
| Excluding one quotient class as if it excluded all solutions | Norm-eight inversion initially covered only part of the possible zero-Tate-class carriers, and even the completed zero-class exclusion says nothing about nonzero classes. | Derive the full minimum-norm quotient spectrum, close every carrier of the named class, and list all surviving classes and existence questions explicitly. |
| Universal cubic gradedness and similar overstrong extrapolations | Explicit counterexamples show that useful minimal-boundary classifications do not extend to every cubic representative. | Retain the valid boundary-minimal theorem and record the universal claim as falsified rather than silently narrowing its statement. |
| Reopening finite cubic tensor sweeps after an all-orders closure | Early coordinate-axis and parameter-plane audits left the six singular squarefree quartic families apparently open, but `KDSQ6` closes every quartic nongauge complement and `SSADPALL` preserves cotangent saturation and the six-generated non-Cartier different through every compatible formal tail. | Treat the squarefree formal-tail queue as closed.  Work instead on the missing boundary-geometric `S_2`/local-CM or Cartier implication, the non-squarefree leading symbols, and global Keller compatibility. |
| Treating moment--nullcone equality and SIC safety as the same question | The full-rank bidegree-`(3,3)` Rodrigues point is semistable with every pure moment zero, so it falsifies moment--nullcone equality, yet its all-order integration-by-parts cutoff proves that entire orbit SIC-safe. | Record the equality as falsified and keep only the complete bidegree-`(3,3)` SIC classification open.  A semistable pure-zero point is not automatically a mixed-moment counterexample. |
| Certifying a support census from totals or a set complement | A duplicated processed support and an omitted unprocessed support can preserve aggregate totals; calling the complement of the last batch “earlier closed” does not prove that it equals the union of earlier theorem families. | Reconstruct the exact ordered support universe, require every record key exactly once, recompute outcome counts, and compare the explicit union of predecessor families with the final complement.  Offer this as a cheap audit independent of the expensive elimination replay. |

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
| Elliptic curves over `Q` | Curve 302 has unconditional rank at least 31 with exact local arithmetic; curves 273 and 356 have rank at least 30 and 29; the complete pinned norm-twelve atlas and seven exact native quotient audits supply published-R17 and alternate-Q80 calibration fibres | Exact rank 31 or rank at least 32; completed residual descent for the high-rank and low-conductor targets; saturated section transports for the 57 still-unknown recognized-fibre quotients; sharper exact conductor/rank records. |
| Elkies--K3 reconstruction and foundry | Published R17, direct alternate Q80, and the noncyclic `4A1/MW13` bridge are equation-explicit over `QQ` with saturated bases.  Smooth and singular-genus-one character maps are exhausted on the two direct norm-twelve charts, the pinned 43-chart/474-curve atlas gives 69 native controls, and seven priority fibres have exact quotient/visibility audits. | The next construction milestone is arithmetic MW17 on a different NS.  The exact classifier leaves 1 existing control possible, excludes NS0024 and NS0031, and keeps 63 candidates `UNKNOWN`; no equation handoff exists until one of those unknown rows passes the full rational-marking gate. Rank `>=32`, 57 native quotient rows, nonzero product quotient classes, and global genus-two injectivity remain open. |

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
- [`elkies-k3/DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md`](elkies-k3/DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md)
  for the current different-NS arithmetic construction gate;
- [`elkies-k3/ELKIES_K3_PROCESS_ATLAS.md`](elkies-k3/ELKIES_K3_PROCESS_ATLAS.md)
  for the detailed K3 chronology and method ledger.
