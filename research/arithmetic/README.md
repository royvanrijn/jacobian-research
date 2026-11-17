# Arithmetic fibres and discriminant geometry

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Constructive finite etale fibres, arithmetic failures and discriminant geometry.

**Authored assessment:** The fixed quintic has exact arithmetic decision criteria, full leading counts on its two Hasse lines and for distinct original-line algebras, and a complete three-line fixed-quadratic classification. Ambient rational local/global counts have sharp order B^3; exact ambient Hasse growth and full ambient constants remain open. All new assurance is internal written proof and scoped checks.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| Over Characteristic Zero, Every Finite Etale Algebra of Rank at Least Three Is a Full Keller Fiber | Unpublished; independent review pending | proved | [Source](../../unpublished/common-arithmetic-fibers/README.md) |
| Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map | Unpublished; independent review pending | proved | [Source](../../unpublished/fixed-map-hasse-failures/README.md) |
| The Exact Nonproperness Locus of the Quadratic-Gauge Keller Maps | Unpublished; independent review pending | open | [Source](../../unpublished/quadratic-gauge-nonproperness/README.md) |
| Generic Discriminants of Polynomial Tangent Pencils | Unpublished; independent review pending | proved | [Source](../../unpublished/discriminant-pencils/README.md) |
| Exact Real Chamber Spectra of Explicit Keller Maps | Unpublished; independent review pending | proved | [Source](../../unpublished/exact-real-chamber-spectra/README.md) |
| Fixed-map local classification and stronger arithmetic counts | Internal result | proved | [Source](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md) |

**Computed obligations:** 25 unresolved; 3 resolved within their exact scope.


Scope limits:

- The characteristic-zero headline is an ordinary A3 polynomial Keller realization, not merely a localized chart or a positive-characteristic model. Its internal JC-failure consequence is displayed explicitly alongside the cubic route.
- The supplied char!=2 proposition retains its nonzero first/third Hasse-coefficient hypotheses and is not formalized end to end. Characteristic-zero formal scope is separate from analytic counting and classical rank-two exclusion.
- The latest fixed-map manuscript claims complete rational Hasse-failure classification and Cphi B^2/log B counting on one specified target line, while the earlier B/sqrt(log B) theorem counts a sufficient integer subfamily. Neither is an ambient A3-target classification. Specialist analytic/source review and publication workflow are separate obligations.
- The manuscript retains its original-line scope. Canonical research now supplies both formerly open leading constants, the full fixed-quadratic line classification, and ambient local/global order B^3. Modern projective thin-cover and classical larger-sieve inputs have exact retained source interfaces. No formal or external-review assurance propagates along these arrows.

Tracked obligations (computed states):

- **verification / open:** [Review the manuscript against its formal theorem and external inputs.](../../unpublished/common-arithmetic-fibers/README.md) — GAP-ARITHMETIC-1
- **workflow / open:** [Complete editorial and priority review before submission.](../../unpublished/common-arithmetic-fibers/README.md) — GAP-ARITHMETIC-2
- **verification / open:** [Complete the priority and source review recorded in LITERATURE_AUDIT.md.](../../unpublished/fixed-map-hasse-failures/README.md) — GAP-ARITHMETIC-3
- **verification / open:** [Obtain specialist review of the complete rational count before preparing a submission bundle.](../../unpublished/fixed-map-hasse-failures/README.md) — GAP-ARITHMETIC-4
- **verification / open:** [Review the geometric proof and exact symbolic audit.](../../unpublished/quadratic-gauge-nonproperness/README.md) — GAP-ARITHMETIC-5
- **workflow / open:** [Decide whether to develop this companion into a submission.](../../unpublished/quadratic-gauge-nonproperness/README.md) — GAP-ARITHMETIC-6
- **verification / open:** [Reassess the scope and complete the proof review before restarting publication work.](../../unpublished/discriminant-pencils/README.md) — GAP-ARITHMETIC-7
- **verification / open:** [Reassess the manuscript against the canonical fibre-spectrum proof and current literature.](../../unpublished/exact-real-chamber-spectra/README.md) — GAP-ARITHMETIC-8
- **verification / open:** [Obtain specialist review of the exact rational target-line classification/count and its Selberg-Delange interface.](../../unpublished/fixed-map-hasse-failures/metadata.json) — VERIFY-RATIONAL-LINE-COUNT
- **mathematical / open:** [Determine the sharp order or leading asymptotic of the ambient rational Hasse-target count A(B). Current proven bounds are Cphi B^2/logB from below and O(B^3) from above; all-target local solubility has sharp order B^3.](FIXED_MAP_AMBIENT_COUNT_BOUNDS.md) — GAP-ARITH-AMBIENT-TARGET-COUNT
- **mathematical / open:** [Refine the sharp-order results G(B),L(B)=Theta(B^3) in actual rational target height to leading asymptotics; the boundary-plane contribution alone need not equal either full count.](FIXED_MAP_AMBIENT_COUNT_BOUNDS.md) — GAP-ARITH-AMBIENT-LOCAL-GLOBAL-CONSTANTS
- **mathematical / open:** [Construct an integral target of fixed Phi with a Hasse-failing fibre or prove none exists. All three fixed-quadratic affine lines avoid Z^3, and the larger-sieve bound is only an upper bound.](FIXED_QUADRATIC_LINE_RIGIDITY.md) — GAP-ARITH-INTEGER-HASSE-EXISTENCE

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
