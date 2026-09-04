# Elliptic K3 / high-rank programme — ACTIVE

The programme is open for theorem-directed breakthrough work. Expensive K3,
neighbour, descent, and specialization calculations require explicit gates and
reproducible, checkpointed outputs. `../MATH_STATUS.json` remains the authority
for exact claim status.

## Current milestone

- **Determinant-948 closure:** the explicit noncyclic chain `published R17 ->
  4A1/MW13 -> published R17` now has thirteen saturated rational sections and
  target-free reverse selection.  Further equations on this same surface are
  supporting work, not the next foundry milestone.
- **Next foundry milestone (OPEN):** construct the first arithmetic MW17
  fibration on a different Neron--Severi lattice that is selected target-free
  by the planner. Determinant-950 `NS0024` and determinant-1184 `NS0031` are
  both excluded because neither can carry the required full rational rank-19
  marking over `QQ`. The planner must rerank the remaining frames through the
  arithmetic marking gate before more equation work; determinant `720` is the
  strongest lattice/corridor control, but its known rational model saturates
  to determinant `20`. See
  [`DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md`](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md).
- **Published R17 over `QQ`:** rootless `24 I1`, Picard rank 19, saturated determinant-948 Mordell–Weil lattice of rank 17.
- **Alternate Q80 over `QQ`:** the canonical equation route is now the direct degree-two hop `norm12-orbit-11952` from published R17. It gives a polynomial K3 model with `(deg A,deg B,deg Delta)=(8,12,24)`, `24 I1`, the alternate determinant-948 rootless frame, and 17 saturated rational sections.
- **Complete norm-twelve public-curve atlas:** all 43 shared-zero degree-two charts form six rational-`PGL2` `j`-classes, and all 474 equations in the pinned ICARM snapshot have exact preimage decisions. There are 69 rational hits and 2,775 misses; all 376 native chart/fibre comparisons are untwisted. Every wgxli component is recognized. Curve 12 is the first native alternate-Q80 rank-at-least-29 control, with exact displayed quotient `Z^12`; curves 273, 302, and 398 miss all six classes.
- **Native ICARM calibration:** seven priority fibres now have exact native displayed quotients and exhaustive fixed-cover visibility spans: curve 12 has `0/12`, curve 395 has `2/11`, curves 363/364/378 have `2/10`, `1/11`, and `6/7`, and curves 393/404 have `2/9` and `1/10`. Curves 12, 395, 363, 364, and 378 also have exact fitted norm-eight genus-one signatures for every quotient-basis direction (51 directions total). All 69 recognized fibres have pinned local/Nagao feature rows; twelve currently have exact quotient labels, while the other 57 remain `UNKNOWN` pending saturated chart transports. See [`R17_NATIVE_ICARM_CALIBRATION_AUDIT_2026-09-04.md`](R17_NATIVE_ICARM_CALIBRATION_AUDIT_2026-09-04.md).
<!-- status-consumer: EC-K3-R17-NORM12-NATIVE-ICARM-CALIBRATION-AUDIT 1b09c81c025e5fc3 -->
- **Noncyclic `4A1/MW13` over `QQ`:** the relative-`U` witness with maximal bridge `Z/4+Z/8` now compiles directly from published R17 to an explicit `4 I2 + 16 I1` equation.  Thirteen rational sections form a saturated basis, and the target-free reverse degree-two hop recovers the literal published R17 equation.
- The alternate pencil has arithmetic generic rank 17 over `QQ`.
- The historical degree-11511 Q80 transport and million-bit third-`q12` reconstruction are superseded operationally and retained only for provenance.
- The four published R17 rank-25--28 controls are **not** rational fibres of the alternate chart: their exact alternate `j`-preimage polynomials have no rational roots.  Native calibration is now supplied independently by curve 12 (rank at least 29) and curves 363, 364, 378, and 395 (ranks at least 24--28).
- **Alternate arithmetic laboratory:** all 121 inherited quadratic covers and their 7,260 pair products are certified.  The cheapest 1,024 native bisections are compiled exactly; together with the inherited covers they give 1,143 distinct rational conics, pairwise-disjoint branch divisors, and no catalogued three-character `V4` closure.  Exact curve intersections give 10,362 intersection-one pairs in the native prefix and 64 promoted rational genus-one `V4` bases; seventeen of their Jacobians have certified rank one.  Product-character inversion exhausts all 63,917 norm-eight traces and the 49 exact norm-twelve residual parities, with no hit in any of the final 833 trace/target cases.  Thus a compatible height-eight section cannot have zero Tate class.  The quotients, their nonzero classes, and the existence of a height-eight product section remain `UNKNOWN`.
<!-- status-consumer: EC-K3-R17-NORM12-11952-PRODUCT-ZERO-TATE-CLASS-EXCLUSION 9e1c09d47fcf0bde -->
- **Rational-normalization boundary:** singular arithmetic-genus-one pencils are now exhausted on both alternate Q80 (63,917 norm-eight classes) and hidden `103b2` (63,925), with no nonsplit rational normalization.  Thus the branch-character map remains injective through that genus-one row.  The two-node genus-two multi-prime and CRT searches also miss on both charts, but remain bounded; global genus-two injectivity is `UNKNOWN`.
- Raw visibility degree has been replaced by a leakage-aware visibility-complexity profile. Four of the 38 named exceptional generators have individual rigid rational-bisection witnesses; all 38 degree-two genus-one witnesses are target-fitted incidence constructions, not predictors.
- Rank `>=32` remains open.
- The target-directed fibration-hopping foundry has one curated planner-ready end-to-end positive control. The 936 bulk source/target rows remain unready and must not be launched as blind searches.

## Canonical current notes

- [`DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md`](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md) — open milestone, arithmetic-first reranking gate, and target-free acceptance contract.
- [`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md) — Fricke-quotient theorem excluding a full rational NS0024 marking and hence arithmetic NS0024/MW17 over `QQ(t)`.
- [`NS0024_DIRECT_QQ_INOSE_OBSTRUCTION_2026-09-04.md`](NS0024_DIRECT_QQ_INOSE_OBSTRUCTION_2026-09-04.md) — narrower direct degree-475 Inose-source obstruction retained as a supporting corollary.
- [`NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md) — split-Clifford and `X_0(37)` theorem excluding a full rational NS0031 marking.
- [`NS0031_MARKED_FORMAL_BRANCH_2026-09-04.md`](NS0031_MARKED_FORMAL_BRANCH_2026-09-04.md) — retained exact one-parameter formally smooth `ZZ_7` marked source branch.
<!-- status-consumer: EC-K3-NS0024-DIRECT-QQ-INOSE-OBSTRUCTION e87afc1b3529a07f -->
<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-SOURCE-PRECURSOR 2e115b35c30a8cea -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->
- [`R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md`](R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md) — direct noncyclic `4A1/MW13` equation, saturated arithmetic MW13 basis, and target-free reverse hop.
- [`R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md`](R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md) — all 43 equations, six exact `j`-classes, record-curve misses, and the exact five-fibre R17 lineage.
- [`R17_NATIVE_ICARM_CALIBRATION_AUDIT_2026-09-04.md`](R17_NATIVE_ICARM_CALIBRATION_AUDIT_2026-09-04.md) — seven native quotient/cover audits, 51 norm-eight signatures on five fibres, and the fail-closed 69-fibre calibration table.
- [`R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md`](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md) — direct alternate-Q80 equation and saturated MW17 basis.
- [`R17_ALTERNATE_Q80_ARITHMETIC_LABORATORY_2026-09-03.md`](R17_ALTERNATE_Q80_ARITHMETIC_LABORATORY_2026-09-03.md) — 121 inherited covers, 7,260 products, and the exact cheapest-1,024 native branch-incidence laboratory.
- [`R17_ALTERNATE_Q80_V4_PRODUCT_TWIST_LABORATORY_2026-09-03.md`](R17_ALTERNATE_Q80_V4_PRODUCT_TWIST_LABORATORY_2026-09-03.md) — 64 exact rational genus-one `V4` bases and seventeen exact rank-one base Jacobians.
- [`R17_ALTERNATE_Q80_PRODUCT_BISECTION_INVERSION_2026-09-03.md`](R17_ALTERNATE_Q80_PRODUCT_BISECTION_INVERSION_2026-09-03.md) — precise integral dictionary and complete 63,917-class norm-eight product-character inversion.
- [`R17_PRODUCT_TATE_COHOMOLOGY_REDUCTION_2026-09-04.md`](R17_PRODUCT_TATE_COHOMOLOGY_REDUCTION_2026-09-04.md) — exact character-glue/Kummer quotient, complete zero-class exclusion, class-sliced equations, and the remaining full-lattice/2-Selmer gate.
- [`R17_NORM12_RATIONAL_NORMALIZATION_BOUNDARY_2026-09-04.md`](R17_NORM12_RATIONAL_NORMALIZATION_BOUNDARY_2026-09-04.md) — complete singular genus-one exclusion on both direct charts and the bounded two-node genus-two miss.
- [`R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md`](R17_ALTERNATE_Q80_ARITHMETIC_RANK_2026-09-03.md) — arithmetic generic rank 17.
- [`R17_VISIBILITY_COMPLEXITY_2026-09-03.md`](R17_VISIBILITY_COMPLEXITY_2026-09-03.md) — rigid, predeclared-pencil, and post-hoc visibility at the rank-25--28 controls.
- [`J2_GEOMETRIC_ACCESSIBILITY_2026-09-03.md`](J2_GEOMETRIC_ACCESSIBILITY_2026-09-03.md) — degree-two accessibility between the two determinant-948 rootless classes.
- [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md) — canonical proof source for target-directed fibration hopping / integral-rank-transfer machinery.
- [`LITERATURE_AND_NOVELTY_MAP_2026-09-03.md`](LITERATURE_AND_NOVELTY_MAP_2026-09-03.md) — literature provenance and conservative novelty map.
- [`INVERSE_ADE_TARGET_PLANNER_2026-09-03.md`](INVERSE_ADE_TARGET_PLANNER_2026-09-03.md) and [`SINGLE_PLANNER_READY_FOUNDRY_ROUTE_2026-09-03.md`](SINGLE_PLANNER_READY_FOUNDRY_ROUTE_2026-09-03.md) — inverse planner and one end-to-end control.
- [`INVERSE_ADE_ADAPTIVE_BACKEND_2026-09-03.md`](INVERSE_ADE_ADAPTIVE_BACKEND_2026-09-03.md) — target-free expanded/orbit/lazy birth backend, terminal controls, and one new NS0024 rootless completion.
- [`MARKED_U_REALIZATION_PLANNER_2026-09-03.md`](MARKED_U_REALIZATION_PLANNER_2026-09-03.md) — separate literal-`U` route planner, ordered R17 controls, and fail-closed foundry boundary.
- [`DEFECT_GRAPH_SMALL_GENUS_DYNAMICS_2026-09-03.md`](DEFECT_GRAPH_SMALL_GENUS_DYNAMICS_2026-09-03.md) — finite-prime reachability/trap controls.
- [`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md) — chronology and historical route context.
- [`../archive/elkies-k3/`](../archive/elkies-k3/) — superseded handoffs and experiments.

## Active fronts

The useful fronts are:

1. identify the stable arithmetic moduli curve for the determinant-720
   marking and decide whether a rational noncuspidal point has exact saturated
   determinant `720`; if not, continue the same arithmetic-first screen down
   the determinant-aware planner queue;
2. solve the 34 missing published-R17 rational visibility directions by a target-directed inverse problem, without enumerating the ambient trisection or quadrisection cosets;
3. compute the nonzero 2-primary quotient classes for the seventeen alternate-Q80 product twists from a full involution lattice or complete two-Selmer calculation, or obtain a complete degree-28 finite-field Frobenius rank bound; the zero class is closed, and the unsliced eight-variable `msolve` campaign remains superseded;
4. study multi-prime defect reachability and possible finite-prime sufficient sets;
5. add independent publication-grade replays or `J1` classification only when needed.

Do **not** restart the old giant Q80 reconstruction, broad foundry enumeration, or ungated rank searches just because those scripts remain in the repository.
Do not count another determinant-948 equation as the different-NS foundry
milestone.
Do not restart an arithmetic NS0024 source search over `QQ`; only geometric
or larger-field NS0024 work remains open.
Do not restart the NS0031 model-157 rational-point or algebraization search
over `QQ`; the full rational marking is now obstructed.

## Reproduction

<!-- status-consumer: EC-K3-R17-NONCYCLIC-4A1-DIRECT-EQUATION f657620e07f8f3f0 -->
<!-- status-consumer: EC-K3-R17-NORM12-SINGULAR-GENUS1-RATIONAL-NORMALIZATION-EXHAUSTION bf05d9b06ccc1502 -->
<!-- status-consumer: OP-K3-DIFFERENT-NS-ARITHMETIC-MW17 c384abf4d95dae7d -->
<!-- status-consumer: OP-EC-NEXT e135b23ef9910845 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION 077c6409d76cbe63 -->
<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-CONTROL-J-PREIMAGES 1ef38474a0d7f629 -->
<!-- status-consumer: EC-K3-ARITHMETIC-RANK-TRANSFER 3031dd2365a29cd5 -->
<!-- status-consumer: EC-K3-UNIVERSAL-DEGREE2-FIBRATION-COMPILER fd4b5d71c9497eaf -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-ALTERNATE-LAB-1024 c2f6309f8d6cc06d -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-V4-RATIONAL-PAIR-SHORTLIST-64 e14368b602eebedb -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-V4-BASE-RANK-SCREEN-64 f706a4396a0b13af -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-PRODUCT-BISECTION-INVERSION 6cfef74eb08601a6 -->
<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 291a539d07b842b9 -->

Use [`../REPRODUCE.md`](../REPRODUCE.md) and the exact checker paths recorded in `../MATH_STATUS.json`. `STATUS.md` is generated and is not an editing surface.
