# Elliptic-curve computation map

Use the [programme map](../README.md), [current inventory](../INVENTORY.md)
and selected claim in [MATH_STATUS.json](../../MATH_STATUS.json). Dated cohort
counts and replay stage totals are evidence for those experiments, not the
current repository inventory.

| Need | Entry point |
|---|---|
| Load Curve302's recovered MW17 parent | [Parent proof](../notes/CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md) · [loader](load_curve302_recovered_parent.sage) |
| Run shared arithmetic or point search | [Runtime contracts and measured improvements](../notes/SHARED_RESEARCH_RUNTIME.md) · [runtime modules](research_runtime/README.md) |
| Understand calibrated V3 | [Exact rule and calibration](../notes/ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md) · [seed and transfer outcomes](../notes/CURVE302_SEEDED_V3_RESULTS_2026-09-08.md) |
| Reuse performance improvements | [Lean maps, cached continuation and box deduplication](../notes/V3_FUTURE_SEARCH_PERFORMANCE_2026-09-08.md) |
| Certify or interpret ranks | [Stable commands](../scripts/README.md) · [replay reference](../REPRODUCE.md) |
| Recover Curve398's source | [Recovery and equivalent MW16 presentations](../notes/ICARM_CURVE398_RANK30_AND_CONSTRUCTION.md) |
| Use the corrected Mestre/Fermigier labels | [Coherent-label audit and historical failure](../notes/PARENT_PORTFOLIO_AND_SECTION_LABEL_AUDIT_2026-09-07.md) |
| Reuse the broad completed search | [Curve ledger](../notes/BROAD_RANK_CURVE_LEDGER_2026-09-12.md) · [rank22 additions](../notes/BROAD_RANK22_CURVE_LEDGER_2026-09-12.md) |
| Investigate the rank-jump mechanism | [Structural reassessment](../notes/RANK_JUMP_REASSESSMENT_2026-09-05.md) · [current mechanism notes](../rank-jump/README.md) |
| Reproduce the six-root family | [Family workflow and exact-rank14 control](newfamily/README.md) |

The shared production entry points are [run_arithmetic_pipeline.py](run_arithmetic_pipeline.py),
[run_mw_search.py](run_mw_search.py), [run_pointed_quartic_search.py](run_pointed_quartic_search.py)
and [run_surface_proof.py](run_surface_proof.py). Read their declared request
formats and canonical proof before scheduling work.

Reuse exact arithmetic contexts, labelled two-torsion algebras, finite
reductions, immutable subgroup state and compatible search receipts. Keep
raw/minimal and coordinate-map identities explicit. A frozen experiment's
historical runtime is not interchangeable with a current amended worker.

Exactly verified independent points prove a lower bound. An unconditional
matching upper bound proves exact rank; incomplete descent, scores and bounded
misses remain unknown. The legacy two-section generic-rank13 entry point is
rejected; use the coherent-label replay and its rank11 generic lower bound.

Find an individual checker through `research.py search QUERY` or `show CLAIM-ID`.
The [full claim catalogue](../../index/elliptic-curves.md) links every canonical
source. The [original 592-line CAS map](../../archive/repository-cleanup-2026-09-12/research__elliptic-curves__cas__README.md.txt)
retains earlier versioned inventories and exact replay descriptions. Its old
descent-before-search instruction has been superseded by the programme's
unconditional rank gates.

Completed and superseded drivers remain in the [EC archive](../../archive/elliptic-curves/README.md).
Missing local outputs are not an instruction to regenerate them; broad searches
and long replays require a separately scoped mathematical reason.

<!-- status-consumer: EC-K3-ICARM-MW16-BLIND-LADDER acfa3bdcebb18137 -->
<!-- status-consumer: EC-FIXED-CUBIC-VARYING-CURVE-LOCAL-KUMMER 46ca45db3e702eb6 -->
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-RANK1 7e488a894d136732 -->
<!-- status-consumer: EC-K3-MW17-JUMP-V2-ZERO-GAIN-RESCUE e5320f1f3bf33148 -->
<!-- status-consumer: EC-HALF-LATTICE-HEIGHT-COMPRESSION 3baeaf370aec751c -->
<!-- status-consumer: EC-CURVE385-HEIGHT-COMPRESSION-BUILDER-PILOT c0d6f2d67018def4 -->
<!-- status-consumer: EC-K3-R17-CURVE385-INDEPENDENT-RESTART-BUDGETS 39cfce110e3e494f -->
<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-RANK-ESCAPE-DETECTOR-V2 eda7a0053b31b7c9 -->
<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER a2d7034fb8977c18 -->
<!-- status-consumer: EC-K3-ICARM-MW16-POINTED-SIEVE cb83c1afae1d0141 -->
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-CASSELS-TATE df45391a84f0e3c9 -->
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-MINIMAL-MODELS 90216b8c456edd20 -->
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-SEARCH-GEOMETRY 678f7beb805a4530 -->
<!-- status-consumer: EC-FIXED-CUBIC-TANGENT-CONIC-GATE 26a49e30ff3128d3 -->
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-SOLVER-COMPARISON 6a178bc3a4ada43b -->
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-LONG-SEARCH 825fb4cd6ed84cb1 -->
<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-VS-SHA-COMPARISON f37417a9fda3ee3f -->
<!-- status-consumer: EC-K3-ICARM-MW16-SENSITIVITY f88886c066d6cb45 -->
<!-- status-consumer: EC-FIXED-FIELD-COMPARISON 02c49a8120aeb7bd -->
