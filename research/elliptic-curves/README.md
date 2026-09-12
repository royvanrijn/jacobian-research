# Elliptic curves over Q

[MATH_STATUS.json](../MATH_STATUS.json) is the authority.
[All elliptic claims](../index/elliptic-curves.md) ·
[Algorithmic lessons](../knowledge/ALGORITHMS.md) ·
[Unknowns and work](../knowledge/work/elliptic-curves.md) ·
[Current curve inventory](INVENTORY.md).

## Current milestone

Curve302 has a certified rank lower bound31. Its
[alternative determinant1092 MW17 parent](notes/CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md)
is explicit with a full saturated generic basis; construction recovery is complete.
[Calibrated V3](notes/ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md) recovers31 from
generic17 on that curve. These results do not prove exact rank31 or rank at least32.

The active task is [finding the next independent point](notes/NEXT_DIRECTION_RANK32_BENCHMARK_2026-09-12.md):
compare search representations on the retained27→28,29→30 and30→31 controls,
then apply a frozen, validated policy to existing rank27 R17/MW16 candidates.
Measure complete CPU to a certified gain beyond the starting subgroup.
The [constructor pilot is closed after failed positive calibration](rank-jump/FRESH_CONSTRUCTOR_TRANSFER_2026-09-12.md);
no fresh fibre ran. The theory results remain retained, but further constructor,
carrier or fixed-word work is not a prerequisite for the rank32 search.

Current curve and conductor totals are generated in the [inventory](INVENTORY.md)
and its [JSON export](data/research_curves/database.json). Dated notes retain
their historical cohort counts.

## Canonical entry points

| Question | Source |
|---|---|
| Load the recovered302 family and basis | [Parent proof](notes/CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md) · [loader](cas/load_curve302_recovered_parent.sage) |
| Understand adaptive point recovery | [V3 rule and calibration](notes/ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md) · [seed and transfer results](notes/CURVE302_SEEDED_V3_RESULTS_2026-09-08.md) |
| Recover points from constructed strict classes | [Both fixed covers solved blindly; exact transport and V3 comparison](rank-jump/BLIND_CONSTRUCTED_CLASS_RECOVERY_2026-09-12.md) |
| Find the next direction toward32 | [Frozen representation benchmark and follow-up protocol](notes/NEXT_DIRECTION_RANK32_BENCHMARK_2026-09-12.md) |
| Inspect the closed constructor pilot | [Failed positive calibration; no fresh fibres](rank-jump/FRESH_CONSTRUCTOR_TRANSFER_2026-09-12.md) |
| Propagate the marked two-class block | [The fixed column-6 word has finitely many Selmer specializations; fresh-dependency construction remains open](rank-jump/FIXED_WORD_HAS_FINITE_SELMER_SPECIALIZATIONS_2026-09-12.md) |
| Reuse the implemented search improvements | [Shared runtime](notes/SHARED_RESEARCH_RUNTIME.md) · [lean maps, cached continuation and box deduplication](notes/V3_FUTURE_SEARCH_PERFORMANCE_2026-09-08.md) |
| Understand what a rank-jump experiment measures | [Structural reassessment](notes/RANK_JUMP_REASSESSMENT_2026-09-05.md) |
| Compare exceptional ancestry and descent | [Principal28 control](notes/PRINCIPAL28_EXCEPTIONAL_ANCESTRY_2026-09-12.md) · [earlier bounded panel](notes/RANK_TRIANGLE_ANCESTRY_AND_DESCENT_2026-09-12.md) |
| Reuse the broad completed search | [Broad ledger](notes/BROAD_RANK_CURVE_LEDGER_2026-09-12.md) · [rank22 additions](notes/BROAD_RANK22_CURVE_LEDGER_2026-09-12.md) |
| Reuse the completed record-scale control | [48 fibres,2352 certified boxes,zero gains](notes/DET1092_RECORD_SCALE_CAMPAIGN_2026-09-07.md) |
| Interpret incomplete arithmetic | [Wide census](notes/WIDE_ARITHMETIC_PROFILE_2026-09-12.md) · [censoring and relation pilot](notes/TWO_CLASS_RELATION_PILOT_2026-09-12.md) |
| Recover Curve398's source | [MW16 recovery and duplicate-presentation proof](notes/ICARM_CURVE398_RANK30_AND_CONSTRUCTION.md) |

## Active fronts

The [generated open queue](../STATUS.md#active-open-problems) records
`OP-EC-NEXT` as the primary objective. `OP-EC-RANK-JUMP-MECHANISM-20260910`
is parked; reopening the constructor route needs a new mathematical reason
and separately scoped authorization.
Construction paths include direct R17/MW17, the deduplicated A1/MW16 family,
Curve302's recovered parent and the [different-NS K3 foundry](../elkies-k3/README.md).
A dated runbook or “ACTIVE” filename does not establish that a worker is running
or authorize a campaign restart.

## Proof and compute gates

Exactly verified independent points prove a lower bound without a completed
descent. Exact rank needs matching unconditional bounds. Generic saturation
does not assert specialized saturation. Scores, incomplete Selmer data,
bounded misses and timeouts cannot exclude a fibre mathematically.

Before another search, check the [method memory](../KNOWLEDGE_BASE.md), retained
exposure and [programme instructions](AGENTS.md). Keep known-record information
out of prospective selection where required by the frozen protocol.

## Reproduction

Use the [short replay guide](../REPRODUCE.md), [elliptic command reference](REPRODUCE.md)
and [complete source catalogue](../index/elliptic-curves.md).
The [previous long README](../archive/repository-cleanup-2026-09-12/research__elliptic-curves__README.md.txt)
is preserved as a dated navigation snapshot.

<!-- status-consumer: EC-PARENT-FOUNDRY-NEW-A1-CONSTRUCTIONS-20260910 5be8375fcf2cdfbe -->
<!-- status-consumer: EC-EUCLIDEAN-FOUR-SPLIT-ADMISSIONS-20260910 7a27319978e6a5fb -->
<!-- status-consumer: EC-DET1092-BLIND-MW16-RECONSTRUCTION-20260908 e799c490f343785a -->
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-SOLVER-COMPARISON 6a178bc3a4ada43b -->
<!-- status-consumer: EC-FIXED-CUBIC-CONIC-LONG-SEARCH 825fb4cd6ed84cb1 -->
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-MINIMAL-MODELS 90216b8c456edd20 -->
<!-- status-consumer: EC-FIXED-CUBIC-RADICAL-SEARCH-GEOMETRY 678f7beb805a4530 -->
<!-- status-consumer: EC-FIXED-CUBIC-TANGENT-CONIC-GATE 26a49e30ff3128d3 -->
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-CASSELS-TATE df45391a84f0e3c9 -->
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-RANK1 7e488a894d136732 -->
<!-- status-consumer: EC-K3-CURVE398-A1-MW16-RECOVERY 75978a18cc26690f -->
<!-- status-consumer: EC-K3-CURVE398-TWO-PARENT-COLLISION 626a440519ff77f3 -->
<!-- status-consumer: EC-K3-ICARM-MW16-BLIND-LADDER acfa3bdcebb18137 -->
<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER a2d7034fb8977c18 -->
<!-- status-consumer: EC-K3-MW17-JUMP-V2-ZERO-GAIN-RESCUE e5320f1f3bf33148 -->
<!-- status-consumer: EC-K3-R17-074D9-HALF-LATTICE-PROMOTION-GATE 9a1f080523c9ecae -->
<!-- status-consumer: EC-HALF-LATTICE-HEIGHT-COMPRESSION 3baeaf370aec751c -->
<!-- status-consumer: EC-CURVE385-HEIGHT-COMPRESSION-BUILDER-PILOT c0d6f2d67018def4 -->
<!-- status-consumer: EC-K3-R17-CURVE385-INDEPENDENT-RESTART-BUDGETS 39cfce110e3e494f -->
<!-- status-consumer: EC-FIXED-CUBIC-VARYING-CURVE-LOCAL-KUMMER 46ca45db3e702eb6 -->
<!-- status-consumer: EC-ICARM-CURVE302-POINT-CLOUD 1e1eb37dd6d4350f -->
<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 8a4c932153e2bb2d -->
<!-- status-consumer: EC-K3-R17-NORM12-ICARM-573-REFRESH a93ce35de34fde21 -->
<!-- status-consumer: EC-CF-NEARMISS-DESCENT-INPUTS 25c9f212e5162216 -->
<!-- status-consumer: OP-EC-NEXT b86e37cc3775f627 -->
<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-SELMER-PANEL 539bd8ec36b36c44 -->
<!-- status-consumer: EC-K3-ICARM-MW16-POINTED-SIEVE cb83c1afae1d0141 -->
<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-VS-SHA-COMPARISON f37417a9fda3ee3f -->
<!-- status-consumer: EC-K3-ICARM-MW16-SENSITIVITY f88886c066d6cb45 -->
<!-- status-consumer: EC-FIXED-FIELD-COMPARISON 02c49a8120aeb7bd -->
