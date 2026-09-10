# X1092 seed and carrier experiment ledger

Consolidated 2026-09-11. This is an experiment-history ledger, not mathematical-status authority. `research/MATH_STATUS.json` and linked proof certificates remain authoritative.

## Executive conclusion

The completed Elkies-style `M19 -> V3` experiment is **not a wholly new idea**. It is the first end-to-end prospective completion on X1092 of a route whose pieces had already been explored separately:

1. X1092 already had explicit rational bisections giving generic `M18` covers and many `M18 -> V3` experiments.
2. X1092 already had two distinct bisections whose common genus-one base carried a generic rank-at-least-19 subgroup, but rational points on that base were not known.
3. Retrospective X948 rank-jump work already exhibited a positive-rank genus-one carrier supporting two marked exceptional directions.
4. A small control family already showed that one quadratic cover can generically add two independent rational directions.

What is genuinely new in the 2026-09-10/11 Elkies-pair lane is the **completed prospective chain on X1092**:

`generic M17 -> two generic bisections -> positive-rank genus-one carrier -> rational M19 fibres -> independently certified M19 packets -> unchanged V3`.

The first fully completed detached V2 run then tested three certified rank-at-least-19 fibres for 100 V3 calls each, with all independent replays passing and **no further rank gain**. This is new negative evidence about amplification from prospectively constructed M19 fibres.

## Completed experiments

| Date | Surface/family | Construction | Starting subgroup | Search/amplification exposure | Result | Interpretation / repeat status |
|---|---|---|---:|---:|---|---|
| 2026-09-08 | X1092, orbit8044 conic | One rational bisection / quadratic base change | M18 | Six fresh fibres, 114 V3 charts each | 18 -> 18 on all six | **Completed bounded no-gain. Do not repeat the same single-conic factory as a generic amplification test without a new selector.** |
| 2026-09-08 | X1092, orbit8044 conic | Same M18 construction, fibre `s=1926/2699` | M18 | 568 V3 charts | 18 -> 21 | Positive proof that constructed M18 can amplify, but does not predict which fibres do so. |
| 2026-09-08/09 | X1092 / curve302 | Known exceptional 302 seeds supplied retrospectively | M18 | 855, 637, 385 charts in sealed examples; larger seed-universality panel later | 18 -> 31 repeatedly | Strong amplification calibration, but **retrospective/oracle-seeded**, not a prospective production mechanism. |
| 2026-09-08 | X1092, orbits 8044 + 127449 | Two independent rational bisections; degree-four common cover | generic M19 over genus-one base | No rational specialization campaign | rank >=19 over `Q(B)`; rational point on `B` UNKNOWN | Direct predecessor of Elkies-pair lane. Geometry solved; arithmetic point source was missing. Superseded as an open gate by the later positive-rank-carrier construction. |
| 2026-09-07-era rank-jump study | X948 fibre `08234-003` | Minimal genus-one carrier for two known native directions | marked two-direction block / M19 over carrier | Retrospective structural analysis | Carrier has rational origin and Jacobian rank 3; chosen fibre realizes two directions | Important mechanism control, but **retrospective and tied to a known +7 fibre**; not a prospective X1092 selector. |
| 2026-09-07-era mechanism control | Small rational elliptic-surface family | One shared squareclass makes two sections rational on one double cover | rank 1 -> generic rank 3 | Symbolic + fixed anchors | Two independent generic directions proved | Shows same-cover rank+2 is real in principle. Not a production X1092 construction. |
| 2026-09-10 | X1092 Elkies-pair V1 | Generic bisection-pair search; first positive-rank carrier | generic M19 carrier | Default scoring/height policy | Carrier succeeded, but only two eligible fibres survived and two reserved controls starved ranked search | **Invalid search exposure, not a mathematical null.** Superseded by V2 starvation-safe scheduler. |
| 2026-09-10 CI smoke | X1092 Elkies-pair | Positive-rank genus-one carrier -> specialization | M19 | One control reached M19; V3 launch attempted | M19 certification passed; V3 backend packaging failed | Engineering smoke only. Superseded by detached V2 end-to-end run. |
| 2026-09-11 | X1092 Elkies-pair V2 | Carrier-quality selection, starvation-safe controls, deep Mestre scoring | Three independently certified M19 fibres: one ranked + two controls | **100 V3 calls each, 300 total; no timeouts; all independent replays passed** | **19 -> 19 on all three** | **Completed bounded no-gain. The end-to-end M19-to-V3 pipeline works; generic M19 alone is not an observed amplification trigger in this panel.** |

Sources for the historical rows:

- `DET1092_SEED_DENSITY_AND_LIMITING_LATTICE_2026-09-08.md`
- `M18_LANDSCAPE_COMPARISON_2026-09-08.md`
- `DET1092_TWO_COVER_BRANCH_GATE_2026-09-08.md`
- `rank-jump/MINIMAL_CARRIER_AND_RATIONAL_SPLITTING_OF_A_TWO_DIRECTION_BLOCK.md`
- `rank-jump/ONE_SQUARE_CONDITION_TWO_RATIONAL_DIRECTIONS.md`
- `X1092_ELKIES_PAIR_SEARCH_2026-09-10.md`

The final V2 row is a user-reported completed detached run on 2026-09-11. Its local report/certificates were not inspected or committed by this ledger update; preserve the local run directory as evidence.

## What has actually been learned

### 1. `M18 -> V3` is old and well tested

The project already had a substantial prospective M18 experiment before the Elkies-pair work. Fresh conic seeds often stayed at 18; one reached 21. Therefore "construct one generic extra section and let V3 amplify it" is **not** a new strategy and should not be proposed again as such.

### 2. `M19 over a genus-one compositum` is also old

The two-cover branch gate already proved a generic rank-at-least-19 subgroup over a genus-one common base. The old blocker was not independence or genus: it was the absence of a rational-point source on that base.

The Elkies-style search is new because it selects a pair whose genus-one base has a certified nontorsion rational point, hence infinitely many rational points, and then carries actual specializations through finite-rank certification into V3.

### 3. Positive-rank two-direction carriers existed before, but retrospectively

The X948 `08234-003` work already found a genus-one carrier of two exceptional directions with auxiliary Jacobian rank 3. That construction used directions from an already successful rank-jump fibre, so it cannot be counted as a prospective production result. It did, however, already establish the geometry/arithmetic pattern that the Elkies-pair lane later made prospective on X1092.

### 4. The new null is specific and useful

The V2 result does **not** prove that M19 fibres never amplify. It does show that, under the current prospective X1092 carrier construction, carrier-quality selection, two-stage Mestre scoring, and unchanged V3 policy, three certified M19 fibres consumed 300 clean calls with no gain.

Therefore the hypothesis

> "getting two constructed generic extra directions at once is by itself enough to unlock 302-like amplification"

has no positive support from this first prospective panel.

## Anti-repeat gate

**Do not launch another ad-hoc X1092 `generic M19 carrier -> ordinary V3` experiment that differs only by more addresses, a larger height cap, another small carrier pair window, or modest Mestre-threshold retuning.** The route is now implemented and has a clean completed no-gain panel.

A new experiment in this family should state which genuinely new axis it tests. Acceptable reopen conditions include at least one of:

- a **same-quadratic-cover** construction with two independent X1092 sections, giving M19 over a rational base rather than a genus-one compositum;
- a carrier giving **M20 or higher** before specialization;
- a different X1092 maximal-rank fibration/frame with a mathematically distinct carrier structure;
- a prospective arithmetic selector tied to the strict/global Kummer or ideal-class structure rather than carrier height/Mestre score alone;
- a new V3/closure theorem or search operator specifically motivated by the 302 seed-universality phenomenon;
- a predeclared statistically meaningful large M19 panel whose purpose is explicitly to estimate amplification frequency, not another exploratory rerun.

Merely finding another positive-rank genus-one M19 carrier is **not** by itself sufficient reason to restart the same amplification experiment.

## Current operational interpretation

- Keep the broad production search independent and running according to its own frozen protocol.
- Treat the X1092 Elkies-pair M19 lane as **implemented, end-to-end validated, and currently bounded-negative for amplification**.
- Preserve it as infrastructure and as a control lane, not as the default next search.
- Theory work should now focus on what distinguishes the exceptional 302 closure from ordinary constructed M18/M19 seeds, rather than on acquiring the first or second generic extra direction again.
