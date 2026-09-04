# AGENTS.md — elliptic K3 / high-rank programme

This directory inherits the repository rules. `../MATH_STATUS.json` is the only authority for what is proved.

## Programme state

**ACTIVE.** The programme is open for theorem-directed breakthrough work.

Large specialization, Selmer, point, neighbour, Q80 reconstruction, or foundry
searches must have an explicit mathematical gate, declared limits, checkpoints,
and a reproducible certificate plan.

## Start here

1. [`README.md`](README.md) — two primary lanes and current certified position.
2. [`DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md`](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md) — open different-NS objective and fail-closed source/endpoint gates.
3. [`R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md`](R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md) — determinant-948 noncyclic closure that motivates the pivot.
4. [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md) — reusable theorem/algorithm layer.
5. [`LITERATURE_AND_NOVELTY_MAP_2026-09-03.md`](LITERATURE_AND_NOVELTY_MAP_2026-09-03.md) — terminology and prior-art boundaries.
6. `../MATH_STATUS.json` and `../REPRODUCE.md` — exact status and replay.

## Current conclusions

- Published R17 and alternate Q80 are both explicit rootless rank-17 fibrations over `QQ` on the pinned determinant-948 K3.
- The alternate chart's canonical route is the direct degree-two `norm12-orbit-11952` hop, not the historical degree-11511 Q80 transport.
- The maximal noncyclic determinant-948 bridge is equation-explicit in both
  directions: `R17 -> 4A1/MW13 -> R17`, with thirteen saturated rational
  sections and target-free reverse selection.
- Lane B is arithmetic-first globally: `T` and its full stable marked curve
  precede `NS=T^perp`, rootlessness, and equation work. Determinants `720`,
  `950`/`NS0024`, and `1184`/`NS0031` are arithmetically excluded. The
  determinant-720 stable curve is `X_0(60)` and has only rational cusps; its
  known rational `3A5` point instead saturates to determinant `20`.
- The construction target is stronger than a plain different-NS MW17:
  require a certified positive-rank low-genus carrier and an independent
  pullback section. The stretch target is an integral `V4`-stable MW lattice
  with character ranks `17+1+1+1` and exact 2-primary graph glue.
- The four published rank-25--28 R17 controls do not transfer to rational alternate-Q80 parameters, but the complete ICARM sweep now supplies native controls: curve 12 has rank at least 29 in class `11952`, while curves 363, 364, 378, and 395 give further rank-at-least-24--28 fibres in classes `08f72` and `11952`.
- Rank `>=32` is open.
- One curated inverse-ADE/foundry route is end-to-end planner-ready; the 936 bulk routes are not.

## Do not reopen by default

- historical degree-11511 alternate-Q80 transport;
- million-bit third-`q12` reconstruction and associated long PRS/Hensel work;
- broad Q80 suffix, q323, or changed-zero route searches;
- ungated Nagao/point/Selmer sweeps;
- bulk foundry route enumeration without complete marked planner inputs.
- treating another determinant-948 equation as the different-NS foundry
  milestone.
- restarting a full-rational-marking NS0024 search over `QQ`.
- restarting the NS0031 model-157 rational-point or algebraization search over
  `QQ`.
- reopening determinant 720 as a full rational rank-19 source over `QQ`.
- inspecting rootless frames or launching coefficient searches before their
  full marked `T` curve has a certified rational non-CM point.

These remain useful provenance/regression material and are indexed from `../archive/elkies-k3/`.

## Claim discipline

- Keep classical infrastructure under established terminology: Shioda–Tate, Nikulin gluing, Kneser–Nishiyama, Kneser neighbours, and fibration hopping.
- Reserve novelty language for narrow inverse/target-directed pieces supported by the provenance map.
- Preserve `UNKNOWN` and bounded-negative qualifiers.
- Do not modify `STATUS.md` manually. If mathematical status genuinely changes, update the canonical proof and `MATH_STATUS.json`, then regenerate status.
- Preserve scripts and generated certificates even when an operational route is archived.

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->
<!-- status-consumer: EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION 972f591d2885f9ba -->
<!-- status-consumer: EC-K3-DIFFERENT-NS-ARITHMETIC-GATE-RERANK 252991e141c42e55 -->
<!-- status-consumer: EC-K3-ARITHMETIC-FIRST-MARKED-T-FOUNDRY 6b9d34ae8d722280 -->
<!-- status-consumer: OP-K3-DIFFERENT-NS-ARITHMETIC-MW17 71f43dc9ef3af620 -->
