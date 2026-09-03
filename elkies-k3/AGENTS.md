# AGENTS.md — elliptic K3 / high-rank programme

This directory inherits the repository rules. `../MATH_STATUS.json` is the only authority for what is proved.

## Programme state

**ACTIVE.** The programme is open for theorem-directed breakthrough work.

Large specialization, Selmer, point, neighbour, Q80 reconstruction, or foundry
searches must have an explicit mathematical gate, declared limits, checkpoints,
and a reproducible certificate plan.

## Start here

1. [`README.md`](README.md) — current milestone and active fronts.
2. [`R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md`](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md) — canonical alternate-Q80 equation route.
3. [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md) — reusable theorem/algorithm layer.
4. [`LITERATURE_AND_NOVELTY_MAP_2026-09-03.md`](LITERATURE_AND_NOVELTY_MAP_2026-09-03.md) — terminology and prior-art boundaries.
5. `../MATH_STATUS.json` and `../REPRODUCE.md` — exact status and replay.

## Current conclusions

- Published R17 and alternate Q80 are both explicit rootless rank-17 fibrations over `QQ` on the pinned determinant-948 K3.
- The alternate chart's canonical route is the direct degree-two `norm12-orbit-11952` hop, not the historical degree-11511 Q80 transport.
- The four published rank-25--28 R17 control fibres do not transfer to rational alternate-Q80 parameters; future alternate work needs native controls.
- Rank `>=32` is open.
- One curated inverse-ADE/foundry route is end-to-end planner-ready; the 936 bulk routes are not.

## Do not reopen by default

- historical degree-11511 alternate-Q80 transport;
- million-bit third-`q12` reconstruction and associated long PRS/Hensel work;
- broad Q80 suffix, q323, or changed-zero route searches;
- ungated Nagao/point/Selmer sweeps;
- bulk foundry route enumeration without complete marked planner inputs.

These remain useful provenance/regression material and are indexed from `../archive/elkies-k3/`.

## Claim discipline

- Keep classical infrastructure under established terminology: Shioda–Tate, Nikulin gluing, Kneser–Nishiyama, Kneser neighbours, and fibration hopping.
- Reserve novelty language for narrow inverse/target-directed pieces supported by the provenance map.
- Preserve `UNKNOWN` and bounded-negative qualifiers.
- Do not modify `STATUS.md` manually. If mathematical status genuinely changes, update the canonical proof and `MATH_STATUS.json`, then regenerate status.
- Preserve scripts and generated certificates even when an operational route is archived.
