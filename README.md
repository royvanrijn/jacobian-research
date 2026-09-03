# Jacobian Research

Exact and computational mathematics with machine-replayable proof and status records.

## Current focus

**Prime-gap calculations are the active compute priority.**

The elliptic-K3 / high-rank elliptic-curve programme is **paused as of 3 September 2026**. Do not launch large K3, descent, neighbour, or specialization searches unless that programme is explicitly resumed.

## Canonical repository state

- [`MATH_STATUS.json`](MATH_STATUS.json) — sole machine-readable authority for proved, conditional, falsified, and open claims.
- [`REPRODUCE.md`](REPRODUCE.md) — replay commands and pinned certificate entry points.
- `STATUS.md` — generated view of `MATH_STATUS.json`; do not edit by hand.
- [`verified/`](verified/) — durable verified theorem notes.
- [`artifacts/generated-results/`](artifacts/generated-results/) — generated certificates.
- [`archive/`](archive/) — superseded or historical research; not a current claim source.
- [`elkies-k3/README.md`](elkies-k3/README.md) — frozen K3/high-rank milestone and resume point.
- [`elliptic-curves/README.md`](elliptic-curves/README.md) — elliptic-curve programme notes.
- [`KNOWLEDGE_BASE.md`](KNOWLEDGE_BASE.md) and [`RESEARCH_TIMELINE.md`](RESEARCH_TIMELINE.md) — broader synthesis and chronology.

## K3 / high-rank milestone at pause

The two determinant-948 rootless rank-17 charts on the pinned K3 are explicit over `QQ`.

- **Published Elkies R17:** certified `24 I1`, geometric and arithmetic generic Mordell–Weil rank 17.
- **Alternate Q80:** canonical direct degree-two hop from `norm12-orbit-11952`; polynomial K3 model with `(deg A, deg B, deg Delta)=(8,12,24)`, `24 I1`, determinant-948 rootless frame, and 17 saturated rational sections.
- The alternate chart has generic arithmetic rank 17 over `QQ`.
- The four published rank-25--28 R17 controls have no rational preimage under the alternate-Q80 `j`-map, so future alternate work needs native calibration fibres.
- Rank `>=32` remains open.
- Target-directed fibration hopping is frozen after one curated planner-ready end-to-end control; the 936 bulk foundry routes are not planner-ready.

Historical Q80 transports, large-coordinate reconstruction routes, and expensive third-`q12` work are provenance only. See [`archive/elkies-k3/`](archive/elkies-k3/).

## Working rules

1. Treat `MATH_STATUS.json` as truth; prose only summarizes it.
2. Keep `UNKNOWN` as `UNKNOWN`.
3. Prefer exact certificates and small independent checks over large speculative searches.
4. Archive superseded research instead of leaving stale instructions on active navigation surfaces.
5. Never edit generated `STATUS.md` manually.

Other research branches remain in the repository with their canonical notes and certificates; this README intentionally stays lean and does not duplicate them.
