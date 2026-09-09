# Autonomous R17 rank hunter — 2026-09-09

## Purpose

This is an additive production controller for a detached, non-AI-in-the-loop search for new high-rank elliptic curves over `Q`. It leaves the frozen historical seed, V3, panel and unattended runners untouched and reuses their exact replay/certificate machinery.

Primary target: certified lower bound **M32**. M31/M30/M29/M28 are automatically retained as significant results. Exact rank is never inferred from a bounded miss.

## Why this search policy

The controller encodes the main prospective evidence from the September 9 programme rather than using current rank as the queue key:

- the 60 fresh R17 fibres beat deep grinding of four old M27/M28 states at similar point-search cost;
- full returned clouds are always reconciled before scheduling decisions;
- productive fresh cascades earn more search only while they retain momentum;
- a stale fixed parent bank is not allowed to absorb thousands of calls merely because its current certified rank is high;
- parent classes are rotated instead of restricting later work to a previously productive span;
- low parameter height receives more exploration, but high-height fibres are never excluded;
- score is a scheduling signal, never a mathematical rank classifier.

## Architecture

`select_autonomous_r17.py` freezes the full 6,144-row retained H4096 population and a repository-known exclusion snapshot. Fresh selection cycles through score-frontier, deterministic low-discrepancy and deterministic hash-spread choices. The twelve-step lane schedule gives low height a 2:1 preference while retaining all three exploration modes at high height. Families are balanced with only a weak empirical prior.

`autonomous_rank_hunter_arithmetic.py` constructs/replays native generic M17 packets, reconciles the complete seed cloud, constructs exact parent banks, preserves an enlarged certified subgroup when a cloud reconciliation changes the basis, and adapts bounded V3 segments. The first banks partition unused exact maximum generic parity classes; later bank generations can move into deterministic broad parity samples rather than staying in the productive-mask span.

`autonomous_rank_hunter_policy.py` makes the explore/exploit decision. With four workers, normally at least two remain on fresh-fibre acquisition. Only a currently productive M29+ cascade may temporarily reduce exploration to one worker. Exploitation is issued in 100-new-call tranches. Rank is a soft proximity term; recent independent directions, calls since the last gain, gain efficiency, seed-cloud multiplicity and late gains dominate the decision.

The fixed-bank no-gain limits are deliberately much smaller than the historical failed unattended grind: M24 100, M25 125, M26 150, M27 200, M28 300, M29 450, M30 600, M31 750 calls. Once a bank is stale, an empirically productive curve may receive a new parent/representation bank; otherwise it is retired. The number of justified banks also grows only with observed cascade evidence and target proximity.

`run_autonomous_rank_hunter.py` is the detached controller. Every expensive arithmetic phase runs under the repository supervisor with wall/RSS limits. Search/replay seals are required before a result can influence scheduling. It supports `start`, `status`, `stop` and `resume`, keeps a process start token and controller lock, preserves per-case ledgers, and stops dispatch on target M32, configured time/call/case bounds, low disk, explicit STOP, or an integrity failure.

## Conductor side objective

Every certified result is given a cheap integral-discriminant screen against exact conductors already present in the research database. This does **not** factor a discriminant and does not claim a conductor record. It only identifies arithmetically small candidates without allowing conductor work to consume point-search workers.

## Evidence retained

Each case keeps the frozen candidate row, generic M17 packet and independent replay, seed search and replay, full seed-cloud reconciliation, exact parent-bank derivation/CVP evidence, every V3 search/replay packet, full-cloud reconciliations, scheduler summary and bounded worker supervision records. M28+ packets are additionally copied to `research/artifacts/generated-results/elliptic-curves/autonomous_rank_hunter_v1/` with their exact points and finite-reduction proof.

Freshness means nonisomorphic to the campaign's frozen repository exclusion snapshot. The autonomous controller does not make an automatic live-world-catalogue novelty claim.

## Launch

From the repository root, after these sources are on the checked-out branch:

```sh
python3 research/elliptic-curves/cas/run_autonomous_rank_hunter.py start \
  --workers 4 --max-hours 168 --max-point-calls 250000 --max-cases 1000 --min-free-gib 10
```

The `start` command launches the controller with `start_new_session=True`; no AI/model/API call is made by the production search.

```sh
python3 research/elliptic-curves/cas/run_autonomous_rank_hunter.py status
python3 research/elliptic-curves/cas/run_autonomous_rank_hunter.py stop
python3 research/elliptic-curves/cas/run_autonomous_rank_hunter.py resume
```

Runtime root:

```text
research/artifacts/local/elliptic-curves/autonomous-rank-hunter-v1/
```

Controller log:

```text
research/artifacts/local/elliptic-curves/autonomous-rank-hunter-v1/controller.log
```

## Validation performed for this change

The pure scheduling regression suite covers the known productive-M27, stalled-M26 and stale-M27 pivot cases, breadth reservation, M29 hot-cascade exploitation and explicit high-height exploration. It passed 6/6. The four new Python modules also pass `py_compile` in the implementation environment.

A Sage end-to-end production preflight and detached launch still has to execute in the actual repository checkout, because the implementation environment used for this change did not expose `/home/royvanrijn/src/jacobian-research` or its local retained-score artifacts. That is an execution-environment limitation, not mathematical evidence.