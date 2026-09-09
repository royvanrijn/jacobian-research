# Frozen sixty-fibre R17 seed/complement panel

This prospective experiment tests whether inexpensive first-seed acquisition
and a short adaptive complementary pass identify fresh fibres worth further
computation. It is motivated by the [fresh-cohort gains](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md),
[second cohort](SECOND_FRESH6_SEED_COHORT_2026-09-09.md), and
[lower-height cloud reconciliation](LOWHEIGHT_FRESH6_SEED_COHORT_2026-09-09.md).
The [302 panel](CURVE302_SEED_UNIVERSALITY_RESULTS_2026-09-09.md) is a separate
retrospective calibration: its 13/14 closures do not predict this panel's rate.

## Frozen selection

Protocol SHA256:
`edd883ef31ee13a475897a10e1174a495f9e80c9677bc1ab0cbf29038b886379`.

Exactly ten candidates come from each of 074d9, 07ca9, 08234, 08f72, 103b2,
and11952. Each family contributes five candidates of parameter height at most
1024 and five of height1025–4096. Within each stratum, the existing selection
score, good-prime count, denominator and signed numerator determine order.
No validation-band score or point/rank outcome enters this ordering.

The population is the same6144 retained score rows. A frozen projection of
named prior result ledgers and point-run rosters excludes equations up to
exact rational isomorphism and reserved family/parameter addresses. All60
selected equations are pairwise nonisomorphic. The three previous six-fibre
cohort reservations participate in exclusions. Freshness is relative to
these snapshots, not every external search or a literature-wide novelty test.

- [Frozen roster](../../artifacts/generated-results/elliptic-curves/r17_60_panel_roster_v1.json).
- [Protocol](../../artifacts/generated-results/elliptic-curves/r17_60_panel_protocol_v1.json).
- [Reservation ledger](../../artifacts/generated-results/elliptic-curves/r17_60_reserved_results_v1.json).
- [Selector](../cas/select_r17_60_panel.py).

## Search and proof gates

Before any point search, all60 exact generic M17 inputs must pass native
specialization and independent finite-certificate replay. Selection is not
refilled using seed-search outcomes.

1. The existing first-M18 constructor tests the43/49 exact maximum generic
   classes: at most86/98 point calls at height125000, ten seconds per call,
   with separate five-second/1-GiB map bounds. It stops at the first standalone
   M18 certificate. All seed attempts receive independent geometry/map/point
   replay, including bounded misses.
2. Replay and reconcile the complete retained seed cloud. Any further
   independent points already returned count without new searches. This can
   make the starting complementary subgroup larger than M18.
3. Select sixteen exact maximum generic classes outside the original winning
   mask's span. In ascending mask order first extend the span, then fill to16.
   Both the integer and rational CVP solvers certify all candidate minima.
4. Run the existing complementary V3 selector for **at most100 total point
   invocations**, summed across all adaptive epochs and reconciliation
   restarts. Every larger cloud must be independently replayed and reconciled
   before a restart. A restart receives only the remaining call allowance.
5. Stop at the budget, target lower bound32, or finite-policy exhaustion.
   There are no deep suffix passes, new candidate refills or outcome-dependent
   budget increases. Resource limits retain an unresolved outcome.

Each expensive search/replay phase has1800 seconds and3GiB; preparation
phases have600 seconds and3GiB. Four independent case processes run through
the shared process-group supervisor. Every phase and map retains logs and
resource receipts. Other active research jobs are outside this controller.

## What is compared

Report final certified subgroup rank, added directions, distinct gaining calls,
exact calls per added direction, last gain position, calls since last gain,
and the complete gain timeline. Arithmetic-only reconciliation gains are
attributed to their original point calls. Keep point and map censoring visible.

The primary table orders by rank descending, number of gaining complementary
calls descending, calls per added direction ascending, then later last gain.
A separate late-gain list retains fibres gaining at complementary call75 or
later, so a lower-rank ongoing cascade is visible beside higher-rank stalls.
Rank and total added directions are redundant here because all start at17;
the number of gaining calls measures a different aspect of the trajectory.

Report height strata and fibrations separately. This selected finite panel
does not establish a general rank distribution, causal benefit of lower
parameter height, exact ranks, conductor improvements or public novelty.
No conductor factorization is part of its point-search budget.

## Execution and evidence

All sixty generic M17 packets passed native specialization and independent
finite-certificate replay before launch:

```text
R17_60_PREFLIGHT_COMPLETE|cases=60|rank=17|charts=0|status=PASS
```

The four-worker detached controller was launched on2026-09-09. Point-search
outcomes are pending; this preflight certifies the inputs, not any new rank
gain. Its receipt is retained as `r17-60-panel-v1/preflight.json` under the
local evidence root below.

Raw checkpoints and per-phase supervisor logs are retained under
`artifacts/local/elliptic-curves/r17-60-panel-v1/`. The sources, projected
selection inputs, roster, atlas and exact maximum-class catalogue are bound
before execution. Completed case packets include full independent point
certificates and evidence hashes. The final report will be exported as
`artifacts/generated-results/elliptic-curves/r17_60_panel_results_v1.json`.

From `research/`:

```sh
python3 elliptic-curves/cas/run_r17_60_panel.py status
```

The launcher requires the all60 zero-chart preflight. An interrupted panel
can be reviewed and resumed with `run_r17_60_panel.py resume`; exclusive
controller/case locks prevent duplicate ownership. It reuses completed
searches and replays. Source/binding or proof failures stop further dispatch.

Ten regression tests cover stratification, outcome-independent ordering,
isomorphism/address exclusions, call attribution, shared budget across
reconciliation restarts, locks, resource classification and exact CVP.
A retained-input zero-search Sage regression reconstructs an M22 seed cloud
twice and independently checks its sixteen exact complementary parents.
The new adapters explicitly normalize JSON values before immutable comparisons;
the initial serialization failures were corrected before this protocol freeze.

A separate [post-execution audit](../cas/audit_r17_60_panel.py) was added after
freeze; it is not imported by selection or execution. It checks the roster,
source/evidence bindings, full epoch replay receipts, cumulative call budgets
and reported gain timing. Five corruption/accounting regressions bring the
targeted test set to15 passing tests. Its optional Sage mode freshly checks
the native generic prefix and finite independence of each completed packet.
The complete gate rejects partial panels and resource-unresolved cases:

```sh
~/.local/bin/sage -python elliptic-curves/cas/audit_r17_60_panel.py --proofs --complete
```

The [post-run analysis](../cas/analyze_r17_60_panel.py) separately reports
seed-cloud and complementary gains, conditional seed costs, height/fibration
groups, finite-policy exhaustion and late gains. It does not alter dispatch.
Two further metric tests bring the targeted test set to17. Partial summaries
are read-only; publication and replay require the complete audit:

```sh
python3 elliptic-curves/cas/analyze_r17_60_panel.py --partial
python3 elliptic-curves/cas/analyze_r17_60_panel.py --write
python3 elliptic-curves/cas/analyze_r17_60_panel.py --check
```

The [portable certificate checker](../cas/verify_r17_60_certificates.py) uses
only repository sources and generated artifacts. It re-specializes the native
M17 basis, checks every final finite independence certificate, and verifies
pairwise rational nonisomorphism. It does not establish the search history;
that remains the separate full local audit above. After the complete export:

```sh
python3 elliptic-curves/cas/verify_r17_60_certificates.py
```

Implementation: [controller](../cas/run_r17_60_panel.py),
[arithmetic adapters](../cas/r17_60_arithmetic.py), and
[retained-input regression](../cas/check_r17_60_adapters.py).
