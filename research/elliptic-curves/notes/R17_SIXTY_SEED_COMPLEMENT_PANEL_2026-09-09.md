# Frozen sixty-fibre R17 seed/complement panel

This prospective experiment tests whether inexpensive first-seed acquisition
and a short adaptive complementary pass identify fresh fibres worth further
computation. It is motivated by the [fresh-cohort gains](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md),
[second cohort](SECOND_FRESH6_SEED_COHORT_2026-09-09.md), and
[lower-height cloud reconciliation](LOWHEIGHT_FRESH6_SEED_COHORT_2026-09-09.md).
The [302 panel](CURVE302_SEED_UNIVERSALITY_RESULTS_2026-09-09.md) is a separate
retrospective calibration: its 13/14 closures do not predict this panel's rate.

## Completed findings

All 60 cases completed on 2026-09-09 in 4014.98 seconds (66.92 minutes),
with four detached workers. Every case passed independent search replay and
saved-cloud reconciliation. The final full-history audit, fresh native-prefix
and finite-independence checks, and portable exported-certificate verification
all passed. All 60 selected equations are pairwise nonisomorphic over Q.

There were **48 first seeds**, including 26 found on the first point call.
Of these 48 seeded fibres, **38 gained further directions** in the complementary
pass; ten did not. The experiment used 1349 seed calls and 4741 complementary
calls, **6090 calls total**, with no point/map censoring or phase-level resource
stops. Across the 60 separate subgroups, seed-chart clouds contributed 85
directions beyond their respective M17 inputs, and complementary passes added
153 more. These are counts across different curves, not the rank of one curve.

The final certified subgroup lower bounds are:

| Lower bound | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cases | 12 | 7 | 4 | 5 | 5 | 4 | 6 | 9 | 4 | 3 | 1 |

Twelve cases stopped after bounded seed misses, 41 exhausted their 100-call
complementary allowance, and seven exhausted the finite complementary policy
earlier. None received a deep suffix or an increased budget. No rank-at-least-28
result was obtained in this panel; this is not an upper bound on any fibre.

### Highest-ranked outcomes

Rows follow the frozen reporting order. All eight used one seed call and 100
complementary calls. "Gain calls" counts distinct complementary calls that
supplied admitted directions, including directions recovered from saved clouds.
"Calls/direction" includes seed calls and divides by the final increase over 17.

| Family, parameter | Seed-cloud bound | Final bound | Gain calls | Calls/direction | Last complementary gain |
| --- | ---: | ---: | ---: | ---: | ---: |
| 103b2, 877/781 | 19 | 27 | 5 | 101/10 | 67 |
| 074d9, -1013/476 | 19 | 26 | 7 | 101/9 | 40 |
| 103b2, -621/281 | 19 | 26 | 6 | 101/9 | 42 |
| 11952, 287/324 | 20 | 26 | 3 | 101/9 | 7 |
| 07ca9, -3267/2257 | 18 | 25 | 7 | 101/8 | 62 |
| 103b2, -550/521 | 19 | 25 | 4 | 101/8 | 26 |
| 08f72, -2318/403 | 18 | 25 | 4 | 101/8 | 18 |
| 07ca9, 49/353 | 20 | 25 | 2 | 101/8 | 11 |

For **103b2 at 877/781**, seed call 1 exposes M19 after full-cloud replay.
Complementary calls 9, 20, 21, 44 and 67 give bounds 20, 21, 24, 26 and 27.
At calls 21 and 44, independently replayed reconciliation adds two and one
directions respectively before rebuilding the landscape. Its three search
segments use 21 + 23 + 56 = 100 calls. The remaining 33 calls after the gain to
27 add nothing certified. This is a prospective bounded cascade without
supplied exceptional points or retrospective target labels, not a rank-27
upper bound or a public-novelty claim.

### Late-gain watchlist

This is the separately frozen last-gain-at-least-75 criterion, not a new
selection rule. These results justify retaining gain timing alongside rank;
they do not prove that extending the budget would succeed.
The table uses the primary ranking order; the raw controller report lists
the same six IDs in roster order.

| Family, parameter | Final bound | Seed calls | Last complementary gain | Subsequent no-gain calls |
| --- | ---: | ---: | ---: | ---: |
| 08234, 2695/232 | 24 | 9 | 92 | 8 |
| 08f72, 1245/2519 | 23 | 5 | 79 | 21 |
| 11952, 443/2529 | 21 | 4 | 89 | 11 |
| 074d9, 569/530 | 20 | 15 | 100 | 0 |
| 074d9, -3921/797 | 20 | 24 | 96 | 4 |
| 11952, -137/2327 | 19 | 62 | 81 | 19 |

### Height and fibration comparisons

| Parameter-height stratum | Seeds | First-call seeds | Mean final lower bound | Seed calls | Complement calls | Seed-cloud / complementary directions |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| At most 1024 | 27/30 | 18 | 332/15 = 22.133… | 368 | 2684 | 56 / 98 |
| 1025–4096 | 21/30 | 8 | 99/5 = 19.8 | 981 | 2057 | 29 / 55 |

| Fibration | Seeds | Best lower bound | Seed-cloud directions | Complementary directions |
| --- | ---: | ---: | ---: | ---: |
| 074d9 | 9/10 | 26 | 10 | 29 |
| 07ca9 | 10/10 | 25 | 15 | 36 |
| 08234 | 5/10 | 24 | 16 | 6 |
| 08f72 | 9/10 | 25 | 13 | 28 |
| 103b2 | 5/10 | 27 | 9 | 25 |
| 11952 | 10/10 | 26 | 22 | 29 |

All four bounds at least 26 occur in the lower-height stratum, while five of
the six late-gain cases occur in the higher-height stratum. In 103b2, all five
lower-height candidates seeded and all five higher-height candidates missed
within the seed allowance. These are descriptive facts about a score-selected
finite panel, not a randomized estimate or a causal effect of parameter height.

The operational lessons are distinct: inexpensive seeds can expose several
directions at once; fresh fibres can exhibit a substantial adaptive cascade;
easy seed incidence does not by itself identify the largest final lower bound;
and some modest-rank fibres are still gaining at the cutoff. For example,
08234 at -3248/79 reaches 24 entirely from its seed cloud, then gains nothing
in 100 complementary calls. A follow-up would therefore compare the leading
27/26 states with the separate late-gain cases. **No such follow-up was launched.**
Exact ranks, conductor improvements, rank 28 or higher on these fibres, and
novelty beyond the frozen internal exclusion snapshots remain open.

- [Complete results and self-contained point packets](../../artifacts/generated-results/elliptic-curves/r17_60_panel_results_v1.json), SHA256 `7744dde0e5b5fb1ec27a1f96ec3c65901b73e18f626f9c203877ab8422bdce99`.
- [Replayed analysis, all 60 trajectories and rankings](../../artifacts/generated-results/elliptic-curves/r17_60_panel_analysis_v1.json), SHA256 `a395bcb72fd339ef768fa025569df011ebe901072b1d78e0c609d83959777a11`.

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

The four-worker detached controller launched and completed on 2026-09-09.
Its terminal state is `COMPLETE_BOUNDED_R17_60_PANEL`, with no active workers.
The zero-chart receipt certifies the inputs, not the later gains; it remains
at `r17-60-panel-v1/preflight.json` under the local evidence root below.

Raw checkpoints and per-phase supervisor logs are retained under
`artifacts/local/elliptic-curves/r17-60-panel-v1/`. The sources, projected
selection inputs, roster, atlas and exact maximum-class catalogue are bound
before execution. Completed case packets include full independent point
certificates and evidence hashes. The final report is exported as
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

Final checks passed: `PASS_COMPLETE_PANEL_AUDIT` with fresh Sage proofs for
all 60 packets; `R17_60_EXPORTED_CERTIFICATES_PASS|cases=60|charts=0`;
analysis `--write` immediately followed by `--check`; and all 17 targeted
regression tests. None of these checks launched point searches.

Implementation: [controller](../cas/run_r17_60_panel.py),
[arithmetic adapters](../cas/r17_60_arithmetic.py), and
[retained-input regression](../cas/check_r17_60_adapters.py).
