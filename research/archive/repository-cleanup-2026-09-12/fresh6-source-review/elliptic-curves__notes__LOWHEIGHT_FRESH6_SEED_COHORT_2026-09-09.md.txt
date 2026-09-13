# Lower-height six-fibre seed cohort

Five of six selected native R17 fibres first supplied certified M18 seeds.
Winning-parent amplification, complementary-parent passes and saved-cloud
reconciliation now certify lower bounds 22, 24, 22, 25 and 23 on those five
fibres (family order 074d9, 07ca9, 08f72, 103b2, 11952).
One bounded seed miss retains M17.
No exact rank, public novelty or record is claimed.

Selection restricts the saved6144-row score population to parameter height
at most1024, retaining the unchanged per-family score order and frozen
rational-isomorphism/address exclusions.366 eligible addresses remain in
this stratum; one per family is selected. Selection and replay take1.109
and1.108seconds. Generic17 preparation/replay take18.896 and13.847seconds.

| Family | Parameter | Initial seed lower bound | Seed calls |
| --- | --- | ---: | ---: |
| 074d9 | 25/3 | 18 | 1 |
| 07ca9 | 142/561 | 18 | 1 |
| 08234 | 383/388 | 17 | 86 |
| 08f72 | 254/17 | 18 | 1 |
| 103b2 | -95/728 | 18 | 1 |
| 11952 | 238/27 | 18 | 1 |

The unchanged exact-maximum-class first-M18 constructor completes91 total
point calls without point timeouts in82.647seconds. Each search has an1800-second/
2GiB cap, each map5seconds/1GiB and each point call10seconds at height125000.
All six full geometry/map/point replays pass. A separate packet checker
reconstructs the17 native sections, finite independence certificates and
pairwise rational nonisomorphism in7.563seconds. All five M18 packets have
completed their winning-parent amplification pass. Full winning witnesses
are retained.

This selected six-case experiment does not establish a general advantage
for lower parameter height or a universal first-chart seed rule.

- [Selector](../cas/select_fresh6_lowheight_r17.py).
- [Exact seed packets](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_first_m18_v1/result.json).
- [Five-seed queue](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_first_m18_v1/amplification-queue.json).
- [Independent packet checker](../cas/verify_fresh6_seed_packets.py).
- [Reservation ledger](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_reserved_results_v1.json).

Raw evidence is in `artifacts/local/elliptic-curves/fresh6-lowheight-selection-v1/`,
`fresh6-lowheight-generic-seeds-v1/` and `fresh6-lowheight-seed-confirmation-v1/`.

## Amplification and preparation for future passes

All five winning-parent runs pass independent geometry/map/point replay.
Their full-cloud audits require reconciliation before further searching.
Replaying the retained clouds and checking standalone finite independence
certificates gives the following subgroup lower bounds:

| Family | Parameter | Retained V3 basis | Reconciled basis | New point calls |
| --- | --- | ---: | ---: | ---: |
| 074d9 | 25/3 | 19 | 22 | 1 |
| 07ca9 | 142/561 | 19 | 20 | 1 |
| 08f72 | 254/17 | 19 | 22 | 1 |
| 103b2 | -95/728 | 21 | 22 | 3 |
| 11952 | 238/27 | 19 | 23 | 1 |

The five searches take 64.648 seconds in total; independent V3 replay takes
28.065 seconds. Reconciliation takes 18.996 seconds and its second exact
replay takes 18.996 seconds. It adds twelve independent directions beyond
the retained terminal bases, with zero additional point searches. These
measurements describe this cohort, not a general speedup factor.

Complementary banks contain at most sixteen exact maximum generic parity
classes outside the original winning-parent span. Their exact CVP values
are checked with both the integer and rational solvers. New preparations
retain these banks and use the reconciled bases, verifying each finite
certificate again. Preparation performs no point search.

The first bounded complementary-parent passes cover family 11952 at
238/27 and family 074d9 at 25/3. All five first complementary-parent passes have now been replayed.
This is bounded exposure, not a record prediction.

Future execution uses the existing optimized workers: cached finite-prime
admission, integer exact CVP, separately bounded lean map construction,
completed-box equivalence deduplication, and checkpointed suffix
continuation. After a full-cloud audit requests reconciliation, replay and
certify that cloud before constructing another search stage. This preserves
the unchanged V3 policy while avoiding searches from an incomplete retained
basis. Use bounded 100-call passes and independent replay; a miss remains
UNKNOWN. No old search campaign needs to be repaired or restarted.

- [Prepared continuation queue](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/continuation-queue.json).
- [Saved-cloud certificate checker](../cas/reconcile_verified_v3_cloud.py).
- [Complementary-parent preparation](../cas/prepare_exact_maximum_complement_parents.py).
- [Reconciled-basis preparation](../cas/prepare_reconciled_v3_continuation.py).

Raw amplification evidence is retained under
`artifacts/local/elliptic-curves/fresh6-lowheight-amplification-v1/`;
compact result and reconciliation packets are in
`artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/`.
No unconditional exact rank, conductor record or public novelty is claimed.

A bounded conductor audit of 11952 at 238/27 reaches its 30-second
factorization limit (31.059 seconds including termination). Certified local
data give lower bound 4573532994030 and retain an unresolved cofactor.
The exact conductor and its record comparison remain UNKNOWN. The
[bounds and timeout receipt](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/11952-conductor-bounds.json)
are retained; no factorization escalation is scheduled.

## First complementary-parent exposure

The first two 100-call passes both complete and pass independent full replay:

| Family | Parameter | Subgroup before/after | Search seconds | Replay seconds |
| --- | --- | --- | ---: | ---: |
| 11952 | 238/27 | 23 → 23 | 107.807 | 66.144 |
| 074d9 | 25/3 | 22 → 22 | 122.651 | 49.774 |

Neither pass has a point or map timeout. These are bounded misses, not
rank upper bounds. The 11952 suffix resumes at centre 50 of 636, leaving
586 centres; the 074d9 suffix also resumes at centre 50. Exact cursor
receipts preserve both completed prefixes. The remaining three first passes are recorded below.

- [11952 result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/11952-M23-complement/result.json).
- [074d9 result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/074d9-M22-complement/result.json).

## Remaining first passes and two gains

All three remaining searches and independent geometry/map/point replays pass.
There are no point or map timeouts.

| Family | Parameter | Certified subgroup before/after | Calls | Search seconds | Replay seconds |
| --- | --- | --- | ---: | ---: | ---: |
| 08f72 | 254/17 | 22 → 22 | 100 | 99.690 | 55.037 |
| 103b2 | -95/728 | 22 → 25 | 100 | 182.604 | 249.277 |
| 07ca9 | 142/561 | 20 → 24 after reconciliation | 1 | 12.563 | 10.460 |

103b2 gains on cumulative calls3,14,23, then completes77 calls at M25.
Its682-centre final stage resumes at centre38, factor_free policy, with
that centre partially unsearched. The standalone M25 finite independence
certificate and all four adaptive stages replay exactly. The 08f72 suffix
resumes at centre50 of486.

07ca9 retains M21 from its first call, then correctly stops because the
complete-cloud audit sees additional directions. Saved-witness reconciliation
certifies M24, adding three points without another point search. Construction
and second reconciliation replay take3.511 and3.507 seconds. A new sealed
preparation uses this M24 basis with the same complementary-parent bank.

- [103b2 M25 result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/103b2-M22-complement/result.json).
- [07ca9 replayed search](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/07ca9-M20-complement/result.json).
- [07ca9 M24 reconciliation](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/07ca9-M20-complement/reconciled.json).
- [08f72 bounded result](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/08f72-M22-complement/result.json).

The 103b2 conductor audit also reaches its30-second factorization cap
(31.042 seconds including termination), leaving exact conductor UNKNOWN.
Its [local bounds and timeout](../../artifacts/generated-results/elliptic-curves/fresh6_lowheight_amplification_v1/103b2-conductor-bounds.json)
are retained. Neither positive subgroup result proves an exact rank,
conductor record or public novelty.

The next higher-rank candidate is curve48's certified M27 with its newly
prepared exact norm12 parents, documented in the
[parent reassessment](PRODUCTIVE_PARENT_SPAN_REASSESSMENT_2026-09-09.md).
The five low-height continuations remain checkpointed for later scheduling.
