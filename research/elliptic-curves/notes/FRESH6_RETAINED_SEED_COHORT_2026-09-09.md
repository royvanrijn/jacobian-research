# Fresh seeds from the retained compact population

The fixed six-fibre cohort now supplies **four independently certified M18
seeds**. The other two retain lower bound17 after bounded misses. Selection
used the frozen equation/roster snapshots below; no global public novelty,
exact rank or record is claimed.

| Family | Parameter | Position in unchanged family score order |
| --- | --- | ---: |
| 074d9 | -2707/3437 | 48 |
| 07ca9 | 587/1074 | 48 |
| 08234 | 2071/4060 | 47 |
| 08f72 | 4032/1663 | 46 |
| 103b2 | 3889/2172 | 46 |
| 11952 | -436/3149 | 46 |

## Selection boundary

The [selector](../cas/select_fresh6_retained_r17_v2.py) uses the6144 saved
H4096 rows from the original retained, discarded-shard and retention512
extended-prime score files. It keeps descending selection-band S1 units,
descending good-prime count, denominator and signed numerator. It chooses
one curve per family. Validation-band columns are removed from selected rows;
measured ranks and points do not enter ordering or specialization.

The frozen exclusions contain1764 distinct equation tuples from compact
result ledgers and previous-equation lists, endpoints and point-run rosters.
The rosters also supply1146 scheduled family/parameter addresses, including
addresses not in this compact pool. Exact rational isomorphism tests exclude
matching equations, and the six selected equations are pairwise nonisomorphic.
At selection time5828 address rows remain eligible; this is an address count,
not a count of independent curves or a predicted high-rank density.

Coverage is explicitly relative to the named snapshots. This is not a scan of
every raw file or every external search. No new public catalogue is loaded.
The first selector stopped on a singular endpoint in the exclusion list;
v2 skips that nonelliptic equation. Its failed v1 evidence remains preserved.

Selection completes in2.048 seconds and exact ordering/exclusion replay in
2.101 seconds, each under120 seconds/1GiB. Raw evidence is in
`artifacts/local/elliptic-curves/fresh6-retained-selection-v2/`. Replay checks
the fixed data and selection; it does not recompute all retained prime traces.

## Certified generic inputs

[Preparation](../cas/prepare_fresh6_generic_seeds.py) specializes exactly the17
atlas sections, checks equality with each selected equation, proves finite
quotient independence, and verifies a rational2-torsion exclusion witness.
All six packets pass an independent exact replay. Preparation takes20.257
seconds and replay14.526 seconds, each under120 seconds/1GiB, without Sage
startup or a point search. Each lower bound is17; no additional direction is
included or inferred.

- [Six replayed packets](../../artifacts/generated-results/elliptic-curves/fresh6_generic17_packets_v1/result.json).
- [Equation-only reservation ledger](../../artifacts/generated-results/elliptic-curves/fresh6_reserved_seed_results_v1.json), so later selections exclude this queued cohort.
- Raw packets: `artifacts/local/elliptic-curves/fresh6-generic-seeds-v1/`.

## First-extra-point confirmation

The bounded constructor uses all exact43/49 generic maximum classes, specialized
exact CVP representatives, lean bounded map workers, finite-group admission
caches and exact signed-permutation box deduplication. It stops at the first
standalone M18 certificate and preserves the full winning chart witness.

| Family | Certified lower bound | Point calls | Search seconds |
| --- | ---: | ---: | ---: |
| 074d9 | 18 | 47 | 38.271 |
| 07ca9 | 18 | 1 | 4.442 |
| 08234 | 17 | 86 | 74.747 |
| 08f72 | 18 | 35 | 27.326 |
| 103b2 | 18 | 2 | 5.016 |
| 11952 | 17 | 98 | 84.951 |

All269 point calls completed without timeout; total search time was234.753
seconds. Each search had a1800-second/2GiB supervisor, each map5seconds/1GiB,
and each point call10seconds at height125000. A bounded miss excludes nothing
beyond this experiment.

All six geometry/map/point replays passed. The first replay encountered a
JSON tuple/list comparison mismatch at the final terminal comparison; v2
normalizes that comparison. The original source, failure and search receipts
remain preserved, and the same search replayed successfully without rerunning
point construction. The other five used v2 throughout.

A separate pure-Python packet checker independently reconstructs all native17
sections, verifies the finite independence certificates and pairwise rational
nonisomorphism. It passed in8.024seconds under120seconds/1GiB.

- [Self-contained exact packets](../../artifacts/generated-results/elliptic-curves/fresh6_first_m18_v1/result.json).
- [Independent checker](../cas/verify_fresh6_seed_packets.py).
- [Four-seed amplification queue](../../artifacts/generated-results/elliptic-curves/fresh6_first_m18_v1/amplification-queue.json).
- Raw search/replay receipts: `artifacts/local/elliptic-curves/fresh6-seed-confirmation-v1/`.

## First productive-parent V3 amplification

All four M18 seeds completed a bounded adaptive pass. The preparation independently
rechecks the native basis, seed proof and winning point witness, then verifies
its generic parent CVP with integer and rational solvers. Each winning parent
has doubled generic norm24. V3 retains its extension scoring and immediate
rebuild rule, with lean separately bounded maps, cached finite admission,
100 total point calls maximum,1800seconds/3GiB per search and replay.

| Family | Final certified lower bound | Point calls | Search seconds | Replay seconds |
| --- | ---: | ---: | ---: | ---: |
| 074d9 | 18 | 6 | 12.624 | 3.658 |
| 07ca9 | 19 | 9 | 23.622 | 6.786 |
| 08f72 | 18 | 6 | 12.142 | 3.718 |
| 103b2 | 18 | 6 | 12.043 | 3.765 |

The07ca9 fibre at587/1074 gains its19th certified direction on chart3;
after rebuilding, six further charts yield no new certified direction.
All four complete their finite productive-parent policy and pass independent
geometry, map, point, gain and complete-cloud replay. These are subgroup
lower bounds; exhausting these small parent subsets proves no rank upper bound.
No record or global novelty is asserted.

- [Preparation source](../cas/prepare_fresh6_productive_seed.py).
- [M19 run evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-fresh-001/result.json).
- [All four results](../../artifacts/generated-results/elliptic-curves/fresh6_first_m18_v1/amplification-queue.json).

## Exact maximum complementary parents: M19 to M23

For07ca9 at587/1074, the next preparation rechecks all42 catalogue maximum
classes outside the original winning-parent span. Both exact CVP solvers
agree. Ascending-mask selection first extends that span and then fills to16
parents, all of doubled norm24; their union with the original parent spans
16 generic parity directions. Preparation takes5.638seconds and uses no
specialized higher-point labels. The same preparation is complete for the
other three fresh fibres, whose complementary searches have not started.

A bounded V3 pass from M19 finds retained gains19→20,20→21 and21→22
in1,5 and1 point calls respectively. It stops after7 calls because the complete
cloud audit certifies23 directions. Search takes38.598seconds; full independent
geometry/map/point/cloud replay passes in35.395seconds. The frozen M22 terminal
and reconciliation stop remain intact.

A separate arithmetic-only reconciliation replays the seven saved chart
witnesses, keeps the terminal22-point basis and admits one additional point
from the final chart. Its standalone finite certificate proves **M23**.
Construction takes6.050seconds; a second complete reconciliation agrees in
5.998seconds, both under120seconds/1GiB. No further point search was needed.
This is a subgroup lower bound, with no exact rank or record claim.

- [Exact complementary-parent preparation](../cas/prepare_exact_maximum_complement_parents.py).
- [Full bounded-run evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-exact-complement/result.json).
- [Reconciled M23 packet](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-exact-complement/reconciled.json).
- [Reconciliation checker](../cas/reconcile_verified_v3_cloud.py).
- [M23 continuation queue](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-exact-complement/continuation-queue.json).

## Rebuilt M23 pass

A fresh sealed preparation rechecks the reconciled23-point certificate and
retains the16 exact complementary parents. The next V3 pass completes100
point calls with no new certified direction and no point or map timeouts.
Search takes117.606seconds; full independent replay takes72.360seconds.
The terminal is `CHART_BUDGET_EXHAUSTED`, preserving the unvisited centre
suffix. This bounded miss does not establish an upper bound.

- [Reconciled continuation preparation](../cas/prepare_reconciled_v3_continuation.py).
- [Replayed M23 pass](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-M23/result.json).
- [Exact suffix cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/07ca9-M23/suffix-queue.json).

The next portfolio action is the prepared complementary-parent pass on another
fresh fibre; the M23 suffix remains available for a cached continuation.
The last gaining chart demonstrates why future continuations must retain all
returned points even when admission stops after the first gain. Existing
curve48 and curve92 suffix queues remain available.

## Fresh074d9 complementary-parent amplification: M18 to M24

The next prepared fibre,074d9 at-2707/3437, completes its100-call pass
with six additional certified directions. Successive epochs raise
18→19→20→21→22→23→24, with immediate rebuilding after each gain.
All point calls and bounded maps complete without timeout. Search takes
180.879seconds; full independent rational-CVP, exact map/point, gain-provenance
and complete-cloud replay passes in171.694seconds.

The terminal is `CHART_BUDGET_EXHAUSTED`, not a rank upper bound. The exact
unvisited suffix is queued for a cached continuation. The fibre has a certified
rank lower bound24; no exact rank, public novelty or conductor record is claimed.

- [M24 evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/074d9-exact-complement/result.json).
- [Continuation cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/074d9-exact-complement/suffix-queue.json).

The cohort now has certified lower bounds24,23,17,18,18,17 in family order
074d9,07ca9,08234,08f72,103b2,11952. Complementary-parent preparations for
08f72 and103b2 remain unsearched.

### Cached M24 continuation

The first cached continuation completes100 additional point calls,200 cumulative,
without a new certified direction or point/map timeout. Search takes97.675seconds
and independent cached replay17.337seconds. The retained lower bound remains24.
Previously verified landscapes and receipts are reused, with no repeated point
search on inherited charts. The exact next cursor is centre58,
`preconditioned_full`, leaving621 centres in the final679-centre landscape.
The earlier queue points to this verified descendant, preventing accidental
repetition of that suffix prefix. This bounded miss gives no rank upper bound.

- [Cached M24 evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/074d9-exact-complement-cached-v1/result.json).
- [Updated continuation cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/074d9-exact-complement-cached-v1/suffix-queue.json).

Next test the still-unsearched complementary parents on08f72 before adding
another suffix budget to this fibre.

## Fresh08f72 complementary-parent amplification: M18 to M24

At08f72 parameter4032/1663, the prepared complementary-parent pass finds six
additional certified directions, raising the native M18 basis to M24 within
100 point calls. Every gain triggers immediate rebuilding. No point or map
timeouts occur. Search takes182.502seconds; full independent rational-CVP,
exact map/point, gain-provenance and complete-cloud replay passes in181.846seconds.
The terminal is `CHART_BUDGET_EXHAUSTED`; the unvisited suffix is preserved.
No exact rank, rank upper bound, public novelty or conductor record is claimed.

- [M24 evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/08f72-exact-complement/result.json).
- [Exact continuation cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/08f72-exact-complement/suffix-queue.json).

Current certified lower bounds in the fixed cohort are24,23,17,24,18,17
in family order074d9,07ca9,08234,08f72,103b2,11952. The103b2
complementary-parent preparation is the last of the four M18-derived branches
yet to receive its first broad-parent pass.

## Fresh103b2 complementary-parent amplification: M18 to M25

At103b2 parameter3889/2172, seven exact rational point gains extend the
native M18 basis to M25 within100 point calls. Immediate rebuilding follows
each gain. No point or map timeouts occur. Search takes234.161seconds;
full independent rational-CVP, exact map/point, gain-provenance and
complete-cloud replay passes in341.550seconds. The terminal is
`CHART_BUDGET_EXHAUSTED`, preserving its unvisited suffix. No exact rank,
rank upper bound, public novelty or conductor record is claimed.

- [M25 evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement/result.json).
- [Exact continuation cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement/suffix-queue.json).

All four first complementary-parent branches are now replayed. The fixed
six-fibre cohort has certified lower bounds24,23,17,24,25,17 in family order
074d9,07ca9,08234,08f72,103b2,11952. The two M17 cases remain bounded
seed misses. The best next amplification candidate is the new M25 branch;
its next pass can reuse the independently verified landscapes and search
previously unvisited boxes. This small cohort does not establish a general
success rate or a new rank/conductor record.

### Cached M25 continuation

The first cached103b2 continuation completes100 additional point calls,
200 cumulative, with no new certified direction and no point/map timeouts.
Search takes98.878seconds; independent cached replay takes19.420seconds.
The lower bound remains25. The exact next cursor is centre76,
`preconditioned_full`, leaving626 of the702 final-epoch centres unvisited.
The prior queue points to this verified descendant. This bounded miss
provides no rank upper bound.

- [Cached M25 evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement-cached-v1/result.json).
- [Updated continuation cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement-cached-v1/suffix-queue.json).

## Bounded post-discovery conductor audit

Four30-second/1GiB audits compute exact local reduction at bad primes up to1000,
checkpoint certified conductor bounds and attempt proof-enabled factorization
of the remaining discriminant. This does not affect candidate selection or rank.

For08f72 at4032/1663, factorization and all remaining local data complete.
A second run reproduces the exact conductor:

```
5010188826156411226384182822992530591002057751473896628642100179406951935900788324249280619427023331755172627345571322858912720052190
```

This exceeds the pinned catalogue's rank-at-least24 benchmark (curve548,
conductor361660950250291954113011326751617907627414070443777168347285410734060239220929310535800). It is therefore not
an improvement over that benchmark. No current external catalogue claim is made.

The074d9,07ca9 and103b2 audits hit their strict wall caps with residual
cofactors of124,104 and113 digits respectively. Their saved bounds do not
establish a record or exclude one; exact conductors remain unknown. No larger
factorization campaign is launched. This evidence favours continued rank
amplification of M25 over conductor-directed effort on08f72.

- [All bounds, exact result and pinned comparison](../../artifacts/generated-results/elliptic-curves/fresh6_conductor_audit_v1/result.json).
- [Proof-enabled local conductor audit](../cas/audit_fresh6_conductors.sage).

### Second cached M25 continuation

The next103b2 suffix pass completes100 additional point calls,300 cumulative,
without a new certified direction or point/map timeout. Search takes101.463seconds;
independent cached replay takes20.098seconds. The lower bound remains25.
The exact next cursor is centre126, `preconditioned_full`, leaving576
of the702 final-epoch centres unvisited. The previous queue points to this
verified descendant. No upper bound follows from these bounded misses.

- [Second cached M25 evidence](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement-cached-v2/result.json).
- [Updated continuation cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-exact-complement-cached-v2/suffix-queue.json).

## Disjoint maximum-parent pass from M25

After two no-gain cached suffix passes, a different parent bank is prepared
from the26 remaining exact maximum generic classes, excluding every mask in
the previous16-parent bank and the original productive class. Both exact CVP
solvers recheck all26 minima; deterministic selection retains16 parents,
extending the original productive span to dimension16. Preparation takes
4.755seconds. These are different generic classes, not a claim of independent
new generic directions beyond the previous bank's entire span.

Starting from the certified M25 basis, the new-parent V3 pass completes100
point calls without a new certified direction or point/map timeout. Search
takes115.838seconds; full independent replay passes in124.490seconds.
The terminal is `CHART_BUDGET_EXHAUSTED`. Its unvisited suffix remains queued,
and all earlier branch receipts remain intact. No rank upper bound follows.

- [Disjoint-parent preparation](../cas/prepare_disjoint_maximum_parents.py).
- [Replayed new-parent pass](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-disjoint-max/result.json).
- [New-parent suffix cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-disjoint-max/suffix-queue.json).

## Lower generic shell coverage from M25

A fixed128-parity sample outside the original productive span is evaluated
with both exact CVP solvers. Its largest sampled doubled minimum is20,
below the maximum-class doubled minimum24. Deterministic selection retains
16 norm10 parents and, with the original productive class, spans all17 generic
parity directions. Preparation takes4.286seconds with no point search.

An exact finite linear-algebra audit shows that all43 maximum classes of103b2
span a16-dimensional subspace, annihilated by parity mask45903. Thirteen of
the16 selected norm10 parents lie outside that hyperplane. This is a statement
about generic parity coverage; it is neither a specialized rank obstruction
nor evidence of rational points in the newly covered classes.

The resulting M25 V3 pass completes100 point calls without a new certified
direction or point/map timeout. Search takes119.167seconds; full independent
replay passes in114.541seconds. Its unvisited suffix is queued. The lower
bound remains25; no rank upper bound or record is established.

- [Sampled continuation preparation](../cas/prepare_sampled_continuation_parents.py).
- [Exact parity-span certificate](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2_sampled_parent_span.json).
- [Replayed lower-shell pass](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-sampled-shell/result.json).
- [Lower-shell suffix cursor](../../artifacts/generated-results/elliptic-curves/fresh6_amplification_v1/103b2-sampled-shell/suffix-queue.json).

## Outside-maximum seed attempt on08234

The two seed misses have different maximum-class parity coverage:08234's
classes span16 generic parity directions, while11952's span all17. A fixed
128-parity sample outside08234's maximum-class span is evaluated with both
exact CVP solvers. The16 selected classes have doubled norm20. This covers
parities absent from the maximum-class seed search; it predicts no point.

At08234 parameter2071/4060, a separately frozen first-M18 constructor tests
all32 bounded map/point boxes without a new certified direction or point
timeout. Preparation takes1.842seconds, search26.694seconds and independent
geometry/map/point replay3.984seconds. The lower bound remains17; the original
maximum-class miss and the new lower-shell miss are preserved separately.
No absence or rank upper bound follows.

- [Frozen generic sample and exact minima](../../artifacts/generated-results/elliptic-curves/fresh6_outside_maximum_seed_masks_v1.json).
- [Sample construction](../cas/prepare_outside_maximum_seed_masks.sage).
- [Bounded constructor](../cas/run_outside_maximum_seed_confirmation.py).
- [Replayed result](../../artifacts/generated-results/elliptic-curves/08234_outside_maximum_seed_v1/result.json).

The [second fresh cohort](SECOND_FRESH6_SEED_COHORT_2026-09-09.md) now supplies five further certified M18 seeds for amplification.
